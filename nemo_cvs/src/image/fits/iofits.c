
/*---------------------------------------------------------------------------
   
   File name	:	iofits.c
   Author	:	N. Devillard / P. Teuben (NEMO adaptation)
   Created on	:	July 1998    / April 2003
   Description	:	convert FITS from a pixel depth to another
			no rescaling is done
   Notice	:	now needs NEMO to compile

   About this program:

   Yes, it is a hack. However, it is a real boost compared to any
   "intelligent" FITS handling software because it does not even try to
   understand what it finds in the FITS header (except for BITPIX and
   NAXIS?, of course). The worst part is pixel to pixel conversion,
   because it required to study it case by case. Nice feature is that it
   handles Intel/Motorola endian-ness dynamically, no need to specify
   anything on the command-line nor at compilation time.

   Conversions almost always use the fastest bitwise method to construct
   the outputs. Heavy use of mmap() allows to go on without allocating
   much memory, I do believe no method could do it faster than direct
   transfer from disk to disk with a shortcut through CPU registers.

   Unsupported are:

   * Machines on which the size of a char is not 8 bits. It would mean a
     complete rewrite of all pixel to pixel conversions. No way.
   * Machines for which a float or a double are not stored as described
     in the IEEE reference (cannot find the exact number... well, the
     usual floats and doubles as described in the FITS format for
     example). It would also mean re-writing most of the conversion
	 functions, forget it.

	Some assumptions taken along:

	* FITS lines are 80 chars wide.
	* A FITS block size is 2880 bytes wide.
	* FITS data types are 8 16 32 -32 and -64 ONLY.
	* The size of an IEEE float is 4 chars.
	* The size of an IEEE double is 8 chars.

	* The size of a short is 16 bits
	* The size of a long is 32 bits

	The two latter are ensured by redefining the following int types:

	short16 is a signed integer on 16 bits
	long32 is a signed integer on 32 bits

 *--------------------------------------------------------------------------*/

/*
	$Id: iofits.c,v 1.2 2003/04/30 04:57:45 pteuben Exp $
	$Author: pteuben $
	$Date: 2003/04/30 04:57:45 $
	$Revision: 1.2 $
*/

/*---------------------------------------------------------------------------
								Includes
 ---------------------------------------------------------------------------*/
#include <stdinc.h>    /* NEMO */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h>


/*---------------------------------------------------------------------------
								Defines
 ---------------------------------------------------------------------------*/

#define NM_SIZ			512
#define FITS_BLSZ		2880
#define ONE_MEGABYTE	(1024 * 1024)

#ifndef SEEK_SET
#define SEEK_SET	0
#endif

#define ENDIAN_UNKNOWN      -1
#define ENDIAN_MOTOROLA      0
#define ENDIAN_INTEL         1


/*---------------------------------------------------------------------------
								Macros
 ---------------------------------------------------------------------------*/

#define ENDIAN_NESS(f)          ((f)&1)
#define IEEE_FLOAT_COMPAT(f)    ((f)&2)
#define IEEE_DOUBLE_COMPAT(f)   ((f)&4)


/*---------------------------------------------------------------------------
								typedefs
 ---------------------------------------------------------------------------*/

typedef short			short16 ;
typedef unsigned short	ushort16 ;
typedef int				long32 ;
/* typedef unsigned char	byte ;                      /* same as NEMO */

/*---------------------------------------------------------------------------
							Function prototypes
 ---------------------------------------------------------------------------*/

/*
 * All functions are private
 */

static int bits_per_byte(void);
static int convert_fits_pixel_depth(char*, char*, int) ;
static int get_FITS_header_size(char*);
static int filesize(char*) ;
static int get_bitpix(char*) ;
static int set_bitpix(char*, int) ;
static int get_npix(char*) ;
static int write_pixbuf(char*, char*, int, int, int) ;
static int test_machine_formats(void) ;
static void swap_bytes(void*, int);
static void check_all_types(int) ;


/*
 * Global variable: numerical compatibility flags
 */

static unsigned compat_flags ;



/*---------------------------------------------------------------------------
								main()	
 ---------------------------------------------------------------------------*/

string defv[] = {
  "in=???\n           Input FITS file",
  "out=???\n          Output FITS file",
  "bitpix=???\n       requested output BITPIX value (8|16|32|-32|-64)",
  "VERSION=1.0\n      29-apr-03 PJT",
  NULL,
};

string usage = "convert FITS from a pixel depth to another without rescaling";

int nemo_main()
{
  int	p = getiparam("bitpix");
  char *i = (char *)getparam("in");
  char *o = (char* )getparam("out");
  
  check_all_types(p);
  if (streq(i,o))
    error("Cannot convert file to itself, change in= and/or out= ");
  return convert_fits_pixel_depth(i,o,p);
}


/*---------------------------------------------------------------------------
   Function	:	convert_fits_pixel_depth()
   In		:	name in, name out, pixel type out
   Out		:	int 0 if Ok, anything else otherwise
   Job		:	convert FITS files from one pixel depth to another
   Notice	:	heavy use of mmap() to speed up the process
 ---------------------------------------------------------------------------*/

static int
convert_fits_pixel_depth(
	char	*	name_in,
	char	*	name_out,
	int			ptype_out
)
{
	int			fd_in,
				fd_out ;
	char	*	buf_in ;
	char	*	buf_out ;
	int			fsize_in,
				fsize_out,
				header_size,
				padd_size ;
	int			ptype_in ;
	int			bpp_out ;
	int			npix ;
	char	*	zero ;
	int			yet_to_dump ;
	int			buf_dump ;
	
	/*
	 * Open input file and get information we need:
	 * - pixel depth
	 * - number of pixels
	 */

	fsize_in = filesize(name_in) ;
	header_size = get_FITS_header_size(name_in) ;
	if ((fd_in = open(name_in, O_RDONLY)) == -1) {
		fprintf(stderr, "cannot open file %s: aborting\n", name_in) ;
		return -1 ;
	}
	buf_in = (char*)mmap(0, fsize_in, PROT_READ, MAP_SHARED, fd_in, 0) ;
	if (buf_in == (char*)-1) {
		perror("mmap") ;
		fprintf(stderr, "cannot mmap file: %s\n", name_in) ;
		close(fd_in) ;
		return -1 ;
	}
	ptype_in = get_bitpix(buf_in) ;
	if (ptype_in == 0) {
		close(fd_in) ; munmap(buf_in, fsize_in) ;
		return -1 ;
	}
	if (ptype_in == -32) {
		if (!IEEE_FLOAT_COMPAT(compat_flags)) {
			fprintf(stderr,
				"this machine does not support IEEE floating point values\n");
			fprintf(stderr,
				"cannot convert from -32 type\n");
			exit(-1) ;
		}
	}
	if (ptype_in == -64) {
		if (!IEEE_DOUBLE_COMPAT(compat_flags)) {
			fprintf(stderr,
				"this machine does not support IEEE double values\n");
			fprintf(stderr,
				"cannot convert from -64 type\n");
			exit(-1) ;
		}
	}

	/*
	 * Compute the size of the output file
	 * same header size, + pixel area + blank padding
	 */
	npix = get_npix(buf_in) ;
	if (npix < 1) {
		close(fd_in) ;
		munmap(buf_in, fsize_in) ;
		return -1 ;
	}
	bpp_out = ptype_out / 8 ;
	if (bpp_out < 0) {
		bpp_out = - bpp_out ;
	}
	fsize_out = header_size + npix * bpp_out ;
	if ((fsize_out % FITS_BLSZ) == 0) {
		padd_size = 0 ;
	} else {
		padd_size = FITS_BLSZ - fsize_out % FITS_BLSZ ;
	}
	fsize_out += padd_size ;

	/*
	 * Now create the output file and fill it with zeros, then mmap it.
	 * The permissions are rw-rw-r-- by default.
	 */
	if ((fd_out=creat(name_out,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH))==-1){
		perror("creat") ;
		fprintf(stderr, "cannot create file %s: aborting\n", name_out) ;
		close(fd_in) ; munmap(buf_in, fsize_in) ;
		return -1 ;
	}
	zero = calloc(ONE_MEGABYTE, 1) ; 
	yet_to_dump = fsize_out ;
	if (fsize_out < ONE_MEGABYTE) {
		buf_dump = fsize_out ;
	} else {
		buf_dump = ONE_MEGABYTE ;
	}
	while (yet_to_dump>0) {
		if (write(fd_out, zero, buf_dump) != buf_dump) {
			perror("write") ;
			fprintf(stderr, "error writing to the output file: aborting\n") ;
			close(fd_in) ; close(fd_out) ; munmap(buf_in, fsize_in) ;
			free(zero) ;
			return -1 ;
		}
		yet_to_dump -= buf_dump ;
		if (yet_to_dump > ONE_MEGABYTE) {
			buf_dump = ONE_MEGABYTE ;
		} else {
			buf_dump = yet_to_dump ;
		}
	}
	free(zero) ;
	close(fd_out) ;
	if ((fd_out = open(name_out, O_RDWR)) == -1) {
		perror("open") ;
		fprintf(stderr, "cannot reopen file %s for writing\n", name_out);
		close(fd_in) ; munmap(buf_in, fsize_in) ;
		return -1 ;
	}

	buf_out = (char*)mmap(0,
						 fsize_out,
						 PROT_READ | PROT_WRITE,
						 MAP_SHARED,
						 fd_out,
						 0) ;
	if (buf_out == (char*)-1) {
		perror("mmap") ;
		fprintf(stderr, "cannot mmap file: %s\n", name_out) ;
		munmap(buf_in, fsize_in) ; close(fd_in) ; close(fd_out) ;
		return -1 ;
	}

	/*
	 * Copy FITS header from input to output, modify BITPIX
	 */
	memcpy(buf_out, buf_in, header_size) ;
	if (set_bitpix(buf_out, ptype_out)!=0) {
		fprintf(stderr, "error writing BITPIX in output: aborting\n");
		munmap(buf_in, fsize_in) ; munmap(buf_out, fsize_out) ;
		close(fd_in) ; close(fd_out) ;
		return -1 ;
	}

	/*
	 * Now write pixels
	 */
	write_pixbuf(buf_in + header_size,
				 buf_out + header_size,
				 npix,
				 ptype_in,
				 ptype_out) ;

	/*
	 * Blank-pad the output file if needed
	 */
	if (padd_size) {
		if (lseek(fd_out, fsize_out-padd_size, SEEK_SET)==-1) {
			perror("lseek");
			fprintf(stderr, "error seeking output file: aborting\n");
		} else {
			zero = calloc(padd_size, 1) ;
			if (write(fd_out, zero, 1)==-1) {
				perror("write") ;
				fprintf(stderr, "error writing into output file: aborting\n");
			}
			free(zero) ;
		}
	}

	/*
	 * Close, unmap, goodbye
	 */
	close(fd_in) ;
	close(fd_out) ;
	munmap(buf_in, fsize_in) ;
	munmap(buf_out, fsize_out) ;
	return 0 ;
}




/*---------------------------------------------------------------------------
   Function :   get_FITS_header_size()
   In       :   FITS file name
   Out      :	int 
   Job      :   compute the size (in bytes) of a FITS header
   Notice   :   should always be a multiple of 2880. This implementation
               assumes only that 80 characters are found per line.
 ---------------------------------------------------------------------------*/


static int get_FITS_header_size(char * name)
{
    FILE    *   in ;
    char        line[81] ;
    int         found = 0 ;
    int         count ;
    int			hs ;

    if ((in = fopen(name, "r")) == NULL) {
        fprintf(stderr, "cannot open %s: aborting\n", name) ;
        return 0 ;
    }
    count = 0 ;
    while (!found) {
        if (fread(line, 1, 80, in)!=80) {
            break ;
        }
        count ++ ;
        if (!strncmp(line, "END ", 4)) {
            found = 1 ;
        }
    }
    fclose(in);

    if (!found) return 0 ;
    /*
     * The size is the number of found cards times 80,
     * rounded to the closest higher multiple of 2880.
     */
    hs = count * 80 ;
    if ((hs % FITS_BLSZ) != 0) {
        hs = (1+(hs/FITS_BLSZ)) * FITS_BLSZ ;
    }
    return hs ;
}


/*---------------------------------------------------------------------------
 * Function :   filesize()
 * In       :   filename
 * Out      :   size of the file in bytes
 * Job      :   strongly non portable. Only on Unix systems!
 * Notice   :   Looks strange, but there is no portable way to answer
 *              the question: how many bytes can I read from this file?
 *--------------------------------------------------------------------------*/

static int filesize(char *filename)
{
    int		size ;
    struct stat	fileinfo ;

    /* POSIX compliant  */
    if (stat(filename, &fileinfo) != 0) {
        size = (int)0 ;
    } else {
        size = (int)fileinfo.st_size ;
    }
    return size ;
}


/*---------------------------------------------------------------------------
   Function	:	get_bitpix()
   In		:	allocated char buffer containing the whole FITS header
   Out		:	int 8 16 32 -32 or -64
   Job		:	gets the value of BITPIX from a FITS header
   Notice	:	returns 0 if cannot find it
 ---------------------------------------------------------------------------*/


static
int get_bitpix(char * buf)
{
	int		bitpix ;
	char *	where ;

	where = strstr(buf, "BITPIX") ;
	if (where == NULL) {
		fprintf(stderr, "cannot find BITPIX in header: aborting\n") ;
		return 0 ;
	}
	sscanf(where, "%*[^=] = %d", &bitpix) ;
	/*
	 * Check the returned value makes sense
	 */
	if ((bitpix != 8) &&
		(bitpix != 16) &&
		(bitpix != 32) &&
		(bitpix != -32) &&
		(bitpix != -64)) {
		bitpix = 0 ;
	}
	return bitpix ;
}



/*---------------------------------------------------------------------------
   Function	:	set_bitpix()
   In		:	allocated char buffer containing the complete FITS header,
				pixel type to write
   Out		:	int 0 if Ok, anything else otherwise
   Job		:	overwrites the value of BITPIX in a FITS header
   Notice	:
 ---------------------------------------------------------------------------*/


static
int set_bitpix(char * buf, int ptype)
{
	char *	where ;
	char	replace[81] ;

	where = strstr(buf, "BITPIX") ;
	if (where == NULL) {
		fprintf(stderr, "cannot find BITPIX in header: aborting\n") ;
		return -1 ;
	}
	sprintf(replace,
			"BITPIX  =                  % 3d / Bits per pixel                                 ",
			ptype) ;
	memcpy(where, replace, 80) ;
	return 0 ;
}


/*---------------------------------------------------------------------------
   Function	:	get_npix()
   In		:	allocated char buffer containing the complete FITS header
   Out		:	int # of pixels in the file
   Job		:	retrieves how many pixels in a FITS file from the header
   Notice	:	does not support extensions!
 ---------------------------------------------------------------------------*/


static
int get_npix(char * buf)
{
	int		naxes ;
	int		npix ;
	int		naxis ;
	char  * where ;
	char	lookfor[80] ;
	int		i ;

	where = strstr(buf, "NAXIS") ;
	if (where == NULL) {
		fprintf(stderr, "cannot find NAXIS in header: aborting\n") ;
		return 0 ;
	}
	sscanf(where, "%*[^=] = %d", &naxes) ;
	if ((naxes<1) || (naxes>999)) {
		fprintf(stderr, "illegal value for %s: %d\n", lookfor, naxes) ;
		return 0 ;
	}
	npix = 1 ;
	for (i=1 ; i<=naxes ; i++) {
		sprintf(lookfor, "NAXIS%d", i) ;
		where = strstr(buf, lookfor) ;
		if (where == NULL) {
			fprintf(stderr, "cannot find %s in header: aborting\n",
					lookfor) ;
			return 0 ;
		}
		sscanf(where, "%*[^=] = %d", &naxis) ;
		if (naxis<1) {
			fprintf(stderr, "error: found %s=%d\n", lookfor, naxis);
			return 0 ;
		}
		npix *= naxis ;
	}
	return npix ;
}


/*---------------------------------------------------------------------------
   Function	:	write_pixbuf()
   In		:	char buffers containing input and output pixel areas,
				number of pixels in input area, input pixel type,
				output pixel type.
   Out		:	int 0 if Ok, anything else otherwise
   Job		:	transfers a pixel buffer from one area to another, doing
				pixel format conversion on the fly if needed.
   Notice	:	optimized code, very limited readability!!!
 ---------------------------------------------------------------------------*/



static
int write_pixbuf(
	char *		bi,
	char *		bo,
	int			npix,
	int			ptype_in,
	int			ptype_out
)
{
	int					bpp ;
	register int		i ;
	register char *		ip ;
	register char *		op ;
	double				d ;
	float				f ;
	char				a2[2] ;
	char				a4[4] ;
	char				a8[8] ;

	/*
	 * Trivial case of identical copies
	 *    8 to  8
	 *   16 to 16
	 *   32 to 32
	 *  -32 to-32
	 *  -64 to-64
	 */

	if (ptype_in == ptype_out) {
		/* How many bytes per pixel? */
		bpp = ptype_in / 8 ; if (bpp < 0) bpp = - bpp ;
		/* Copy the buffer as it is and return */
		memcpy(bo, bi, npix * bpp) ; 
		return 0 ;
	}

	/*
	 * Use registers to boost speed
	 */
	ip = bi ;
	op = bo ;

	/*
	 * ---------- 8 bit to something 
	 */

	/*
	 * 8 to 16
	 * one byte in input becomes the LSB in the 2 output bytes
	 * we transfer it directly in the LSB, this is independent from
	 * machine endian-ness.
	 */

	if ((ptype_in == 8) && (ptype_out == 16)) {
		for (i=0 ; i<npix ; i++) {
			op++ ;
			*op++ |= *ip++ ;
		}
		return 0 ;
	}

	/*
	 * 8 to 32
	 * One byte in input becomes the LSB in the 4 output bytes
	 * we transfer it directly in the LSB, this is independent from
	 * machine endian-ness.
	 */

	if ((ptype_in == 8) && (ptype_out == 32)) {
		for (i=0 ; i<npix ; i++) {
			op++ ;
			op++ ;
			op++ ;
			*op++ |= *ip++ ;
		}
		return 0 ;
	}

	/*
	 * 8 to -32
	 * Conversion achieved through cast of the input byte into a float
	 * this means we use the local float type, need to know endian-ness
	 */

	if ((ptype_in == 8) && (ptype_out == -32)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				f = (float)(*(byte*)ip) ;
				memcpy(op, &f, 4) ;
				ip ++ ;
				op+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				f = (float)(*(byte*)ip) ;
				swap_bytes(&f, sizeof(float)) ;
				memcpy(op, &f, 4) ;
				ip ++ ;
				op+=4 ;
			}
		}
		return 0 ;
	}

	/*
	 * 8 to -64
	 * Conversion achieved through cast of the input byte into a double
	 * this means we use the local float type, need to know endian-ness
	 */
	if ((ptype_in == 8) && (ptype_out == -64)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = (double)(*(byte*)ip) ;
				memcpy(op, &d, 8) ;
				ip++ ;
				op+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				d = (double)(*(byte*)ip) ;
				swap_bytes(&d, sizeof(double)) ;
				memcpy(op, &d, 8) ;
				ip++ ;
				op+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * ---------- 16 bit to something
	 */

	/*
	 * 16 to 8
	 * We only keep the 8 LSBits. If the MSB is set, it means the input
	 * number is negative (FITS 16 bit is signed), so we clip to 0.  If
	 * the MSB is cleared but any other bit in the above 8 is set, it
	 * means the input is greater than 0xff, so we clip to 0xff.
	 */

	if ((ptype_in == 16) && (ptype_out == 8)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		for (i=0 ; i<npix ; i++) {
			if (!*ip) {
				ip++ ;
				*op++ = *ip++ ;
			} else {
				if (*ip & 0x80) {
					*op++ = 0 ;
					ip ++ ;
					ip ++ ;
				} else {
					*op++ = (byte)0xff ;
					ip ++ ;
					ip ++ ;
				}
			}
		}
		return 0 ;
	}

	/*
	 * 16 to 32
	 * We only need to transfer the 16 input bits as LSB for the output
	 * 32 bits. No need to know local endian-ness.
	 */
	if ((ptype_in == 16) && (ptype_out == 32)) {
		for (i=0 ; i<npix ; i++) {
			op++ ;
			op++ ;
			*op++ = *ip++ ;
			*op++ = *ip++ ;
		}
		return 0 ;
	}

	/*
	 * 16 to -32
	 * Conversion to IEEE float is done through a cast of the input
	 * number. Need to know local endian-ness.
	 */
	if ((ptype_in == 16) && (ptype_out == -32)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				/* signed short to float conversion */
				f = (float)(*(short16*)ip) ;
				memcpy(op, &f, 4) ;
				ip ++ ;
				ip ++ ;
				op+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				/* swap to get a local signed short */
				memcpy(a2, ip, 2) ;	
				swap_bytes(a2, 2) ;
				/* signed short to float conversion */
				f = (float)(*(short16*)&a2[0]) ;
				/* swap to write as Motorola data */
				swap_bytes(&f, sizeof(float)) ;
				memcpy(op, &f, 4) ;
				ip ++ ;
				ip ++ ;
				op+=4 ;
			}
		}
		return 0 ;
	}

	/*
	 * 16 to -64
	 * Conversion to IEEE double is done through a cast of the input
	 * number. Need to know local endian-ness.
	 */
	if ((ptype_in == 16) && (ptype_out == -64)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				/* signed short to double conversion */
				d = (double)(*(short16*)ip) ;
				memcpy(op, &d, 8) ;
				ip++ ;
				ip++ ;
				op+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				/* swap to get a local signed short */
				memcpy(a2, ip, 2) ;
				swap_bytes(a2, 2) ;
				/* signed short to double conversion */
				d = (double)(*(short16*)&a2[0]) ;
				/* swap back to Motorola format */
				swap_bytes(&d, sizeof(double)) ;
				memcpy(op, &d, 8) ;
				ip++ ;
				ip++ ;
				op+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * ---------- 32 bit to something
	 */

	/*
	 * 32 to 8
	 * Loss of the 24 upper bits: as for 16->8 we try to find if the
	 * input is negative by looking at the MSB and in that case clip to
	 * zero. If the MSB is cleared but any other of the upper 24 bits is
	 * set, the number is greater than 0xff and is clipped to 0xff.
	 */
	if ((ptype_in == 32) && (ptype_out == 8)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		for (i=0 ; i<npix ; i++) {
			/*
			 * if first bit is set, the number is negative on 32 bits,
			 * should be clipped to 0 on 8 bits
			 */
			if (*ip & 0x80) {
				*op++ = 0 ;
			} else if ((byte)*ip |
					   (byte)*(ip+1) |
					   (byte)*(ip+2)) {
			/*
			 * if any bit is set in the first 24 bits, the number is
			 * greater than 0xff, should be clipped to 0xff
			 */
				*op++ = (byte)0xff ;
			} else {
				*op++ = *(ip+3) ;
			}
			ip += 4 ;
		}
		return 0 ;
	}

	/*
	 * 32 to 16
	 * Loss of the 16 upper bits.
	 * We convert the input to a local signed long to see if it
	 * overflows a signed short. If it does, it is clipped.
	 */
	if ((ptype_in == 32) && (ptype_out == 16)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				if ((*(long32*)ip) > 32767) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					ip += 4 ;
				} else if ((*(long32*)ip) < -32768) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					ip += 4 ;
				} else {
					ip+=2 ;
					*op++ = *ip++ ;
					*op++ = *ip++ ;
				}
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				/* construct a local signed long to check its value */
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				if ((*(long32*)&a4[0]) > 32767) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					ip += 4 ;
				} else if ((*(long32*)&a4[0]) < -32768) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					ip += 4 ;
				} else {
					ip+=2 ;
					*op++ = *ip++ ;
					*op++ = *ip++ ;
				}
			}
		}
		return 0 ;
	}

	/*
	 * 32 to -32
	 * We cast the input to a local float, then output it to disk.
	 */
	if ((ptype_in == 32) && (ptype_out == -32)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				/* conversion through cast to local float */
				f = (float)(*(long32*)ip) ;
				memcpy(op, &f, 4) ;
				ip+=4 ;
				op+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				/* conversion through swap and cast to a local float */
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				f = (float)(*(long32*)&a4[0]) ;
				/* construct a motorola float in the output */
				swap_bytes(&f, sizeof(float)) ;
				memcpy(op, &f, 4) ;
				ip+=4 ;
				op+=4 ;
			}
		}
		return 0 ;
	}

	/*
	 * 32 to -64
	 * Conversion through cast to a local double, then output to disk.
	 */
	if ((ptype_in == 32) && (ptype_out == -64)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = (double)(*(long32*)ip) ;
				memcpy(op, &d, 8) ;
				ip+=4 ;
				op+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				/* build a local double though swap and cast */
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				d = (double)(*(long32*)&a4[0]) ;
				/* construct a Motorola double for output */
				swap_bytes(&d, sizeof(double)) ;
				memcpy(op, &d, 8) ;
				ip+=4 ;
				op+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * ---------- -32 to something
	 */

	/*
	 * -32 to 8
	 * Construct the input float through cast (need to know endian-ness).
	 * Clip it to [0..255], round the value to the smaller integer.
	 */
	if ((ptype_in == -32) && (ptype_out == 8)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				f = *((float*)ip) ;
				if (f>255.0) {
					*op++ = (byte)0xff ;
				} else if (f<0.0) {
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (byte)(f+0.5) ;
				}
				ip += 4;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				f = *((float*)&a4[0]) ;
				if (f>255.0) {
					*op++ = (byte)0xff ;
				} else if (f<0.0) {
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (byte)(f+0.5) ;
				}
				ip += 4;
			}
		}

		return 0 ;
	}

	/*
	 * -32 to 16
	 * Construct the input float through cast (need to know endian-ness)
	 * Clip it to [-32768..32767], round the value to the smaller
	 * integer.
	 */
	if ((ptype_in == -32) && (ptype_out == 16)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				f = *((float*)ip) ;
				if (f>32767.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
				} else if (f<-32768.0) {
					*op++ = (byte)0x80 ;
					*op++ = 0x00 ;
				} else {
					*op++ = (short16)f >> 8 ;
					*op++ = (short16)f & (byte)0xff ;
				}
				ip+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				f = *((float*)&a4[0]) ;
				if (f>32767.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
				} else if (f<-32768.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (short16)f >> 8 ;
					*op++ = (short16)f & (byte)0xff ;
				}
				ip+=4 ;
			}
		}
		return 0 ;
	}

	/*
	 * -32 to 32
	 * Construct the input float through cast (need to know endian-ness)
	 * Clip it to [-2147483648 .. 2147483647], then round to the smaller
	 * integer.
	 */
	if ((ptype_in == -32) && (ptype_out == 32)) {
		fprintf(stderr, "warning: probable loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				f = *((float*)ip) ;
				if (f>2147483647.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
				} else if (f<-2147483648.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (long32)f >> 24 ;
					*op++ = ((long32)f >> 16) & 0xff ;
					*op++ = ((long32)f >> 8)  & 0xff ;
					*op++ = (long32)f & 0xff ;
				}
				ip+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				f = *((float*)&a4[0]) ;
				if (f>2147483647.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
				} else if (f<-2147483648.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (long32)f >> 24 ;
					*op++ = ((long32)f >> 16) & 0xff ;
					*op++ = ((long32)f >> 8)  & 0xff ;
					*op++ = (long32)f & 0xff ;
				}
				ip+=4 ;
			}
		}
		return 0 ;
	}

	/*
	 * -32 to -64
	 * Conversion achieved through cast. Need to know endian-ness.
	 */
	if ((ptype_in == -32) && (ptype_out == -64)) {
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = (double)(*(float*)ip) ;
				memcpy(op, &d, 8) ;
				ip+=4 ;
				op+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a4, ip, 4) ;
				swap_bytes(a4, 4) ;
				d = (double)(*(float*)&a4[0]) ;
				swap_bytes(&d, sizeof(double)) ;
				memcpy(op, &d, 8) ;
				ip+=4 ;
				op+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * ---------- -64 to something
	 */


	/*
	 * -64 to 8
	 * Heavy loss of input precision!
	 * The input double is constructed through swap and cast. It is then
	 * clipped to [0..255] and rounded to the smaller integer.
	 */
	if ((ptype_in == -64) && (ptype_out == 8)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = *((double*)ip) ;
				if (d>255.0) {
					*op++ = (byte)0xff ;
				} else if (d<0.0) {
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (byte)(d+0.5) ;
				}
				ip+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a8, ip, 8) ;
				swap_bytes(a8, 8) ;
				d = *((double*)&a8[0]) ;
				if (d>255.0) {
					*op++ = (byte)0xff ;
				} else if (d<0.0) {
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (byte)(d+0.5) ;
				}
				ip+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * -64 to 16
	 * Heavy loss of input precision!
	 * The input double is constructed through swap and cast. It is then
	 * clipped to [-32768..32767] and rounded to the smaller integer.
	 */
	if ((ptype_in == -64) && (ptype_out == 16)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = *((double*)ip) ;
				if (d>32767.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
				} else if (d<-32768.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (short16)d >> 8 ;
					*op++ = (short16)d & (byte)0xff ;
				}
				ip+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a8, ip, 8) ;
				swap_bytes(a8, 8) ;
				d = *((double*)&a8[0]) ;
				if (d>32767.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
				} else if (d<-32768.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (short16)d >> 8 ;
					*op++ = (short16)d & (byte)0xff ;
				}
				ip+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * -64 to 32
	 * Heavy loss of input precision!
	 * The input double is constructed through swap and cast. It is then
	 * clipped to [-2147483648..2147483647] and rounded to the smaller
	 * integer.
	 */
	if ((ptype_in == -64) && (ptype_out == 32)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				d = *((double*)ip) ;
				if (d>2147483647.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
				} else if (d<-2147483648.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (long32)d >> 24 ;
					*op++ = ((long32)d >> 16) & (byte)0xff ;
					*op++ = ((long32)d >> 8)  & (byte)0xff ;
					*op++ = (long32)d & (byte)0xff ;
				}
				ip+=8 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a8, ip, 8) ;
				swap_bytes(a8, 8) ;
				d = *((double*)&a8[0]) ;
				if (d>2147483647.0) {
					*op++ = (byte)0x7f ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
					*op++ = (byte)0xff ;
				} else if (d<-2147483648.0) {
					*op++ = (byte)0x80 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
					*op++ = (byte)0x00 ;
				} else {
					*op++ = (long32)d >> 24 ;
					*op++ = ((long32)d >> 16) & 0xff ;
					*op++ = ((long32)d >> 8)  & 0xff ;
					*op++ = (long32)d & 0xff ;
				}
				ip+=8 ;
			}
		}
		return 0 ;
	}

	/*
	 * -64 to -32
	 * Heavy loss of input precision!
	 * The input double is constructed through swap and cast. It is then
	 * simply cast to a float using the local algorithm for this
	 * conversion.
	 */
	if ((ptype_in == -64) && (ptype_out == -32)) {
		fprintf(stderr, "warning: loss of precision in output\n") ;
		if (ENDIAN_NESS(compat_flags) == ENDIAN_MOTOROLA) {
			for (i=0 ; i<npix ; i++) {
				f = (float)(*(double*)ip) ;
				memcpy(op, &f, 4) ;
				ip+=8 ;
				op+=4 ;
			}
		} else {
			for (i=0 ; i<npix ; i++) {
				memcpy(a8, ip, 8) ;
				swap_bytes(a8, 8) ;
				f = (float)(*(double*)&a8[0]) ;
				swap_bytes(&f, sizeof(float)) ;
				memcpy(op, &f, 4) ;
				ip+=8 ;
				op+=4 ;
			}
		}
		return 0 ;
	}
	return 0 ;
}



static int bits_per_byte(void)
{
    unsigned char c=1 ;
    int i=0 ;
    while (c) { c <<=1 ; i++ ; }
    return i ;
}



/*---------------------------------------------------------------------------
   Function	:	test_machine_formats()
   In		:	void
   Out		:	flags (1 unsigned int)
				bit 0 stands for endian-ness: 0=Motorola, 1=Intel
				bit 1 stands for IEEE float compatible (1=yes, 0=no)
				bit 2 stands for IEEE double compatible (1=yes, 0=no)
   Job		:	test out endian-ness and IEEE compatibility of the local
				machine
   Notice	:	this means the machine which compiles the software
				should preferrably be the same as the only running it.
 ---------------------------------------------------------------------------*/


static int test_machine_formats(void)
{
    ushort16  ps = 0xff ;
    int       endian = ENDIAN_UNKNOWN ;
    unsigned  flags = 0;
    float     f ;
    double    d ;
    byte	* b ;

    /*
     * Test endian-ness for this machine
     */
    endian = ((*((char*)(&ps))) ? ENDIAN_INTEL : ENDIAN_MOTOROLA ) ;
    flags |= endian ;

    /*
     * Test whether the internal float format is IEEE 32bit floating
     */

    if (sizeof(float)==4) {
        f = (float)2.5e-16 ;
        b = (byte*)(&f) ;
        if (endian == ENDIAN_INTEL) {
            swap_bytes(&f, sizeof(float)) ;
        }
        if ((b[0] == (byte)0x25) &&
            (b[1] == (byte)0x90) &&
            (b[2] == (byte)0x1d) &&
            (b[3] == (byte)0x7d)) {
            flags |= 2 ;
        }
    }

    if (sizeof(double)==8) {
        d = (double)3.14e32 ;
        b = (byte*)&d ;
        if (endian == ENDIAN_INTEL) {
            swap_bytes(&d, sizeof(double)) ;
        }
        if ((b[0] == (byte)0x46) &&
            (b[1] == (byte)0xae) &&
            (b[2] == (byte)0xf6) &&
            (b[3] == (byte)0x79) &&
            (b[4] == (byte)0x70) &&
            (b[5] == (byte)0xae) &&
            (b[6] == (byte)0xec) &&
            (b[7] == (byte)0xd7)) {
            flags |= 4 ;
        }
    }
    return flags ;
}


/*---------------------------------------------------------------------------
   Function	:	swap_bytes()
   In		:	void *, int
   Out		:	void
   Job		:	swap bytes in a even-sized memory area
   Notice	:	's' shall always be even
				the area pointed to by p is modified (byte swapped) as:
				AB       -> BA
				ABCD     -> DCBA
				ABCDEFGH -> HGFEDCBA
 ---------------------------------------------------------------------------*/

static void swap_bytes(void * p, int s)
{
    byte tmp, *a, *b ;

    a = (byte*)p ;
    b = a + s ;

    while (a<b) {
        tmp = *a ;
        *a++ = *--b ;
        *b = tmp ;
    }
}

/*---------------------------------------------------------------------------
   Function	:	check_all_types()
   In		:	int (requested pixel type for output)
   Out		:	void
   Job		:	test out if the requested pixel type corresponds to a
				valid FITS type. Also tests out the sizes of all basic
				types (short16, long32, byte) to
				ensure they are the correct size. Also checks out that a
				float is IEEE compatible and a double is what it
				pretends to be.
   Notice	:	If one the of the tests goes wrong, this routine exits
				the main program.
 ---------------------------------------------------------------------------*/

static void
check_all_types(int ptype)
{
	/*
	 * Check the requested pixel type is correct
	 */
	if ((ptype!=8)&&(ptype!=16)&&(ptype!=32)&&(ptype!=-32)&&(ptype!=-64)) {
		fprintf(stderr, "illegal pixel type: %d\n", ptype) ;
		fprintf(stderr, "should be 8 16 32 -32 or -64\n") ;
		exit(-1) ;
	}

	/*
	 * Check the following:
	 * bitsize of byte is 8 bits
	 * bitsize of short is 16 bits
	 * bitsize of long is 32 bits
	 */

	if (bits_per_byte() != 8) {
		fprintf(stderr, "this machine does not support 8-bit chars\n") ;
		fprintf(stderr, "cannot go any further -- sorry.\n") ;
		exit(-1) ;
	}

	if (sizeof(short16)!=2) {
		fprintf(stderr, "internal: wrong definition for short16\n") ;
		fprintf(stderr, "edit definition and recompile.\n") ;
		exit(-1) ;
	}
	
	if (sizeof(long32)!=4) {
		fprintf(stderr, "internal: wrong definition for long32\n") ;
		fprintf(stderr, "edit definition and recompile.\n") ;
		exit(-1) ;
	}
	
	/*
	 * If -32 or -64 are requested in input check out that
	 * the local CPU works with float and double in IEEE format,
	 * otherwise, we cannot do much (arithmetic and bitwise
	 * interpretations needed, not implemented here...).
	 */

	compat_flags = test_machine_formats() ;
	if (ptype == -32) {
		if (!IEEE_FLOAT_COMPAT(compat_flags)) {
			fprintf(stderr,
				"this machine does not support IEEE floating point values\n");
			fprintf(stderr,
				"cannot convert from -32 type\n");
			exit(-1) ;
		}
	}
	if (ptype == -64) {
		if (!IEEE_DOUBLE_COMPAT(compat_flags)) {
			fprintf(stderr,
				"this machine does not support IEEE double values\n");
			fprintf(stderr,
				"cannot convert from -64 type\n");
			exit(-1) ;
		}
	}
	return ;
}

/*------------------------------ end of file -------------------------------*/
