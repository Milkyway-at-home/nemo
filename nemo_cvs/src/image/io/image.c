/* write_image(), read_image(), free_image(), create_image(), create_cube()   */
/* map2_image(), map3_image() */
/* TESTBED: main(), ini_matrix() */
/*
 * IMAGE-like format for 2 & 3D images using binary filestructure
 *	Simple NEMO example 
 *
 *  25-jun-87   V1.0 created 			P.J. Teuben
 *  30-jun-87	V2.0 using structures			PJT
 *  31-jul-87   V2.0a output modification  		PJT
 *  30-oct-87        some output to stderr 		PJT
 *  10-mar-88	V2.1 added data history 		PJT
 *   1-Jun-88	V2.2 new filestruct			PJT
 *  22-Dec-88   V2.3 added velocity into header		PJT
 *  10-jan-89   V2.4 added create_image			PJT
 *  22-jan-89   V2.5 added 3D support + create_cube	PJT
 *  30-jan-89   V3.1 deleted Vel from header - now in Z PJT
 *   1-feb-89   V4.0 FORDEF/CDEF of matrices used       PJT
 *  27-jun-89	V4.1 added free_image			PJT
 *  28-sep-90   V4.2 use get_data_coerced		PJT
 *   7-mar-92   V4.2a  happy gcc2.0			PJT
 *   5-mar-93   V4.2b fixed free() declaration for alpha    PJT
 *  30-jul-93   V4.3  Added Unit()/Time() member of image   PJT
 *  21-feb-94   V4.4  ANSI
 *  13-apr-96   V5.0 support for (orthogonal) non-linear axes	PJT
 *                   in header image.h
 *  21-feb-00   V6.1 mapX_image() routines, cleaned up code a bit PJT
 *   9-sep-02   V6.2 added copy_image()
 *  13-nov-02   V??? added boolean images for masking operations  PJT
 *   7-may-04   V7.0 (optional) proper astronomical axes          PJT
 *  16-mar-05   V7.1 protection to re-use a pointer for different images PJT
 *			
 *
 *  Note: bug in TESTBED section; new items (Unit) not filled in
 *
 *	  Example of usage: see SNAPCCD.C	for writing
 *			        CCDSMOOTH.C  	for reading
 * 
 *  Note on VO usage:  the WCS in SIAP is a simple, but common, way
 *   to express a WCS of an astronomical image:
 *   naxis, cframe, equinox, crpix, crval, cdelt, rotang, proj
 */

#include <stdinc.h>
#include <filestruct.h>
#include <history.h>
#include <image.h>

#define DLEV   5		/* local default debug output level */

local char *mystrcpy(char *);

/*	storage of matrices can be done in several ways: */

local char *matdef[] = {"ForDef", "CDef", NULL };

/* one of the above XXXDEF's must be defined below, and
   set the appropriate Tag value in the Image descriptor
   If it is not done, the compiled versions will (should)
   complain about missing integer 'idef'
 */
#if defined(FORDEF)
local int idef = 0;
#endif
#if defined(CDEF)
local int idef = 1;
#endif


/*
 *  WRITE_IMAGE: writes out a matrix, including header, in binary format
 *
 */


int write_image (stream outstr, imageptr iptr)
{
  put_history(outstr);
  put_set (outstr,ImageTag);
    put_set (outstr,ParametersTag);
      put_data (outstr,NxTag,  IntType,  &(Nx(iptr)),   0);
      put_data (outstr,NyTag,  IntType,  &(Ny(iptr)),   0);
      put_data (outstr,NzTag,  IntType,  &(Nz(iptr)),   0);
      put_data (outstr,XminTag,RealType, &(Xmin(iptr)), 0);
      put_data (outstr,YminTag,RealType, &(Ymin(iptr)), 0);
      put_data (outstr,ZminTag,RealType, &(Zmin(iptr)), 0);
      put_data (outstr,DxTag,  RealType, &(Dx(iptr)),   0);
      put_data (outstr,DyTag,  RealType, &(Dy(iptr)),   0);
      put_data (outstr,DzTag,  RealType, &(Dz(iptr)),   0);
      put_data (outstr,XrefTag,RealType, &(Xref(iptr)), 0);
      put_data (outstr,YrefTag,RealType, &(Yref(iptr)), 0);
      put_data (outstr,ZrefTag,RealType, &(Zref(iptr)), 0);
      put_data (outstr,MapMinTag, RealType, &(MapMin(iptr)), 0);
      put_data (outstr,MapMaxTag, RealType, &(MapMax(iptr)), 0);
      put_data (outstr,BeamTypeTag, IntType, &(BeamType(iptr)), 0);
      put_data (outstr,BeamxTag, RealType, &(Beamx(iptr)), 0);
      put_data (outstr,BeamyTag, RealType, &(Beamy(iptr)), 0);
      put_data (outstr,BeamzTag, RealType, &(Beamz(iptr)), 0);
      if (Namex(iptr))
         put_string (outstr,NamexTag,Namex(iptr));
      if (Namey(iptr))
         put_string (outstr,NameyTag,Namey(iptr));
      if (Namez(iptr))
         put_string (outstr,NamezTag,Namez(iptr));
      if (Unit(iptr))
         put_string (outstr,UnitTag,Unit(iptr));
      put_data (outstr,TimeTag,  RealType, &(Time(iptr)), 0);
      put_string(outstr,StorageTag,matdef[idef]);
      put_data (outstr,AxisTag,  IntType, &(Axis(iptr)), 0);
    put_tes (outstr, ParametersTag);
         
    put_set (outstr,MapTag);
		if (Nz(iptr)==1)
	                put_data (outstr,MapValuesTag,RealType,
        	        	Frame(iptr),Nx(iptr),Ny(iptr),0);
		else
	                put_data (outstr,MapValuesTag,RealType,
        	        	Frame(iptr),Nx(iptr),Ny(iptr),Nz(iptr),0);
    put_tes (outstr, MapTag);
  put_tes (outstr, ImageTag);
  return 1;
}
 	

/*
 * READ_IMAGE: read a matrix from a stream
 *	      returns 0 on error
 *	              1 if read seems OK
 *	To improve:	One cannot use a external buffer for iptr, iptr
 *			is always assumed to be initialized by one of
 *			the XXX_image() routines here. 
 *                    It will refuse to read an image that doesn't have
 *                    the same shape [nx,ny,nz] as the old one
 */
 
int read_image (stream instr, imageptr *iptr)
{
    string read_matdef;
    int nx=0, ny=0, nz=0;

    get_history(instr);         /* accumulate history */

    if (!get_tag_ok (instr,ImageTag))
        return 0;			/* not an image available */
        
    if (*iptr==NULL) {		/* allocate image if neccessary */
    	*iptr = (imageptr ) allocate(sizeof(image));
	dprintf (DLEV,"Allocated image @ %d ",*iptr);
	if (*iptr == NULL) return 0;	/* no memory available */
	Frame(*iptr) = NULL;
    } else {
        nx = Nx(*iptr);
        ny = Ny(*iptr);
        nz = Nz(*iptr);
    	dprintf (DLEV,"Image %dx%dx%d already allocated @ %d\n",
		 nx,ny,nz,*iptr);
    }
    	
    get_set (instr,ImageTag);
        get_set (instr,ParametersTag);
            get_data (instr,NxTag,IntType, &(Nx(*iptr)), 0);
            get_data (instr,NyTag,IntType, &(Ny(*iptr)), 0);
            get_data (instr,NzTag,IntType, &(Nz(*iptr)), 0);
	    if ((nx>0 || ny>0 || nz>0) &&
		(nx != Nx(*iptr) || ny != Ny(*iptr) || nz != Nz(*iptr)))
	      error("Cannot read different sized images in old pointer yet");
	    if (get_tag_ok(instr,AxisTag))
	      get_data (instr,AxisTag,IntType, &(Axis(*iptr)), 0);
	    else
	      Axis(*iptr) = 0;
	    if (Axis(*iptr) == 1) {
	      get_data_coerced (instr,XrefTag,RealType, &(Xref(*iptr)), 0);
	      get_data_coerced (instr,YrefTag,RealType, &(Yref(*iptr)), 0);
	      get_data_coerced (instr,ZrefTag,RealType, &(Zref(*iptr)), 0);
	    } else {
	      Xref(*iptr) = 0.0;
	      Yref(*iptr) = 0.0;
	      Zref(*iptr) = 0.0;
	    }

            get_data_coerced (instr,XminTag,RealType, &(Xmin(*iptr)), 0);
            get_data_coerced (instr,YminTag,RealType, &(Ymin(*iptr)), 0);
            get_data_coerced (instr,ZminTag,RealType, &(Zmin(*iptr)), 0);
            get_data_coerced (instr,DxTag,RealType, &(Dx(*iptr)), 0);
            get_data_coerced (instr,DyTag,RealType, &(Dy(*iptr)), 0);
            get_data_coerced (instr,DzTag,RealType, &(Dz(*iptr)), 0);
	    get_data_coerced (instr,MapMinTag, RealType, &(MapMin(*iptr)), 0);
	    get_data_coerced (instr,MapMaxTag, RealType, &(MapMax(*iptr)), 0);
	    get_data (instr,BeamTypeTag, IntType, &(BeamType(*iptr)), 0);
	    get_data_coerced (instr,BeamxTag, RealType, &(Beamx(*iptr)), 0);
	    get_data_coerced (instr,BeamyTag, RealType, &(Beamy(*iptr)), 0);
	    get_data_coerced (instr,BeamzTag, RealType, &(Beamz(*iptr)), 0);
            if (get_tag_ok(instr,NamexTag))             /* X-axis name */
                Namex(*iptr) = get_string(instr,NamexTag);
            else
                Namex(*iptr) = NULL;
            if (get_tag_ok(instr,NameyTag))             /* Y-axis name */
                Namey(*iptr) = get_string(instr,NameyTag);
            else
                Namey(*iptr) = NULL;
            if (get_tag_ok(instr,NamezTag))             /* Z-axis name */
                Namez(*iptr) = get_string(instr,NamezTag);
            else
                Namez(*iptr) = NULL;
            if (get_tag_ok(instr,UnitTag))             /* units  */
                Unit(*iptr) = get_string(instr,UnitTag);
            else
                Unit(*iptr) = NULL;
            if (get_tag_ok(instr,TimeTag))             /* time  */
   	    	get_data_coerced (instr,TimeTag, RealType, &(Time(*iptr)), 0);
   	    else
   	    	Time(*iptr) = 0.0;
            read_matdef = get_string(instr,StorageTag);
	    if (!streq(read_matdef,matdef[idef]))
                dprintf(0,"read_image: StorageTag = %s, compiled with %s\n",
		        read_matdef, matdef[idef]);
         get_tes (instr,ParametersTag);

         get_set (instr,MapTag);
            if (Frame(*iptr)==NULL) {        /* check if allocated */
                Frame(*iptr) = (real *) 
			allocate(Nx(*iptr)*Ny(*iptr)*Nz(*iptr)*sizeof(real));   
                dprintf (DLEV,"Frame allocated @ %d ",Frame(*iptr));
                if (Frame(*iptr)==NULL) {
                    printf ("READ_IMAGE: Not enough memory to allocate image\n");
                    return 0;
                }
            } else
                dprintf (DLEV,"Frame already allocated @ %d\n",Frame(*iptr));
	    if (Nz(*iptr)==1)
                get_data_coerced (instr,MapValuesTag,RealType, Frame(*iptr), 
                                Nx(*iptr), Ny(*iptr), 0);
            else
                get_data_coerced (instr,MapValuesTag,RealType, Frame(*iptr),
                                Nx(*iptr), Ny(*iptr), Nz(*iptr), 0);
         get_tes (instr,MapTag);
      get_tes (instr,ImageTag);

      dprintf (DLEV,"Frame size %d * %d \n",Nx(*iptr), Ny(*iptr));
      
      return 1;		/* succes return code  */
}

/*
 * FREE_IMAGE: free an image, previously allocated by one of the image(3NEMO)
 *	       routines
 *             Only free's up the big data, doesn't free string space of
 *             axis names, units etc. since they are frequently in private
 *             space - sloppy programming
 */
 
int free_image (imageptr iptr)
{
  free ((char *) Frame(iptr));
  free ((char *) iptr);
  return  0;
}

int free_image_mask (image_maskptr mptr)
{
  free ((char *) Frame(mptr));
  free ((char *) mptr);
  return  0;
}

/*
 * CREATE_IMAGE: create a blank image Nx by Ny in size ; forced 2D
 *		see: create_cube for a full 3D version
 */
 
int create_image (imageptr *iptr, int nx, int ny)
{
    *iptr = (imageptr ) allocate(sizeof(image));
    dprintf (DLEV,"create_image:Allocated image @ %d size=%d * %d",*iptr,nx,ny);
    if (*iptr == NULL) return 0;	/* no memory available */
    	
    Frame(*iptr) = (real *) allocate(nx*ny*sizeof(real));	
    dprintf (DLEV,"Frame allocated @ %d ",Frame(*iptr));
    if (Frame(*iptr)==NULL) {
        printf ("CREATE_IMAGE: Not enough memory to allocate image\n");
        return 0;
    }
    Axis(*iptr) = 0;            /* old style axis with no reference pixel */
    Nx(*iptr) = nx;             /* old style ONE map, no cube */
    Ny(*iptr) = ny;
    Nz(*iptr) = 1;
    Xmin(*iptr) = 0.0;          /* start lower left corner at 0.0 */
    Ymin(*iptr) = 0.0;
    Zmin(*iptr) = 0.0;
    Dx(*iptr) = 1.0;            /* unity pixels */
    Dy(*iptr) = 1.0;
    Dz(*iptr) = 1.0;
    Xref(*iptr) = 0.0;
    Yref(*iptr) = 0.0;
    Zref(*iptr) = 0.0;
    MapMin(*iptr) = 0.0;
    MapMax(*iptr) = 0.0;
    BeamType(*iptr) = 0;
    Beamx(*iptr) = 0.0;         /* name beams */
    Beamy(*iptr) = 0.0;
    Beamz(*iptr) = 0.0;
    Namex(*iptr) = NULL;        /* no axis names */
    Namey(*iptr) = NULL;
    Namez(*iptr) = NULL;
    Unit(*iptr)  = NULL;        /* no units */
    Time(*iptr)  = 0.0;
    Storage(*iptr) = matdef[idef];
    Axis(*iptr) = 0;
    Mask(*iptr) = NULL;

    return 1;		/* succes return code  */
}

int create_image_mask(imageptr iptr, image_maskptr *mptr)
{
  int nx = Nx(iptr);
  int ny = Ny(iptr);
  int nz = Nz(iptr);
  int np = nx*ny*nz;

  *mptr = (image_maskptr ) allocate(sizeof(image_mask));
  dprintf (DLEV,"create_image_mask:Allocated image_mask @ %d size=%d * %d * %d",*mptr,nx,ny,nz);
  if (*mptr == NULL) return 0;	/* no memory available */
    	
  Frame(*mptr) = (bool *) allocate(np*sizeof(bool));	
  dprintf (DLEV,"Frame allocated @ %d ",Frame(*mptr));
  if (Frame(*mptr)==NULL) {
    printf ("CREATE_IMAGE_MASK: Not enough memory to allocate image\n");
    return 0;
  }
  Nx(*mptr) = nx;
  Ny(*mptr) = ny;
  Nz(*mptr) = nz;
  return 1;
}

int copy_image (imageptr iptr, imageptr *optr)
{
  int nx,ny,nz;
  nx = Nx(iptr);
  ny = Ny(iptr);
  nz = Nz(iptr);

  *optr = (imageptr ) allocate(sizeof(image));
  dprintf (DLEV,"copy_image:Allocated image @ %d size=%d * %d * %d",*optr,nx,ny,nz);
    	
  Frame(*optr) = (real *) allocate(nx*ny*nz*sizeof(real));	
  dprintf (DLEV,"Frame allocated @ %d ",Frame(*optr));
  if (Frame(*optr)==NULL) {
    warning ("CREATE_IMAGE: Not enough memory to allocate image\n");
    return 0;
  }
  Nx(*optr) = nx;
  Ny(*optr) = ny;
  Nz(*optr) = nz;
  Xmin(*optr) = Xmin(iptr);
  Ymin(*optr) = Ymin(iptr);
  Zmin(*optr) = Zmin(iptr);
  Dx(*optr) = Dx(iptr);
  Dy(*optr) = Dy(iptr);
  Dz(*optr) = Dz(iptr);
  Namex(*optr) = mystrcpy(Namex(iptr));
  Namey(*optr) = mystrcpy(Namey(iptr));
  Namez(*optr) = mystrcpy(Namez(iptr));
  Xref(*optr) = Xref(iptr);
  Yref(*optr) = Yref(iptr);
  Zref(*optr) = Zref(iptr);
  Storage(*optr) = matdef[idef];
  Axis(*optr) = Axis(iptr);
  
  return 1;		/* succes return code  */
}


/*
 * CREATE_CUBE: create a blank cube Nx by Ny by Nz in size
 */
int create_cube (imageptr *iptr, int nx, int ny, int nz)
{
    *iptr = (imageptr ) allocate(sizeof(image));
    dprintf (DLEV,"CREATE_CUBE: Allocated image @ cube %d size=%d * %d * %d",
		*iptr,nx,ny,nz);
    if (*iptr == NULL)
	return 0;	/* no memory available */
    	
    Frame(*iptr) = (real *) allocate(nx*ny*nz*sizeof(real));	
    dprintf (DLEV,"Frame allocated @ %d ",Frame(*iptr));
    if (Frame(*iptr)==NULL) {
        printf ("CREATE_CUBE: Not enough memory to allocate image cube\n");
        return 0;
    }
    Nx(*iptr) = nx;             /* cube dimension */
    Ny(*iptr) = ny;
    Nz(*iptr) = nz;
    Xmin(*iptr) = 0.0;          /* start lower left corner at 0.0 */
    Ymin(*iptr) = 0.0;
    Zmin(*iptr) = 0.0;
    Xref(*iptr) = 0.0;
    Yref(*iptr) = 0.0;
    Zref(*iptr) = 0.0;
    Dx(*iptr) = 1.0;            /* unity pixels */
    Dy(*iptr) = 1.0;
    Dz(*iptr) = 1.0;
    MapMin(*iptr) = 0.0;
    MapMax(*iptr) = 0.0;
    BeamType(*iptr) = 0;
    Beamx(*iptr) = 0.0;         /* name beams */
    Beamy(*iptr) = 0.0;
    Beamz(*iptr) = 0.0;
    Namex(*iptr) = NULL;        /* no axis names yet */
    Namey(*iptr) = NULL;
    Namez(*iptr) = NULL;
    Unit(*iptr)  = NULL;        /* no units */
    Time(*iptr)  = 0.0;
    Storage(*iptr) = matdef[idef];
    Axis(*iptr) = 0;
    Mask(*iptr) = NULL;
    
    return 1;		/* succes return code  */
}

/*
 * MAPx_IMAGE: return more user friendly pointer thingos
 *             a little modeled after the NumRec routines
 */
real **map2_image (imageptr iptr)
{
    real *base = Frame(iptr);
    real **map;
    int nx, ny, ix;

    nx = Nx(iptr);
    ny = Ny(iptr);

#if defined(CDEF)
    map = (real **) allocate(sizeof(real *) * ny);
    for (ix=0; ix<nx; ix++) {               /* CDEF */
        map[ix] = base;
        base += ny;
    }
    return map;
#else
    error("map2_image: not implemented for !CDEF");    
#endif    
}

real ***map3_image (imageptr iptr)
{
    real *base = Frame(iptr);
    real ***cube;
    int nx, ny, nz, ix, iy;

    nx = Nx(iptr);
    ny = Ny(iptr);
    nz = Nz(iptr);

#if defined(CDEF)
    cube = (real ***) allocate(sizeof(real **) * nz);
    for (iy=0; iy<ny; iy++) {
    for (ix=0; ix<nx; ix++) {               /* CDEF */
        cube[ix][iy] = base;
        base += nz;
    }
    }
    return cube;
#else
    error("map3_image: not implemented for !CDEF");    
#endif    
}

char *mystrcpy(char *a)
{
  char *b;
  if (a==NULL || *a == 0) return NULL;
  b = allocate(strlen(a)+1);
  strcpy(b,a);
  return b;
}


#ifdef TESTBED

#include <getparam.h>

#define N 10
/*  #define AR 1  */

string defv[] = {
	"mode=w\n	Read (r) or Write (w)",
	"VERSION=4.4\n	21-feb-94 pjt",
	NULL
};

string usage = "testbed for image I/O";

void ini_matrix(imageptr *, int, int);
	
nemo_main()
{
    real frame[N][N];
    image f1;
    imageptr fp1, fp2;	/* or image *fp1 */
    stream instr,outstr;

    if (strcmp(getparam("mode"),"w")==0) {		/* write test */
	printf ("WRITING test (mode=w) foo.dat\n");
	outstr = stropen ("foo.dat","w");

	f1.frame = &frame[0][0];  /* compiler complained when = frame used ???? */
				  /* would be same as fp2->frame = ... */
	fp1 = &f1;		  /* to initialize structures, we need pointer to them */
	ini_matrix (&fp1,N,N);    /* ini_matrix MUST have pointer to pointer to image */
	
	fp2 = NULL;		  /* This is to check dynamic allocation of whole thing */
	ini_matrix (&fp2,N,N);

	write_image (outstr,fp2);
	write_image (outstr,&f1);		/* or fp1 */
	strclose(outstr);
	exit(0);
    } else {
	printf ("READING test (mode<>w) foo.dat\n");
	fp2=NULL;					/* read test */
	instr = stropen ("foo.dat","r");
	while (read_image(instr,&fp2)) {
            printf ("Read image\n");
            printf ("with MapValue(5,5)=%f\n",MapValue(fp2,5,5));
	}
	strclose(instr);
    }
}

void ini_matrix(iptr, nx, ny)
imageptr *iptr;
int nx,ny;
{
	int ix,iy;
	
	printf ("iptr @ %d    (sizeof = %d)\n",*iptr,sizeof(image));
	
	if (*iptr==0) {
		*iptr = (imageptr ) allocate((sizeof(image)));
		printf ("allocated image 'iptr' @ %d\n",*iptr);
		Frame(*iptr) = (real *) allocate (nx*ny*sizeof(real));
		printf ("Allocated frame @ %d\n",Frame(*iptr));
	} else {
		printf ("Image already allocated @ %d\n",*iptr);
		printf ("with Frame @ %d\n",Frame(*iptr));
	}
	Nx(*iptr) = nx;
	Ny(*iptr) = ny;
	Xmin(*iptr) = 1.0;
	Ymin(*iptr) = 2.0;
	Dx(*iptr) = 0.1;
	Dy(*iptr) = 0.1;
	for (ix=0; ix<nx; ix++)
	for (iy=0; iy<ny; iy++)
		MapValue(*iptr,ix,iy)  = (ix*ny+iy)*0.1;
}
#endif

