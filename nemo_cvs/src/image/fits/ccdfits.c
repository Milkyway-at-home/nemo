/* 
 *      CCDFITS: convert CCD image to a fits file
 *
 *	Also contains a controversial 'long strings' implementation
 *
 *      29-apr-88  V1.0 created (again) PJT
 *       2-jun-88  V1.1 renamed 'wfits' to 'ccdfits'  
 *                      new filestruct: code is almost same
 *      23-dec-88  V1.2 added velocity information to header
 *      30-jan-89  V1.3 velocity information is now in Z
 *       6-mar-89  V1.3a 2D cube written when NZ=1
 *      28-jul-90  V2.0 Bob's new fitsio.c routines, finally
 *       1-oct-90  V2.1 default scale is now 1, not 1/60   PJT
 *                      crpix is now real !!
 *	11-oct-90  V2.2 blocking factor introduced 	   PJT
 *      14-feb-91  V2.3 write axes names for AIPS          PJT
 *	 7-mar-92  V2.4 gcc2.0 happy - new style	   PJT
 *	23-feb-93  V2.5  axisnames (CTYPE's) now ok	   PJT
 *      18-may-93  V2.6 history trial following sci.astro.fits proposal  PJT
 *      21-feb-94  V2.6a ansi
 *      22-nov-94      b ???
 *	20-mar-95      c long histories now written multiline		PJT
 *	 8-oct-95      d force CTYPEs for some WCS freaks		pjt
 *      28-mar-98      e write out Headline[] as COMMENT                PJT
 *      11-mar-99      e fixed bug in the above change                  pjt
 *	21-mar-99  V2.7  fixed scaling problems in bitpix=8		pjt
 *      18-may-99  V2.8  cdmatrix optional output                       pjt
 *      15-oct-99  V2.9  force a more silly RA---SIN/DEC--SIN axis      pjt
 *      23-mar-01  V3.0  allow fits reference image to inherit WCS from PJT
 *       8-apr-01      a fixed SINGLEPREC operation
 *      23-may-01  V3.1  added dummy= to be able to not write dummy axes  PJT
 *                        NOT WORKING YET
 *	 7-aug-01  
 *      18-dec-01  V4.0  work with new fitsio_nemo.h
 *       6-may-02      b implemented dummy=f
 *       6-jul-02   5.0  changed ref*= to crval/cdelt/crpix=
 *      10-jul-02   5.1  better handling of long COMMENT fields
 *       4-feb-04   5.2  also listen to changed crval/cdelt/crpix= without refmap
 *       8-may-05   5.3  deal with the new  axis type 1 images          PJT
 *
 *  TODO:
 *      reference mapping has not been well tested, especially for 2D
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>
#include <fitsio_nemo.h>

string defv[] = {
        "in=???\n        Input image filename",
        "out=???\n       Output fits filename",
	"bitpix=-32\n	 FITS bitpix value {16,32,-32}",
        "scale=1,1,1\n   Extra scalefactor for cdelt&crval",
        "iscale=1,0\n    Scale and Offset Intensity conversion factor",
        "object=\n       Object name",
        "comment=\n      Extra comments",
	"cdmatrix=f\n    Use standard CD matrix instead of CDELT?",
	"blocking=1\n	 Blocking factor for output (blocksize/2880)",
	"refmap=\n       reference map to inherit WCS from",
	"crpix=\n        reference pixel, if different from default",
	"crval=\n        reference value, if different from default",
	"cdelt=\n        pixel value increment, if different from default",
	"radecvel=f\n    Enforce reasonable RA/DEC/VEL axis descriptor",
	"dummy=t\n       Write dummy axes also ?",
	"nfill=0\n	 Add some dummy comment cards to test fitsio",
	"ndim=\n         Testing if only that many dimensions need to be written",
	"select=1\n      Which image (if more than 1 present) to select",
        "VERSION=5.4\n   20-jun-09 PJT",
        NULL,
};

string usage = "convert image to a fits file";

string cvsid = "$Id: ccdfits.c,v 1.14 2009/06/21 02:59:37 pteuben Exp $";

stream  instr, outstr;                         /* file streams */

imageptr iptr=NULL;                     /* image, allocated dynamically */
int  isel = 0;

double scale[3];        /* scale conversion for FITS (CDELT) */
double iscale[2];	/* intensity rescale */
string object;           /* name of object in FITS header */
string comment;          /* extra comments */
string headline;         /* optional NEMO headline, added as COMMENT */
bool Qcdmatrix;         /* writing out new-style cdmatrix ? */
bool Qradecvel;         /* fake astronomy WCS header */
bool Qrefmap;
bool Qcrval, Qcdelt, Qcrpix;
bool Qdummy;            /* write dummy axes ? */

int   nref = 0, nfill = 0;
FLOAT ref_crval[3], ref_crpix[3], ref_cdelt[3];
char  ref_ctype[3][80], ref_cunit[3][80];

void setparams(void);
void write_fits(string,imageptr);
void stuffit(FITS *, string, string);
void set_refmap(string);
void permute(int *x ,int *idx, int n);

void nemo_main()
{
  int i;

  setparams();                               /* set cmdln par's */
  instr = stropen (getparam("in"), "r");     /* open image file */
  for (i=0; i<isel; i++)
    if (read_image (instr,&iptr)==0)
      error("Cannot process image select=%d",i);
  headline = ask_headline();                 /* possible headline */
  strclose(instr);                           /* close image file */
  write_fits(getparam("out"),iptr);          /* write fits file */
  free_image(iptr);
}

void setparams(void)
{
  int i,n;
  real tmpr[3];

  isel = getiparam("select");

  switch (nemoinpd(getparam("scale"),scale,3)) {
  case 1:
    scale[1] = scale[0];
    scale[2] = 1.0;
    break;
  case 2:
    scale[2] = 1.0;
    break;
  case 3:
    break;
  case 0:
    scale[0] = scale[1] = scale[2] = 1.0;
    break;
  default:
    error("parsing error scale=%s",getparam("scale"));
  }
  if (nemoinpd(getparam("iscale"),iscale,2) != 2)
    error("parsing error scale=%s - must be 2 numbers",getparam("scale"));
  object = getparam("object");
  comment = getparam("comment");
  Qcdmatrix = getbparam("cdmatrix");
  Qradecvel = getbparam("radecvel");
  Qrefmap = hasvalue("refmap");
  if (Qrefmap)
    set_refmap(getparam("refmap"));

  Qcrpix = hasvalue("crpix");
  if (Qcrpix) {
    nref =  nemoinpr(getparam("crpix"),tmpr,3);
    for (i=0; i<nref; i++)
      ref_crpix[i] = tmpr[i];
  }
  Qcrval = hasvalue("crval");
  if (Qcrval) {
    nref =  nemoinpr(getparam("crval"),tmpr,3);
    for (i=0; i<nref; i++)
      ref_crval[i] = tmpr[i];
  }
  Qcdelt = hasvalue("cdelt");
  if (Qcdelt) {
    nref =  nemoinpr(getparam("cdelt"),tmpr,3);
    for (i=0; i<nref; i++)
      ref_cdelt[i] = tmpr[i];
  }
  Qdummy = getbparam("dummy");
  nfill = getiparam("nfill");
}

static string ctypes[3] = { "CTYPE1",   "CTYPE2",   "CTYPE3" };
static string cdelts[3] = { "CDELT1",   "CDELT2",   "CDELT3" };
static string crvals[3] = { "CRVAL1",   "CRVAL2",   "CRVAL3" };
static string crpixs[3] = { "CRPIX1",   "CRPIX2",   "CRPIX3" };
static string radeve[3] = { "RA---SIN", "DEC--SIN", "VELO-LSR" };
static string xyz[3]    = { "X",        "Y",        "Z" };

void write_fits(string name,imageptr iptr)
{
    FLOAT tmpr,xmin[3],xref[3],dx[3],mapmin,mapmax;   /* fitsio FLOAT !!! */
    FITS *fitsfile;
    char *cp, origin[80];
    string *hitem, axname[3];
    float *buffer, *bp;
    int i, j, k, axistype, bitpix, keepaxis[3], nx[3], p[3], nx_out[3], ndim=3;
    double bscale, bzero;
    
    if (hasvalue("ndim")) ndim = getiparam("ndim");
    nx[0] = Nx(iptr);
    nx[1] = Ny(iptr);
    nx[2] = Nz(iptr);   if (nx[2] <= 0) nx[2] = 1;
    xmin[0] = Xmin(iptr)*scale[0];
    xmin[1] = Ymin(iptr)*scale[1];
    xmin[2] = Zmin(iptr)*scale[2];
    dx[0] = Dx(iptr)*scale[0];
    dx[1] = Dy(iptr)*scale[1];
    dx[2] = Dz(iptr)*scale[2];
    xref[0] = Xref(iptr)+1.0;
    xref[1] = Yref(iptr)+1.0;
    xref[2] = Zref(iptr)+1.0;
    axistype = Axis(iptr);
    axname[0] = (Namex(iptr) ? Namex(iptr) : xyz[0]);
    axname[1] = (Namey(iptr) ? Namey(iptr) : xyz[1]);
    axname[2] = (Namez(iptr) ? Namez(iptr) : xyz[2]);
    mapmin = MapMin(iptr);
    mapmax = MapMax(iptr);
    if (Qdummy) 
      for (i=0; i<3; i++) p[i] = i;
    else {
      if (Qrefmap) warning("dummy=f and usage of refmap will result in bad headers");
      permute(nx,p,3);
      dprintf(0,"Reordering axes: %d %d %d\n",p[0],p[1],p[2]);
    }
#if 1
    for (i=0; i<3; i++)  nx_out[i] = nx[p[i]];
    /* fix this so CubeValue works */
    Nx(iptr) = nx_out[0];
    Ny(iptr) = nx_out[1];
    Nz(iptr) = nx_out[2];
#else
    for (i=0; i<3; i++)  nx_out[i] = nx[i];
#endif
    sprintf(origin,"NEMO ccdfits %s",getparam("VERSION"));

    dprintf(1,"NEMO Image file written to FITS disk file\n");
    dprintf(1,"%d %d %d   %f %f %f   %f %f %f  %f %f %f   %f %f \n",
	    nx[0],nx[1],nx[2],xmin[0],xmin[1],xmin[2],dx[0],dx[1],dx[2],xref[0],xref[1],xref[2],
	    mapmin,mapmax);
    dprintf(1,"keepaxis(%d,%d,%d)\n",keepaxis[0],keepaxis[1],keepaxis[2]);
    
    fit_setblocksize(2880*getiparam("blocking"));
    bitpix = getiparam("bitpix");
    fit_setbitpix(bitpix);
    if (bitpix == 16) {      /* scale from -2^(bitpix-1) .. 2^(bitpix-1)-1 */
        bscale = (mapmax - mapmin) / (2.0*32768.0 - 1.0);
        bzero = mapmax - bscale*32767.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 32) {
        bscale = (mapmax - mapmin) / (2.0*2147483648.0 - 1.0);
        bzero = mapmax - bscale*2147483647.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 8) {
    	bscale = (mapmax - mapmin) / (2.0*128.0 - 1.0);
    	bzero = mapmin;
    	fit_setscale(bscale,bzero);
    }
    dprintf(1,"bscale,bzero=%g %g\n",bscale,bzero);

    fitsfile = fitopen(name,"new",ndim,nx_out);
    if (fitsfile==NULL) error("Could not open fitsfile %s for writing\n",name);

    if (Qrefmap || Qcrpix) {
      fitwrhdr(fitsfile,"CRPIX1",ref_crpix[0]);       
      fitwrhdr(fitsfile,"CRPIX2",ref_crpix[1]);       
      if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",ref_crpix[2]);
    } else {
      if (axistype==1) {
	fitwrhdr(fitsfile,"CRPIX1",xref[0]);      
	fitwrhdr(fitsfile,"CRPIX2",xref[1]);
	if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",xref[2]);
      } else {
	fitwrhdr(fitsfile,"CRPIX1",1.0);        /* CRPIX = 1 by Nemo definition */
	fitwrhdr(fitsfile,"CRPIX2",1.0);
	if (ndim>2) fitwrhdr(fitsfile,"CRPIX3",1.0);
      }
    }
    if (Qrefmap || Qcrval) {
      fitwrhdr(fitsfile,"CRVAL1",ref_crval[0]);
      fitwrhdr(fitsfile,"CRVAL2",ref_crval[1]);
      if (ndim>2) fitwrhdr(fitsfile,"CRVAL3",ref_crval[2]);
    } else {
      fitwrhdr(fitsfile,"CRVAL1",xmin[p[0]]);
      fitwrhdr(fitsfile,"CRVAL2",xmin[p[1]]);
      if (ndim>2) fitwrhdr(fitsfile,"CRVAL3",xmin[p[2]]);
    }

    if (Qcdmatrix) {
      fitwrhdr(fitsfile,"CD1_1",dx[p[0]]);    
      fitwrhdr(fitsfile,"CD2_2",dx[p[1]]);    
      if (ndim>2) fitwrhdr(fitsfile,"CD3_3",dx[p[2]]);    
    } else {
      if (Qrefmap || Qcdelt) {
	fitwrhdr(fitsfile,"CDELT1",ref_cdelt[0]*scale[0]);
	fitwrhdr(fitsfile,"CDELT2",ref_cdelt[1]*scale[1]);
	if (ndim>2) fitwrhdr(fitsfile,"CDELT3",ref_cdelt[2]*scale[2]);
      } else {
	fitwrhdr(fitsfile,"CDELT1",dx[p[0]]);    
	fitwrhdr(fitsfile,"CDELT2",dx[p[1]]);    
	if (ndim>2) fitwrhdr(fitsfile,"CDELT3",dx[p[2]]);
      }
    }

    if (Qradecvel) {
      dprintf(0,"[Axes names written as %s, %s, %s\n",
	      radeve[p[0]],radeve[p[1]],radeve[p[2]]);
      fitwrhda(fitsfile,"CTYPE1",radeve[p[0]]);
      fitwrhda(fitsfile,"CTYPE2",radeve[p[1]]);
      if (ndim>2) fitwrhda(fitsfile,"CTYPE3",radeve[p[2]]);
    } else {
      if (Qrefmap) {
        fitwrhda(fitsfile,"CTYPE1",ref_ctype[0]);
        fitwrhda(fitsfile,"CTYPE2",ref_ctype[1]);
        if (ndim>2) fitwrhda(fitsfile,"CTYPE3",ref_ctype[2]);
      } else {
	fitwrhda(fitsfile,"CTYPE1",axname[p[0]]);
	fitwrhda(fitsfile,"CTYPE2",axname[p[1]]);
	if (ndim>2) fitwrhda(fitsfile,"CTYPE3",axname[p[2]]);
      }
    }

    fitwrhdr(fitsfile,"DATAMIN",mapmin);
    fitwrhdr(fitsfile,"DATAMAX",mapmax);
    fitwrhda(fitsfile,"ORIGIN",origin);

    cp = getenv("USER");                                /* AUTHOR */
    if (cp)
        fitwrhda(fitsfile,"AUTHOR",cp);
    else
        fitwrhda(fitsfile,"AUTHOR","NEMO");

    if (object)                                        /* OBJECT */
        fitwrhda(fitsfile,"OBJECT",object);

    if (comment)                                       /* COMMENT */
        stuffit(fitsfile,"COMMENT",comment);
    if (headline)
        stuffit(fitsfile,"COMMENT",headline);

    hitem = ask_history();                              /* HISTORY */
    fitwra(fitsfile,"HISTORY","NEMO: History in reversed order");
    for (i=0, cp=hitem[0]; cp != NULL; i++) {
    	stuffit(fitsfile,"HISTORY",cp);
        cp = hitem[i+1];
    }

    for(i=0; i<nfill; i++)   /* debugging header I/O */
        fitwra(fitsfile,"COMMENT","Dummy filler space");

    buffer = (float *) allocate(nx[p[0]]*sizeof(float));

    for (k=0; k<nx_out[2]; k++) {          /* loop over all planes */
        fitsetpl(fitsfile,1,&k);
        for (j=0; j<nx_out[1]; j++) {      /* loop over all rows */
            for (i=0, bp=buffer; i<nx_out[0]; i++, bp++)
                *bp =  iscale[0] * CubeValue(iptr,i,j,k) + iscale[1];
            fitwrite(fitsfile,j,buffer);
        }
    }
    free(buffer);
    fitclose(fitsfile);
}


/*

 stuff a character string accross the 80-line boundary
 NOTE: the exact implementation of this routines is still controversial

 From: pence@tetra.gsfc.nasa.gov Sat May 15, 1993 10:37 "Continuation Keywords"

 Option 1.   Blank Continuation Keyword Convention

 A continued string value is recognized by a backslash (\)
 as the last character in the first keyword string, followed
 by a keyword with a blank name followed by the continuation
 of the keyword value enclosed in quotes.  Example:

 KEYNAME = 'This is a very long keyword string value which \'
         'continues over several lines\'
                 ' of the FITS header.'

 */
                 
void stuffit(FITS *fitsfile, string fkey, string cp)
{
  char line[81], *hp;
  int i;
  int maxlen = 70;	/* comment field , minus 1 */
  
  hp = cp;
  strncpy(line,hp,maxlen);
  line[maxlen] = 0;
  fitwra(fitsfile,fkey,line);
  while ((int)strlen(hp) > maxlen) {
    hp += maxlen;
    strncpy(line,hp,maxlen);
    line[maxlen] = 0;
    fitwra(fitsfile," ",line);
  }
}

/* refmap stuff */

void set_refmap(string name)
{
  FITS *fitsfile;
  FLOAT tmpr, defval;
  int ndim = 3;
  int naxis[3], tmpi;
  int wcsaxes = -1;

  fitsfile = fitopen(name,"old",ndim,naxis);
  dprintf(0,"[Reading reference map %s [%d,%d,%d]\n",name,naxis[0],naxis[1],naxis[2]);

  /* set defaults according to Greisen & Calabretta 2002 WCS paper-I */
  defval = 1.0;
  fitrdhdr(fitsfile,"CDELT1",&ref_cdelt[0], defval);
  fitrdhdr(fitsfile,"CDELT2",&ref_cdelt[1], defval);
  fitrdhdr(fitsfile,"CDELT3",&ref_cdelt[2], defval);
  /* ieck; what if no CDELT's present, but all in CD matrix */

  defval = 0.0;
  fitrdhdr(fitsfile,"CRPIX1",&ref_crpix[0], defval);
  fitrdhdr(fitsfile,"CRPIX2",&ref_crpix[1], defval);
  fitrdhdr(fitsfile,"CRPIX3",&ref_crpix[2], defval);

  defval = 0.0;
  fitrdhdr(fitsfile,"CRVAL1",&ref_crval[0], defval);
  fitrdhdr(fitsfile,"CRVAL2",&ref_crval[1], defval);
  fitrdhdr(fitsfile,"CRVAL3",&ref_crval[2], defval);

  fitrdhda(fitsfile,"CTYPE1",ref_ctype[0],"");
  fitrdhda(fitsfile,"CTYPE2",ref_ctype[1],"");
  fitrdhda(fitsfile,"CTYPE3",ref_ctype[2],"");

  fitrdhda(fitsfile,"CUNIT1",ref_cunit[0],"");
  fitrdhda(fitsfile,"CUNIT2",ref_cunit[1],"");
  fitrdhda(fitsfile,"CUNUT3",ref_cunit[2],"");

  fitrdhdi(fitsfile,"WCSAXES",&wcsaxes, -1);
  if (wcsaxes != -1) warning("WCSAXES = %d\n",wcsaxes);

  fitclose(fitsfile);


}


/* return a reverse order index array  */
/* if done to an array 1 2 3, it will result in buggy data */
/* this routine should just shift the 1's to the end , not sort */

void permute (int *x ,int *p, int n)
{
  int    i, j, tmp;

  for (i=0; i<n; i++)
    p[i]=i;               /*  one-to-one */
                
  for (j=0; j<n; j++) {
    for (i=1; i<n; i++) {
      if (x[p[i-1]]==1 && x[p[i]] > 1) {
	tmp = p[i];
	p[i]   = p[i-1];
	p[i-1] = tmp;
      }
    }
  }
}


/* right now this is a silly bubble sort !!! */

void permute_old (int *x ,int *idx, int n)
{
    int    gap, i, j;
    int    tmp;
    for (i=0; i<n; i++)
        idx[i]=i;               /*  one-to-one */
                
    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[idx[j]] > x[idx[j+gap]])
                    break;          /* in order */
                tmp = idx[j];
                idx[j]    = idx[j+gap];
                idx[j+gap]= tmp;
            }
}
