/* xyio: exerimental MIRIAD-like line-access to 2D and 3D NEMO images 
 *
 *	Although the format can work on full N-dimensional data, 
 *	only 2- and 3D is currently supported.
 *      Data is accessed in slices using the random access I/O in
 *	filestruct.
 *	All ideas are based on the image I/O within miriad (xyio.c)
 *	for which Bob Sault deserves all the credit.
 *
 *  Note: xyio in NEMO is  0 based, not 1 based, as is in the Miriad
 *        package
 *
 *	may 1992 - created				pjt
 *	jul 1993 - added Unit()				pjt
 *	feb 1994 - ansi					pjt
 *	apr 1995 - no more ARGS, fix up prototypes
 *    6-jan-1998 - fixed up TESTBED bit more            pjt
 */

#include <stdinc.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>

#define DLEV   5		/* local default debug output level */


/*	storage of matrices can be done in several ways: */

local char *matdef[] = {"ForDef", "CDef", NULL };

/* 
 * one of the above XXXDEF's must be defined below, and
 * set the appropriate Tag value in the Image descriptor
 * If it is not done, the compiled versions will (should)
 * complain about missing integer 'idef'
 */

#if defined(FORDEF)
local int idef = 0;
#endif
#if defined(CDEF)
local int idef = 1;
#endif

/* The xyio (MIRIAD) interface handles images in the following way:
 * [this discussion ignores the mask data, they are ignored in NEMO now]
 *
 *	xyopen (handle, name, status, naxis, axes[])
 *	xyclose(handle) 
 *	xyread (handle,index,array)
 *	xywrite(handle,index,array)
 *	xysetpl(handle,naxis,axes[])
 *	
 */

#define PUT 1
#define GET 2
#define MAXNAXIS 3   /* for now !!! */
#define MAXOPEN 10

local struct { 
    stream str;
    image *iptr;
    int naxis;
    int axes[MAXNAXIS];
    int offset;
    int access;     /* PUT or GET ; maybe both in future for inplace mods */
} images[MAXOPEN];

local bool first=TRUE;


local void xy_init(void);
local void xy_handle(int, string);


image *xyopen(int *handle, string name, string status, int naxis, int *axes)
{
    int i, access;
    string read_matdef;
    image *iptr;
    stream str;

    if(naxis>MAXNAXIS) 
        error("naxis=%d not supported: MAXNAXIS=%d",naxis, MAXNAXIS);

    if (first) xy_init();

    str = stropen(name,status);                       /* open file */
    switch (*status) {
        case 'r':               /* "r", "old" */
        case 'o':
            access = GET;
            break;
        case 'w':               /* "w", "new" */
        case 'n':
            access = PUT;
            break;
        default:
            error("xyopen: Unsupported mode %s",status);
    }
    iptr = (imageptr )allocate(sizeof(image));  /* new image */

    if(access == GET){
      get_history(str);
      get_set (str,ImageTag);
        get_set (str,ParametersTag);
          get_data (str,NxTag,IntType, &(Nx(iptr)), 0);
          get_data (str,NyTag,IntType, &(Ny(iptr)), 0);
          get_data (str,NzTag,IntType, &(Nz(iptr)), 0);
          get_data_coerced (str,XminTag,RealType, &(Xmin(iptr)), 0);
          get_data_coerced (str,YminTag,RealType, &(Ymin(iptr)), 0);
          get_data_coerced (str,ZminTag,RealType, &(Zmin(iptr)), 0);
          get_data_coerced (str,DxTag,RealType, &(Dx(iptr)), 0);
          get_data_coerced (str,DyTag,RealType, &(Dy(iptr)), 0);
          get_data_coerced (str,DzTag,RealType, &(Dz(iptr)), 0);
	  get_data_coerced (str,MapMinTag, RealType, &(MapMin(iptr)), 0);
	  get_data_coerced (str,MapMaxTag, RealType, &(MapMax(iptr)), 0);
	  get_data (str,BeamTypeTag, IntType, &(BeamType(iptr)), 0);
	  get_data_coerced (str,BeamxTag, RealType, &(Beamx(iptr)), 0);
	  get_data_coerced (str,BeamyTag, RealType, &(Beamy(iptr)), 0);
	  get_data_coerced (str,BeamzTag, RealType, &(Beamz(iptr)), 0);
          if (get_tag_ok(str,NamexTag))             /* X-axis name */
            Namex(iptr) = get_string(str,NamexTag);
          else
            Namex(iptr) = NULL;
          if (get_tag_ok(str,NameyTag))             /* Y-axis name */
            Namey(iptr) = get_string(str,NameyTag);
          else
            Namey(iptr) = NULL;
          if (get_tag_ok(str,NamezTag))             /* Z-axis name */
            Namez(iptr) = get_string(str,NamezTag);
          else
            Namez(iptr) = NULL;
          if (get_tag_ok(str,UnitTag))             /* units  */
            Unit(iptr) = get_string(str,UnitTag);
          else
            Unit(iptr) = NULL;
          read_matdef = get_string(str,StorageTag);
	  if (!streq(read_matdef,matdef[idef]))
             dprintf(0,"read_image: StorageTag = %s, compiled with %s\n",
		        read_matdef, matdef[idef]);
        get_tes(str,ParametersTag);
        get_set(str,MapTag);
        if(Nz(iptr)<=1)
            get_data_set(str,MapValuesTag,RealType,Nx(iptr),Ny(iptr),0);
        else
            get_data_set(str,MapValuesTag,RealType,Nx(iptr),Ny(iptr),Nz(iptr),0);
        for (i=0; i<naxis; i++) axes[i] = 1;
        axes[0] = Nx(iptr); axes[1] = Ny(iptr); axes[2] = Nz(iptr);
    } else { /* PUT */
      Nx(iptr) = Ny(iptr) = Nz(iptr) = 1;
      if (naxis>0) Nx(iptr)  = axes[0];
      if (naxis>1) Ny(iptr)  = axes[1];
      if (naxis>2) Nz(iptr)  = axes[2];

      put_history(str);
      put_set (str,ImageTag);
        put_set (str,ParametersTag);
          put_data (str,NxTag,  IntType,  &(Nx(iptr)),   0);
          put_data (str,NyTag,  IntType,  &(Ny(iptr)),   0);
          put_data (str,NzTag,  IntType,  &(Nz(iptr)),   0);
          put_data (str,XminTag,RealType, &(Xmin(iptr)), 0);
          put_data (str,YminTag,RealType, &(Ymin(iptr)), 0);
          put_data (str,ZminTag,RealType, &(Zmin(iptr)), 0);
          put_data (str,DxTag,  RealType, &(Dx(iptr)),   0);
          put_data (str,DyTag,  RealType, &(Dy(iptr)),   0);
          put_data (str,DzTag,  RealType, &(Dz(iptr)),   0);
          put_data (str,MapMinTag, RealType, &(MapMin(iptr)), 0);
          put_data (str,MapMaxTag, RealType, &(MapMax(iptr)), 0);
          put_data (str,BeamTypeTag, IntType, &(BeamType(iptr)), 0);
          put_data (str,BeamxTag, RealType, &(Beamx(iptr)), 0);
          put_data (str,BeamyTag, RealType, &(Beamy(iptr)), 0);
          put_data (str,BeamzTag, RealType, &(Beamz(iptr)), 0);
          if (Namex(iptr))
            put_string (str,NamexTag,Namex(iptr));
          if (Namey(iptr))
            put_string (str,NameyTag,Namey(iptr));
          if (Namez(iptr))
            put_string (str,NamezTag,Namez(iptr));
      	  if (Unit(iptr))
            put_string (str,UnitTag,Unit(iptr));
          put_string(str,StorageTag,matdef[idef]);
        put_tes(str, ParametersTag);
        put_set(str, MapTag);
        if(Nz(iptr)<=1)
          put_data_set(str,MapValuesTag,RealType,Nx(iptr),Ny(iptr),0);
        else
          put_data_set(str,MapValuesTag,RealType,Nx(iptr),Ny(iptr),Nz(iptr),0);
    }

    *handle = -1;
    for(i=0; i<MAXOPEN; i++) {        /* look for a new table entry */
      if(images[i].str == NULL) *handle = i;
    }
    if(*handle < 0) 
        error("xyopen: No more free slots; too many open images");
    for (i=0; i<MAXNAXIS; i++)
        images[*handle].axes[i] = 1;

    images[*handle].str     = str;
    images[*handle].iptr    = iptr;
    images[*handle].offset  = 0;
    images[*handle].access  = access;
    images[*handle].naxis   = (Nz(iptr)<=1 ? 2 : 3);
    images[*handle].axes[0] = Nx(iptr);
    images[*handle].axes[1] = Ny(iptr);
    images[*handle].axes[2] = Nz(iptr);

    return iptr;
}

void xyclose(int handle)
{
    xy_handle(handle,"xyclose");

    switch(images[handle].access) {
      case PUT:
        put_data_tes(images[handle].str,MapValuesTag);
        put_tes(images[handle].str,MapTag);
        put_tes(images[handle].str,ImageTag);
        break;
      case GET:
        get_data_tes(images[handle].str,MapValuesTag);
        get_tes(images[handle].str,MapTag);
        get_tes(images[handle].str,ImageTag);
        break;
      default:
        error("xyclose: Illegal access mode %d",images[handle].access);
    }
    strclose(images[handle].str);
    images[handle].str=NULL;
    
}

void xyread(int handle, int index, real *array)
{
  int offset,length;

  xy_handle(handle,"xyread");
  length = images[handle].axes[0];
  offset = images[handle].offset + index * length;

  get_data_ran(images[handle].str,MapValuesTag,array,offset,length);
}

void xywrite(int handle, int index, real *array)
{
  int offset,length;

  xy_handle(handle,"xywrite");
  length = images[handle].axes[0];
  offset = images[handle].offset + index * length;
  put_data_ran(images[handle].str,MapValuesTag,array,offset,length);
}

void xysetpl(int handle, int naxis, int *axes)
{
  int size,i;

  xy_handle(handle,"xysetpl");
  if(naxis+2 > MAXNAXIS)
     error("xysetpl: Too many dimensions: MAXNAXIS=%d naxis+2=%d",
                MAXNAXIS, naxis+2);
  size = 0;
  for(i=naxis-1; i >= 0; i--){
    if(axes[i] < 0 || axes[i] >= images[handle].axes[i+2])
	error("xysetpl: Dimension error naxis[%d]=%d",i,axes[i]);
    size = (size + axes[i]) * images[handle].axes[i+1];
  }
  images[handle].offset = size * images[handle].axes[0];
}

/*
 *      initializing local data structures
 */
local void xy_init()
{
    int i;

    if (!first) error("xy_init: cannot initialize twice");

    first = FALSE;
    for (i=0; i<MAXOPEN; i++)       /* clear table */
        images[i].str = NULL;
}

/*
 *      double check validity and error out if a bad handle
 */
local void xy_handle(int h, string msg)
{
    if (h<0 || h>=MAXOPEN) 
        error("%s: illegal handle %d",msg, h);
}


#ifdef TESTBED

#include <getparam.h>

#define N 10

string defv[] = {
        "in=\n      Input image file to be listed",
        "out=\n     Output image file",
	"size=4,4\n Size of output map",
	"fmt=%g\n   Format used in output",
	"VERSION=1.2\n	5-jan-99 PJT",
	NULL
};

string usage = "testbed for xyio image I/O";
	
nemo_main()
{
  int i, j, k, naxes[3], handle, naxis;
  real *data;
  char fmt[20];

  strcpy(fmt, getparam("fmt"));
  strcat(fmt, " ");
  if (hasvalue("in")) {
    xyopen(&handle,getparam("in"),"r",3,naxes);
    dprintf(0,"Image size: %d * %d * %d\n",naxes[0],naxes[1],naxes[2]);
    data = (real *) allocate(sizeof(real)*naxes[0]);
    for(k=0; k<naxes[2];k++) {
        xysetpl(handle,1,&k);
        printf("plane: %d\n",k);
        for(j=naxes[1]-1;j>=0; j--) {
            xyread(handle,j,data);
            printf("%d:",j);
            for (i=0; i<naxes[0]; i++)
                printf(fmt,data[i]);
            printf("\n");
        }
    }
    xyclose(handle);
  } else if (hasvalue("out")) { 
    naxis = nemoinpi(getparam("size"),naxes,3);
    if (naxis<1) error("Bad format n=%s\n",getparam("size"));
    for (i=naxis; i<3; i++) naxes[i] = 1;
    xyopen(&handle,getparam("out"),"w",naxis,naxes);
    data = (real *) allocate(sizeof(real)*naxes[0]);
    for(k=0; k<naxes[2];k++) {
      xysetpl(handle,1,&k);
      for(j=0; j<naxes[1]; j++) {
        for (i=0; i<naxes[0]; i++) data[i] = i+10*j+100*k;
        xywrite(handle,j,data);
      }
    }
    xyclose(handle);
  } else
    warning("Either in= or out= needs to used");
}

#endif

