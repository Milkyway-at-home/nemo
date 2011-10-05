/*
 * YAPP: Yet Another Plotting Package.
 *
 *	5-dec-2003	cloned off yapp_pgplot.c
 */

#include <stdinc.h>
#include <yapp.h>
#include <pgplot/cpgplot.h>     /* we keep a version in NEMO */

extern string yapp_string;	/* a kludge, see: getparam.c */
extern int debug_level;         /* see dprintf.c   from DEBUG env.var. */

#define COLOR

local real   dxymax;    /* size of user window */
local int    iterm;     /* terminal number */

#ifdef COLOR
# define MAXCOLOR 256
  local int ncolors=0;
  local float red[MAXCOLOR];		/* RGB color tables		    */
  local float blue[MAXCOLOR];
  local float green[MAXCOLOR];
  local void cms_rgbsetup();
#else
# define MAXCOLOR 0
#endif

/*
 * PLINIT: initalize the plotting package.
 */

int plinit(string pltdev, real xmin, real xmax, real ymin, real ymax)
{
    float width, height, x1,x2,y1,y2, zero, one;
    int   nx, ny, ask, units, dummy = 0;
    
    iterm = 1;
    if (yapp_string == NULL || *yapp_string == 0) {
	if (pltdev == NULL || *pltdev == 0)
        	iterm = 0;
	else
	    yapp_string = pltdev;
    }

    nx = ny = 1;        /* only one window on the page */
    if (iterm) {
      iterm = cpgbeg(dummy,yapp_string, nx, ny)==1;
    } else {
      iterm = cpgbeg(dummy,"?", nx, ny)==1;
    }
    if (iterm==0) return 0;
#if 1
    ask = 0;        /* Logical = int in cpgplot */
    cpgask(ask);    /* here we want ask=FALSE */
#endif

    units = 2;
    cpgqvsz(units, &x1, &x2, &y1, &y2);
    dprintf(1,"PGQVSZ: X= %g - %g Y= %g - %g\n",x1,x2,y1,y2);
    zero = 0.0;  one = 1.0;    /* zero is good for "max size of device" */
    cpgsvp(zero,one,zero,one); /* xleft,xright,ybot,ytop */
#if 0

	/* it's better to use the PGPLOT_PS_ * environment variables */
    zero = 7.874;	       /* 20.0 cm is good for most YAPP devices */
    zero = 9.874;
#endif    
    cpgpap(zero,one);    /* (width,aspect) set size and make it come out square  */

    x1=xmin;  x2=xmax;
    y1=ymin;  y2=ymax;
    cpgswin(x1,x2,y1,y2);         /* set the viewport */

    if (ymax - ymin < xmax - xmin) {		/* not used for now */
        dxymax = xmax - xmin;
        width = 1.0;
        height = (ymax - ymin) / dxymax;
    } else {
        dxymax = ymax - ymin;
        width = (xmax - xmin) / dxymax;
        height = 1.0;
    }
#if defined(COLOR)
    cms_rgbsetup();
    plcolor(1);     /* set initial color to the forground color */
#endif
    return 0;
}

/*
 * PLSWAP: does nothing.
 */

int plswap() { return 0; }

/*
 * PLXSCALE, PLYSCALE: transform from usr to plotting coordinates.
 * At present, these do nothing (identity transformation).
 */

real plxscale(real x, real y)
{
    return x;
}

real plyscale(real x, real y)
{
    return y;
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

int plltype(int lwid, int lpat)
{
    int lw, ls;

    if (iterm==0) return 0;       /* no graphics output requested */

    if (lwid > 0) {
        lw = lwid;
	cpgslw(lw);			/* set line width */
    }
    if (lpat > 0) {
        ls = lpat;
	cpgsls(ls);			/* set line style */
    }
    return 0;
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

int plline(real x, real y)
{
    float xp,yp;

    if (iterm==0) return 0;       /* no graphics output requested */

    xp=x; yp=y;         /* RECALC !! */
    cpgdraw(xp,yp);
    return 0;
}

int plmove(real x, real y)
{
   float xp,yp;

   if (iterm==0) return 0;       /* no graphics output requested */

   xp=x; yp=y;          /* RECALC !! */
   cpgmove(xp,yp);
   return 0;
}

int plpoint(real x, real y)
{
    int ipoint=-1, npoint=1;
    float xp,yp;
    
    if (iterm==0) return 0;       /* no graphics output requested */

    xp=x; yp=y;                             /* RECALC !! */
    cpgpt(npoint, &xp, &yp, ipoint);       /* draw 1 dot */
    return 0;
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

int plcircle(real x, real y, real r)
{
    int npnts, i;
    real theta;

    if (iterm==0) return 0;       /* no graphics output requested */

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
        theta = TWO_PI * ((real) i) / ((real) npnts);
        plline(x + r * cos(theta), y + r * sin(theta));
    }
    return 0;
}

int plcross(real x, real y, real s)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    if (s > 0.0) {
        plmove(x - s, y);
        plline(x + s, y);
        plmove(x, y - s);
        plline(x, y + s);
    } else {
        s = s / 1.4142;
        plmove(x - s, y - s);
        plline(x + s, y + s);
        plmove(x - s, y + s);
        plline(x + s, y - s);
    }
    return 0;
}

int plbox(real x, real y, real s)
{
    if (iterm==0) return 0;       /* no graphics output requested */

    if (s > 0.0) {
        plmove(x - s, y - s);
        plline(x + s, y - s);
        plline(x + s, y + s);
        plline(x - s, y + s);
        plline(x - s, y - s);
    } else {
        s = s * 1.4142;
        plmove(x - s, y);
        plline(x, y - s);
        plline(x + s, y);
        plline(x, y + s);
        plline(x - s, y);
    }
    return 0;
}

/*
 * PLJUST: specify justification of strings and numbers.
 * Imports: jus: 
 */

static float fjust = 0.0;   /* pgplot default: left justified */

int pljust(int jus)       /* -1, 0, 1 for left, mid, right just */
{
    if (iterm==0) return 0;       /* no graphics output requested */

    fjust = (jus < -1 ? 0.0 : (jus > 1 ? 1.0 : 0.5));
    return 0;
}

/*
 * PLTEXT: plot a text string.
 */

int pltext(string msg, real x, real y, real hgt, real ang)
{
    float xp, yp, ap;
    float newsize;

    if (iterm==0) return 0;       /* no graphics output requested */

    xp=x; yp=y; ap=ang;         /* copy into local variables */

    newsize = 2 * hgt;          /* pgplot scaling factor */        
    cpgsch(newsize);           /* set height of char */

    cpgptxt(xp, yp, ap, fjust, msg);        /* plot it */
    return 0;
}

/*
 * PLFLUSH: output any pending graphics.
 */

int plflush() 
{ 
    if (iterm==0) return 0;

    cpgupdt();
    return 0;
}

/*
 * PLFRAME: advance to next frame.
 */

int plframe()
{
    if (iterm==0) return 0;       /* no graphics output requested */

    cpgpage();
    return 0;
}

/*
 * PLSTOP: finish up a display.
 */

#define ONEMIN (60 * 1000 * 1000)

int plstop()
{
    if (iterm==0) return 0;       /* no graphics output requested */


    if (debug_level > 0)
        cpgiden();    

    cpgend();
    return 0;
}

int pl_matrix(real *frame,int nx,int ny,real xmin,real ymin,
	      real cell,real fmin,real fmax,real findex, real blankval)
{
    int ix,iy,ix0,ix1,iy0,iy1;
    float gray0, gray1, tr[6], *data, *dp;
    real dval,dfac;
    
    if (iterm==0) return 0;       /* no graphics output requested */

    dprintf(0,"PL_MATRIX development version for YAPP_PGPLOT\n");
    if (findex < 0.0) {
        dprintf(0,"Swapping min and max to : %g %g\n",fmax,fmin);
        findex = -findex;
        dval = fmin;
        fmin = fmax;
        fmax = dval;
    }
    if (fmin==fmax) {
        warning("Cannot greyscale a uniform image yet");
        return 1;
    } else {
        dfac = 1.0/(fmax-fmin);
    }
    ix0 = 1;
    iy0 = 1;
    ix1 = nx;
    iy1 = ny;
    gray0 = 0.0;
    gray1 = 1.0;
    tr[2] = tr[4] = 0.0;    /* off axis matrix scaling */
    tr[1] = 16.0/nx;        /* linear scaling to make it fit in 16x16cm */
    tr[5] = 16.0/ny;
    tr[0] = 2.0 - 8.0/nx;   /* offset to make pgplot fit the ll corner */
    tr[3] = 2.0 - 8.0/ny;
    data = (float *) allocate(sizeof(float)*nx*ny);
    for(iy=0, dp=data; iy<ny; iy++)
        for(ix=0; ix<nx; ix++) {
            dval =  *(frame + ix*ny + iy);
	    if (dval == blankval) {
	      *dp++ = pow(fmin,findex)-1;   /* could it be +1 if fmin>fmax ?? */
            } else {
	      if (fmin < fmax) {
                if (dval < fmin) dval=fmin;
                if (dval > fmax) dval=fmax;
                dval = (dval-fmin)*dfac;
	      } else {
                if (dval < fmax) dval=fmax;
                if (dval > fmin) dval=fmin;
                dval = (dval-fmax)*dfac + 1.0;
	      }
	      *dp++ = pow(dval,findex);
	    }
        }
    cpggray(data,nx,ny,ix0,ix1,iy0,iy1,gray1,gray0,tr);
    free(data);
    return 0;
}

/* not functional yet */

int pl_contour(real *frame,int nx,int ny, int nc, real *c)
{
    int ix,iy,ix0,ix1,iy0,iy1;
    float tr[6], *data, *dp, *cnt;
    
    if (iterm==0) return 0;       /* no graphics output requested */

    dprintf(0,"PL_CONTOUR development version for YAPP_PGPLOT\n");
    ix0 = 1;
    iy0 = 1;
    ix1 = nx;
    iy1 = ny;
    tr[2] = tr[4] = 0.0;    
    tr[1] = 16.0/nx;
    tr[5] = 16.0/ny;
    tr[0] = 2.0 - 8.0/nx;
    tr[3] = 2.0 - 8.0/ny;
    data = (float *) allocate(sizeof(float)*nx*ny);
    for(iy=0, dp=data; iy<ny; iy++)
        for(ix=0; ix<nx; ix++)
            *dp++ = *(frame + ix*ny + iy);
    cnt = (float *) allocate(sizeof(float)*nc);
    for(ix=0; ix<nc; ix++)
        cnt[ix] = c[ix];
    cpgcont(data,nx,ny,ix0,ix1,iy0,iy1,cnt,nc,tr);
    free(data);
    free(cnt);
    return 0;
}


int pl_screendump(string fname)
{
  printf("pl_screendump(%s): Not implemented for yapp_pgplot\n",fname);
  return 0;
}

int local bell()
{
    int ring=7;
    putchar(ring);
    putchar('\n');		/* send a line feed to flush the buffer */
    return 0;
}

int pl_getpoly(float *x, float *y, int n)
{
    int nn, k, symbol;
    float xold,yold, xnew,ynew;
    char ch[10];

    bell();

    printf("Define a polygon:\n");
    printf("   LEFT   = define first/next vertex of polygon\n");
    printf("   MIDDLE = delete previous vertex point from polygon\n");
    printf("	    (this option cannot remove already plotted vertex lines\n");
    printf("   RIGHT  = close polygon\n");

    nn = 0;                           /* count points in polygon */
    for(;;) {
        k = cpgcurs(&xnew, &ynew, ch);   /* check if length of 'ch' ok */
        if(k==0) {
            warning("Device has no cursor...");
            return 0;
        }
        if (ch[0]=='A') {                     /* left button = DEFINE POLYGON */
            xold = xnew;
            yold = ynew;
            dprintf (2,"Button A at x=%f y=%f\n",xold,yold);
            if (nn==0) {
                cpgmove(xold,yold);
                k=1;
                symbol=2;
                cpgpt(k,&xold,&yold,symbol);
#if 0
                set_marker_symbol(43);      /* a 'plus' symbol */
                marker_abs_2(xold,yold);
        
#endif
            } else
                cpgdraw(xold,yold);
            if (nn<n) {
               x[nn] = xold;
               y[nn] = yold;
               nn++;
            } else
                warning("No more space to store this data");
        } else if (ch[0]=='D') {                  /* middle  = DELETE PREVIOUS POINT */
            if (nn>0)
                nn--;
            if (nn!=0)
               cpgmove(x[nn-1],y[nn-1]);    /* reset */
        } else if (ch[0]=='X') {                          /* right = QUIT */
            break;
        } else
            warning("Unknown return char %c from PGCURSE (expected A,D,X)",ch[0]);
    }   /* for(;;) */
    dprintf (2,"Button C pressed to finish up polygon (%d)\n",nn);
    if (nn>0) {
        cpgdraw(x[0],y[0]);       /* close polygon */
    }
    return(nn<3 ? 0 : nn);        
}

int pl_cursor(real *x, real *y, char *c)
{
    char inf[8], ans[8];
    int len, inf_len, ans_len;
    permanent float xsave, ysave;
                
    strcpy(inf,"CURSOR");
    inf_len = strlen(inf);
    ans_len = 1;
    cpgqinf(inf, ans, &len);
    if (ans[0]=='n' || ans[0]=='N') return 0;
    cpgcurs(&xsave, &ysave, c);
    *x = xsave;
    *y = ysave;
    return 1;
}


/*  The rest of this file is for optional color support */
#if defined(COLOR)
/*
 * CMS_RGBSETUP: default color table initialization.
 */

#define	BLACK		0
#define RED		1
#define GREEN		2
#define BLUE		3
#define CYAN 		4	/* (Green + Blue) */	
#define MAGENTA 	5	/* (Red + Blue) */	
#define YELLOW		6	/*   (Red + Green) */	
#define ORANGE		7
#define GREEN_YELLOW	8
#define GREEN_CYAN	9
#define BLUE_CYAN	10
#define BLUE_MAGENTA	11
#define RED_MAGENTA	12
#define DARK_GRAY	13
#define LIGHT_GRAY	14
#define	WHITE   	15

local void cms_rgbsetup()
{
  ncolors = WHITE+1;
  /* default PGPLOT  colors: although, defined here, plpalette is not called */
  red[BLACK]=0.00;         green[BLACK]=0.00;          blue[BLACK]=0.00;
  red[RED]=1.00;           green[RED]=0.00;            blue[RED]=0.00;
  red[GREEN]=0.00;         green[GREEN]=1.00;          blue[GREEN]=0.00;
  red[BLUE]=0.00;          green[BLUE]=0.00;           blue[BLUE]=1.00;
  red[CYAN]=0.00;          green[CYAN]=1.00;           blue[CYAN]=1.00;
  red[MAGENTA]=1.00;       green[MAGENTA]=0.00;        blue[MAGENTA]=1.00;
  red[YELLOW]=1.00;        green[YELLOW]=1.00;         blue[YELLOW]=0.00;
  red[ORANGE]=1.00;        green[ORANGE]=0.50;         blue[ORANGE]=0.00;
  red[GREEN_YELLOW]=0.50;  green[GREEN_YELLOW]=1.00;   blue[GREEN_YELLOW]=0.00;
  red[GREEN_CYAN]=0.00;    green[GREEN_CYAN]=1.00;     blue[GREEN_CYAN]=0.50;
  red[BLUE_CYAN]=0.00;     green[BLUE_CYAN]=0.50;      blue[BLUE_CYAN]=1.00;
  red[BLUE_MAGENTA]=0.50;  green[BLUE_MAGENTA]=0.00;   blue[BLUE_MAGENTA]=1.00;
  red[RED_MAGENTA]=1.00;   green[RED_MAGENTA]=0.00;    blue[RED_MAGENTA]=0.50;
  red[DARK_GRAY] =0.33;    green[DARK_GRAY]=0.33;      blue[DARK_GRAY]=0.33;
  red[LIGHT_GRAY]=0.66;    green[LIGHT_GRAY]=0.66;     blue[LIGHT_GRAY]=0.66;
  red[WHITE]=1.00;         green[WHITE]=1.00;          blue[WHITE]=1.00;
}

/*
 * PLCOLOR: specify new plotting color as an integer between 0 and ncolors-1;
 * values outside this range are mapped to the nearest endpoint.
 */

void plcolor(int color)
{
    if (color < 0)
	color = 0;
    else if (color > ncolors - 1)
	color = ncolors - 1;
    cpgsci(color);
}

/*
 * PLNCOLORS: return current value of local variable ncolors.
 */

int plncolors()
{
    return ncolors;
}

/*
 * PLPALETTE: re-initialize color table from user-supplied values.
 */

void plpalette(real *r, real *g, real *b, int nc)
{
    int i;

    if (nc > MAXCOLOR)
	error("plpalette: cannot define more than %d colors", MAXCOLOR);
    ncolors = nc;
    for (i = 0; i < ncolors; i++) {
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
	dprintf(1,"->PGSCR_(%d,%g,%g,%g) \n",i,red[i],green[i],blue[i]);
        cpgscr(i,red[i],green[i],blue[i]);
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;    /* terminate */
    plcolor(1);			/* reset default color to foreground */
    dprintf(0,"Setting default color to: %g %g %g\n",
            red[1], green[1], blue[1]);
}

/*
 * PLLUT:  read a new RGB palette from a LUT color table
 *	   The lut must be a simple ascii file, with 1st, 2nd and 3rd column being
 *	   the R, G and B response, a real number between 0.0 and 1.0
 */
#define NOCOLOR(x)  ((x)<0||(x)>1)
void pllut(string fname, bool compress)
{
    stream cstr;
    int ncolors, nskip=0;
    real red[MAXCOLOR], green[MAXCOLOR], blue[MAXCOLOR];
    float r, g, b, r_old, g_old, b_old;
    char line[256];

    if (fname==NULL || *fname == 0) return;
    cstr = stropen(fname, "r");
    ncolors = 0;
    while (fgets(line,256,cstr)) {
       if(ncolors>=MAXCOLOR) error("(%d/%d): Too many colors in %s",
              ncolors+1, MAXCOLOR, fname);
       sscanf(line,"%f %f %f",&r, &g, &b);
       if (NOCOLOR(r) || NOCOLOR(g) || NOCOLOR(b)) {
          warning("Skipping RGB=%f %f %f",r,g,b);
	  continue;
       }
       if (ncolors>0 && r==r_old && g==g_old && b==b_old && compress) {
          nskip++;
       } else {
          dprintf(1,"LUT(%d): %f %f %f\n",ncolors,r,g,b);
          red[ncolors]   = r;   r_old = r;
          green[ncolors] = g;   g_old = g;
          blue[ncolors]  = b;   b_old = b;
          ncolors++;
       }
    }
    strclose(cstr);
    dprintf(0,"Colortable %s: %d entries; %d redundant entries skipped\n",
            fname, ncolors,nskip);
    plpalette(red, green, blue, ncolors);
}

#endif
