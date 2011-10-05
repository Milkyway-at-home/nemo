/*
 * YAPP_CG_NEW.C: Yet Another Plotting Package, w/ Color Graphics.
 * Peter Wisnovsky	   May   1987  I. A. S.  Princeton, NJ.
 * Joshua Barnes	 7 March 1989  I. A. S.  Princeton, NJ.
 *
 * This version works with the SunCore graphics library,
 * and drives a color display through a simple interface.
 *
 *  2-jul-87  added dummy pl_matrix()		[PJT].
 * 15-jul-87  added color-table loading		[JEB].
 *  7-mar-89  added support for rasters		[JEB].
 * 14-mar-90  make GCC happy			[PJT].
 *  9-dec-90  pl_matrix
 */

#include "stdinc.h"
#include <usercore.h>

local struct vwsurf vwsurf;		/* SunCore viewing surface object   */

#define MAXCOLORS  128			/* number of colors available       */

local int ncolors;			/* number of colors defined	    */

local float red[256];			/* RGB color tables		    */
local float blue[256];
local float green[256];
  
local double dxymax;			/* size of user window		    */

local float width, height;		/* size of NDC window		    */

local struct suncore_raster raster;	/* raster object filling window	    */

local bool raster_initialized = FALSE;	/* true once raster is initialized  */

local int  cms_rgbsetup();		/* color setup */

/*
 * PLINIT: initalize the plotting package.
 */

plinit(pltdev, xmin, xmax, ymin, ymax)
char *pltdev;		/* output device name (ignored) */
double xmin, xmax;	/* user plotting area */
double ymin, ymax;
{
    initialize_core(BASIC, SYNCHRONOUS, TWOD);
    my_get_view_surface(&vwsurf, NULL, MAXCOLORS);
    initialize_view_surface(&vwsurf, FALSE);
    select_view_surface(&vwsurf);
    cms_rgbsetup();
    define_color_indices(&vwsurf, 0, ncolors, red, green, blue);
    plcolor(ncolors);				/* init display to "white"  */
    if (ymax - ymin < xmax - xmin) {
	dxymax = xmax - xmin;
	width = 1.0;
	height = (ymax - ymin) / dxymax;
    } else {
	dxymax = ymax - ymin;
	width = (xmax - xmin) / dxymax;
	height = 1.0;
    }
    set_ndc_space_2(width, height);
    set_viewport_2(0.0, width, 0.0, height);
    set_window(xmin, xmax, ymin, ymax);
    set_font(ROMAN);
    set_charprecision(CHARACTER);
    initialize_device(KEYBOARD, 1);
    set_echo(KEYBOARD, 1, 0);		/* no echo */
    create_temporary_segment();
}

/*
 * MY_GET_VIEW_SURFACE: routine to determine proper graphics driver.
 */

#include "my_get_view_surface.c"

/*
 * CMS_RGBSETUP: default color table initialization.
 */

#define	BLACK	0
#define	RED	1
#define	YELLOW	2
#define	GREEN	3
#define	CYAN	4
#define	BLUE	5
#define	MAGENTA	6
#define	WHITE	7

local cms_rgbsetup()
{
    ncolors = WHITE+1;
    red[BLACK]   = 0.0;  green[BLACK]   = 0.0;  blue[BLACK]   = 0.0;
    red[RED]     = 1.0;  green[RED]     = 0.0;  blue[RED]     = 0.0;
    red[YELLOW]  = 1.0;  green[YELLOW]  = 1.0;  blue[YELLOW]  = 0.0;
    red[GREEN]   = 0.0;  green[GREEN]   = 1.0;  blue[GREEN]   = 0.0;
    red[CYAN]    = 0.0;  green[CYAN]    = 0.8;  blue[CYAN]    = 0.8;
    red[BLUE]    = 0.2;  green[BLUE]    = 0.6;  blue[BLUE]    = 1.0;
    red[MAGENTA] = 1.0;  green[MAGENTA] = 0.4;  blue[MAGENTA] = 1.0;
    red[WHITE]   = 1.0;  green[WHITE]   = 1.0;  blue[WHITE]   = 1.0;
}

/*
 * PLCOLOR: specify plotting color as integer between 0 and ncolors-1;
 * values outside this range are maped to the nearest endpoint.
 */

plcolor(color)
int color;
{
    if (color < 0)
	color = 0;
    else if (color > ncolors - 1)
	color = ncolors - 1;
    set_line_index(color);
    set_text_index(color);
}

/*
 * PLNCOLORS: return current value of local variable ncolors.
 */

int plncolors()
{
    return (ncolors);
}

/*
 * PLPALETTE: re-initialize color table from user-supplied values.
 */

plpalette(r, g, b, nc)
real r[], g[], b[];
int nc;
{
    int i;

    if (nc > MAXCOLORS)
	error("plpalette: cannot define more than %d colors\n", MAXCOLORS);
    ncolors = nc;
    for (i = 0; i < ncolors; i++) {
	red[i] = r[i];
	green[i] = g[i];
	blue[i] = b[i];
    }
    red[ncolors] = green[ncolors] = blue[ncolors] = 0.0;
    define_color_indices(&vwsurf, 0, ncolors, red, green, blue);
    plcolor(ncolors - 1);			/* init display to "white"  */
}

/*
 * PLGETRASTER: get raster image of display.
 */

plgetraster(bitsptr, widthptr, heightptr, depthptr)
byte **bitsptr;
int *widthptr;
int *heightptr;
int *depthptr;
{
    if (! raster_initialized) {
	size_raster(&vwsurf, 0.0, width, 0.0, height, &raster);
	allocate_raster(&raster);
	if (raster.bits == NULL)
	    error("plgetraster: allocate_raster [%d][%d][%d] failed\n",
		  raster.width, raster.height, raster.depth);
	raster_initialized = TRUE;
    }
    get_raster(&vwsurf, 0.0, width, 0.0, height, 0, 0, &raster);
    *bitsptr = (byte *) raster.bits;
    *widthptr = raster.width;
    *heightptr = raster.height;
    *depthptr = raster.depth;
}

/*
 * PLPUTRASTER: (re)display raster image from plgetraster.
 */

plputraster()
{
    float winx, winy;

    map_ndc_to_world_2(0.0, 0.0, &winx, &winy);
    move_abs_2(winx, winy);
    put_raster(&raster);
}

/*
 * PLLTYPE: select line width and dot-dash pattern.
 */

plltype(lwid, lpat)
int lwid;		/* line width */
int lpat;		/* line pattern */
{
    if (lwid > 0)
	set_linewidth(0.1 * (lwid - 1));
    if (lpat > 0)
	set_linestyle(lpat - 1);
}

/*
 * PLLINE, PLMOVE, PLPOINT: plot lines, moves, and points.
 */

plline(x, y)
double x, y;		/* user coordinates */
{
    line_abs_2(x, y);
}

plmove(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x, y);
}

plpoint(x, y)
double x, y;		/* user coordinates */
{
    move_abs_2(x, y);
    line_rel_2(0.0, 0.0);
}

/*
 * PLCIRCLE, PLCROSS, PLBOX: plot various symbols.
 */

plcircle(x, y, r)
double x, y;
double r;
{
    int npnts, i;
    double theta, sin(), cos();

    npnts = MAX(2400 * r / dxymax, 6.0);
    plmove(x + r, y);
    for (i = 1; i <= npnts; i++) {
	theta = TWO_PI * ((double) i) / ((double) npnts);
	plline(x + r * cos(theta), y + r * sin(theta));
    }
}

plcross(x, y, s)
double x, y;
double s;
{
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
}

plbox(x, y, s)
double x, y;
double s;
{
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
}

/*
 * PLJUST: specify justification of strings and numbers.
 * Imports: jus: 
 */

static int textjust = -1;

pljust(jus)
int jus;		/* -1, 0, 1 for left, mid, right just */
{
    textjust = (jus < -1 ? -1 : (jus > 1 ? 1 : jus));
}

/*
 * PLTEXT: plot a text string.
 */

pltext(msg, x, y, hgt, ang)
char *msg;		/* message to plot, with NULL termination */
double x, y;		/* user coordinates (modified by justification) */
double hgt;		/* height of characters in user coordinates */
double ang;		/* angle of message, counterclockwise in degrees */
{
    double c, cos(), s, sin();
    float dx, dy;

    c = cos(ang / 57.296);
    s = sin(ang / 57.296);
    set_charpath_2(c, s);
    set_charup_2(- s, c);
    set_charsize(0.75 * hgt, hgt);
    inquire_text_extent_2(msg, &dx, &dy);
    move_abs_2(x - 0.58 * (1 + textjust) * dx,
               y - 0.55 * (1 + textjust) * dy);
				/* note fudge factors in offsets */
    text(msg);
}

/*
 * PL_MATRIX: does nothing, and I dont even know what that aint!
 */

pl_matrix (frame,nx,ny,xmin,ymin,cell,fmin,fmax,findex)
double *frame, xmin, ymin, cell, fmin, fmax, findex;
int nx, ny;
{
    double x,y,f,grayscale,ds,pow();
    int ix,iy;

#if TRUE
    printf("PL_MATRIX is unsupported in suncore version of YAPP_CG\n");
    return (0);
#endif
}

/*
 * PLFLUSH: output any pending graphics.
 */

plflush() { }

/*
 * PLFRAME: advance to next frame.
 */

plframe()
{
    close_temporary_segment();
    new_frame();
    create_temporary_segment();
}

/*
 * PLSTOP: finish up a display.
 */

#define FIVEMIN (5 * 60 * 1000 * 1000)

plstop()
{
    char msg[80];
    int len;

    close_temporary_segment();
    await_keyboard(FIVEMIN, 1, msg, &len);
    deselect_view_surface(&vwsurf);
    terminate_core();
}

#ifdef TESTBED

main(argc, argv)
int argc;
string argv[];
{
    int i, j, nc, rwidth, rheight, rdepth;
    byte *rbits;

    plinit("***", 0.0, 20.0, 0.0, 20.0);
    plmove(0.0, 0.0);
    plline(20.0, 0.0);
    plline(20.0, 20.0);
    plline(0.0, 20.0);
    plline(0.0, 0.0);
    plline(20.0, 20.0);
    plmove(20.0, 0.0);
    plline(0.0, 20.0);
    plltype(12, 0);
    plmove(4.0, 18.0);
    plline(16.0, 18.0);
/*  plltype(-6, 0);  */
    plltype(6, 0);
    plcolor(0);
    plmove(6.0, 18.0);
    plline(14.0, 18.0);
    plcolor(32767);
    plltype(10, 0);
    nc = plncolors();
    for (i = 1; i < nc; i++) {
	plcolor(i);
	plmove(6.0 + 2 * (i - 1.0) / nc, 16.0);
	plline(6.5 + 2 * (i - 1.0) / nc, 17.0);
    }
    for (i = 1; i <= 4; i++) {
	plltype(i, 1);
        plmove(1.0, 13.0 - i);
        plline(3.0, 13.0 - i);
        plpoint(3.5, 13.0 - i);
	plltype(1, i);
	for (j = 1; j <= 4; j++) {
	    plmove(1.5, 13.0 - i - 0.2*j);
	    plline(1.5 + j, 13.0 - i - 0.2*j);
	}
    }
    plltype(1, 1);
    plcircle(15.0, 9.0, 0.5);
    plcircle(16.0, 9.0, 0.25);
    plcircle(17.0, 9.0, 0.125);
    plcircle(18.0, 9.0, 0.0625);
    plbox(16.0, 8.0, 0.4);
    plbox(17.0, 8.0, 0.2);
    plbox(18.0, 8.0, -0.2);
    plcross(16.0, 7.0, 0.4);
    plcross(17.0, 7.0, 0.2);
    plcross(18.0, 7.0, -0.2);
    pltext("Foo Bar!", 8.0, 5.0, 0.5, 0.0);
    pltext("Fum Bar!", 8.0, 3.0, 0.25, 0.0);
    for (i = 0; i <= 4; i++)
	pltext(" testing angles", 16.0, 10.0, 0.2, 45.0*i);
    for (i = 0; i < 32; i++ ) {
	plmove(12.0, 2.0 + i / 31.0);
	plline(13.0, 2.0 + i / 31.0);
    }
    plmove(10.0, 8.5);
    plline(10.0, 11.5);
    pljust(-1);
    pltext("left justified",  10.0,  9.0, 0.25, 0.0);
    pljust(0);
    pltext("centered",        10.0, 10.0, 0.25, 0.0);
    pljust(1);
    pltext("right justified", 10.0, 11.0, 0.25, 0.0);
    sleep(10);
    plgetraster(&rbits, &rwidth, &rheight, &rdepth);
    printf("width: %d  height: %d  depth: %d\n", rwidth, rheight, rdepth);
    for (i = 0; i < 512; i++)
	rbits[2*i] = (byte) i;
    plputraster();
    plstop();
}

#endif
