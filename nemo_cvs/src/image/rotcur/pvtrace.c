/* 
 * PVTRACE:   PV diagram envelope tracing
 *
 *	5-may-01	created			PJT
 *      6-may-01    V1.1    added clip=, center= and linear interpolation
 *                          to get the trace velocity
 *                    a     fixed RPD !!
 *     29-dec-01      b     added extra security
 *
 *     13-jul-02    V1.2    If input file is a cube, reduce the cube to a velocity field
 *     18-mar-09    V1.3    Added rotcur=
 *      8-jun-09    V1.4    Gentile et al's. MET - not very robust yet
 *
 * Shane, W.W. \& Bieger-Smith, G.P. 1966 B.A.N. 18, 263
 *      v_x = v_m + 1/Y_m \int_(v_M)^(V+T) { Y(v) dv }
 *  or
 *      v_x = v_M - DV/2 + 1/YM sum { Y(v) DV }
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
  "in=???\n       Input image file - must be an PV diagram",
  "eta=0.2\n	  [0..1] weight factor for I_max adding to I_lc",
  "ilc=0\n        Lowest Contour value",
  "sign=1\n       1: Rotation largest for positive coordinates, -1: reverse",
  "sigma=0\n	  Velocity dispersion correction factor",
  "clip=\n        One (for -clip:clip) or Two clipping values, else includes all",
  "vsys=0\n       Systemic velocity",
  "inc=90\n       Inclination of disk (ignored for cubes)",
  "center=0\n     Center of galaxy along position axis",
  "rotcur=f\n     Fold all numbers positive so two rotation curves can be overplotted",
  "out=\n         Output velocity field if input was cube",
  "mode=2\n       Tracing mode to output for cube->map  (2=ET, 3=moment, 4=peak etc.)",
  "VERSION=1.4a\n 10-jun-09 PJT",
  NULL,
};

string usage="PV diagram envelope tracing ";

string cvsid="$Id: pvtrace.c,v 1.9 2009/06/10 17:46:40 pteuben Exp $";


#define RPD (PI/180.0)
#ifndef HUGE
#define HUGE 1.0e20
#endif

local void pv_trace(imageptr iptr, int vsign, 
		    real eta, real ilc, real sigma, real vsys, real sini, real psys,
                    real *clip, bool Qrotcur);
local void xyv_trace(imageptr iptr, int vsign, 
		    real eta, real ilc, real sigma, real xsys, real ysys, real vsys,
		    real *clip, int mode);
local int get_vels(int n, real *s, real v, real dv, real it, real *clip,
		   real *vels);

local real gfit1(int n, real *s, real *clip, int ipeak, real sig, int mode);

local real qabs(real x, bool Q);


nemo_main()
{
    stream  instr;
    imageptr iptr = NULL;
    real vel, eta, ilc, sigma, sini, vsys, psys, clip[2], center[2], xsys, ysys;
    int nc, vsign = getiparam("sign");
    bool Qrotcur = getbparam("rotcur");
    int mode = getiparam("mode");

    instr = stropen(getparam("in"), "r");     /* get file name and open file */
    read_image( instr, &iptr);                /* read image */
    strclose(instr);                          /* close file */

    eta = getdparam("eta");
    ilc = getdparam("ilc");
    sigma = getdparam("sigma");
    vsys = getdparam("vsys");
    sini = sin(RPD * getdparam("inc"));
    nc = nemoinpd(getparam("clip"),clip,2);
    if (nc == 1) {
      clip[1] = clip[0];
      clip[0] = -clip[1];
    } else if (nc == 0) {
      clip[0] = clip[1] = 0.0;
    } else if (nc != 2)
      error("Syntax error; need 0, 1 or 2 numbers for clip=%s",
	    getparam("clip"));

    if (Nz(iptr)==1) {
      psys = getdparam("center");
      pv_trace(iptr, vsign, eta, ilc, sigma, vsys, sini, psys, clip, Qrotcur);
    } else {
      warning("cube mode");
      nc = nemoinpd(getparam("center"),center,2);
      if (nc == 2) {
	xsys = center[0];
	ysys = center[1];
      } else if (nc == 1) {
	xsys = ysys = center[0];
      } else if (nc == 0) {
	xsys = ysys = 0.0;
      } else if (nc != 2)
	error("Syntax error; need 0, 1 or 2 numbers for center=%s",
	      getparam("center"));      
      xyv_trace(iptr, vsign, eta, ilc, sigma, xsys, ysys, vsys, clip, mode);
    }
}

/*
 * TRACE:  trace around in PV diagram
 *
 */

#define MAXV  32

local void pv_trace(imageptr iptr, int vsign, 
		    real eta, real ilc, real sigma, real vsys, real sini, real psys,
		    real *clip, bool Qrotcur)
{
  int  ix, iy, j, nx, ny, nv;
  real pos, v0, dv, it, imax, vel_corr;
  real *spec, *vel, vels[MAXV];

    
  nx = Nx(iptr);     /* assumed to be position for now */
  ny = Ny(iptr);     /* assumed to be velocity for now */
  spec = (real *) allocate(ny*sizeof(real));
  vel  = (real *) allocate(ny*sizeof(real));

  imax = MapMax(iptr);

  it = sqrt(sqr(eta*imax)+sqr(ilc));
  dprintf(0,"Map [%d POS x %d VEL] I_t=%g vsys=%g sini=%g\n",
	  nx,ny,it,vsys,sini);
  if (it<=0.0) error("I_t = %g too small",it);

  for (ix=0; ix<nx; ix++) {                 /* loop over all positions in P's */
    pos = ix*Dx(iptr) + Xmin(iptr) - psys;
    v0 = Ymin(iptr);
    dv = Dy(iptr);
    if (pos*vsign > 0.0) {                 /* extract spectrum that always runs 'up' */
      for (iy=0; iy<ny; iy++) {
	vel[iy] = v0 + iy*dv;
	spec[iy] = MapValue(iptr,ix,iy);
      }
    } else {
      for (iy=0; iy<ny; iy++) {
	vel[ny-iy-1] = v0 + iy*dv;
	spec[ny-iy-1] = MapValue(iptr,ix,iy);
      }
    }
    if (pos*vsign < 0.0) {           /* fix up our mini WCS header */
      v0 = v0 + (ny-1)*dv;     
      dv = -dv;
    }
    nv = get_vels(ny,spec,v0,dv,it,clip,vels);    /* get the vel's */
    printf("%g ",qabs(pos,Qrotcur));
    vel_corr = (vels[0]-vsys)/sini;
    if (vel_corr > 0)
      vel_corr -= sigma;
    else
      vel_corr += sigma;
    printf(" %g", qabs(vel_corr,Qrotcur));

    for (j=1; j<nv; j++) {
      printf(" %g",qabs((vels[j]-vsys)/sini,Qrotcur));
    }
    printf("\n");
  }
  
  free(spec);
}


local void xyv_trace(imageptr iptr, int vsign, 
		    real eta, real ilc, real sigma, 
		    real xsys, real ysys, real vsys, 
		    real *clip, int mode)
{
  int  ix, iy, iz, j, nx, ny, nz, nv;
  real xpos, ypos, v0, dv, it, imax, vel_corr, vmin,vmax, dmin, dmax;
  real *spec, *vel, vels[MAXV];
  imageptr optr = NULL;
  stream ostr = stropen(getparam("out"),"w");
    
  nx = Nx(iptr);     /* assumed to be x-position for now */
  ny = Ny(iptr);     /* assumed to be y-position for now */
  nz = Nz(iptr);     /* assumed to be velocity for now */
  spec = (real *) allocate(nz*sizeof(real));
  vel  = (real *) allocate(nz*sizeof(real));

  create_cube(&optr,nx,ny,1);    /* trace 'velocity' field */
  dmin = dmax = 0.0;

  imax = MapMax(iptr);

  it = sqrt(sqr(eta*imax)+sqr(ilc));
  vmin = Zmin(iptr);
  vmax = vmin + (nz-1)*Dz(iptr);
  dprintf(0,"Map [%d x %d POS x %d VEL] I_t=%g vsys=%g (vrange: %g %g)\n",
	  nx,ny,nz,it,vsys,vmin,vmax);
  if (it<=0.0) error("I_t = %g too small",it);

  v0 = Zmin(iptr);
  dv = Dz(iptr);
  for (iz=0; iz<nz; iz++)
    vel[iz] = v0 + iz*dv;

  for (ix=0; ix<nx; ix++) {                 /* loop over all positions in P's */
    xpos = ix*Dx(iptr) + Xmin(iptr) - xsys;
    for (iy=0; iy<ny; iy++) {
      ypos = iy*Dy(iptr) + Ymin(iptr) - ysys;
      for (iz=0; iz<nz; iz++)
	spec[iz] = CubeValue(iptr,ix,iy,iz);
      nv = get_vels(nz,spec,v0,dv,it,clip,vels);    /* get the vel's */
      if (mode<2 || mode-1>nv) error("Illegal mode=%d found %d vel values",mode,nv);
      vel_corr = (vels[mode-2]-vsys);
      if (vel_corr > 0)
	vel_corr -= sigma;
      else
	vel_corr += sigma;
      CubeValue(optr,ix,iy,0) = vel_corr;
      dmin = MIN(dmin, vel_corr);
      dmax = MAX(dmax, vel_corr);
    }
  }
  free(spec);
  MapMin(optr) = dmin;
  MapMax(optr) = dmax;
  dprintf(0,"MinMax in map: %g %g\n",dmin,dmax);
  write_image(ostr,optr);
}


/* get_vels:
 *    this spectrum is always sorted such that we're looking for
 *    a peak at the upper end of the array, the sign of 'dv' will
 *    say where this is in velocity, but the array is always running
 *    "up"
 */

local int get_vels(int n, real *s, real v0, real dv, real smin, real *clip, real *vels)
{
  real sum1=0.0, sum0=0.0, sum2=0.0, v, speak = s[0], vpeak=0, vfit;
  real v1, v2, v3, sig;
  int i, nret = 0, ipeak = 0, itrace = -1;

  for (i=0,v=v0; i<n; i++,v+=dv) {     /* loop over spectrum and get basics: mom1 and peak */
    if (s[i] <= clip[0] || s[i] >= clip[1]) { /* but only outside clipping interval */
      sum0 += s[i];
      sum1 += s[i]*v;
      sum2 += s[i]*v*v;
      if (s[i] >= speak) {  /* find maximum */
	ipeak = i;
	vpeak = v;
	speak = s[i];
      }
    }
  }

  for (i=n-1; i>=0; i--) {      /* find the highest pixels above smin */
    if (itrace < 0 && s[i] > smin) {
      itrace = i;
      break;
    }
  }

  /* envelope trace (Sofue method) */
  vfit = v0 + dv*itrace;      /* pixel at which s[i] > smin */
  if (itrace > 0) {
    vfit += dv * (smin-s[itrace])/(s[itrace+1]-s[itrace]); /* small linear correction */
  } else if (itrace == 0 || itrace == n-1) {
    vfit += 0.0;    /* no correction can be made, vfit outside V range !! */
  } else
    vfit = 0.0;     /* bad itrace, should not happen */
  vels[nret++] = vfit;

  /* simple 1st order moment */
  vels[nret++] = (sum0 == 0.0 ? 0.0 : sum1/sum0);

  /* estimste sigma to aid in gaussian fitting */
  if (sum0 == 0.0)
    sig = 1.0;
  else {
    sig = sqrt(sum2/sum0-sum1*sum1/(sum0*sum0))/ABS(dv);
  }


  /* peak location */
  vels[nret++] = vpeak;

  /* peak location, from fitting paraboloid */
  if (ipeak == 0 || ipeak == n-1)
    vfit = vpeak;
  else {
    v1 = s[ipeak-1];
    v2 = s[ipeak];
    v3 = s[ipeak+1];
    if (v1+v3 == 2*v2) 
      vfit = 0.0;
    else {
      vfit = 0.5*(v1-v3)/(v1+v3-2*v2);
      vfit = dv*vfit + vpeak;
    }
  }
  vels[nret++] = vfit;

  /* simple gaussian fit */
  vfit = gfit1(n,s,clip,ipeak,sig,0);
  if (vfit > 0) {
    vels[nret++] = v0 + dv*vfit;
  } else
    vels[nret++] = 0.0;
  
  /* half gaussian fit */
  vfit = gfit1(n,s,clip,ipeak,sig,1);
  if (vfit > 0) {
    vels[nret++] = v0 + dv*vfit;
  } else
    vels[nret++] = 0.0;

  /* half gaussian fit minus corrections (MET) */
  vfit = gfit1(n,s,clip,ipeak,sig,2);
  if (vfit > 0) {
    vels[nret++] = v0 + dv*vfit;
  } else
    vels[nret++] = 0.0;

  /* AIPS++ window method */

  return nret;
}



local real qabs(real x, bool Q)
{
  if (Q && x < 0) return -x;
  return x;
}


/*
 * fit GAUSS1d:       y = a + b * exp( - (x-c)^2/(2*d^2) )
 *
 */


static real func_gauss1d(real *x, real *p, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  return p[0] + p[1] * exp(-arg);
}

static void derv_gauss1d(real *x, real *p, real *e, int np)
{
  real a,b,arg;
  a = p[2]-x[0];
  b = p[3];
  arg = a*a/(2*b*b);
  e[0] = 1.0;
  e[1] = exp(-arg);
  e[2] = -p[1]*e[1] * a   / (b*b);
  e[3] =  p[1]*e[1] * a*a / (b*b*b);
}

#define MAXP  1024

/*
 *  gfit fits a gauss, returns the central velocity
 *  in integer coordinates, which must be 0..n-1
 *  returns -1 if bad fit, no data etc.
 *  mode = 0: full gaussian
 *  mode = 1: half gaussian
 *  mode = 2: half gaussian + 0.5*FWHM (needed for MET)
 *
 */

local real gfit1(int n, real *s, real *clip, int ipeak, real sig, int mode)
{
  real x[MAXP], y[MAXP], fpar[4], epar[4], tol, lab;
  int i, npt, mpar[4], nrt, itmax;

  if (ipeak+1 >= n) return -1.0;

  /* accumulate the points; for half-gaussian, take one point before the peak
   * and then all after, as long as it's not inside the clip range 
   */

  if (mode == 0) {         /* gaussian */
    for (i=0, npt=0; i<n;  i++) {
      if (s[i] <= clip[0] || s[i] >= clip[1]) {
	if (npt==MAXP) error("Too many points, MAXP=%d",MAXP);
	x[npt] = i;
	y[npt] = s[i];
	npt++;
      }
    }
  } else {                /* half gaussian */
    for (i=ipeak-1, npt=0; i<n;  i++) {
      if (s[i] <= clip[0] || s[i] >= clip[1]) {
	if (npt==MAXP) error("Too many points, MAXP=%d",MAXP);
	x[npt] = i;
	y[npt] = s[i];
	npt++;
      }
    }
  } 

  if (npt < 3) return -1.0;    /* not enuf points for fit */


  while (sig > 0.5) {
    /* determine some reasonable initial estimates */
    fpar[0] = 0.0;      /* we'll fix this, there better be no continuum here */
    fpar[1] = s[ipeak]; /* peak value */
    fpar[2] = ipeak;    /* fitting is done in integer channel space */
    fpar[3] = sig;      /* this is the tricky one : we iterate on it */

    mpar[0] = 0;
    mpar[1] = mpar[2] = mpar[3] = 1;

    tol = 0.001;
    lab = 0.001;
    itmax = 50;
    
    nrt = nllsqfit(x,1,y,NULL,NULL,npt,fpar,epar,mpar,4,tol,itmax,lab, 
		 func_gauss1d, derv_gauss1d);
    if (nrt < 0) {
      sig = sig - 1.0;  /* lower sig (units are channel #s) for next iter */
      continue;
    }

    dprintf(1,"@ %d/%d %d %d %g: %g %g    %g %g  %d\n",ipeak,n,npt,mode,sig,
	    fpar[2],epar[2], fpar[3],epar[3],nrt);

    if (nrt > 0) {
      if (mode==2) return fpar[2] + 1.1774*fpar[3];   /* MET ;  1.1774 = sqrt(2ln2) */
      return fpar[2];
    } else {
#if 0
      printf("# %g %g %d %g\n",0.0,s[ipeak],ipeak,sig);
      for(i=0; i<npt; i++)
	printf("%g %g\n",x[i],y[i]);
#endif
      return -1;
    }
  }
  return -1;
}
