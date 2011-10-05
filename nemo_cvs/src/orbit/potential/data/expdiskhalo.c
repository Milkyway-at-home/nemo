/* expdiskhalo.c

Uses the bulge and halo model from mpl4.c, but includes an exponential disk instead of a Miyamoto-Nagai disk

parameters: 0,diska,diskmass,bulger,bulgemass,vhalo,haloq,halod
*/

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double a = 1.0;
local double mass = 1.0;
local double plu_rc = 1.0;
local double plu_mass = 1.0;
local double vhalo = 1.0;
local double q = 1.0;
local double d = 1.0;

/* 20070427 bwillett commented this out
#if !defined(Y) && !defined(Z)
                                          default: XZY setup (re-oriented)
#define X 0
#define Y 2
#define Z 1
#else
                      the normal axisymmetric case, as defined in e.g. B&T */
#define X 0
#define Y 1
#define Z 2
// #endif

extern double bessi0(double), bessk0(double), bessi1(double), bessk1(double);

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) a = par[1];
    if (n>2) mass = par[2];
    if (n>3) plu_rc = par[3];
    if (n>4) plu_mass = par[4];
    if (n>5) vhalo = par[5];
    if (n>6) q = par[6];
    if (n>7) d = par[7];
    if (n>8) warning("mpl: only first 8 parameters recognized");

}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double apar, qpar, spar, ppar, lpar, rpar, rcyl;
    int i;

// 20070312 bwillett added ppar - the plummer denominator: r + rc = sqrt(x^2+y^2+z^2) + rc
// 20070312 bwillett added lpar - the logarithmic argument: R^2 + (z/q)^2 + d^2
// 20070427 bwillett added rpar - the spherical radius: r + rc - rc = ppar - plu_rc
// 20070501 bwillett added apar - a + qpar
// 20070507 bwillett took out pow statements
// 20070507 bwillett used hypot from math.h

    rcyl = hypot(pos[X],pos[Y]);
    ppar = sqrt ((pos[X]*pos[X])+(pos[Y]*pos[Y])+(pos[Z]*pos[Z])) + plu_rc;
    rpar = ppar - plu_rc;
    lpar = (rcyl*rcyl) + ((pos[Z]/q)*(pos[Z]/q)) + (d*d);

// This is only valid for 3 dimensions, and is in (x,y,z)
// Recall F_mu = -grad_mu U
// So a_mu = -grad_mu Phi
// I did these derivatives in Mathematica, and will try to keep it consistent with the conventions written above

	acc[X] = - ( ( (2.0*vhalo*vhalo*pos[X])/(lpar) ) + ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) );
	acc[Y] = - ( ( (2.0*vhalo*vhalo*pos[Y])/(lpar) ) + ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) );
	acc[Z] = - ( ( (2.0*vhalo*vhalo*pos[Z])/(lpar) ) + ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) );

// Copied from expdisk.c
    double r2, r, arg, i0, k0, i1, k1, f;

	double alpha;
	alpha = 1.0/a;
	r = rpar;
	r2 = r*r;
    arg = 0.5*alpha*r;

	//printf("%f %f %f %f %f\n", a, mass, r, r2, x);
        i0=bessi0(arg);
        k0=bessk0(arg);
        i1=bessi1(arg);
        k1=bessk1(arg);

	// 20080928 - willeb added exponential disk to acceleration field 
        *pot = -mass*arg*(i0*k1-i1*k0);
        f = -0.5*alpha*alpha*alpha*mass*(i0*k0-i1*k1);
        acc[X] += f*pos[X];
        acc[Y] += f*pos[Y];
        acc[Z] += f*pos[Z];

// 20080928 - willeb added bulge and halo to potential 
    *pot += (-(plu_mass)/ppar);
*pot += (vhalo*vhalo*log(lpar));
}
