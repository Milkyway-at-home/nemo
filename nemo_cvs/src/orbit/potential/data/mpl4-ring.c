/*
 * mpl.c: procedures for intializing and calculating the forces and
 *             	potential of a combination of 3 potentials:
 *			Miyamoto Nagai Disk
 *			Plummer Sphere
 *			Logarithmic Halo
 *          Refs: BT pp. 43-44; Miyamoto and Nagai PASJ 27, 533 (1975)
 *	Miyamoto Potential Phi_m (R, z) = -GMdisk / sqrt (R^2+(a+sqrt(z^2+b^2))^2)
 * 	       Parameters: a, b (shape parameters), M (mass); G=1
 *             Names used: miya_ascal, miya_bscal, miya_mass
 *	Plummer Potential Phi_p (r) = -GMsph / (r + rc) = -GMsph/ (sqrt(R^2 + z^2) + rc)
 *		Parameters: rc, Msph; G=1
 *		Names used: plu_rc, plu_mass
 * 	Logarithmic Halo Phi_l (R, z) = vhalo^2 ln(R^2 + (z^2/q^2) + d^2)
 *		Parameters: vhalo, q, d
*		Names used: vhalo, q, d
 *  March 90 Stefano Casertano, University of Pittsburgh
 * 10-Nov-90 inserted omega as first parameter for standard Nemo  PJT
 *  6-oct-91 fixed bug - and made code accept both XYZ and XZY versions (pjt)
 *  7-mar-92 merged sun and 3b1 versions once more			 pjt
 *    oct-93 get_pattern
 *  12-mar-07 bwillett changed from miyamoto.c to mpl.c
 *  27-apr-07 bwillett modifying mpl3.c to change the acceleration fields
 *  1-may-07 bwillett created mpl4.c - took out gravitational constant
 *	all masses scaled as 1 M.U. = 222288.47 Ms 
 *	length unit: 1 kpc time unit: 1 Gyr
* 21-march-10 added ring in x-y plane
*
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>

#define PI 3.141592654

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 0.0;
local double miya_bscal = 1.0;
local double miya_mass = 1.0;
local double plu_rc = 1.0;
local double plu_mass = 1.0;
local double vhalo = 1.0;
local double q = 1.0;
local double d = 1.0;
local double ringmass = 1.0;
local double a = 1.0;

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

double f1(double mu);
double f2(double mu);
double fenergy(double x, double y, double z, double a);

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) miya_ascal = par[1];
    if (n>2) miya_bscal = par[2];
    if (n>3) miya_mass  = par[3];
    if (n>4) plu_rc = par[4];
    if (n>5) plu_mass = par[5];
    if (n>6) vhalo = par[6];
    if (n>7) q = par[7];
    if (n>8) d = par[8];
    if (n>9) ringmass = par[9];
    if (n>10) a = par[10];
    if (n>11) warning("mpl: only first 11 parameters recognized");

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
    qpar = hypot(pos[Z],miya_bscal);
    apar = miya_ascal+qpar;
    spar = (pos[X]*pos[X]) + (pos[Y]*pos[Y]) + ((miya_ascal+qpar)*(miya_ascal+qpar));
    ppar = sqrt ((pos[X]*pos[X])+(pos[Y]*pos[Y])+(pos[Z]*pos[Z])) + plu_rc;
    rpar = ppar - plu_rc;
    lpar = (rcyl*rcyl) + ((pos[Z]/q)*(pos[Z]/q)) + (d*d);
	
// 20070312 bwillett added plummer sphere potential and logarithmic potential
    *pot = (-(miya_mass)/sqrt(spar)) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(lpar));

// 20070312 bwillett rest left unchanged
// 20070427 bwillett changed the acceleration field

// This is only valid for 3 dimensions, and is in (x,y,z)
// Recall F_mu = -grad_mu U
// So a_mu = -grad_mu Phi
// I did these derivatives in Mathematica, and will try to keep it consistent with the conventions written above

double mplaccx, mplaccy, mplaccz;

	mplaccx = - ( ( (2.0*vhalo*vhalo*pos[X])/(lpar) ) + ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[X])/(pow(spar,1.5)) ) );
	mplaccy = - ( ( (2.0*vhalo*vhalo*pos[Y])/(lpar) ) + ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Y])/(pow(spar,1.5)) ) );
	mplaccz = - ( ( (2.0*vhalo*vhalo*pos[Z])/(q*q*lpar) ) + ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Z]*apar) / (qpar*pow(spar,1.5)) ) );

        double theta, sintheta, costheta, xi, mu, f1val, f2val, quant, psi;

        costheta = pos[Z] / sqrt(pos[X]*pos[X] + pos[Y]*pos[Y] + pos[Z]*pos[Z]);

        if(pos[Z] >= 0) {
                theta = acos(costheta);
                sintheta = cos(PI/2.0 - theta);
        } else {
                theta = -acos(costheta);
                sintheta = cos(PI/2.0 - theta);
        }

        psi = atan2(pos[Y],pos[X]);
        xi = sqrt(pos[X]*pos[X] + pos[Y]*pos[Y] + pos[Z]*pos[Z])/a;

        quant = 1.0 + xi*xi;

        mu = (2.0*xi*sintheta)/quant;

        // Calculate the integrals

        f1val = f1(mu);
        f2val = f2(mu);

        // The energy (FIX THIS, it doesn't matter unless you want to separate out parts of the tail)

    *pot = *pot + ((-G*ringmass)/(2*PI)) * fenergy(pos[X],pos[Y],pos[Z],a);

        // The accelerations

        double accmag;
        accmag = (-G*ringmass / (PI*a*a)) * (xi*sintheta*pow(quant, -1.5)*f1val - pow(quant, -1.5)*f2val);

        acc[X] = mplaccx + accmag*cos(psi);
        acc[Y] = mplaccy + accmag*sin(psi);
        acc[Z] = mplaccz + (-G*ringmass / (PI*a*a)) * (xi*costheta*pow(quant, -1.5)*f1val);

}

// The integrals
double f1(double mu) {
	double dalpha = PI / 10.0;
	double alpha = 0.0;
	double quant = 0.0;
	double total = 0.0;

	for(alpha=0.0;alpha<=PI;alpha+=dalpha) {
		quant = 1.0 - mu*cos(alpha);
		total += pow(quant, -1.5)*dalpha;
	}

	return total;
}

double f2(double mu) {
	double dalpha = PI / 10.0;
	double alpha = 0.0;
	double quant = 0.0;
	double total = 0.0;

	for(alpha=0.0;alpha<=PI;alpha+=dalpha) {
		quant = 1.0 - mu*cos(alpha);
		total += cos(alpha)*pow(quant, -1.5)*dalpha;
	}

	return total;
}

double fenergy(double x, double y, double z, double a) {
	double dalpha = 2*PI / 20.0;
	double alpha = 0.0;
	double quant = 0.0;
	double total = 0.0;

	for(alpha=0.0;alpha<=2*PI;alpha+=dalpha) {
		quant = alpha / pow((x-a*cos(alpha))*(x-a*cos(alpha)) + (y-a*sin(alpha))*(y-a*sin(alpha)) + z*z, 0.5);
		total += quant*dalpha;
	}

	return total;
}

