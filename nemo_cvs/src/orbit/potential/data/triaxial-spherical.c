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
 *  1-feb-09 bwillett created mpl4-triaxial.c, implemented Law09's triaxial halo
* 21-may-10 bwillett created triaxial-spherical.c, creating inner triaxial halo and outer spherical halo in the same potential
*
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>
#define PI2 3.141592654

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 0.0;
local double miya_bscal = 1.0;
local double miya_mass = 1.0;
local double plu_rc = 1.0;
local double plu_mass = 1.0;
local double vhalo = 1.0;
local double phi = 0.0; 
local double C1=1.0, C2=1.0, C3=1.0;
local double q1 = 1.0, q2 = 1.0, qz = 1.0;
local double d = 1.0;

local double vhalo2 = 1.0;
local double qhalo2 = 1.0;
local double halod2 = 1.0;

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
    if (n>7) q1 = par[7];
    if (n>8) q2 = par[8];
    if (n>9) qz = par[9];
    if (n>10) phi = par[10];
    if (n>11) d = par[11];
    if (n>12) vhalo2 = par[12];
    if (n>13) qhalo2 = par[13];
    if (n>14) halod2 = par[14];
    if (n>15) warning("mpl: only first 12 parameters recognized");

}

void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double apar, qpar, spar, ppar, rpar, rcyl,lpar;
    double arg;
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
    lpar = (rcyl*rcyl) + ((pos[Z]/qhalo2)*(pos[Z]/qhalo2)) + (halod2*halod2);
	
//    *pot = 1;

// 20070312 bwillett rest left unchanged
// 20070427 bwillett changed the acceleration field

// This is only valid for 3 dimensions, and is in (x,y,z)
// Recall F_mu = -grad_mu U
// So a_mu = -grad_mu Phi
// Added the Law09 triaxial halo 

// Define C1, C2, C3
// phi is in degrees

C1 = ((1.0/(q1*q1)) * cos(phi*PI2/180.0) * cos(phi*PI2/180.0)) + ((1.0/(q2*q2)) * sin(phi*PI2/180.0) * sin(phi*PI2/180.0));
C2 = ((1.0/(q2*q2)) * cos(phi*PI2/180.0) * cos(phi*PI2/180.0)) + ((1.0/(q1*q1)) * sin(phi*PI2/180.0) * sin(phi*PI2/180.0));
C3 = 2*sin(phi*PI2/180.0)*cos(phi*PI2/180.0)*( (1.0/(q1*q1)) - (1.0/(q2*q2)));

arg = C1*pos[X]*pos[X] + C2*pos[Y]*pos[Y] + C3*pos[X]*pos[Y] + ((pos[Z]*pos[Z])/(qz*qz)) + d*d;

    *pot = (-(miya_mass)/sqrt(spar)) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(arg)) + (vhalo2*vhalo2*log(lpar));

	acc[X] = - ( ( (vhalo*vhalo*(2.0*C1*pos[X] + C3*pos[Y]))/(arg) ) + ( (2.0*vhalo2*vhalo2*pos[X])/(lpar) ) + ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[X])/(pow(spar,1.5)) ) );
	acc[Y] = - ( ( (vhalo*vhalo*(2.0*C2*pos[Y] + C3*pos[X]))/(arg) ) + ( (2.0*vhalo2*vhalo2*pos[Y])/(lpar) ) + ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Y])/(pow(spar,1.5)) ) );
	acc[Z] = - ( ( (2.0*vhalo*vhalo*pos[Z])/(qz*qz*arg) ) + ( (2.0*vhalo2*vhalo2*pos[Z])/(qhalo2*qhalo2*lpar) ) + ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Z]*apar) / (qpar*pow(spar,1.5)) ) );

}