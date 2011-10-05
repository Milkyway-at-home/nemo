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
*
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

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
local double ring_mass = 1.0;
local double ring_radius = 1.0;

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
    if (n>7) q = par[7];
    if (n>8) d = par[8];
    if (n>9) ring_mass = par[9];
    if (n>10) ring_radius = par[10];
    if (n>11) warning("mpl: only first 11 parameters recognized");

}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double apar, qpar, spar, ppar, lpar, rpar, rcyl;
    double i1,i2,i3;

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
    //*pot = (-(miya_mass)/spar) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(lpar));

// 20070312 bwillett rest left unchanged
// 20070427 bwillett changed the acceleration field

// This is only valid for 3 dimensions, and is in (x,y,z)
// Recall F_mu = -grad_mu U
// So a_mu = -grad_mu Phi
// I did these derivatives in Mathematica, and will try to keep it consistent with the conventions written above

	acc[X] = - ( ( (2.0*vhalo*vhalo*pos[X])/(lpar) ) + ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[X])/(pow(spar,1.5)) ) );
	acc[Y] = - ( ( (2.0*vhalo*vhalo*pos[Y])/(lpar) ) + ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Y])/(pow(spar,1.5)) ) );
	acc[Z] = - ( ( (2.0*vhalo*vhalo*pos[Z])/(lpar) ) + ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Z]*apar) / (qpar*pow(spar,1.5)) ) );



// willeb 20080911 - do the integrals required for the ring potential
// X-integral

        double x,y,z;
        double R = ring_radius;
        double theta1 = 0, theta2 = 0, theta3 = 0, dtheta = 0;
        double N = 100;
        double result1 = 0, result2 = 0, result3 = 0;

        x = pos[X];
        y = pos[Y];
        z = pos[Z];
        dtheta = (2*PI)/N;

// X-integral
        for(i1=0;i1<=N;i1++) {
                theta1 = ((2*PI)/N) * i1;

                result1 += ((x-R*cos(theta1)) / pow(z*z + pow(x-R*cos(theta1), 2.0) + pow(y-R*sin(theta1), 2.0), 1.5))*dtheta;
        }
	//printf("%f %f %f: %f\n", x,y,z,result1);

	// Get the units right
	result1 = -(ring_mass/(2*PI)) * result1;

	acc[X] += result1;

//  Y-integral
        for(i2=0;i2<=N;i2++) {
                theta2 = ((2*PI)/N) * i2;

                result2 += ((y-R*sin(theta2)) / pow(z*z + pow(x-R*cos(theta2), 2.0) + pow(y-R*sin(theta2), 2.0), 1.5) )*dtheta;
        }

	//printf("Y: %f\n", result);

	// Get the units right
	result2 = -(ring_mass/(2*PI)) * result2;

	acc[Y] += result2;

// Z-integral
        for(i3=0;i3<=N;i3++) {
                theta3 = ((2*PI)/N) * i3;

                result3 += ((z) / pow(z*z + pow(x-R*cos(theta3), 2.0) + pow(y-R*sin(theta3), 2.0), 1.5))*dtheta;
        }
	//printf("Z: %f\n", result);

	// Get the units right
	result3 = -(ring_mass/(2*PI)) * result3;

	acc[Z] += result3;

}
