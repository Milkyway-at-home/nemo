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
*
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>
#define PI2 3.141592654

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 1.0;
local double miya_bscal = 1.0;
local double miya_mass = 1.0;
local double plu_rc = 1.0;
local double plu_mass = 1.0;
local double vhalo = 1.0;
local double eath = 90.0;
local double eaph = 90.0; 
local double eaps = 90.0;
local double q1 = 1.0, q2 = 1.0, qz = 1.0;
local double d = 1.0;
local double rot[3][3];
local double rotinv[3][3];  

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
    if (n>8) qz = par[8];
    if (n>9) eath = par[9];
    if (n>10) eaph = par[10];	
    if (n>11) eaps = par[11];	
    if (n>12) d = par[12];	
    if (n>13) warning("mpl: only first 13 parameters recognized");
// convert angles to radians
eath *= 0.017453278;
eaph *= 0.017453278;
eaps *= 0.017453278;
}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    local double apar, qpar, spar, ppar, rpar, rcyl;
    local double arg;
    local int i;

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

// Define the elements of the rotation matrix (using the xyz pitch-roll-yaw convention e.g. Wolfram)

rot[0][0] = cos(eath)*cos(eaph);
rot[0][1] = cos(eath)*sin(eaph);
rot[0][2] = -sin(eath);

rot[1][0] = sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph);
rot[1][1] = sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph);
rot[1][2] = cos(eath)*sin(eaps);

rot[2][0] = cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph);
rot[2][1] = cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph);
rot[2][2] = cos(eath)*cos(eaps);

rotinv[0][0] = ((sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[0][1] = (-(cos(eath)*sin(eaph)*cos(eath)*cos(eaps)+sin(eath)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[0][2] = (cos(eath)*sin(eaph)*cos(eath)*sin(eaps)+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[1][0] = (cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))-(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[1][1] = (cos(eath)*cos(eaph)*cos(eath)*cos(eaps)+sin(eath)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[1][2] = (-(sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))+cos(eath)*cos(eaph)*cos(eath)*sin(eaps)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[2][0] = ((sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[2][1] = (cos(eath)*sin(eaph)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))-cos(eath)*cos(eaph)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

rotinv[2][2] = (cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph)))/(cos(eath)*cos(eaph)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*cos(eath)*cos(eaps)-cos(eath)*cos(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-sin(eath)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*(cos(eaps)*sin(eath)*sin(eaph)-sin(eaps)*cos(eaph))-cos(eath)*sin(eaph)*(sin(eaps)*sin(eath)*cos(eaph)-cos(eaps)*sin(eaph))*cos(eath)*cos(eaps)+cos(eath)*sin(eaph)*cos(eath)*sin(eaps)*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph))+sin(eath)*(sin(eaps)*sin(eath)*sin(eaph)+cos(eaps)*cos(eaph))*(cos(eaps)*sin(eath)*cos(eaph)+sin(eaps)*sin(eaph)));

local double xp;
local double yp;
local double zp;

q2 = 1.0;

// Construct the potential
//xp = (rotinv[0][0]*(pos[X]/q1) + rotinv[0][1]*(pos[Y]/q2) + rotinv[0][2]*(pos[Z]/qz));
//yp = (rotinv[1][0]*(pos[X]/q1) + rotinv[1][1]*(pos[Y]/q2) + rotinv[1][2]*(pos[Z]/qz));
//zp = (rotinv[2][0]*(pos[X]/q1) + rotinv[2][1]*(pos[Y]/q2) + rotinv[2][2]*(pos[Z]/qz));
xp = (rot[0][0]*(pos[X]/q1) + rot[0][1]*(pos[Y]/q1) + rot[0][2]*(pos[Z]/q1));
yp = (rot[1][0]*(pos[X]/q2) + rot[1][1]*(pos[Y]/q2) + rot[1][2]*(pos[Z]/q2));
zp = (rot[2][0]*(pos[X]/qz) + rot[2][1]*(pos[Y]/qz) + rot[2][2]*(pos[Z]/qz));
arg = xp*xp + yp*yp + zp*zp + d*d;

*pot = (-(miya_mass)/sqrt(spar)) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(arg));

// Calculate the accelerations

// Do the disk and bulge first (from mpl4)
acc[X] = - ( ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[X])/(pow(spar,1.5)) ) );
acc[Y] = - ( ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Y])/(pow(spar,1.5)) ) );
acc[Z] = - ( ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Z]*apar) / (qpar*pow(spar,1.5)) ) ); 

// Add the contributions from the triaxial halo (the minus sign is because a = - grad Phi) 

//acc[X] += -(2.0*vhalo*vhalo/(arg*q1))*(xp*rotinv[0][0] + yp*rotinv[1][0] + zp*rotinv[2][0]);
//acc[Y] += -(2.0*vhalo*vhalo/(arg*q2))*(xp*rotinv[0][1] + yp*rotinv[1][1] + zp*rotinv[2][1]);
//acc[Z] += -(2.0*vhalo*vhalo/(arg*qz))*(xp*rotinv[0][2] + yp*rotinv[1][2] + zp*rotinv[2][2]);

acc[X] += -(2.0*vhalo*vhalo/(arg))*(xp*(rot[0][0]/q1) + yp*(rot[1][0]/q2) + zp*(rot[2][0]/qz));
acc[Y] += -(2.0*vhalo*vhalo/(arg))*(xp*(rot[0][1]/q1) + yp*(rot[1][1]/q2) + zp*(rot[2][1]/qz));
acc[Z] += -(2.0*vhalo*vhalo/(arg))*(xp*(rot[0][2]/q1) + yp*(rot[1][2]/q2) + zp*(rot[2][2]/qz));
}
