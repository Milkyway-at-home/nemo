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

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double miya_ascal = 0.0;
local double miya_bscal = 1.0;
local double miya_mass = 1.0;
local double plu_rc = 1.0;
local double plu_mass = 1.0;
local double vhalo0 = 1.0;
local double vhaloa1 =1.0; 
local double vhaloa2 =1.0; 
local double vhaloa3 =1.0; 
local double vhaloa4 =1.0; 
local double vhaloa5 =1.0; 
local double vhaloa6 =1.0; 
local double vhaloa7 =1.0; 
local double vhaloa8 =1.0; 
local double vhalob1 =1.0; 
local double vhalob2 =1.0; 
local double vhalob3 =1.0; 
local double vhalob4 =1.0; 
local double vhalob5 =1.0; 
local double vhalob6 =1.0; 
local double vhalob7 =1.0; 
local double vhalob8 =1.0; 
local double q0 = 1.0;
local double qa1 =1.0; 
local double qa2 =1.0; 
local double qa3 =1.0; 
local double qa4 =1.0; 
local double qa5 =1.0; 
local double qa6 =1.0; 
local double qa7 =1.0; 
local double qa8 =1.0; 
local double qb1 =1.0; 
local double qb2 =1.0; 
local double qb3 =1.0; 
local double qb4 =1.0; 
local double qb5 =1.0; 
local double qb6 =1.0; 
local double qb7 =1.0; 
local double qb8 =1.0; 
local double d = 1.0;
local double Lvhalo = 1.0;
local double Lq = 1.0;

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
    int num,i;

    num = *npar;
    if (num>0) omega = par[0];
    if (num>1) miya_ascal = par[1];
    if (num>2) miya_bscal = par[2];
    if (num>3) miya_mass  = par[3];
    if (num>4) plu_rc = par[4];
    if (num>5) plu_mass = par[5];
    if (num>6) vhalo0 = par[6];
    if (num>7) vhaloa1 = par[7];
    if (num>8) vhaloa2 = par[8];
    if (num>9) vhaloa3 = par[9];
    if (num>10) vhaloa4 = par[10];
    if (num>11) vhaloa5 = par[11];
    if (num>12) vhaloa6 = par[12];
    if (num>13) vhaloa7 = par[13];
    if (num>14) vhaloa8 = par[14];
    if (num>15) vhalob1 = par[15];
    if (num>16) vhalob2 = par[16];
    if (num>17) vhalob3 = par[17];
    if (num>18) vhalob4 = par[18];
    if (num>19) vhalob5 = par[19];
    if (num>20) vhalob6 = par[20];
    if (num>21) vhalob7 = par[21];
    if (num>22) vhalob8 = par[22];
    if (num>23) q0 = par[23];
    if (num>24) qa1 = par[24];
    if (num>25) qa2 = par[25];
    if (num>26) qa3 = par[26];
    if (num>27) qa4 = par[27];
    if (num>28) qa5 = par[28];
    if (num>29) qa6 = par[29];
    if (num>30) qa7 = par[30];
    if (num>31) qa8 = par[31];
    if (num>32) qb1 = par[32];
    if (num>33) qb2 = par[33];
    if (num>34) qb3 = par[34];
    if (num>35) qb4 = par[35];
    if (num>36) qb5 = par[36];
    if (num>37) qb6 = par[37];
    if (num>38) qb7 = par[38];
    if (num>39) qb8 = par[39];
    if (num>40) d = par[40];
    if (num>41) Lvhalo = par[41];
    if (num>42) Lq = par[42];
    if (num>43) warning("mpl: only first 43 parameters recognized");
}

void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
    double apar, qpar, spar, ppar, lpar, rpar, rcyl;
    int i;
    double vhalo, q;
    double t = time[0];
    double pi = 3.141592654;

    vhalo = (vhalo0/2.0) + vhaloa1*cos(1.0*pi*t/Lvhalo) + vhaloa2*cos(2.0*pi*t/Lvhalo) + vhaloa3*cos(3.0*pi*t/Lvhalo) + vhaloa4*cos(4.0*pi*t/Lvhalo) + vhaloa5*cos(5.0*pi*t/Lvhalo) + vhaloa6*cos(6.0*pi*t/Lvhalo) + vhaloa7*cos(7.0*pi*t/Lvhalo) + vhaloa8*cos(8.0*pi*t/Lvhalo) + vhalob1*sin(1.0*pi*t/Lvhalo) + vhalob2*sin(2.0*pi*t/Lvhalo) + vhalob3*sin(3.0*pi*t/Lvhalo) + vhalob4*sin(4.0*pi*t/Lvhalo) + vhalob5*sin(5.0*pi*t/Lvhalo) + vhalob6*sin(6.0*pi*t/Lvhalo) + vhalob7*sin(7.0*pi*t/Lvhalo) + vhalob8*sin(8.0*pi*t/Lvhalo);

    q = (q0/2.0) + qa1*cos(1.0*pi*t/Lq) + qa2*cos(2.0*pi*t/Lq) + qa3*cos(3.0*pi*t/Lq) + qa4*cos(4.0*pi*t/Lq) + qa5*cos(5.0*pi*t/Lq) + qa6*cos(6.0*pi*t/Lq) + qa7*cos(7.0*pi*t/Lq) + qa8*cos(8.0*pi*t/Lq) + qb1*sin(1.0*pi*t/Lq) + qb2*sin(2.0*pi*t/Lq) + qb3*sin(3.0*pi*t/Lq) + qb4*sin(4.0*pi*t/Lq) + qb5*sin(5.0*pi*t/Lq) + qb6*sin(6.0*pi*t/Lq) + qb7*sin(7.0*pi*t/Lq) + qb8*sin(8.0*pi*t/Lq);

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

	acc[X] = - ( ( (2.0*vhalo*vhalo*pos[X])/(lpar) ) + ( (plu_mass*pos[X])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[X])/(pow(spar,1.5)) ) );
	acc[Y] = - ( ( (2.0*vhalo*vhalo*pos[Y])/(lpar) ) + ( (plu_mass*pos[Y])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Y])/(pow(spar,1.5)) ) );
	acc[Z] = - ( ( (2.0*vhalo*vhalo*pos[Z])/(q*q*lpar) ) + ( (plu_mass*pos[Z])/(rpar*ppar*ppar) ) + ( (miya_mass*pos[Z]*apar) / (qpar*pow(spar,1.5)) ) );

}
