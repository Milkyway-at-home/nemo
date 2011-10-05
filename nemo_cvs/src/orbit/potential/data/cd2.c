/*
 * cd.c: procedures for intializing and calculating the forces and
 *             	of a combination of 2 potentials:
 *			Triaxial luminous disk	
 *			Triaxial Dark Matter Halo
 *          Refs: Capuzzo-Dolcetta 2007
 * 	       Parameters: a, b, c, M Luminous 
 * 	       Parameters: a, b, c Dark Matter 
 *  13-august-08 bwillett created cd.c - took out gravitational constant
 *	all masses scaled as 1 M.U. = 222288.47 Ms 
 *	length unit: 1 kpc time unit: 1 Gyr
*
 */

#include <stdinc.h>
#include <potential_float.h>
#include <math.h>
#include <omp.h>

local double G = 1;
local double omega = 0.0;		/* pattern speed */
local double al = 1.0;
local double bl = 1.0;
local double cl = 1.0;
local double M = 1.0;
local double adm = 1.0;
local double bdm = 1.0;
local double cdm = 1.0;

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
    if (n>1) al = par[1];
    if (n>2) bl = par[2];
    if (n>3) cl  = par[3];
    if (n>4) M = par[4];
    if (n>5) adm = par[5];
    if (n>6) bdm = par[6];
    if (n>7) cdm = par[7];
    if (n>8) warning("mpl: only first 8 parameters recognized");

}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{
// I get the expressions for the accelerations from Capuzzo-Dolcetta, 2007
// We do the integrals numerically using the trapezoid rule
double upperlimit1,n1,h1,total1=0.0,oldtotal1=-1.0,i1,prefactor1,func1;
double upperlimit2,n2,h2,total2=0.0,oldtotal2=-1.0,i2,prefactor2,func2;
double upperlimit3,n3,h3,total3=0.0,oldtotal3=-1.0,i3,prefactor3,func3;
double m, mp;
double accxl,accyl,acczl,accxdm,accydm,acczdm;

// Define the coordinates
m = (pow(pos[X],2.0)/pow(al,2.0)) +  (pow(pos[Y],2.0)/pow(bl,2.0)) + (pow(pos[Z],2.0)/pow(cl,2.0));
mp = (pow(pos[X],2.0)/pow(adm,2.0)) +  (pow(pos[Y],2.0)/pow(bdm,2.0)) + (pow(pos[Z],2.0)/pow(cdm,2.0));

#pragma omp parallel sections num_threads(3) 
{

#pragma omp section 
{
// The x-luminous integral
upperlimit1 = 5000.0;
n1 = 1000000.0;

h1 = (upperlimit1)/n1;

prefactor1 = h1/2.0;

for(i1=0;i1<=upperlimit1;i1+=h1) {
	if(i1==0 || i1==upperlimit1) {
		func1 = 1.0 / ( pow((i1+pow(al,2.0))*(i1+pow(bl,2.0))*(i1+pow(cl,2.0)),0.5) * (pow(al,2.0) + i1) );
	} else {
		func1 = 2.0 / ( pow((i1+pow(al,2.0))*(i1+pow(bl,2.0))*(i1+pow(cl,2.0)),0.5) * (pow(al,2.0) + i1) );
	}

	total1 += prefactor1 * func1;

	if(total1 == oldtotal1) {
		break;
	} else {
		oldtotal1 = total1;
	}
}

// The x-luminous acceleration
accxl = -((G * M * pos[X])/( pow(1.0+m,3.0)*m)) * total1;
}

#pragma omp section 
{ 
// The y-luminous integral
upperlimit2 = 5000.0;
n2 = 1000000.0;

h2 = (upperlimit2 - 0.0)/n2;

prefactor2 = h2/2.0;

for(i2=0;i2<=upperlimit2;i2+=h2) {
	if(i2==0 || i2==upperlimit2) {
		func2 = 1.0 / ( pow((i2+pow(al,2.0))*(i2+pow(bl,2.0))*(i2+pow(cl,2.0)),0.5) * (pow(bl,2.0) + i2) );
	} else {
		func2 = 2.0 / ( pow((i2+pow(al,2.0))*(i2+pow(bl,2.0))*(i2+pow(cl,2.0)),0.5) * (pow(bl,2.0) + i2) );
	}

	total2 += prefactor2 * func2;

	if(total2 == oldtotal2) {
		break;
	} else {
		oldtotal2 = total2;
	}
}

// The y-luminous acceleration
accyl = -((G * M * pos[Y])/( pow(1.0+m,3.0)*m)) * total2;
}

#pragma omp section 
{
// The y-luminous integral
upperlimit3 = 5000.0;
n3 = 1000000.0;

h3 = (upperlimit3 - 0.0)/n3;

prefactor3 = h3/2.0;

for(i3=0;i3<=upperlimit3;i3+=h3) {
	if(i3==0 || i3==upperlimit3) {
		func3 = 1.0 / ( pow((i3+pow(al,2.0))*(i3+pow(bl,2.0))*(i3+pow(cl,2.0)),0.5) * (pow(cl,2.0) + i3) );
	} else {
		func3 = 2.0 / ( pow((i3+pow(al,2.0))*(i3+pow(bl,2.0))*(i3+pow(cl,2.0)),0.5) * (pow(cl,2.0) + i3) );
	}

	total3 += prefactor3 * func3;

	if(total3 == oldtotal3) {
		break;
	} else {
		oldtotal3 = total3;
	}
}

// The z-luminous acceleration
acczl = -((G * M * pos[Z])/( pow(1.0+m,3.0)*m)) * total3;
}

}
	acc[X] = accxl; 
	acc[Y] = accyl; 
	acc[Z] = acczl; 

}
