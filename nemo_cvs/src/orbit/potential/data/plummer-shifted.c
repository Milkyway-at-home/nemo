/*
 * plummer.c:  (spherical) plummer potential
 *
 *	sep-2001	provide both a _double and _float version for the new potproc interface
 *
 */

// modified by willeb to shift the center to (xc,yc,zc)
// 20100420 willeb fixed shift problem
/*CTEX
 *	{\bf potname=plummer
 *       potpars={\it $\Omega,M,R$}}
 *
 *  Plummer potential (BT, pp.42, eq. 2.47, see also MNRAS 71, 460 (1911))
 *
 * $$
 *    \Phi = -  {  M  \over
 *                    {   {(r_c^2 + (r-r0)^2)}^{1/2} }  }
 * $$
 */                     
  
#include <stdinc.h>
#include <math.h>
#include <vectmath.h>	/* define DIMensionality */
 
local double omega = 0.0;
local double plummer_mass = 1.0;
local double plummer_radius = 1.0;	/* radius=0 gives point mass */
local double xc = 0.0;
local double yc = 0.0;
local double zc = 0.0;
local double r2;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) plummer_mass = par[1];
    if (n>2) plummer_radius = par[2];
    if (n>3) xc = par[3];
    if (n>4) yc = par[4];
    if (n>5) zc = par[5];
    if (n>6) warning("plummer: npar=%d only 6 parameters accepted",n);

#if defined(TWODIM)
    dprintf (1,"INIPOTENTIAL Plummer: [2d opt]\n");
#else
    dprintf (1,"INIPOTENTIAL Plummer: [3d version]\n");
#endif
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  mass, radius = %f %f \n",plummer_mass,plummer_radius);
    dprintf (1,"  center point (xc, yc, zc) = (%f, %f, %f)\n", xc, yc, zc);
	
    r2 = sqrt(plummer_radius);
    par[0] = omega;
}


void potential_double (int *ndim,double *pos,double *acc,double *pot,double *time)
{
    register double tmp;

#if defined(TWODIM)
    /* TWODIM: Saves two MUL's - hardly making any difference... */
    tmp = 1.0/(r2 + (pos[0]-xc)*(pos[0]-xc) + (pos[1]-yc)*(pos[1]-yc));
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*(pos[0]-xc);
    acc[1] = tmp*(pos[1]-yc);
    acc[2] = 0.0;
#else
    tmp = 1.0/(r2 + (pos[0]-xc)*(pos[0]-xc) + (pos[1]-yc)*(pos[1]-yc) + (pos[2]-zc)*(pos[2]-zc));
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*(pos[0]-xc);
    acc[1] = tmp*(pos[1]-yc);
    acc[2] = tmp*(pos[2]-zc);

#endif
}

void potential_float (int *ndim,float *pos,float *acc,float *pot,float *time)
{
    register float tmp;

#if defined(TWODIM)
    /* TWODIM: Saves two MUL's - hardly making any difference... */
    tmp = 1.0/(r2 + (pos[0]-xc)*(pos[0]-xc) + (pos[1]-yc)*(pos[1]-yc));
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*(pos[0]-xc);
    acc[1] = tmp*(pos[1]-yc);
    acc[2] = 0.0;
#else
    tmp = 1.0/(r2 + (pos[0]-xc)*(pos[0]-xc) + (pos[1]-yc)*(pos[1]-yc) + (pos[2]-zc)*(pos[2]-zc));
    *pot = -sqrt(tmp);
    tmp *= (*pot) * plummer_mass;
    *pot *= plummer_mass;
    acc[0] = tmp*(pos[0]-xc);
    acc[1] = tmp*(pos[1]-yc);
    acc[2] = tmp*(pos[2]-zc);

#endif
}

