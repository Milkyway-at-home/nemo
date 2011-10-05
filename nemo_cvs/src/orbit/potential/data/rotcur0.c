/*
 * rotcur0.c: potential & forces as defined by a rotation curve
 *            that is linear to (r0,v0) and flat thereafter
 *
 *	29-dec-01	derived from rotcur
 *      19-sep-04       float/double
 */

/*CTEX
 *  {\bf potname=rotcur0
 *       potpars={\it $\Omega,r_0,v_0$}}
 *	 
 * The forces returned are the axisymmetric forces as defined by
 * a linear-flat rotation curve as defined by the turnover point $r_0,v_0$.
 * The potential is not computed, instead the interpolated rotation
 * curve is returned in as the potential value.
 */
 

#include <stdinc.h>
#include <spline.h>
#include <table.h>
#include <potential_float.h>

local double omega = 0.0;
local double r0  = 1.0;
local double v0 = 1.0;

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0];
    if (n>1) r0    = par[1];
    if (n>2) v0    = par[2];
    if (n>3) warning("Rotcur potential: only 3 parameters usable");
    
    dprintf (1,"INIPOTENTIAL Rotcur0 potential %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  R0 = %g  V0 = %g\n", r0, v0);

    par[0] = omega;     /* return pattern speed again */
}
    
void potential_double (int *ndim, double *pos,double *acc,double *pot,double *time)
{
    real r, r2, v, f;
    int    i;

    for (i=0, r2=0.0; i<2; i++)
        r2 += sqr(pos[i]);
    r=sqrt(r2);

    if (r < r0)
        v = (r/r0)*v0;
    else
        v = v0;

    if (r > 0)
	f = sqr(v/r);
    else
    	f = 0;
    dprintf(2,"r=%g v=%g f=%g\n",r,v,f);

    *pot = 0.0;             /* no potentials... for now */
    acc[0] = -f*pos[0]; 
    acc[1] = -f*pos[1]; 
    acc[2] = 0.0;
}
