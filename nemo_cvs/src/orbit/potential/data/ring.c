#include <stdinc.h>
#include <potential_float.h>
#include <math.h>

local double omega = 0; // Pattern speed
local double G = 1;
local double ringmass = 1E5;
local double a = 15.0;

#define X 0
#define Y 1
#define Z 2
#define PI 3.141592654

double f1(double mu);
double f2(double mu);
double fenergy(double x, double y, double z, double a);

void inipotential (int *npar, double *par, string name)
{
    int n;

    n = *npar;
    if (n>0) omega = par[0]; // Pattern speed 
    if (n>1) ringmass = par[1];
    if (n>2) a = par[2];
    if (n>3) warning("mpl: only first 3 parameters recognized");

}
    
void potential_double(int *ndim,double *pos,double *acc,double *pot,double *time)
{

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

    *pot = ((-G*ringmass)/(2*PI)) * fenergy(pos[X],pos[Y],pos[Z],a);

	// The accelerations

	double accmag;	
	accmag = (-G*ringmass / (PI*a*a)) * (xi*sintheta*pow(quant, -1.5)*f1val - pow(quant, -1.5)*f2val);
	
	acc[X] = accmag*cos(psi);
	acc[Y] = accmag*sin(psi);
	acc[Z] = (-G*ringmass / (PI*a*a)) * (xi*costheta*pow(quant, -1.5)*f1val);

}

// The integrals
double f1(double mu) {
	double dalpha = PI / 10000.0;
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
	double dalpha = PI / 10000.0;
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
	double dalpha = 2*PI / 20000.0;
	double alpha = 0.0;
	double quant = 0.0;
	double total = 0.0;

	for(alpha=0.0;alpha<=2*PI;alpha+=dalpha) {
		quant = alpha / pow((x-a*cos(alpha))*(x-a*cos(alpha)) + (y-a*sin(alpha))*(y-a*sin(alpha)) + z*z, 0.5);
		total += quant*dalpha;
	}

	return total;
}

