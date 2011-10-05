#include <omp.h>
#include <public/defman.h>
#include <public/basic.h>
#include <public/tools.h>
#include <public/io.h>
#include <public/bodyfunc.h>
#include <sstream>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979

namespace falcON { namespace Manipulate {
  class corefinder : public manipulator {
  private:
	mutable output      OUT;
  public:
    const char*name    () const { return "corefinder"; }
    const char*describe() const {
      return "Prints out the bodies that are unbound";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::f; }
    fieldset change () const { return fieldset::f; }
    //--------------------------------------------------------------------------
    corefinder(const char*p, const char*f) : OUT (f? f : ".")  {
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*S) const {
	int i=0, j=0, k=0;

	// Sgr properties TODO: make these parameters
	int N=10000;

	double positions[10000][3];
	double mag2[10000], mag1[10000], mag[10000];

	// Gradient search stuff
	double potential, newpotential;
	double gradpx, gradpy, gradpz;
	double x=0,y=0,z=0;
	double xnew,ynew,znew;
	double xold, yold, zold;
	double epsilon = 0.001;
	double h = 1.0;
	double MAXITER = 10000;
	double lambda = 1.0;
	int iter;
	
  	LoopSubsetBodies(S,b) {
		positions[i][0] = pos(b)[0];
		positions[i][1] = pos(b)[1];
		positions[i][2] = pos(b)[2];
		i++;
  	}

	// Gradient search for the minumum potential value
	// This is just the potential of the dwarf, we want to find the location of highest density

	for(iter=0;iter<=MAXITER;iter++) {

		potential = 0;
		gradpx = 0;
		gradpy = 0;
		gradpz = 0;

		// Loop through the particles and calculate the potentials and gradients (we can do this analytically)

		//#pragma omp parallel for
		for(i=0;i<=N-1;i++) {
			mag1[i] = pow(pow(positions[i][0] - x, 2) + pow(positions[i][1] - y, 2) + pow(positions[i][2] - z, 2), 0.5) + epsilon;
			
			potential += - 1.0 / mag1[i];
		}

		//#pragma omp parallel for	
		for(j=0;j<=N-1;j++) {	
			mag[j] = pow(pow(positions[j][0] - x, 2) + pow(positions[j][1] - y, 2) + pow(positions[j][2] - z, 2), 0.5) * pow(pow(pow(positions[j][0] - x, 2) + pow(positions[j][1] - y, 2) + pow(positions[j][2] - z, 2), 0.5) + epsilon, 2);

			gradpx += -(positions[j][0] - x)  / mag[j];
			gradpy += -(positions[j][1] - y)  / mag[j];
			gradpz += -(positions[j][2] - z)  / mag[j];
		}


		// New values
	
		xnew = x - h*lambda*gradpx;
		ynew = y - h*lambda*gradpy;
		znew = z - h*lambda*gradpz;
	
		newpotential = 0.0;

		// Calculate the new potential	
		
		//#pragma omp parallel for	
		for(k=0;k<=N-1;k++) {
			mag2[k] = pow(pow(positions[k][0] - xnew, 2) + pow(positions[k][1] - ynew, 2) + pow(positions[k][2] - znew, 2), 0.5) + epsilon;
			
			newpotential += - 1.0 / mag2[k];
		}

		// Check if the new value is smaller than the old one
		if(newpotential < potential) {
			xold = x;
			yold = y;
			zold = z;
			x = xnew;
			y = ynew;
			z = znew;
			lambda = lambda * 1.03;
		} else {
			xold = x;
			yold = y;
			zold = z;
			lambda = lambda * 0.8;
		}	
		
		//printf("OLDVALS: %f %f %f NEWVALS: %f %f %f\n", xold, yold, zold, xnew, ynew, znew);

		// The stop condition

		if(lambda < 0.0001)
			break;	
	}
		
	printf("Core: %f %f %f\n", xnew, ynew, znew);

	return false;

}
  };
} }

////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::corefinder);
