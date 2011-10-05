#include <public/defman.h>
#include <public/basic.h>
#include <public/tools.h>
#include <public/io.h>
#include <public/bodyfunc.h>
#include <sstream>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979

int compare_doubles (const void *vala, const void *valb);

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
	int i=0, j=0;

	// Sgr properties TODO: make these parameters
	double N=10000;
	double rs = 0.9, mass=4498.66;
	double thresh = 1.0; 

	double positions[10000][3];
	double velocities[10000][3];
	int numberofneighbors[10000][2];
	double mag;
	char line[500];
	int dp;

  	LoopSubsetBodies(S,b) {
		positions[i][0] = pos(b)[0];
		positions[i][1] = pos(b)[1];
		positions[i][2] = pos(b)[2];
		velocities[i][0] = vel(b)[0];
		velocities[i][1] = vel(b)[1];
		velocities[i][2] = vel(b)[2];
		i++;
  	}

	// Loop through the particles
	for(i=0;i<=N-1;i++) {
		for(j=0;j<=N-1;j++) {
			mag = pow(pow(positions[i][0] - positions[j][0], 2) + pow(positions[i][1] - positions[j][1], 2) + pow(positions[i][2] - positions[j][2], 2), 0.5);
			if(mag < thresh) {
				numberofneighbors[i][0]++;
				numberofneighbors[i][1] = i;
			}
		}
	}
	
	// sort the numberofneighbors list

	qsort(numberofneighbors, N-1, sizeof(double), compare_doubles);

	dp = numberofneighbors[0][1];	

	printf("Core: %f %f %f\n", positions[dp][0], positions[dp][1], positions[dp][2]);	

	return false;
      }
  };
} }

int compare_doubles (const void *vala, const void *valb)
{
        double da, db;
        da = *(double *)vala;
        db = *(double *)valb;

        if(da < db) {
                return -1;
        } else if (da > db) {
                return +1;
        }
        return 0;
}

////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::corefinder);
