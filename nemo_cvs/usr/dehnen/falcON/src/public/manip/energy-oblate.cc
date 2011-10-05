#include <public/defman.h>
#include <public/basic.h>
#include <public/tools.h>
#include <public/io.h>
#include <public/bodyfunc.h>
#include <sstream>
#include <stdio.h>
#include <math.h>

#define X 0
#define Y 1
#define Z 2
#define PI 3.14159265358979

FILE *orbitfile;

int compare_doubles (const void *vala, const void *valb);

namespace falcON { namespace Manipulate {
  class bb : public manipulator {
  private:
	mutable output      OUT;
  public:
    const char*name    () const { return "bb"; }
    const char*describe() const {
      return "Prints out the bodies that are unbound";
    }
    //--------------------------------------------------------------------------
    fieldset need   () const { return fieldset::o; }
    fieldset provide() const { return fieldset::f; }
    fieldset change () const { return fieldset::f; }
    //--------------------------------------------------------------------------
    bb(const char*p, const char*f) : OUT (f? f : ".")  {
    //bb(const char*p, const char*f) falcON_THROWING  {
    //bb(const char *file) : OUT (file? file : ".")  {
	//falcON_WarningN("Hi from bb!\n");
	//OUT << "Hi from bb to file!\n";
    }
    //--------------------------------------------------------------------------
    bool manipulate(const snapshot*S) const {
	int i=0;
	int j=0;
	int k=0;	

	// Sgr properties TODO: make these parameters
	double N=10000;
	double rs = 0.9, mass=4498.66;

	double positions[10000][3];
	double velocities[10000][3];

	//falcON_WarningN("Hi from inside the bb manipulator!\n");
	//OUT << "Hi from inside the bb manipulator to the file!\n";
  	LoopSubsetBodies(S,b) {
		positions[i][0] = pos(b)[0];
		positions[i][1] = pos(b)[1];
		positions[i][2] = pos(b)[2];
		velocities[i][0] = vel(b)[0];
		velocities[i][1] = vel(b)[1];
		velocities[i][2] = vel(b)[2];
		i++;
  	}

        // Now, calculate the potential energy from the background 
        // From mpl4.c
        // TODO: Make these parameters


	// Loop through the particles
	for(i=0;i<=N-1;i++) {

                double apar, qpar, spar, ppar, lpar, rpar, rcyl;
                double G = 1;
                double omega = 0.0;                
                double miya_ascal = 6.5;
                double miya_bscal = 0.26;
                double miya_mass = 4.45865888E5;
                double plu_rc = 0.7;
                double plu_mass = 1.52954402E5;
                double vhalo = 116.5866;
                double q = 0.9;
                double d = 12.0;


		double r;
		double KE, PE;	
		double totenergy = 0.0;
		rcyl = sqrt(positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1]);
                qpar = sqrt(positions[i][2]*positions[i][2] + miya_bscal*miya_bscal);
                apar = miya_ascal+qpar;
                spar = (positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1]) + ((miya_ascal+qpar)*(miya_ascal+qpar));
                ppar = sqrt (positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1] + positions[i][2]*positions[i][2]) + plu_rc;
                rpar = ppar - plu_rc;
                lpar = (rcyl*rcyl) + ((positions[i][2]/q)*(positions[i][2]/q)) + (d*d);

		KE = 0.5*(mass/N)*(velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2]);

		PE = (mass/N) * ((-(miya_mass)/sqrt(spar)) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(lpar)));

		totenergy += KE;
		totenergy += PE;

		/*for(j=0;j<=N-1;j++) {
			if(i != j) {
				r = sqrt ((positions[i][0]-positions[j][0])*(positions[i][0]-positions[j][0]) + (positions[i][1]-positions[j][1])*(positions[i][1]-positions[j][1]) + (positions[i][2]-positions[j][2])*(positions[i][2]-positions[j][2]));
			
				totenergy += -((mass/N)*(mass/N))/r;
			}
		}*/
		
		OUT << S->time() << "	" << i << "	" << totenergy << "\n"; 
	}

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
__DEF__MAN__ALT(falcON::Manipulate::bb);
