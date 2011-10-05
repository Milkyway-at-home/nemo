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
	
	// Sgr properties TODO: make these parameters
	double N=10000;
	double rs = 0.9, mass=4498.66;

	double positions[10000][3];
	double velocities[10000][3];
	double BE, KE, BGE;
	double simtime[10000];
	double orbittime;
	double orbitx, orbity, orbitz, orbitvx, orbitvy, orbitvz;

	//falcON_WarningN("Hi from inside the bb manipulator!\n");
	//OUT << "Hi from inside the bb manipulator to the file!\n";
  	LoopSubsetBodies(S,b) {
    		//W += mass(b);
		positions[i][0] = pos(b)[0];
		positions[i][1] = pos(b)[1];
		positions[i][2] = pos(b)[2];
		velocities[i][0] = vel(b)[0];
		velocities[i][1] = vel(b)[1];
		velocities[i][2] = vel(b)[2];
		simtime[i] = S->time();

    		//OUT << S->time() << "	" << i << " X = " << positions[i][0] << " Y = " << positions[i][1] << " Z = " << positions[i][2] << "\n";
		i++;
  	}
	//falcON_WarningN("AFTER READING BODIES: i=%i\n",i);

	// Calculate the total binding energy (per particle) of the plummer sphere (assuming negligible mass loss and size change)
	//BE = (3*PI*mass*mass)/(32*rs*N); // From Sparke & Gallagher, 2nd ed., p. 120
	BE = (mass*(mass/N))/(rs);

	// Loop through the particles
	for(i=0;i<=N-1;i++) {
		
		// Calculate the kinetic energy of the particle	
		KE = (0.5)*(mass/N)*(velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2]);

		/* Now, calculate the potential energy from the background 
		// From mpl4.c
		// TODO: Make these parameters

		double apar, qpar, spar, ppar, lpar, rpar, rcyl;
		double G = 1;
		double omega = 0.0;                
		double miya_ascal = 6.5;
		double miya_bscal = 0.26;
		double miya_mass = 4.45865888E5;
		double plu_rc = 0.7;
		double plu_mass = 1.52954402E5;
		double vhalo = 116.5866;
		double q = 1.0;
		double d = 12.0;

		rcyl = sqrt(positions[i][X]*positions[i][X] + positions[i][Y]*positions[i][Y]);
    		qpar = sqrt(positions[i][Z]*positions[i][Z] + miya_bscal*miya_bscal);
    		apar = miya_ascal+qpar;
    		spar = (positions[i][X]*positions[i][X]) + (positions[i][Y]*positions[i][Y]) + ((miya_ascal+qpar)*(miya_ascal+qpar));
    		ppar = sqrt ((positions[i][X]*positions[i][X])+(positions[i][Y]*positions[i][Y])+(positions[i][Z]*positions[i][Z])) + plu_rc;
    		rpar = ppar - plu_rc;
    		lpar = (rcyl*rcyl) + ((positions[i][Z]/q)*(positions[i][Z]/q)) + (d*d);
		*/

		BGE = (mass/N)*(((miya_mass)/sqrt(spar)) + ((plu_mass)/ppar) - (vhalo*vhalo*log(lpar)));
		
		//falcON_WarningN("BE = %f KE = %f BG = %f\n",BE,KE,BGE);

		// Roche limit
		double roche;
	
		roche = rs*pow(2*(miya_mass + plu_mass + 100.0*vhalo*vhalo)/mass, 0.33);
		if(rpar < roche) {
		// If the KE + BGE > BE, the particle is unbound	
		//if(KE < BE) {
			//falcON_WarningN("GOT ONE!\n");
    			OUT << S->time() << "	" << i << "\n";
		//}
		}
	}

	return false;
      }
  };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::bb);
