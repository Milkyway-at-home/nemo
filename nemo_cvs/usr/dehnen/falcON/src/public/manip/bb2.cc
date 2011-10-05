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
	double BE, KE;
	double simtime;
	double * orbittime = NULL;
	double * step = NULL;
	double * orbitx = NULL, * orbity = NULL, * orbitz = NULL, * orbitvx = NULL, * orbitvy = NULL, * orbitvz = NULL;
	char line[500];

	double worbitx=0.0, worbity=0.0, worbitz=0.0, worbitvx=0.0, worbitvy=0.0, worbitvz=0.0;

	// Read in the orbit data
	orbitfile = fopen("orbit.dat", "rt");
		
	while(fgets(line,500,orbitfile)) {
		k++;
		step = (double*) realloc (step, k * sizeof(double));
		orbittime = (double*) realloc (orbittime, k * sizeof(double));
		orbitx = (double*) realloc (orbitx, k * sizeof(double));
		orbity = (double*) realloc (orbity, k * sizeof(double));
		orbitz = (double*) realloc (orbitz, k * sizeof(double));
		orbitvx = (double*) realloc (orbitvx, k * sizeof(double));
		orbitvy = (double*) realloc (orbitvy, k * sizeof(double));
		orbitvz = (double*) realloc (orbitvz, k * sizeof(double));

		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &step[k-1], &orbittime[k-1], &orbitx[k-1], &orbity[k-1], &orbitz[k-1], &orbitvx[k-1], &orbitvy[k-1], &orbitvz[k-1]);
	}


	fclose(orbitfile);
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
		simtime = S->time();

    		//OUT << S->time() << "	" << i << " X = " << positions[i][0] << " Y = " << positions[i][1] << " Z = " << positions[i][2] << "\n";
		i++;
  	}

	// Find the correct orbit point
	for(j=0;j<=k;j++) {
		if(abs(simtime - orbittime[j]) < 0.0001) {
			worbitx = orbitx[j];				
			worbity = orbity[j];				
			worbitz = orbitz[j];				
			worbitvx = orbitvx[j];				
			worbitvy = orbitvy[j];				
			worbitvz = orbitvz[j];				
		}
	}

	// Check to see if it found the orbit point
	if(worbitx == 0 && worbity == 0 && worbitz == 0 && worbitvx == 0 && worbitvy == 0 && worbitvz == 0) {
		printf("The orbit point was not found!\n");
	}

	//falcON_WarningN("AFTER READING BODIES: i=%i\n",i);

	// Loop through the particles
	for(i=0;i<=N-1;i++) {

		// Calculate the kinetic energy of the particle	
		KE = (0.5)*(mass/N)*((velocities[i][0]-worbitvx)*(velocities[i][0]-worbitvx) + (velocities[i][1]-worbitvy)*(velocities[i][1]-worbitvy) + (velocities[i][2]-worbitvz)*(velocities[i][2]-worbitvz));


		// Potential energy of the particle, assuming negligible mass loss of the dwarf and that the dwarf is precisely where the orbit predicts it would be
		double r;

		r = pow((positions[i][0] - worbitx)*(positions[i][0] - worbitx) + (positions[i][1] - worbity)*(positions[i][1] - worbity) + (positions[i][2]-worbitz)*(positions[i][2]-worbitz), 0.5);

		BE = - (((N-1.0)/N)*mass * (mass/N))/r;

		if(KE + BE > 0) {
			// If the KE + BE > 0, the particle is unbound	
    			OUT << S->time() << "	" << i << "\n";
		}
	}

	return false;
      }
  };
} }
////////////////////////////////////////////////////////////////////////////////
__DEF__MAN__ALT(falcON::Manipulate::bb);
