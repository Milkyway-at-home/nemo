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
	int i=0; int j=0;
	int k=0;	

	// Sgr properties TODO: make these parameters
	double N=10000;
	double rs = 0.9, mass=4498.66;

	double positions[10000][3];
	double velocities[10000][3];
	double simtime;
	double * orbittime = NULL;
	double orbitstep;
	double * orbitx = NULL, * orbity = NULL, * orbitz = NULL, * orbitvx = NULL, * orbitvy = NULL, * orbitvz = NULL;
	char line[500];

	double worbitx=0.0, worbity=0.0, worbitz=0.0, worbitvx=0.0, worbitvy=0.0, worbitvz=0.0;

	// Read in the orbit data
	orbitfile = fopen("orbit.dat", "rt");
		
	while(fgets(line,500,orbitfile)) {
		k++;
		orbittime = (double*) realloc (orbittime, k * sizeof(double));
		orbitx = (double*) realloc (orbitx, k * sizeof(double));
		orbity = (double*) realloc (orbity, k * sizeof(double));
		orbitz = (double*) realloc (orbitz, k * sizeof(double));
		orbitvx = (double*) realloc (orbitvx, k * sizeof(double));
		orbitvy = (double*) realloc (orbitvy, k * sizeof(double));
		orbitvz = (double*) realloc (orbitvz, k * sizeof(double));

		sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf\n", &orbitstep, &orbittime[k-1], &orbitx[k-1], &orbity[k-1], &orbitz[k-1], &orbitvx[k-1], &orbitvy[k-1], &orbitvz[k-1]);
	}

	fclose(orbitfile);

	//falcON_WarningN("Hi from inside the bb manipulator!\n");
	//OUT << "Hi from inside the bb manipulator to the file!\n";
  	LoopSubsetBodies(S,b) {
		positions[i][0] = pos(b)[0];
		positions[i][1] = pos(b)[1];
		positions[i][2] = pos(b)[2];
		velocities[i][0] = vel(b)[0];
		velocities[i][1] = vel(b)[1];
		velocities[i][2] = vel(b)[2];
		simtime = S->time();
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

        // Now, calculate the potential energy from the background 
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

                rcyl = sqrt(worbitx*worbitx + worbity*worbity);
                qpar = sqrt(worbitz*worbitz + miya_bscal*miya_bscal);
                apar = miya_ascal+qpar;
                spar = (worbitx*worbitx) + (worbity*worbity) + ((miya_ascal+qpar)*(miya_ascal+qpar));
                ppar = sqrt ((worbitx*worbitx)+(worbity*worbity)+(worbitz*worbitz)) + plu_rc;
                rpar = ppar - plu_rc;
                lpar = (rcyl*rcyl) + ((worbitz/q)*(worbitz/q)) + (d*d);


	// Find the tidal radius (the place where the galactic force equals the dwarf force)
	double mag; 
	double acc[3], BGforce, Dforce, rfromgal, rfromdwarf, enclosedmass;
	double rt = 0;

	for(mag=0;mag<=1;mag+=0.00001) {
		// Find the background force on a test particle of mass = Mdwarf/N

                rcyl = sqrt(mag*worbitx*mag*worbitx + mag*worbity*mag*worbity);
                qpar = sqrt(mag*worbitz*mag*worbitz + miya_bscal*miya_bscal);
                apar = miya_ascal+qpar;
                spar = (mag*worbitx*mag*worbitx) + (mag*worbity*mag*worbity) + ((miya_ascal+qpar)*(miya_ascal+qpar));
                ppar = sqrt ((mag*worbitx*mag*worbitx)+(mag*worbity*mag*worbity)+(mag*worbitz*mag*worbitz)) + plu_rc;
                rpar = ppar - plu_rc;
                lpar = (rcyl*rcyl) + ((mag*worbitz/q)*(mag*worbitz/q)) + (d*d);

		acc[X] = - ( ( (2.0*vhalo*vhalo*mag*worbitx)/(lpar) ) + ( (plu_mass*mag*worbitx)/(rpar*ppar*ppar) ) + ( (miya_mass*mag*worbitx)/(pow(spar,1.5)) ) );
		acc[Y] = - ( ( (2.0*vhalo*vhalo*mag*worbity)/(lpar) ) + ( (plu_mass*mag*worbity)/(rpar*ppar*ppar) ) + ( (miya_mass*mag*worbity)/(pow(spar,1.5)) ) );
		acc[Z] = - ( ( (2.0*vhalo*vhalo*mag*worbitz)/(q*q*lpar) ) + ( (plu_mass*mag*worbitz)/(rpar*ppar*ppar) ) + ( (miya_mass*mag*worbitz)/(pow(spar,1.5)) ) );

		rfromgal = pow(((mag)*worbitx)*((mag)*worbitx) + ((mag)*worbity)*((mag)*worbity) + ((mag)*worbitz)*((mag)*worbitz), 0.5);
		rfromdwarf = pow(((mag)*worbitx)*((1-mag)*worbitx) + ((1-mag)*worbity)*((1-mag)*worbity) + ((1-mag)*worbitz)*((1-mag)*worbitz), 0.5);

		BGforce = (mass/N) * pow(acc[X]*acc[X] + acc[Y]*acc[Y] + acc[Z]*acc[Z], 0.5);

		enclosedmass = (mass * rfromdwarf * rfromdwarf) / (rs*rs + rfromdwarf*rfromdwarf); // Enclosed plummer mass given by Sparke & Gallagher p 310

		// Find the dwarf force on a test particle of mass = Mdwarf / N
		Dforce = (enclosedmass * (mass/N)) / (rfromdwarf*rfromdwarf); 

		if(abs(BGforce - Dforce) < 10) {
			rt = rfromdwarf;
			break;
		}	
		
	}

	// Check
	if(rt != 0.0) {
		printf("The tidal radius: %lf\n", rt);
	} else {
		printf("There was a problem finding the tidal radius\n");
	}

	// Loop through the particles
	for(i=0;i<=N-1;i++) {
		double r;
		r = pow((positions[i][0] - worbitx)*(positions[i][0] - worbitx) + (positions[i][1] - worbity)*(positions[i][1] - worbity) + (positions[i][2]-worbitz)*(positions[i][2]-worbitz), 0.5);

		rcyl = sqrt(positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1]);
                qpar = sqrt(positions[i][2]*positions[i][2] + miya_bscal*miya_bscal);
                apar = miya_ascal+qpar;
                spar = (positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1]) + ((miya_ascal+qpar)*(miya_ascal+qpar));
                ppar = sqrt (positions[i][0]*positions[i][0] + positions[i][1]*positions[i][1] + positions[i][2]*positions[i][2]) + plu_rc;
                rpar = ppar - plu_rc;
                lpar = (rcyl*rcyl) + ((positions[i][2]/q)*(positions[i][2]/q)) + (d*d);

		double totenergy, KE, PE, theta;	

		KE = 0.5*(mass/N)*(velocities[i][0]*velocities[i][0] + velocities[i][1]*velocities[i][1] + velocities[i][2]*velocities[i][2]);

		PE = (mass/N) * ((-(miya_mass)/sqrt(spar)) + (-(plu_mass)/ppar) + (vhalo*vhalo*log(lpar)));

		totenergy = KE + PE;

		if(r > 5*rt) {
			// If r > 2*rt, the particle is unbound	
    			OUT << S->time() << "	" << i <<  "	" << totenergy << "\n";
		}
	}

	// Memory management
	free(orbittime);
	free(orbitx);
	free(orbity);
	free(orbitz);
	free(orbitvx);
	free(orbitvy);
	free(orbitvz);
		
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
