/*
 *  Convert NBODYx output datasets ('OUT3' files) to snapshots
 *
 *   7-apr-93   1.0  Created
 *   9-apr-93   1.1  NB3 to include NK length of A() array      PJT
 *  29-mar-94   1.2  optional alen= override; also ansi         pjt
 *  23-may-95   1.2a changed routine names in nbody_io.f        pjt
 *  10-jun-95   1.2c fixed default of a blank alen=		pjt
 *   8-aug-95   1.3  added a swap= keyword if to swap bytes     pjt
 *  16-jun-97   1.3a warning swap/FIO
 *  28-feb-06   1.4  with the new nbody4 a mode= is needed      pjt
 *                   and in general header= needs known        
 *   4-mar-06   1.5  header= now blank default                  pjt


nbody1
      WRITE (3)  N, MODEL, NRUN, NK
      WRITE (3)  (A(K),K=1,NK), (BODY(J),J=1,N),
     &           ((XS(K,J),K=1,3),J=1,N), ((XDOT(K,J),K=1,3),J=1,N),
     &           (NAME(J),J=1,N)
nbody2
      WRITE (3)  N, MODEL, NRUN, NK
      WRITE (3)  (A(K),K=1,NK), (BODY(J),J=1,N),
     &           ((XS(K,J),K=1,3),J=1,N), ((XDOT(K,J),K=1,3),J=1,N),
     &           (NAME(J),J=1,N)
nbody4
      WRITE (3)  NTOT, MODEL, NRUN, NK
      WRITE (3)  (AS(K),K=1,NK), (BODYS(J),J=1,NTOT),
     &           ((XS(K,J),K=1,3),J=1,NTOT), ((VS(K,J),K=1,3),J=1,NTOT),
     &           (RHO1(J),J=1,NTOT),(PHI1(J),J=1,NTOT),
     &           (NAME(J),J=1,NTOT),(KSTAR(J),J=1,NTOT)


 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <snapshot/snapshot.h>
#include <unfio.h>
#include "nbody_io.h"

string defv[] = {
    "in=???\n       Input file (in NBODY's OUT3 format)",
    "out=???\n      Output file (snapshot(5NEMO) format)",
    "frame=0\n      Frames to read (0=all)",
    "alen=\n        If given, override length of A array (0=do not read A)",
    "swap=f\n       Swap byte (big vs. little endian only)",
    "mode=1\n       NBODYx mode (valid are 1,2,4)",
    "key=name\n     snapshot Key comes from 'name' or 'key'?",
    "header= \n     if used, force unfio header size (4 or 8)",
    "integer=4\n    Size of integers in dataset (2 or 4) ** 2 is deprecated **",
    "VERSION=1.4\n  1-mar-06 PJT",
    NULL,
};

string usage = "Convert NBODY output to snapshot";

string cvsid="$Id: u3tos.c,v 1.4 2006/03/05 05:05:48 pteuben Exp $";


#define MAXFRAME  512
#define MAXHEADER  32

extern int nemo_file_size(string);


void nemo_main(void)
{
    string fname = getparam("in");
    string keysel = getparam("key");
    int ilen, alen, saved=0, alen_fix=-1;
    int frame[MAXFRAME], nframe, iframe;
    int nbody, model, run, *name, *key, i, j, k, ibody, mode, nwflt, nwint;
    int coordsys = CSCode(Cartesian, 3, 2);
    float *mass, *pos, *vel, *phi, *aux, a[MAXHEADER];
    real *rmass, *rphase, tsnap;
    bool Qswap;
    stream outstr;

    if (hasvalue("header"))
      unfsize(getiparam("header"));   

    nframe = nemoinpi(getparam("frame"),frame,MAXFRAME);
    if(nframe<0) error("Parsing frame=%s",getparam("frame"));
    if(nframe==1 && frame[0]==0) nframe=0;
    if (hasvalue("alen")) alen_fix = getiparam("alen");
    Qswap = getbparam("swap");
    ilen = getiparam("integer");
    mode = getiparam("mode");

    if (nemo_file_size(fname) < 0) error("File %s does not exist",fname);
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);

#ifdef FIO
    nb3open_(fname,strlen(fname));
    if (Qswap) warning("-DFIO compiled cannot swap");
#else
    nb3open_c(fname, ilen, Qswap);
#endif

    for(iframe=0;;) {                         /* Loop over all frames */
      alen = alen_fix;
#ifdef FIO
      nb3header_(&nbody, &model, &run, &alen);
#else
      nb3header_c(&nbody, &model, &run, &alen);
#endif
      if (nbody == 0) break;     /* end of file */
      if (nbody<0 || nbody > 1000000) 
        warning("Strange value for nbody=%d; perhaps need to change to swap=%c",
                nbody, Qswap ? "false" : "true");
      if (alen_fix > 0) alen = alen_fix;       /* override alen */
      nwflt = (mode == 4 ?  9  :  6);
      nwint = (mode == 4 ?  2  :  1);
      dprintf(1,"Header: nbody=%d model=%d run=%d nk=%d Snapshotsize=%d bytes\n",
                nbody,model,run,alen,
                32 + sizeof(float)*(alen+nwflt*nbody) + sizeof(int)*nwint*nbody);
      if (nbody < 0) break;      /* something bad surely */
      if (alen<1 || alen>MAXHEADER)
         error("Bad headerlength nk=%d\n",alen);

      mass  = (float *) allocate(nbody*sizeof(float));
      pos   = (float *) allocate(nbody*sizeof(float)*3);
      vel   = (float *) allocate(nbody*sizeof(float)*3);
      phi   = (float *) allocate(nbody*sizeof(float));    /* only for nbody 4 */
      aux   = (float *) allocate(nbody*sizeof(float));    /* only for nbody 4 */
      name  = (int *)   allocate(nbody*sizeof(int));
      key   = (int *)   allocate(nbody*sizeof(int));      /* only for nbody 4 */
#ifdef FIO
      nb3data_(&nbody,&alen,a,mass,pos,vel,name);
#else
      nb3data_c(&nbody,&alen,&mode,a,mass,pos,vel,phi,aux,name,key);
#endif

      dprintf(1,"Data  : a(%d)=",alen);
      for (i=0; i<alen; i++)
        dprintf(1," %g",a[i]);
      dprintf(1,"\n");

      if (nframe==0 || model==frame[iframe]) {
        tsnap = a[0];
        dprintf(0,"Saving model %d (nbody=%d time=%g)\n",model,nbody,tsnap);
        rmass  = (real *) allocate(nbody*sizeof(real));
        rphase = (real *) allocate(nbody*sizeof(real)*6);
        for (ibody=0, i=0, j=0, k=0; ibody<nbody ; ibody++) {
            rmass[ibody] = mass[ibody];
            rphase[k++] = pos[i++];
            rphase[k++] = pos[i++];
            rphase[k++] = pos[i++];
            rphase[k++] = vel[j++];
            rphase[k++] = vel[j++];
            rphase[k++] = vel[j++];
        }
        put_set(outstr, SnapShotTag);
        put_set(outstr, ParametersTag);
        put_data(outstr, NobjTag, IntType, &nbody, 0);
        put_data(outstr, TimeTag, RealType, &tsnap, 0);
        put_tes(outstr, ParametersTag);
        put_set(outstr, ParticlesTag);
        put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
        put_data(outstr, MassTag, RealType, rmass, nbody, 0);
        put_data(outstr,PhaseSpaceTag,RealType,rphase,nbody,2,3,0);
	if (mode == 4) {
	  for (ibody=0; ibody<nbody ; ibody++)
	    rmass[ibody] = phi[ibody];
	  put_data(outstr, PotentialTag, RealType, rmass, nbody, 0);

	  for (ibody=0; ibody<nbody ; ibody++)
	    rmass[ibody] = aux[ibody];
	  put_data(outstr, AuxTag, RealType, aux, nbody, 0);
	}
	if (*keysel == 'n')
	  put_data(outstr,KeyTag,IntType,name,nbody,0);
	else if (*keysel == 'k')
	  put_data(outstr,KeyTag,IntType,key,nbody,0);
        put_tes(outstr, ParticlesTag);
        put_tes(outstr, SnapShotTag);
        free(rmass);
        free(rphase);
        iframe++;
        saved++;
      }

      free(mass);
      free(pos);
      free(vel);
      free(name);

      if(nframe>0 && iframe>=nframe) break;      /* early stop */
    }

#ifdef FIO
    nb3close_();
#else
    nb3close_c();
#endif
    strclose(outstr);

    dprintf(0,"Saved %d snapshots\n",saved);

}
    
