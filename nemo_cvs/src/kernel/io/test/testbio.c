/*
 * Test of Nemo's blocked i/o stuff
 *      2-jun-05   created       PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>


string defv[] = {
    "in=\n         Input file (if provided)",
    "out=\n        output file (if provided)",
    "real=1.234\n  Real number to write (if out=)",
    "count=1\n     Number of reals/block to write",
    "blocks=1\n    Number of blocks to write",
    "set=t\n       Wrap inside of a set/tes?"
    "VERSION=1.4\n 7-jun-05 PJT",
    NULL,
};

string usage = "testing NEMO's blocked I/O";

string cvsid="$Id: testbio.c,v 1.4 2005/06/07 20:41:49 pteuben Exp $";


extern void get_nanf(float *);
extern void get_nand(double *);

static bool Qset;
static string setName = "TestingBlockIO";


int get_number(string name, real *x, int nb)
{
    stream str;
    real *buf;
    int  i,n,m;

    str = stropen(name,"r");
    get_history(str);
    if (Qset)  get_set(str,setName);

    if (!get_tag_ok(str,"n"))
        error("%s: missing \"n\", Not a valid TESTIO input file",name);
    get_data(str,"n",IntType,&n,0);
    buf = (real *) allocate(n*sizeof(real));

    if (!get_tag_ok(str,"x"))
        error("%s: missing \"x\", Not a valid TESTIO input file",name);
    if (nb==1)
      get_data(str,"x",RealType,buf,n,0);
    else {
      get_data_set(str,"x",RealType,n,0);
      m = n/nb;
      for (i=0; i<nb; i++)
	get_data_blocked(str,"x",buf,m);
      get_data_tes(str,"x");
    }
    if (Qset)  get_tes(str,setName);
    strclose(str);

    *x = *buf;
    dprintf(1,"Read number %f from file %s\n",*x,name);
    for (i=0; i<n; i++)
      dprintf(2,"buf[%d] = %g\n",i,buf[i]);
    free( (char *) buf);
    return n;
}

void put_number(string name, real x, int n, int r)
{
    stream str;
    real *buf;
    int nout,i;

    str = stropen(name,"w");
    put_history(str);
    nout = n*r;
    if (Qset)  put_set(str,setName);
    put_data(str,"n",IntType,&nout,0);
    buf = (real *) allocate(nout*sizeof(real));
    for (i=0; i<nout; i++)
      buf[i] = x+i;
    if (r==1) {
        put_data(str,"x",RealType,buf,n,0);
        put_data(str,"y",RealType,buf,n,0);
    } else {
#if 1
      // the right way
        put_data_set(str,"x",RealType,nout,0);
        for (i=0;i<r;i++)
            put_data_blocked(str,"x",buf,n);
        put_data_tes(str,"x");

        put_data_set(str,"y",RealType,nout,0);
        for (i=r;i>0;i--)
            put_data_blocked(str,"y",buf,n);
        put_data_tes(str,"y");
#else
	// this should fail (for now)
        put_data_set(str,"y",RealType,nout,0);
        put_data_set(str,"x",RealType,nout,0);
        for (i=0;i<r;i++) {
            put_data_ran(str,"x",&buf[i*n],i,n);
            put_data_ran(str,"y",&buf[i*n],i,n);
	}
        put_data_tes(str,"x");
        put_data_tes(str,"y");
#endif
    }
    if (Qset) put_tes(str,setName);
    strclose(str);
    dprintf(1,"Wrote number %f to file %s\n",x,name);
    free( (char *) buf);
}

/* -------------------------------------------------------------------------- */

void nemo_main()
{
    real x,y = 0;
    int count, random;
    string name;
    int work = 0;
    bool Qset = getbparam("set");

    x = getdparam("real");    
    count = getiparam("count");
    random = getiparam("blocks");
    
    name = getparam("in");
    if (*name) {
      work++;
      get_number(name,&y,random);
    } else
      dprintf(1,"No input file specified\n");
    
    name = getparam("out");
    if (*name) {
      work++;
      put_number(name,y+x,count,random);
    } else
      dprintf(1,"No output file specified\n");

    if (work==0) warning("No work done, use in=, out=, or both");
}
