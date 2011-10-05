/*
 * SNAPTIME.C: copy an N-body snapshot.
 *	
 *      15-nov-2006      cloned from snapcopy but never committed    PJT
 *      30-apr-2008      committed
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
  "in=???\n           Input snapshot",
  "out=???\n          Output snapshot",
  "time=0.0\n         Time to set",
  "VERSION=0.1\n      30-apr-08 PJT",
  NULL,
};

string usage="copy an N-body snapshot with a new time";


void nemo_main(void)
{
  stream instr, outstr;
  real   tsnap, tsnap0;
  string times;
  Body   *btab = NULL, *bpi, *bpo;
  int    i, nbody, bitsi, nsnap = 0;

  tsnap0 = getdparam("time");
  instr = stropen(getparam("in"), "r");
  outstr = stropen(getparam("out"), "w");

  get_history(instr);
  put_history(outstr);		

  for (;;) {                /* grab each snapshot */
    get_history(instr);		/* skip over stuff we can forget */
    if (!get_tag_ok(instr, SnapShotTag))
      break;			/* done with work in loop */
    get_snap(instr, &btab, &nbody, &tsnap, &bitsi);
    if ((bitsi & MassBit) == 0 && (bitsi & PhaseSpaceBit) == 0) {
      continue;       /* just skip */
    }
    put_snap(outstr, &btab, &nbody, &tsnap0, &bitsi);
    nsnap++;
  }
  if (nsnap > 1) warning("%d snapshots were given time=%g",nsnap,tsnap0);
}
