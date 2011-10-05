// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file   src/public/exe/s2a.cc
///
/// \author Walter Dehnen
/// \date   2002-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2002-2010 Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
// 
// history:
//
// v 0.0   29/10/2003  WD created
// v 0.1   10/02/2004  WD debugged (body.cc: SPH)
// v 1.0   10/02/2004  WD automized output; report history, run info 
// v 1.1   30/04/2004  WD improved automisation; new body.h; happy icc 8.0
// v 2.0   20/09/2004  WD write defaults to all present data; snapshot;
// v 2.1   24/10/2004  WD parameter header
// v 2.2   11/11/2004  WD format string in out; index= 
// v 2.3   20/05/2005  WD several minor updates
// v 3.0   24/06/2004  WD new falcON 
// v 3.1   27/06/2004  WD deBUGged
// v 4.0   15/11/2005  WD added filter (proprietary only)
// v 4.0.1 16/11/2005  WD improvements in filter (BodyFunc<>)
// v 4.1   04/07/2006  WD made public (along with bodyfunc.h)
// v 4.2   27/02/2007  WD print out "nan" or "inf" if float is nan or inf
// v 4.3   05/03/2007  WD enabled filter in public version
// v 4.4   30/08/2007  WD fivename() for field headers
// v 4.5   27/09/2007  WD Nsink output
// v 4.5.1 13/02/2008  WD minor bug (added "std::" to isnan & std::isinf)
// v 4.5.2 10/09/2008  WD happy gcc 4.3.1
// v 4.5.3 21/08/2009  WD duBUGged (very minor)
// v 4.6   09/02/2010  WD fixed problem with sink bodies
// v 4.7   12/03/2010  WD using snapshot::apply_filter, add keys if wanted
////////////////////////////////////////////////////////////////////////////////
#define falcON_VERSION   "4.7"
#define falcON_VERSION_D "12-mar-2010 Walter Dehnen                          "
//------------------------------------------------------------------------------
#ifndef falcON_NEMO                                // this is a NEMO program
#  error You need NEMO to compile "s2a"
#endif
#define falcON_RepAction 0                         // no action reporting
//------------------------------------------------------------------------------
#include <public/bodyfunc.h>                       // body functions
#include <main.h>                                  // NEMO basics & main 
#include <cstdio>                                  // C std I/O
//------------------------------------------------------------------------------
const char*defv[] = {
  "in=???\n         input snapshot file                                ",
  "out=-\n          file for ascii output (may contain format string)  ",
  "times=all\n      times to process                                   ",
  "write=\n         select data to write out (default: all)            ",
  "iformat=%d\n     format for integers                                ",
  "rformat=%14.6E\n format for real numbers                            ",
  "index=\n         first index if out contains format string          ",
  "header=t\n       write header                                       ",
  "filter=\n        bodyfunc filter (boolean): which bodies to write   ",
  "pars=\n          parameters for filter, if any                      ",
  "zeromissing=f\n  set body properties missing for filter to zero?    ",
  falcON_DEFV, NULL };
//------------------------------------------------------------------------------
const char*usage = "s2a -- Walter's simple snapshot to ascii converter";
//------------------------------------------------------------------------------
namespace {
  //----------------------------------------------------------------------------
  falcON::fieldset OUTPUT = falcON::fieldset::empty;
  FILE*            OUT    = 0;
  char             RFORMAT[32], IFORMAT[32];
  //----------------------------------------------------------------------------
  template<typename T> struct Print {
    static void print(T const&t) {
      fprintf(OUT,IFORMAT,(int)(t));
    } };
  //----------------------------------------------------------------------------
  template<> struct Print<falcON::real> {
    static void print(falcON::real const&t) {
      if     (std::isnan(t)) fprintf(OUT,"nan ");
      else if(std::isinf(t)) fprintf(OUT,"inf ");
      else                   fprintf(OUT,RFORMAT,t);
    } };
  //----------------------------------------------------------------------------
  template<> struct Print<falcON::vect> {
    static void print(falcON::vect const&t) {
      Print<falcON::real>::print(t[0]);
      Print<falcON::real>::print(t[1]);
      Print<falcON::real>::print(t[2]);
    } };
  //----------------------------------------------------------------------------
  template<typename T> 
  void print(T const&t) { Print<T>::print(t); }
  //----------------------------------------------------------------------------
  template<int BIT> struct BodyPrint    {
    static void act(falcON::body const&b) {
      if(OUTPUT.contain(falcON::fieldbit(BIT)))
	print(falcON::const_datum<BIT>(b));
    }
  };
  //----------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
void falcON::main() falcON_THROWING {
  nemo_in        IN   (getparam("in"));
  bool           KEYS (getioparam_z("write").contain(fieldbit::k));
  const fieldset NEED (getioparam_a("write"));
  fieldset       READ;
  snapshot       SHOT;
  int            INDEX(getiparam_z("index"));
  const bool     FOUT (hasvalue("out") && strcmp(getparam("out"),"-"));
  strcpy(IFORMAT,getparam("iformat")); strcat(IFORMAT," ");
  strcpy(RFORMAT,getparam("rformat")); strcat(RFORMAT," ");
  BodyFilter F(getparam_z("filter"),getparam_z("pars"));
  if(!FOUT) OUT = stdout;
  while(IN.has_snapshot()) {
    if(! SHOT.read_nemo(IN,READ,NEED,getparam("times"),0))
      continue;
    if(KEYS) {
      SHOT.add_field(fieldbit::k);
      READ |= fieldset::k;
    }
    SHOT.apply_filter(F,getbparam("zeromissing"));
    if(0 == SHOT.N_bodies()) {
      falcON_Warning("nobody left in snapshot after filtering at time %g\n",
		     SHOT.time());
      continue;
    }
    if(FOUT) {
      char FNAME[256];
      sprintf(FNAME,getparam("out"),INDEX++);
      if(strcmp(getparam("out"), FNAME)) {
	if(OUT) fclose(OUT);
	OUT = fopen(FNAME,"w");
      } else if(OUT==0)
	OUT = fopen(FNAME,"w");
    }
    OUTPUT = READ & NEED;
    if(getbparam("header")) {
      if(RunInfo::cmd_known())
	fprintf(OUT,"#\n# %s\n#\n", RunInfo::cmd());
      fprintf(OUT,"# run %s\n",RunInfo::time());
      if(RunInfo::user_known())
	fprintf(OUT,"#  by user  %s\n",RunInfo::user());
      if(RunInfo::host_known())
	fprintf(OUT,"#  on host  %s\n",RunInfo::host());
      if(RunInfo::pid_known())
	fprintf(OUT,"#  with pid %s\n",RunInfo::pid());
      fprintf(OUT,"#\n# time: %f\n"
	      "# Ntot: %d, Nsink:%d, Ngas: %d, Nstd: %d\n#\n#",
	      SHOT.time(), SHOT.N_bodies(),
	      SHOT.N_bodies(bodytype::sink),
	      SHOT.N_bodies(bodytype::gas),
	      SHOT.N_bodies(bodytype::std) );
      for(fieldbit f; f; ++f)
	if(OUTPUT.contain(f) && SHOT.have(f))
	  fprintf(OUT," '%s'",fivename(f));
      fprintf(OUT,"\n#\n");
    }
    for(bodytype t; t; ++t) {
      fieldset set=t.allows();
      LoopTypedBodies(&SHOT,b,t) {
	LoopFields<BodyPrint>::const_some(b,set);
	fprintf(OUT,"\n");
      }
    }
  }
  if(FOUT && OUT) fclose(OUT);
}
