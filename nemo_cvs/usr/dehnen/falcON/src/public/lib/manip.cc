// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
//
/// \file /src/public/manip.cc
//
// Copyright (C) 2004-2011  Walter Dehnen
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
// version 0.0  17/09/2004 WD  created based on acceleration.cc v3.2
// version 1.0  19/05/2005 WD  allowed for manips given in single file
// version 1.1  20/05/2005 WD  added manpath as 4th argument to Manipulator
// version 1.2  14/06/2005 WD  fixed bug (if empty arguments)
// version 1.3  22/06/2005 WD  no $FALCON, nemo_dprintf -> debug_info
// version 1.4  12/07/2005 WD  if manippath given, do not search elsewhere
// version 1.5  08/11/2005 WD  if no manippath given, try $MANIPPATH
// version 1.6  04/08/2006 WD  if manippath given, put it in top of seach path
// version 1.7  10/06/2008 WD  debug_info -> DebugInfo 
// version 1.8  10/09/2008 WD  nemo++.h to avoid #including nemo headers
// version 1.9  18/03/2009 WD  avoid warnings from -Wshadow
// version 1.10 08/03/2011 WD  debugged.
////////////////////////////////////////////////////////////////////////////////
#include <public/manip.h>              // the header we are implementing
#include <public/nemo++.h>             // the header we are implementing
#include <fstream>                     // C++ file I/O
#include <cstring>                     // C type string manipultions
#ifdef falcON_MPI
# include <parallel/parallel.h>
#endif
////////////////////////////////////////////////////////////////////////////////
int falcON::Manipulator::parse(const char*params, double*pars, int maxp)
{
  if(params == 0 || *params == 0) return 0;
  int npar = nemoinp(params,pars,maxp);
  if(npar > maxp)
    falcON_THROW("Manipulator::parse(): too many parameters (%d > %d)",
		 npar,maxp);
  if(npar < 0)
    falcON_THROW("Manipulator::parse(): parsing error in parameters: \"%s\"",
		 params);
  return npar;
}
////////////////////////////////////////////////////////////////////////////////
int falcON::Manipulator::parse(char*data, char sep, char**list, int nmax)
{
  int  n  = 0;
  list[0] = data;
  char *d = data;
  for(; *d && n!=nmax; ++d)
    if(*d == sep) {
      *d = 0;
      list[++n] = d+1;
    }
  if(*d && n==nmax)
    for(; *d; ++d) if(*d==sep) ++n;
  return n+1;
}      
////////////////////////////////////////////////////////////////////////////////
namespace {
  using namespace falcON;
//   using falcON::mapsys;
//   using falcON::findfn;
//   using falcON::loadobj;
//   using falcON::localsymbols;
//   using falcON::getparam;
//   using falcON::NewArray;
//   using falcON::exception;
//   using falcON::message;
//   using falcON::manipulator;
//   using falcON::Manipulator;
  //////////////////////////////////////////////////////////////////////////////
  //                                                                            
  // single_manipulator()                                                       
  //                                                                            
  // 1. we try to match manname against those already successfully done;        
  //    if match found, we simply call the corresponding inimanip() and return  
  //                                                                            
  // 2. we try to load the routines inimanip() and, if successful, call it.     
  //                                                                            
  //////////////////////////////////////////////////////////////////////////////
  /// type of pointer to inimanip(), see file defman.h
#ifdef MANIP_PARSE_AT_INIMANIP
  typedef void(*iniman_pter) (const manipulator**, const char*, const char*);
#else
  typedef void(*iniman_pter) (const manipulator**, const double*, int,
			      const char*);
#endif
  /// boolean indicating whether we have loaded local symbols
  bool first = true;

  typedef void(*proc)();
  // routine for finding a function  "NAME", "NAME_" or "NAME__"
  inline proc findfunc(const char*func)
  {
    char fname[256];
    strcpy(fname,func);
    mapsys(fname);
    proc f = findfn(fname);
    for(int i=0; f==0 && i!=2; ++i) {
      strcat(fname,"_");
      f= findfn(fname);
    }
    return f;
  }
  //----------------------------------------------------------------------------
  const unsigned IniMnMax=256;
  const unsigned MnNamMax=128;
  unsigned       IniMnInd=0;
  char           MnNames[IniMnMax][MnNamMax];
  iniman_pter    IniMn[IniMnMax] = {0};
  //----------------------------------------------------------------------------
  void single_manipulator(const manipulator*&manip,
			  const char*manname,
			  const char*manpars,
			  const char*manfile,
			  const char*manpath) falcON_THROWING
  {
    //
    // NOTE: the present implementation will NOT try to compile a source code
    //       but abort if no .so file is found.
    //
    DebugInfo(3,"Manipulator: trying to initialize manipulator with\n"
	      "  name=\"%s\",\n  pars=\"%s\",\n  file=\"%s\"\n",
	      manname,manpars,manfile);
    // 1. parse the parameters
#ifndef MANIP_PARSE_AT_INIMANIP
    const int MAXPAR = 256;
    double pars[MAXPAR];
    int    npar = Manipulator::parse(manpars,pars,MAXPAR);
#endif
    // 2. load local symbols
    if(first) {
      localsymbols();
      first = false;
    }
    // 3. try to find manname in list of mannames already done
    for(unsigned i=0; i!=IniMnInd; ++i)
      if(0 == strcmp(manname, MnNames[i])) {
	DebugInfo(3,"Manipulator: name=\"%s\": known already: "
		  "no need to load it again\n",manname);
#ifdef MANIP_PARSE_AT_INIMANIP
	(*IniMn[i])(&manip,manpars,manfile);
#else
	(*IniMn[i])(&manip,pars,npar,manfile);
#endif
	return;
      }

    // 4. seek file manname.so and load it
    // 4.1 put search path together
    char manpaths[2048] = {0};
    if(manpath && *manpath) {                     // get input path
      strcpy(manpaths,manpath);
      strcat(manpaths,"/:");
    }
    strcat(manpaths,".");                         // try "."
    const char *path;
    path = falcON::directory();                   // try falcON/manip/
    if(path) {
      strcat(manpaths,":");
      strcat(manpaths,path);
#ifdef falcON_MPI
      if(falcON::MPI::Initialized())
	strcat(manpaths,"/parallel");             // try $falcON/parallel/manip/
#endif
      strcat(manpaths,"/manip/");
    }
    path = getenv("MANIPPATH");                   // try $MANIPPATH
    if(path) {
      strcat(manpaths,":");
      strcat(manpaths,path);
    }
    // 4.2 seek file in path and load it
    char name[256];
    strcpy(name,manname);
    strcat(name,".so");
    DebugInfo(3,"Manipulator: searching file \"%s\" in path \"%s\" ...\n",
	      name,manpaths);
    const char*fullname = pathfind(manpaths,name);// seek for file in manpaths
    if(fullname == 0)
      falcON_THROW("Manipulator: cannot find file \"%s\" in path \"%s\"",
		   name,manpaths);
    DebugInfo(3,"             found one: \"%s\"; now loading it\n",fullname);
    loadobj(fullname);

    // 5. try to get inimanip()
    //    if found, remember it, call it and return manipulator given by it
    iniman_pter im = (iniman_pter) findfunc("inimanip");
    if(im) {
      if(IniMnInd < IniMnMax && strlen(manname) < MnNamMax) {
	strcpy(MnNames[IniMnInd],manname);
	IniMn[IniMnInd] = im;
	IniMnInd++;
      }
#ifdef MANIP_PARSE_AT_INIMANIP
      im (&manip,manpars,manfile);
#else
      im (&manip,pars,npar,manfile);
#endif
#ifdef falcON_MPI
      if(falcON::MPI::Initialized() && !manip->is_mpi())
	falcON_THROW("Manipulator \"%s\" does not support MPI\n",manip->name());
#endif
    } else {
      manip = 0;
      falcON_THROW("Manipulator: cannot find function \"inimanip\" "
		   "in file \"%s\"",fullname);
    }
  } // single_manipulator()
  //////////////////////////////////////////////////////////////////////////////
  inline bool is_sep(char const&c, const char*seps) {
    for(const char*s=seps; *s; ++s)
      if( c == *s ) return true;
    return false;
  }
  //----------------------------------------------------------------------------
  template<int NMAX>
  int splitstring(char*data, char*list[NMAX], const char*seps)
  {
    int     n = 0;
    list[n]   = data;
    for(char* d=data; *d!=0 && n!=NMAX; ++d)
      if(is_sep(*d,seps)) {
	*d = 0;
	list[++n] = d+1;
      }
    return n+1;
  }
  //////////////////////////////////////////////////////////////////////////////
  char NameSeps[3] = {',','+',0};
  char ParsSeps[3] = {';',' ',0};
  char FileSeps[3] = {';',' ',0};

} // namespace {
////////////////////////////////////////////////////////////////////////////////
//                                                                              
// class Manipulator                                                            
//                                                                              
////////////////////////////////////////////////////////////////////////////////
falcON::Manipulator::Manipulator(const char*mannames,
				 const char*manparss,
				 const char*manfiles,
				 const char*manpaths) falcON_THROWING :
  N     ( 0 ),
  NEED  ( fieldset::empty ),
  CHNG  ( fieldset::empty ),
  PRVD  ( fieldset::empty ),
  IS_MPI( true )
{
  if((mannames == 0 || *mannames == 0) &&
     (manfiles == 0 || *manfiles == 0)) return;
  const static unsigned Lnames=NMAX*16, Lparss=NMAX*32, Lfiles=NMAX*32;
  char *_names=falcON_NEW(char,Lnames), *_name[NMAX]={0};
  char *_parss=falcON_NEW(char,Lparss), *_pars[NMAX]={0};
  char *_files=falcON_NEW(char,Lfiles), *_file[NMAX]={0};
#define CLEANUP					\
  falcON_DEL_A(_names);				\
  falcON_DEL_A(_parss);				\
  falcON_DEL_A(_files);
  // 1 parse input into arrays of _name[], _pars[], _file[]
  if(mannames == 0 || *mannames == 0 || *mannames=='+') { 
    // 1.1 no mannames given: 
    //     read _names, _parss, & _files from _file given in manfiles
    if(manfiles == 0 || manfiles[0] == 0) {
      CLEANUP;
      falcON_THROW("Manipulator: neither names nor files given:"
		   " cannot initialize\n");
    }
    DebugInfo(3,"Manipulator: initializing from file \"%s\"\n",manfiles);
    std::ifstream in(manfiles);
    if(! in.is_open() ) {
      CLEANUP;
      falcON_THROW("Manipulator: cannot open file \"%s\"\n",manfiles);
    }
    N = 0;
    char* n=_name[N]=_names; *_names=0;
    char* p=_pars[N]=_parss; *_parss=0;
    char* f=_file[N]=_files; *_files=0;
    while(! in.eof() && N != NMAX) {
      char  line[1024], *l=line;
      in.getline(line,1023);                           // read line
      while(*l &&  isspace(*l)) l++;                   // skip space @ start
      if(*l ==  0 ) continue;                          // skip empty line
      if(*l == '#') continue;                          // skip line if '#'
      while(*l && !isspace(*l)) *n++ = *l++;           // copy _name
      if(_name[N][0]==0) continue;                     // no _name? next line!
      *n++ = 0;                                        // close _name
      if(*l && *l=='#') continue;                      // skip line if '#'
      while(*l &&  isspace(*l)) l++;                   // skip space after _name
      if(*l && *l=='#') continue;                      // skip line if '#'
      if(!isalpha(*l) && *l !='/')
	while(*l && !isspace(*l)) *p++ = *l++;         // copy _pars
      *p++ = 0;                                        // close _pars
      if(*l && *l=='#') continue;                      // skip line if '#'
      while(*l &&  isspace(*l)) l++;                   // skip space after _pars
      if(*l && *l=='#') continue;                      // skip line if '#'
      while(*l && !isspace(*l)) *f++ = *l++;           // copy _file
      *f++ = 0;                                        // close _file
      N++;                                             // increment N
      _name[N] = n;                                    // make ready for next
      _pars[N] = p;                                    //   line
      _file[N] = f;
    }
    if(N == 0) {
      CLEANUP;
      falcON_THROW("Manipulator: no manipulator name found in file \"%s\"\n",
		   manfiles);
    }
    if(N > NMAX) {
      CLEANUP;
      falcON_THROW("Manipulator: too many manipulators (%d > NMAX=%d)\n",
		   N,NMAX);
    }
  } else {
    // 1.2 mannames given: assume list of manipulators
    DebugInfo(3,"Manipulator: initializing from\n"
	      "  names=\"%s\"\n  parss=\"%s\"\n  files=\"%s\"\n",
	      mannames,manparss,manfiles);
    // 1.2.1 split mannames and count number of manipulators
    if(strlen(mannames) > Lnames) {
      CLEANUP;
      falcON_THROW("Manipulator: names too long (>%d)\n",Lnames);
    }
    strcpy(_names,mannames);
    N = splitstring<NMAX>(_names,_name,NameSeps);
    if(N > NMAX) {
      CLEANUP;
      falcON_THROW("Manipulator: too many manipulators (%d > NMAX=%d)\n",
		   N,NMAX);
    }
    // 1.2.2 split manparss, allow for empty manparss -> all _pars[]=0
    if(manparss) {
      if(strlen(manparss) > Lparss) {
	CLEANUP;
	falcON_THROW("Manipulator: manpars too long (>%d)\n",Lparss);
      }
      strcpy(_parss,manparss);
      int n = splitstring<NMAX>(_parss,_pars,ParsSeps);
      if(N != n) {
	CLEANUP;
	falcON_THROW("Manipulator: %d names but %d parss",N,n);
      }
    }
    // 1.2.3 split manfiles, allow for empty manfiles -> all _file[]=0
    if(manfiles) {
      if(strlen(manfiles) > Lfiles) {
	CLEANUP;
	falcON_THROW("Manipulator: manfile too long (>%d)\n",Lfiles);
      }
      strcpy(_files,manfiles);
      int n = splitstring<NMAX>(_files,_file,FileSeps);
      if(N != n) {
	CLEANUP;
	falcON_THROW("Manipulator: %d names but %d files",N,n);
      }
    }
  }
  if(debug(2)) {
    std::cerr<<"Manipulator: parsed "<<N<<" manipulators:\n";
    for(int n=0; n!=N; ++n) {
      std::cerr<<" \""<<_name[n]<<"\",";
      if(_pars[n]) std::cerr<<" pars=\""<<_pars[n]<<"\",";
      else         std::cerr<<" no pars,";
      if(_file[n]) std::cerr<<" file=\""<<_file[n]<<"\"\n";
      else         std::cerr<<" no file\n";
    }
  }
  // 2 get manipulators and set data fields
  size_t namesize = 0;
  size_t dscrsize = 0;
  for(int i=0; i!=N; ++i) {
    try {
      single_manipulator(MANIP[i], _name[i], _pars[i], _file[i], manpaths);
    } catch (falcON::exception E) {
      CLEANUP;
      falcON_RETHROW(E);
    }
    namesize += strlen(MANIP[i]->name());
    dscrsize += strlen(MANIP[i]->describe());
    NEED     |= MANIP[i]->need() & ~PRVD;
    CHNG     |= MANIP[i]->change();
    PRVD     |= MANIP[i]->provide();
    IS_MPI    = IS_MPI && MANIP[i]->is_mpi();
  }
  CLEANUP;
  NAME = falcON_NEW(char,2+namesize+N);
  DSCR = falcON_NEW(char,2+dscrsize+3*N);
  strcpy(NAME,MANIP[0]->name());
  strcpy(DSCR,MANIP[0]->describe());
  for(int i=1; i!=N; ++i) {
    strcat(NAME,"+");
    strcat(NAME,MANIP[i]->name());
    strcat(DSCR,"\n# ");
    strcat(DSCR,MANIP[i]->describe());
  }
}
////////////////////////////////////////////////////////////////////////////////
