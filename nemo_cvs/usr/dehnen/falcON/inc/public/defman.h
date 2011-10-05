// -*- C++ -*-                                                                 |
//-----------------------------------------------------------------------------+
//                                                                             |
// defman.h                                                                    |
//                                                                             |
// Copyright (C) 2004-2008  Walter Dehnen                                      |
//                                                                             |
// This program is free software; you can redistribute it and/or modify        |
// it under the terms of the GNU General Public License as published by        |
// the Free Software Foundation; either version 2 of the License, or (at       |
// your option) any later version.                                             |
//                                                                             |
// This program is distributed in the hope that it will be useful, but         |
// WITHOUT ANY WARRANTY; without even the implied warranty of                  |
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU           |
// General Public License for more details.                                    |
//                                                                             |
// You should have received a copy of the GNU General Public License           |
// along with this program; if not, write to the Free Software                 |
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                   |
//                                                                             |
//-----------------------------------------------------------------------------+
#ifndef __cplusplus
# error defman.h is a C++ header
# error to be included by C++ source code only. 
#endif

#ifndef falcON_included_defman_h
#define falcON_included_defman_h

#ifndef falcON_included_manip_h
#  include <public/manip.h>                        // definition of manipulator 
#endif
//------------------------------------------------------------------------------
//                                                                              
// define inimanip() as extern "C" (so that it may be loaded at run time)       
//                                                                              
//------------------------------------------------------------------------------
extern "C" {
#ifdef MANIP_PARSE_AT_INIMANIP
  void inimanip(const falcON::manipulator**,       // O: manipulator            
		const char                *,       // I: parameter set          
		const char                *);      // I: data file              
#else
  void inimanip(const falcON::manipulator**,       // O: manipulator            
		const double              *,       // I: parameter              
		int                        ,       // I: # parameters           
		const char                *);      // I: data file              
#endif
}
//------------------------------------------------------------------------------
//                                                                              
// define macros that define inimanip                                           
//                                                                              
// To use them, you must provide a non-abstract class NAME derived from class   
// manipulator and with constructor of the form                                 
//                                                                              
//    NAME::NAME(const double*pars,       // I: parameters of manipulator       
//               int          npar,       // I: # parameters                    
//               const char  *file);      // I: data file                       
//                                                                              
//------------------------------------------------------------------------------
#define __DEF__MAN(NAME)					\
void inimanip(const falcON::manipulator**manip,			\
	      const char                *pars,			\
	      const char                *file)			\
{								\
  const int MAXPAR = 256;					\
  double p[MAXPAR];						\
  int    n = falcON::Manipulator::parse(pars,p,MAXPAR);		\
  *manip = new NAME(p,n,file);					\
  falcON::CheckAgainstLibrary(falcON::CurrentStatus(),		\
  falcON::message("Manipulator \"%s\"",(*manip)->name()));	\
}

#define __DEF__MAN__ALT(NAME)					\
void inimanip(const falcON::manipulator**manip,			\
	      const char                *pars,			\
	      const char                *file)			\
{								\
  *manip = new NAME(pars,file);					\
  falcON::CheckAgainstLibrary(falcON::CurrentStatus(),		\
  falcON::message("Manipulator \"%s\"",(*manip)->name()));	\
}
//------------------------------------------------------------------------------
#endif // falcON_included_defman_h
