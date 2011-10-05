// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/exception.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-2010
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2010  Walter Dehnen
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
#include <exception.h>
#include <cstdlib>
#include <cerrno>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#ifdef _OPENMP
#  include <omp.h>
#endif
extern "C" {
#if defined(unix) || defined(__DARWIN_UNIX03)
#  include <unistd.h>
#  include <sys/time.h>
#endif
#ifdef WIN32
#  include <windows.h>
#endif
}

#ifdef __INTEL_COMPILER
#pragma warning (disable:981) /* operands are evaluated in unspecified order */
#endif

//                                                                              
// class RunInfo                                                                
//                                                                              
WDutils::RunInfo::RunInfo()
  : __host_known(0),
    __user_known(0),
    __pid_known(0),
    __name_known(0),
    __is_mpi_proc(0),
    __debug(0)
{
  try {
    // set wall-clock time
    {
#if defined(unix) || defined(__DARWIN_UNIX03)
      timeval now;
      gettimeofday(&now, NULL);
      __sec = now.tv_sec;
      __usec = now.tv_usec;
#elif defined(WIN32)
      LARGE_INTEGER tmp;
      QueryPerformanceCounter(&tmp);
      QueryPerformanceFrequency(&tmp);
      __timecount = tmp.QuadPart;
      __timetick = 1.0/double(tmp.QuadPart);
#endif
    }
    // set run time
    {
      time_t now = ::time(0);
      SNprintf(__time,100,ctime(&now));
      __time[24] = 0;
    }
#if defined(unix) || defined(__DARWIN_UNIX03)
    // set host name
    {
      gethostname(__host,100);
      __host_known = 1;
    }
    // set user name
    {
      const char*user__ = getenv("USER");
      if(user__) {
	SNprintf(__user,100,user__);
	__user_known = 1;
      } else
	SNprintf(__user,100,"unknown.user");
    }
    // set pid
    {
      __pid_num = getpid();
      SNprintf(__pid,20,"%d",__pid_num);
      __pid_known  = 1;
    }
    // set command and name of executable
    {
      char file[64];
      if(__pid_known) {
	SNprintf(file,64,"/proc/%s/cmdline",__pid);
	std::ifstream in(file);
	if(in) {
	  int i,e=0;
	  for(i=0; i!=1024; ++i) __cmd[i]=0;
	  in.getline(__cmd,1023);
	  for(i=1023; i!=0; --i)
	    if(__cmd[i]==0 || isspace(__cmd[i])) __cmd[i] = ' ';
	    else if(e==0) e=i;
	  __cmd[e+1] = 0;
	  for(i=0; !isspace(__cmd[i]); ++i)
	    __name[i] = __cmd[i];
	  __name[i] = 0;
	  __name[i] = 0;
	  __cmd_known  = 1;
	  __name_known = 1;
	}
      }
    }
#else // unix
    {
      SNprintf(__host,100,"unknown.host");
      SNprintf(__user,100,"unknown.user");
      SNprintf(__name,100,"unknown.name");
    }
#endif
    // set # proc available for openMP
    {
#ifdef _OPENMP
      if(omp_in_parallel())
	WDutils_Error("RunInfo::RunInfo() called inside OMP parallel region\n");
      __omp_proc = omp_get_num_procs();
#else
      __omp_proc = 1;
#endif
      __omp_size = __omp_proc;
    }
  } catch(exception E) {
    WDutils_RETHROW(E);
  }
}
//
void WDutils::RunInfo::set_omp(const char*arg) WDutils_THROWING
{
#ifdef _OPENMP
  if(arg==0 || arg[0]==0 || arg[0]=='t')
    Info.__omp_size = Info.__omp_proc;
  else if(arg[0] == 'f')
    Info.__omp_size = 1;
  else if(arg && arg[0]) {
    Info.__omp_size = strtol(arg,0,10);
    if(errno == EINVAL || errno == ERANGE)
      WDutils_THROW("RunInfo::set_omp('%s')\n",arg);
    if(Info.__omp_size < 1) {
      Info.__omp_size = 1;
      WDutils_WarningN("RunInfo::set_omp('%s') assume '1'\n",arg);
    }
  }
  omp_set_num_threads(Info.__omp_size);
#else
  Info.__omp_size = 1;
#endif
}
//
void WDutils::RunInfo::header(std::ostream&out)
{
  if(out) {
    if(RunInfo::cmd_known())
      out<<"# \""<<RunInfo::cmd()<<"\"\n#\n";
    out<<"# run at  "  <<RunInfo::time()<<"\n";
    if(RunInfo::user_known())  out<<"#     by  \""<<RunInfo::user()<<"\"\n";
    if(RunInfo::host_known())  out<<"#     on  \""<<RunInfo::host()<<"\"\n";
    if(RunInfo::pid_known())   out<<"#     pid  " <<RunInfo::pid() <<"\n";
    if(RunInfo::is_mpi_proc()) out<<"#     mpi  " <<RunInfo::mpi_size()<<"\n";
    out<<"#\n";
  }
}
//
#if defined(unix) || defined(__DARWIN_UNIX03)
double WDutils::RunInfo::WallClock()
{
  timeval now;
  gettimeofday(&now, NULL);
  return (now.tv_sec - Info.__sec) + (now.tv_usec - Info.__usec) * 1.e-6;
}
void WDutils::RunInfo::WallClock(unsigned&sec, unsigned&usec)
{
  timeval now;
  gettimeofday(&now, NULL);
  if(now.tv_usec > Info.__usec) {
    sec  = now.tv_sec  - Info.__sec;
    usec = now.tv_usec - Info.__usec;
  } else {
    sec  = now.tv_sec  - Info.__sec - 1;
    usec = (1000000 + now.tv_usec) - Info.__usec;
  }
}
#elif defined(WIN32)
double WDutils::RunInfo::WallClock()
{
  LARGE_INTEGER now;
  QueryPerformanceCounter(&now);
  return (now.QuadPart - Info.__timecount) * Info.__timetick;
}
#endif
//
WDutils::RunInfo WDutils::RunInfo::Info;

//

namespace {
  using namespace WDutils;

#if defined(__PGI) || defined(__SUNPRO_CC)
  using ::snprintf;
  using ::vsnprintf;
#else
  using std::snprintf;
  using std::vsnprintf;
#endif

  inline void printerr(const char*lib,
		       const char*issue,
		       const char*fmt,
		       va_list   &ap,
		       int        depth= 0,
		       const char*file = 0,
		       int        line = 0,
		       bool       name = true)
  {
    const int size=1024;
    int s=size, w=0;
    char ffmt[size], *t=ffmt;
    char dpth[24] = ":                     ";
    if(depth>20) depth=20;
    dpth[2+depth]=0;
    w=snprintf(t,s,"### %s %s",lib,issue);
    t+=w; s-=w;
    if(name && RunInfo::name_known()) {
      w=snprintf(t,s," [%s]",RunInfo::name());
      t+=w; s-=w;
    }
    if(RunInfo::is_mpi_proc()) {
      w=snprintf(t,s," @%2d",RunInfo::mpi_proc());
      t+=w; s-=w;
#ifdef _OPENMP
    } else if(omp_in_parallel()) {
      w=snprintf(t,s," @%2d",omp_get_thread_num());
      t+=w; s-=w;
#endif
    }
    if(file) {
      w=snprintf(t,s," [%s:%d]",file,line);
      t+=w; s-=w;
    }
    if (fmt[strlen(fmt)-1] != '\n')
      w=snprintf(t,s,"%s%s\n",dpth,fmt);
    else
      w=snprintf(t,s,"%s%s",dpth,fmt);
    t+=w; s-=w;
    vfprintf(stderr,ffmt,ap);
    fflush(stderr);
  }
}
//
void WDutils::Error::operator()(const char* fmt, ...) const
{
  va_list  ap;
  va_start(ap,fmt);
  printerr(lib, "Error", fmt, ap, 0, file, line);
  va_end(ap);
  std::exit(1);
}
//
void WDutils::Warning::operator()(const char* fmt, ...) const
{
  va_list  ap;
  va_start(ap,fmt);
  printerr(lib, "Warning", fmt, ap, 0, file, line);
  va_end(ap);
}
//
void WDutils::DebugInformation::operator()(const char* fmt, ...) const
{
  va_list  ap;
  va_start(ap,fmt);
  printerr(lib, "Debug Info", fmt, ap, 0, file, line, false);
  va_end(ap);
}
//
void WDutils::DebugInformation::operator()(int deb, const char* fmt, ...) const
{
  if(RunInfo::debug(deb)) {
    va_list  ap;
    va_start(ap,fmt);
    printerr(lib, "Debug Info", fmt, ap, deb, file, line, false);
    va_end(ap);
  }
}
//
WDutils::exception WDutils::Thrower::operator()(const char*fmt, ...) const
{
  size_t size = 1024;
  char   buffer[1024], *buf=buffer;
  if(file) {
    int len = SNprintf(buf,size,"in [%s:%d]: ",file,line);
    buf  += len;
    size -= len;
  }
  va_list  ap;
  va_start(ap,fmt);
  vsnprintf(buf, size, fmt, ap);
  va_end(ap);
  return exception(buffer);
}
//
WDutils::exception::exception(const char*fmt, ...)
{
  const int msize=1024;
  char __text[msize];
  va_list  ap;
  va_start(ap,fmt);
  int w = vsnprintf(__text,msize,fmt,ap);
  if(w>=msize) {
    WDutils_Warning("WDutils::exception::exception(): "
		    "string size of %d characters exceeded\n",msize);
    __text[msize-1]=0;
  }
  if(w<0)
    WDutils_Warning("WDutils::exception::exception(): formatting error\n");
  va_end(ap);
  std::string::operator= (__text);
}
//
WDutils::message::message(const char*fmt, ...) throw(exception)
{
  va_list  ap;
  va_start(ap,fmt);
  int w = vsnprintf(__text,size,fmt,ap);
  if(w>=static_cast<int>(size))
    throw exception("WDutils::message::message(): "
		    "string size of %d characters exceeded\n",size);
  if(w < 0)
    throw exception("WDutils::message::message(): "
		    "formatting error\n");
  va_end(ap);
}
//
int WDutils::snprintf(char*str, size_t l, const char* fmt, ...)
  WDutils_THROWING
{
  va_list ap;
  va_start(ap,fmt);
  int w = vsnprintf(str,l,fmt,ap);
  va_end(ap);
  if(w==static_cast<int>(l))
    WDutils_THROW("snprintf(): trailing 0 lost");
  if(w >static_cast<int>(l))
    WDutils_THROW("snprintf(): string size exceeded [%d:%lu]",w,l);
  if(w <0)
    WDutils_THROW("snprintf(): formatting error");
  return w;
}
//
int WDutils::snprintf__::operator()(char*str, size_t l, const char* fmt, ...)
  WDutils_THROWING
{
  va_list ap;
  va_start(ap,fmt);
  int w = vsnprintf(str,l,fmt,ap);
  va_end(ap);
  if(w==static_cast<int>(l))
    WDutils_THROW("[%s:%d]: snprintf(): trailing 0 lost",file,line);
  if(w >static_cast<int>(l))
    WDutils_THROW("[%s:%d]: snprintf(): string size exceeded [%d:%lu]",
		  file,line,w,l);
  if(w <0)
    WDutils_THROW("[%s:%d]: snprintf(): formatting error",file,line);
  return w;
}
////////////////////////////////////////////////////////////////////////////////
