/*
   Routines for obtaining various information about the current
   process (header file).

   Copyright (C) 2004, 2009 Oleksandr Yermolenko
   <oleksandr.yermolenko@gmail.com>

   This file is part of MDTK, the Molecular Dynamics Toolkit.

   MDTK is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   MDTK is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MDTK.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef mdtk_PROCMON_h
#define mdtk_PROCMON_h

#include <iostream>

#ifdef __WIN32__
#include <windows.h>

namespace procmon
{
  FILETIME   getCreationTimeFT(HANDLE);
  FILETIME   getExitTimeFT(HANDLE);
  FILETIME   getKernelTimeFT(HANDLE);
  FILETIME   getUserTimeFT(HANDLE);

  SYSTEMTIME getCreationTimeST(HANDLE);
  SYSTEMTIME getExitTimeST(HANDLE);
  SYSTEMTIME getKernelTimeST(HANDLE);
  SYSTEMTIME getUserTimeST(HANDLE);
  
  void       loadFILETIME(std::istream&,FILETIME&);
  void       saveFILETIME(std::ostream&,const FILETIME&);
  void       printSYSTEMTIME(std::ostream&,const SYSTEMTIME&);
  __int64    FILETIME2SECONDS(const FILETIME&);
  __int64    FILETIME2INT64(const FILETIME&);
  
  class ProcmonTimer
  {
    __int64 startTime;
    __int64 deltaPrevTime;
  public:
    ProcmonTimer()
    :startTime(FILETIME2INT64(getUserTimeFT(GetCurrentProcess()))),
     deltaPrevTime(FILETIME2INT64(getUserTimeFT(GetCurrentProcess())))
    {
    }
    double getTimeInSeconds()
    {
      __int64 currentTime = FILETIME2INT64(getUserTimeFT(GetCurrentProcess()));
      return (currentTime - startTime)/(10000000.0);
    }
    __int64 getTime()
    {
      __int64 currentTime = FILETIME2INT64(getUserTimeFT(GetCurrentProcess()));
      return currentTime - startTime;
    }
    __int64 getDeltaTime()
    {
      __int64 tmp = deltaPrevTime;
      deltaPrevTime = FILETIME2INT64(getUserTimeFT(GetCurrentProcess()));
      return deltaPrevTime - tmp;
    }  
    double getDeltaTimeInSeconds()
    {
      __int64 tmp = deltaPrevTime;
      deltaPrevTime = FILETIME2INT64(getUserTimeFT(GetCurrentProcess()));
      return (deltaPrevTime - tmp)/(10000000.0);
    }  
  };  

void
unUseCPU();
}

#else
namespace procmon
{

typedef long unsigned int	 PMTime;

const double PM_JIFFY = 100.0;

  PMTime   getKernelTime();
  PMTime   getUserTime();

  class ProcmonTimer
  {
    PMTime startTime;
    PMTime deltaPrevTime;
  public:
    ProcmonTimer()
    :startTime(getUserTime()),
     deltaPrevTime(getUserTime())
    {
    }
    double getTimeInSeconds()
    {
      PMTime currentTime = getUserTime();
      return (currentTime - startTime)/(PM_JIFFY);
    }
    PMTime getTime()
    {
      PMTime currentTime = getUserTime();
      return currentTime - startTime;
    }
    PMTime getDeltaTime()
    {
      PMTime tmp = deltaPrevTime;
      deltaPrevTime = getUserTime();
      return deltaPrevTime - tmp;
    }  
    double getDeltaTimeInSeconds()
    {
      PMTime tmp = deltaPrevTime;
      deltaPrevTime = getUserTime();
      return (deltaPrevTime - tmp)/(PM_JIFFY);
    }  
  };  

inline
void
unUseCPU(){}

}

#endif


#endif

  
