/*
   Routines for obtaining various information about the current
   process.

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

#include <iomanip>

#include "procmon.hpp"

#include <stdio.h>

#ifdef __WIN32__

namespace procmon
{
  #define CHECK_ON_SUCCESS(cmd) \
  {\
    if ((cmd)==0)  \
    { \
    std::cerr << "Error requesting process timing info." << std::endl; \
    throw; \
    } \
  }
  
  #define GetProcessTimes_ALL(handle) \
  FILETIME ftCreationTime; \
  FILETIME ftExitTime;     \
  FILETIME ftKernelTime;   \
  FILETIME ftUserTime;     \
 \
  SYSTEMTIME stCreationTime; \
  SYSTEMTIME stExitTime; \
  SYSTEMTIME stKernelTime; \
  SYSTEMTIME stUserTime; \
 \
  CHECK_ON_SUCCESS( \
      GetProcessTimes \
      ( \
        handle, \
        &ftCreationTime, \
        &ftExitTime, \
        &ftKernelTime, \
        &ftUserTime \
      ) \
      ); \
  { \
    CHECK_ON_SUCCESS(FileTimeToSystemTime(&ftCreationTime,&stCreationTime)); \
    CHECK_ON_SUCCESS(FileTimeToSystemTime(&ftExitTime,&stExitTime)); \
    CHECK_ON_SUCCESS(FileTimeToSystemTime(&ftKernelTime,&stKernelTime)); \
    CHECK_ON_SUCCESS(FileTimeToSystemTime(&ftUserTime,&stUserTime)); \
  }

  FILETIME getCreationTimeFT(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return ftCreationTime;
  }  

  FILETIME getExitTimeFT(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return ftExitTime;
  }  

  FILETIME getKernelTimeFT(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return ftKernelTime;
  }  

  FILETIME getUserTimeFT(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return ftUserTime;
  }  

  SYSTEMTIME getCreationTimeST(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return stCreationTime;
  }  

  SYSTEMTIME getExitTimeST(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return stExitTime;
  }  

  SYSTEMTIME getKernelTimeST(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return stKernelTime;
  }  

  SYSTEMTIME getUserTimeST(HANDLE handle)
  {
    GetProcessTimes_ALL(handle);
    return stUserTime;
  }  
  
  void       loadFILETIME(std::istream& is,FILETIME& ftStruct)
  {
    is  >> ftStruct.dwLowDateTime
        >> ftStruct.dwHighDateTime;
  }
    
  void       saveFILETIME(std::ostream& os,const FILETIME& ftStruct)
  {
    os << ftStruct.dwLowDateTime  << "\n"
       << ftStruct.dwHighDateTime << std::endl;
  }  

  void       printSYSTEMTIME(std::ostream& os,const SYSTEMTIME& time)
  {
    os << std::setw(4) << time.wYear  << "."
       << std::setw(2) << time.wMonth << "."
       << std::setw(2) << time.wDay   << " "
       << std::setw(2) << time.wHour   << ":"
       << std::setw(2) << time.wMinute << ":"
       << std::setw(2) << time.wSecond << "\n" ;
  }  

  __int64
  FILETIME2SECONDS(const FILETIME& ft)
  {
    return FILETIME2INT64(ft) / (__int64)(10000000);
  }  

  __int64
  FILETIME2INT64(const FILETIME& ft)
  {
    return 
    (
      (__int64)(ft.dwLowDateTime) | ((__int64)(ft.dwHighDateTime) << (__int64)32)
    );
  }  

#include <windows.h>

void
unUseCPU()
{
  Sleep(1);
}

}

#else
namespace procmon
{

  void
  GetProcessTimes_ALL(PMTime &userTime,PMTime &kernelTime)
  {
    FILE* proc_self_stat = fopen("/proc/self/stat","r");
    if (!proc_self_stat)
    {
    std::cerr << "Error requesting process timing info." << std::endl;
//    throw;
    userTime = 0; kernelTime = 0;
    return;
    }
    int pid;
    char comm[1000];
    char state;
    int  ppid;
    int  pgrp;
    int  session;
    int  tty_nr,tpgid;
    unsigned flags;
    PMTime  minflt,cminflt,majflt,cmajflt,utime,stime;
    fscanf(proc_self_stat,"%d %s %c %d %d %d %d %d %u"
" %lu %lu %lu %lu %lu %lu",
&pid,comm,&state,&ppid,&pgrp,&session,&tty_nr,&tpgid,&flags,
              &minflt,&cminflt,&majflt,&cmajflt,&utime,&stime);

    fclose(proc_self_stat);
    userTime = utime; kernelTime = stime;
  }

  PMTime getUserTime()
  {
    PMTime userTime, kernelTime;
    GetProcessTimes_ALL(userTime,kernelTime);
    return userTime;
  }  

  PMTime getKernelTime()
  {
    PMTime userTime, kernelTime;
    GetProcessTimes_ALL(userTime,kernelTime);
    return kernelTime;
  }  
}
#endif

