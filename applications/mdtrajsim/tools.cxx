/* 
   Some useful functions for mdtrajsim

   Copyright (C) 2015 Oleksandr Yermolenko
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

#ifdef __WIN32__
#include <windows.h>
#include <WinBase.h>
#else
#include <sys/types.h>
#include <unistd.h>
#include <sys/utsname.h>
#endif

#include <vector>
#include <fstream>

#ifdef MPIBATCH
#include "mpi.h"
#endif

#include <yaatk/yaatk.hpp>

#include "tools.hpp"

using namespace yaatk;

#ifdef MPIBATCH
int
comm_rank()
{
  int comm_rank;
  MPI_TEST_SUCCESS(MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank));
  return comm_rank;
}

int
comm_size()
{
  int comm_size;
  MPI_TEST_SUCCESS(MPI_Comm_size(MPI_COMM_WORLD,&comm_size));
  return comm_size;
}

std::string
comm_name()
{
  char comm_name[MPI_MAX_PROCESSOR_NAME+1];
  int len;
  MPI_TEST_SUCCESS(MPI_Get_processor_name(comm_name,&len));
  return comm_name;
}
#endif

std::string
getProcessID()
{
  std::ostringstream ossid;
#ifdef MPIBATCH
  ossid << "-" << comm_rank();
#endif
#ifdef __WIN32__
  ossid << "-" << GetCurrentProcessId();
#else
  ossid << "-" << getpid();
  struct utsname buffer;
  if (uname(&buffer) == 0)
    ossid << "-" << buffer.nodename;
#endif
  return ossid.str();
}

const std::string lockFilenameBase = "mdtrajsim.lock";
std::string myLockFilename = lockFilenameBase;

bool
isLockedByOthers()
{
  std::vector<std::string> cwdfiles
    = yaatk::listFiles(yaatk::getcwd());
  for(size_t i = 0; i < cwdfiles.size(); ++i)
  {
    const std::string& cwdfile = cwdfiles[i];
    if (cwdfile.find(lockFilenameBase) != std::string::npos &&
        cwdfile != myLockFilename)
      return true;
  }

  return false;
}

void
placeMyLock()
{
  std::ifstream test(myLockFilename.c_str());
  if (!test)
  {
    std::ofstream folock(myLockFilename.c_str());
    folock.close();
  }
  else
    throw Exception("Lockfile aready exists");
}

void
removeMyLock()
{
  std::ifstream test(myLockFilename.c_str());
  if (test)
  {
    test.close();
    int retval = yaatk::remove(myLockFilename);
    if (retval)
      throw Exception("Cannot remove lockfile");
  }
}

void
shrinkLogFiles()
{
  int retcode;

  yaatk::DataState ds;

  const std::string dElogfilename="dE.dat";
  if (yaatk::exists(dElogfilename))
    return;

  const std::string logfilename="stdout.txt";
  {
    size_t numberOfLogfiles = 0;
    std::vector<std::string> cwdfiles
      = yaatk::listFiles(yaatk::getcwd());
    for(size_t i = 0; i < cwdfiles.size(); ++i)
    {
      const std::string& cwdfile = cwdfiles[i];
      if (cwdfile.find(logfilename) != std::string::npos)
      {
        ++numberOfLogfiles;
        if (numberOfLogfiles > 1)
          return;
      }
    }
  }

  yaatk::text_ofstream stdout_txt_filtered(logfilename + "-filtered");
  REQUIRE(stdout_txt_filtered.isOpened());
  {
    yaatk::text_ofstream dE_dat(dElogfilename);
    REQUIRE(dE_dat.isOpened());
    {
      std::ifstream stdout_txt(logfilename.c_str());
      REQUIRE(stdout_txt);
      std::vector<std::string> stdout_txt_lines(100);
      size_t stdout_txt_lines_index = 0;
      size_t stdout_txt_lines_read = 0;

      char cline[1000];
      while(stdout_txt.getline(cline, 1000-1, '\n'))
      {
        std::string line = cline;

        stdout_txt_lines[stdout_txt_lines_index++] = line;
        stdout_txt_lines_index %= stdout_txt_lines.size();

        stdout_txt_lines_read++;
        if (stdout_txt_lines_read <= stdout_txt_lines.size())
          stdout_txt_filtered << line << "\n";

        if (line.find("Eo+") != std::string::npos)
        {
          const std::string description = "dE/(Eo+Eb) : ";
          if (line.find(description) == 0)
            line.replace(0, description.size(), "");
          dE_dat << line << "\n";
        }
      }
      stdout_txt.close();

      size_t tail_lines_count = stdout_txt_lines.size();
      if (stdout_txt_lines_read <= stdout_txt_lines.size()*2)
      {
        if (tail_lines_count > stdout_txt_lines.size()*2 - stdout_txt_lines_read)
          tail_lines_count -= stdout_txt_lines.size()*2 - stdout_txt_lines_read;
        else
          tail_lines_count = 0;
      }
      else
      {
        stdout_txt_filtered << "~~~~~~~~~~~~~~ cut ~~~~~~~~~~~~~~\n";
      }

      if (tail_lines_count > stdout_txt_lines_index)
        stdout_txt_lines_index += stdout_txt_lines.size();
      stdout_txt_lines_index -= tail_lines_count;

      while (tail_lines_count > 0)
      {
        stdout_txt_filtered << stdout_txt_lines[stdout_txt_lines_index++] << "\n";
        stdout_txt_lines_index %= stdout_txt_lines.size();

        --tail_lines_count;
      }
    }
    dE_dat.close();
  }
  stdout_txt_filtered.close();

  retcode = yaatk::remove(logfilename);
  REQUIRE(retcode == 0);
  retcode = yaatk::rename(stdout_txt_filtered.getZippedFileName(),
                          logfilename + stdout_txt_filtered.getZippedExt());
  REQUIRE(retcode == 0);
}

void
removeXVAfiles()
{
  std::vector<std::string> cwdfiles
    = yaatk::listFiles(yaatk::getcwd());
  for(size_t i = 0; i < cwdfiles.size(); ++i)
  {
    const std::string& cwdfile = cwdfiles[i];
    if (cwdfile.find("xva") != std::string::npos &&
        cwdfile.find("mde") == 0)
    {
      int retval = yaatk::remove(cwdfile);
      if (retval)
        throw Exception("Cannot remove XVA file");
    }
  }
}

void SleepForSeconds(int seconds)
{
#ifdef __WIN32__
  Sleep(seconds*1000);
#else
  sleep(seconds);
#endif
}
