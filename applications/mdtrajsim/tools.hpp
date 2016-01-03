/* 
   Some useful functions for mdtrajsim (header file)

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

#ifndef mdtk_apps_mdtrajsim_tools_hpp
#define mdtk_apps_mdtrajsim_tools_hpp

#include <string>
#include <sstream>

#ifdef MPIBATCH
int comm_rank();
int comm_size();
std::string comm_name();
#endif

std::string getProcessID();

bool isLockedByOthers();
void placeMyLock();
void removeMyLock();

void shrinkLogFiles();
void removeXVAfiles();

extern std::string myLockFilename;

void SleepForSeconds(int seconds);

void setupSignalHandlers();

#endif