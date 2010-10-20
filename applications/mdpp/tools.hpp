/* 
   tools.hpp (molecular dynamics postprocessor, tools)

   Copyright (C) 2007, 2008, 2009, 2010 Oleksandr Yermolenko
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

#ifndef mdtk_apps_mdpp_tools_hpp
#define mdtk_apps_mdpp_tools_hpp


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>

#include <algorithm>

#include "../common.h"

#include <dirent.h>

namespace mdepp
{

typedef bool (*FProcessTrajectory)(const char* trajDirName);

inline
bool 
trajProcessAll(const char* trajDirName)
{
  return true;
}

inline
bool 
trajProcess_Cu13_at_C60_Only(const char* trajDirName)
{
  return strstr(trajDirName,"Cu13_in_C60_on_Cu_");
}

inline
void
addTrajDirNames(std::vector<std::string> &stateFileNames,const char *trajsetDir_, FProcessTrajectory fpt)
{
  std::cout << "Adding states from " << trajsetDir_ << std::endl;

  char trajsetDir[10000];
  strcpy(trajsetDir,trajsetDir_);

  if (trajsetDir[strlen(trajsetDir)-1] != DIR_DELIMIT_CHAR)
    strcat(trajsetDir,DIR_DELIMIT_STR);
  
  {
    {
      char trajdir_src[10000];
      char stateFileName[10000];

      DIR* trajsetDirHandle = opendir(trajsetDir);
      REQUIRE(trajsetDirHandle != NULL);


      struct dirent* entry = readdir(trajsetDirHandle);
      while (entry != NULL)
      {
        if (entry->d_type == DT_DIR && strcmp(entry->d_name,".") && strcmp(entry->d_name,".."))
        if (entry->d_name[0] == 'C')
        {
          std::sprintf(trajdir_src,"%s%s",trajsetDir,entry->d_name);
          std::sprintf(stateFileName,"%s"DIR_DELIMIT_STR,trajdir_src);
          if (fpt(stateFileName))
	    stateFileNames.push_back(stateFileName);
        }
        entry = readdir(trajsetDirHandle);
      };

      int res_closedir = closedir(trajsetDirHandle);
      REQUIRE(res_closedir == 0);
    }  
  }  
}  

inline
void
findTrajDirs(const char* trajsetFileName, const char* stdTrajsetDir,
	     std::vector<std::string>& trajDirNames, mdepp::FProcessTrajectory fpt = mdepp::trajProcessAll)
{
  char trajsetDir[10000];

  {
    using mdtk::Exception;
    std::ifstream fi(trajsetFileName);
    if (fi)
    {
      while(fi.getline(trajsetDir, 1000-1, '\n'))
      {
        if (strlen(trajsetDir) > 0)
          mdepp::addTrajDirNames(trajDirNames,trajsetDir,fpt);
      };
      fi.close();
    }
    else
    {
      sprintf(trajsetDir,"%s",stdTrajsetDir);
      mdepp::addTrajDirNames(trajDirNames,trajsetDir,fpt);
    }
  }
  std::sort(trajDirNames.begin(),trajDirNames.end());
}

  
};

#endif
