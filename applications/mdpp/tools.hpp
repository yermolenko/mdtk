/* 
   tools.hpp (molecular dynamics postprocessor, tools)

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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
  return strstr(trajDirName,"Cu13_in_C60_on_Cu_") && !strstr(trajDirName,"1000eV");
}

inline
bool 
trajProcess_Custom(const char* trajDirName)
{
  return 
    strstr(trajDirName,"C60_on_Cu_")
    && !strstr(trajDirName,"1000eV");
}

inline
bool 
trajProcess_Custom1(const char* trajDirName)
{
  return 
    strstr(trajDirName,"Cu_by_Cu");
}

inline
bool 
trajProcess_Custom2(const char* trajDirName)
{
  return 
    strstr(trajDirName,"000");
}

inline
bool 
isAlreadyFinished(const char* trajDirName)
{
  {
    std::string finalStateFile 
      = std::string(trajDirName)+DIR_DELIMIT_STR+"mde_final";
    if (yaatk::exists(finalStateFile.c_str())) return true;
  }
  {
    std::string finalStateFile 
      = std::string(trajDirName)+DIR_DELIMIT_STR+"completed.ok";
    if (yaatk::exists(finalStateFile.c_str())) return true;
  }
  return false;
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
        if (entry->d_name[0] == '0')
        {
          std::sprintf(trajdir_src,"%s%s",trajsetDir,entry->d_name);
          std::sprintf(stateFileName,"%s" DIR_DELIMIT_STR,trajdir_src);
          if (fpt(stateFileName) && isAlreadyFinished(stateFileName))
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

inline
void
removeDuplicates(std::vector<std::string>& states)
{
  size_t i;
  std::vector<std::string> states_new;
  
  if (states.size() >= 1) states_new.push_back(states[0]);
  for(i = 1; i < states.size(); i++)
  {
    char it_prev_s[100];
    char it_s[100];
    strcpy(it_prev_s,states[i-1].substr(3,10).c_str());
    strcpy(it_s,     states[i].  substr(3,10).c_str());
    int it_prev; sscanf(it_prev_s,"%d",&it_prev);
    int it;      sscanf(it_s,     "%d",&it);
    if (it-it_prev>5)
      states_new.push_back(states[i]);
  }

  states = states_new;
}

inline
void
findIntermediateStates(std::string trajDir,std::vector<std::string>& states)
{
  {
    {
      DIR* trajsetDirHandle = opendir(trajDir.c_str());
//      REQUIRE(trajsetDirHandle != NULL);


      struct dirent* entry = readdir(trajsetDirHandle);
      while (entry != NULL)
      {
        if (entry->d_type == DT_REG)
        {
          if (entry->d_name[0] == 'm' && entry->d_name[1] == 'd' && entry->d_name[2] == 'e' &&
              entry->d_name[3] == '0')
          {
            states.push_back(entry->d_name);
//            TRACE(entry->d_name);
          }
/*
          std::sprintf(trajdir_src,"%s%s",trajsetDir,entry->d_name);
          std::sprintf(stateFileName,"%s"DIR_DELIMIT_STR,trajdir_src);          stateFileNames.push_back(stateFileName);
          TRACE(stateFileName);
*/
        }
        entry = readdir(trajsetDirHandle);
      };

      /*int res_closedir = */closedir(trajsetDirHandle);
//      REQUIRE(res_closedir != NULL);
    }  
  }  
  sort(states.begin(),states.end());
  removeDuplicates(states);
}

struct _SavedStateSortStruct
{
  std::string fullTrajDirName;
  std::string shortTrajDirName;
  _SavedStateSortStruct(std::string stateFileName)
   :fullTrajDirName(stateFileName), shortTrajDirName(yaatk::extractLastItem(yaatk::extractDir(stateFileName)))
  {
  }
  friend int operator<(const _SavedStateSortStruct& left, const _SavedStateSortStruct& right);
  friend int operator<(_SavedStateSortStruct& left, _SavedStateSortStruct& right);
};
    
inline
int operator<(const _SavedStateSortStruct& left, const _SavedStateSortStruct& right)
{
  return left.shortTrajDirName < right.shortTrajDirName;
}
        
inline
int operator<(_SavedStateSortStruct& left, _SavedStateSortStruct& right)
{
  return left.shortTrajDirName < right.shortTrajDirName;
}

/*
inline
mdtk::Float getMassInAMU(const std::vector<mdtk::Atom>& atoms)
{
  mdtk::Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    moleculeMass += atom.M;
  } 
  return mdtk::academic_round(moleculeMass/mdtk::amu);
}

inline
mdtk::Vector3D getVelocity(const std::vector<mdtk::Atom>& atoms)
{
  using mdtk::Exception;    
    
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  mdtk::Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.V*atom.M;
  };
  return sumOfP/sumOfM;    
}  
*/

/*
inline
void printGlobalIndexes(const std::vector<mdtk::Atom>& atoms,
			std::ostream& fo)
{
  fo << atoms.size() << std::endl;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    fo << atom.globalIndex << std::endl;
  } 
} 
*/
  
}

#endif
