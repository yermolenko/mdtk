/* 
   Molecular dynamics postprocessor, main classes

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

#include "StatPostProcess.hpp"

#include <algorithm>

//using namespace std;
#include <fstream>
#include <mdtk/tools.hpp>


namespace mdepp
{

mdtk::AtomsContainer dummy_ac;

void
depthHist2file(const char* filename, std::vector<Float>& depth, Float scale = 1.0)
{
  using mdtk::Exception;    

/*  
    const Float minDepth = -10.0;
    const Float maxDepth =  20.0;
    const int n = (maxDepth-minDepth)/(2.547);// for polyethylene // prev was (0.50);
*/

/*  // Graphite
    const Float minDepth_desired   = -40.0;
    const Float maxDepth_desired   =  40.0;
    const Float matchPoint         =   0.0;
    const Float c = 6.708;
    const Float histStep           =   c/6.0;//2.547; // for polyethylene // prev was (0.50);
*/
/*  // Cu
    const Float c = 3.61;
    const Float minDepth_desired   = -c/2.0*15.123;
    const Float maxDepth_desired   =  c/2.0*15.123;
    const Float matchPoint         = -c/4.0;
    const Float histStep           =   c/2.0;//2.547; // for polyethylene // prev was (0.50);
*/

  // PE
    const Float c = 2.547;
    const Float minDepth_desired   = -100.0;
    const Float maxDepth_desired   =  100.0;
    const Float matchPoint         =    0.0;
    const Float histStep           =   c/3.0;//2.547; // for polyethylene // prev was (0.50);


    REQUIRE(minDepth_desired < matchPoint);
    REQUIRE(maxDepth_desired > matchPoint);
    const int n_below_matchPoint = int( (matchPoint       - minDepth_desired)/histStep ) +1;
    const int n_above_matchPoint = int( (maxDepth_desired - matchPoint      )/histStep ) +1;

    const int n = n_above_matchPoint + n_below_matchPoint;
    const Float minDepth   = matchPoint - n_below_matchPoint*histStep;
    const Float maxDepth   = matchPoint + n_above_matchPoint*histStep;


  {
    std::ofstream fo((std::string(filename)+".plt").c_str());
    fo << "reset\n#set yrange [*:0]\nset yrange [0:*]\nset xrange [-3:25]\nset format x \"%.1f\"\nset xtics " << c/2.0 << "\nset grid xtics\nplot \'" << filename << "\' with boxes\n";
    fo.close();
  }  
  std::ofstream fo(filename);


  fo << "# min depth = " << minDepth << " Ao" << std::endl
     << "# max depth = " << maxDepth << " Ao" << std::endl
     << "# number of bins = " << n << std::endl;
    gsl_histogram * h = gsl_histogram_alloc (n);
    gsl_histogram_set_ranges_uniform (h, minDepth, maxDepth);
    for(size_t i = 0; i < depth.size(); i++)
      gsl_histogram_increment (h, depth[i]);
  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    fo << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i)*scale << std::endl;
  }  
    gsl_histogram_free (h);
  fo.close();
}  

void saveDepth(std::vector<Float>& depth,const char *filename)
{
  std::ofstream fo(filename);
  fo << depth.size() << std::endl;
  for(size_t i = 0; i < depth.size(); i++)
    fo << depth[i] << std::endl;
  fo.close();
}



void
saveHistogram(gsl_histogram *h, const char *datFileName)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());
  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << std::endl;
  }  
//  gsl_histogram_free (h);
  foByEscapeHist.close();
  foByEscapeHistPlt << "#reset\nset yrange [0:*]\nplot \'" << byEscapeTimeDatHist << "\' with histeps\n";
  foByEscapeHistPlt.close();
}  

void
saveHistogram_new(gsl_histogram *h, const char *datFileName)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());
  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << std::endl;
  }  
//  gsl_histogram_free (h);
  foByEscapeHist.close();
  foByEscapeHistPlt << "reset\nset yrange [0:*]\nset xrange [-180:180]\nplot \'" << byEscapeTimeDatHist << "\' with boxes\n";
  foByEscapeHistPlt << "pause -1 \"Press Enter\"\n";
  foByEscapeHistPlt.close();
}  

void
saveHistogram_polar(gsl_histogram *h, const char *datFileName/*, bool halfshift = false*/)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());

  int n;
  n = gsl_histogram_bins(h);

  for(size_t i = 0; i < gsl_histogram_bins(h)+1; i++)
  {
    double lower, upper;
    double index = i;
    if (i == gsl_histogram_bins(h)) index = 0;
    gsl_histogram_get_range (h, index, &lower, &upper);
    Float ang = (lower+upper)/2.0;//-(halfshift?((360.0/n)/2.0):0.0);
//    if (ang < -180) ang += ;
    foByEscapeHist << ang << " " << gsl_histogram_get(h,index) << std::endl;
  }  
//  gsl_histogram_free (h);
  foByEscapeHist.close();
  foByEscapeHistPlt << "reset\nset style fill pattern 1\nset polar\nset angles degrees\nset size ratio -1\n\
\nset grid polar\n\nplot \'" << byEscapeTimeDatHist << "\' with filledcurves lw 2 notitle,\\\n \'" << byEscapeTimeDatHist << "\' with impulses lw 2 notitle\n";
//  foByEscapeHistPlt << "pause -1 \"Press Enter\"\n";
  foByEscapeHistPlt.close();
}  

void
findIntermediateStates(std::string trajDir,std::vector<std::string>& states);

void
StatPostProcess::execute() 
{
  using mdtk::Exception;

#define MDEPP_MIN_TRANSLATION_ACCOUNTED 1.0*mdtk::Ao
//int totalTransitionCandidates = 0;

  cout << "PostProcess::execute() started." << std::endl;
  cerr << "PostProcess::execute() started." << std::endl;

  for(size_t i = 0; i < /*getStatesSize()*/trajData.size(); i++)
  {
    mdtk::SimLoop* currentState = new mdtk::SimLoop();
    currentState->allowToFreePotentials = true;
    currentState->allowToFreeAtoms = true;
    setupPotentials(*currentState);
    std::string trajFinalName = trajData[i].trajDir+"mde_init"; //was mde_final
    cout << "Loading state " << trajFinalName << std::endl;

    yaatk::text_ifstream fi(trajFinalName.c_str()); 
//    YAATK_FSTREAM_CREATE(ifstream,fi,savedStateNames[i].c_str()); 

currentState->initNLafterLoading = false;
    currentState->loadFromStream(fi);
UPD_NL_REF_POT(currentState);
    setTags(currentState);
//currentState->initNLafterLoading = true;

    fi.close();


    dummy_ac.setPBC(currentState->atoms_.getPBC());
    TRACE(dummy_ac.getPBC()/mdtk::Ao);
    TRACE(&dummy_ac);

{
    std::vector<std::string> interStates;
    findIntermediateStates(trajData[i].trajDir,interStates);
TRACE(interStates.size());
    int stateIndex = interStates.size()-1;
//    for(int stateIndex = interStates.size()-1; stateIndex >= 0; stateIndex--)
    {
      TRACE(interStates[stateIndex]); //exit(1);
      std::string mde_inter_filename = trajData[i].trajDir+interStates[stateIndex];
      //to remove .GZ simply resize
      mde_inter_filename.resize(mde_inter_filename.size()-3);
      TRACE(mde_inter_filename);
/*
      YAATK_IFSTREAM_CREATE_ZIPPED_OPT(std::ifstream,fi,mde_inter_filename.c_str(),std::ios::binary); 
      currentState->loadFromStreamXVA_bin(fi);
*/


mdtk::AtomsContainer atoms_start;
  for(size_t i1 = 0; i1 < currentState->atoms_.size(); i1++)
  {
   mdtk::Atom *new_atom;
   new_atom = new mdtk::Atom();
   *new_atom = *(currentState->atoms_[i1]);
   atoms_start.push_back(new_atom);
  }  

      yaatk::text_ifstream fi(mde_inter_filename.c_str()); 
      currentState->loadFromStreamXVA(fi/*,false*/);
UPD_NL_REF_POT(currentState);
      setTags(currentState);
      fi.close(); 

mdtk::AtomsContainer atoms_end;
  for(size_t i1 = 0; i1 < currentState->atoms_.size(); i1++)
  {
   mdtk::Atom *new_atom;
   new_atom = new mdtk::Atom();
   *new_atom = *(currentState->atoms_[i1]);
   atoms_end.push_back(new_atom);
  }  

  REQUIRE(atoms_start.size() == atoms_end.size());
  REQUIRE(atoms_start.size() == currentState->atoms_.size());

  for(size_t i1 = 0; i1 < atoms_start.size(); i1++)
  {
    Float trans_dist = fabs(atoms_start[i1]->coords.z - atoms_end[i1]->coords.z);
    if (trans_dist > MDEPP_MIN_TRANSLATION_ACCOUNTED)
    {
//      TRACE("adding transition....");
//      TRACE(trans_dist/mdtk::Ao);
//      if (trans_dist > 3.61*mdtk::Ao/2.0) totalTransitionCandidates++;
      trajData[i].trans.translations.push_back(Translation(*(atoms_start[i1]),*(atoms_end[i1])));
    }
  }

//TRACE("ADDING COMPLETED");

  for(size_t i1 = 0; i1 < atoms_start.size(); i1++) delete atoms_start[i1];
  for(size_t i1 = 0; i1 < atoms_end.size(); i1++)   delete atoms_end[i1];

    }
}

    TRACE("mde_final loaded!");
//    YAATK_FSTREAM_CLOSE(fi); 


//    currentState->fpot.NL_init(currentState->atoms_);
//    currentState->fpot.NL_UpdateIfNeeded(currentState->atoms_); /// DO NOT DELETE!!!!!!!!!!!!!
//    NeighbourList_Update(&(currentState->fpot),currentState->atoms_,SPOTTED_DISTANCE);


//    ProcessState(*currentState,savedStateFileNames[i]);
    ProcessState(*currentState,i);
    delete currentState;
  }  

//  TRACE(totalTransitionCandidates);

  cout << "PostProcess::execute() done." << std::endl;
  cerr << "PostProcess::execute() done." << std::endl;
}  


void
StatPostProcess::ProcessState(mdtk::SimLoop& state, size_t i)
{
//  trajData.push_back(TrajData());

//  trajData[i].trajDir = extractDir(trajNameFinal);
//  (trajData.end()-1)->trajNameFinal = trajNameFinal;
  trajData[i].aboveSpottedHeight = getAboveSpottedHeight(state);
  buildSpottedMolecules(state,i);
}  


bool
Molecule::hasAtom(mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
      found = true;
  return found;
}  

void
Molecule::buildFromAtom(mdtk::Atom& a, mdtk::SimLoop& ml,double SPOTTED_DISTANCE)
{
//TRACE("adding atom");
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);atoms[atoms.size()-1].container = &dummy_ac;

//  TRACE(a.fixed);
//  TRACE(atoms[atoms.size()-1].fixed);
  

//TRACE(a.globalIndex);  
//TRACE(a.ID);  
//TRACE("1");  
//TRACE(a.nl(&(ml.fpot)).size());
  for(size_t i = 0; i < REF_POT_OF(ml.fpot)->NL(a).size(); i++)
  {
//TRACE("2");
    mdtk::Atom& nb_a = *(REF_POT_OF(ml.fpot)->NL(a)[i]);
//    nb_a.container = &dummy_ac;
//       a.container = &dummy_ac;
    if (!isHandled(nb_a)) continue;
    Float distance = REF_POT_OF(ml.fpot)->r_vec_module_no_touch(a,nb_a);//sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
    if (Rc(a,nb_a) >= distance)
    {
      if (nb_a.coords.z < SPOTTED_DISTANCE)
      {
//TRACE("Next dive");        
        buildFromAtom(nb_a,ml,SPOTTED_DISTANCE);
      }
      else
      {
//TRACE("atoms.clear");        
        atoms.clear();
      }  
      if (atoms.size() == 0) break;
    }  
  }  
}  


int
StatPostProcess::getAboveSpottedHeightTotal() const
{
  int totalSpotted = 0;
  for(size_t i = 0; i < trajData.size();i++)
  {
    totalSpotted += trajData[i].aboveSpottedHeight;
  }  
  return totalSpotted;
}  

int
StatPostProcess::getAboveSpottedHeight(mdtk::SimLoop& state) const
{
  int spotted;
  {
    spotted = 0;
    for(size_t atomIndex = 0; atomIndex < state.atoms_.size(); atomIndex++)
    {
      mdtk::Atom &atom = *(state.atoms_[atomIndex]);
      if (atom.coords.z < SPOTTED_DISTANCE)
      {
        spotted++;
      }  
    }  
//    cout << "Trajectory " << trajIndex;
//    cout << " has " << spotted << " spotted atoms." << std::endl;
  }  
  return spotted;
}  

int
StatPostProcess::getYieldSum( FProcessMolecule fpm) const
{
  int totalSpotted = 0;
  for(size_t i = 0; i < trajData.size();i++)
  {
    totalSpotted += getYield(i,fpm);
  }  
  return totalSpotted;
}  

int
StatPostProcess::getYield(size_t trajIndex, FProcessMolecule fpm) const
{
  int spotted = 0;
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const Molecule& mol = td.molecules[mi];
    if (!fpm(mol)) continue;
    spotted += mol.atoms.size();
/*
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      mdtk::Atom& atom = td.molecules[mi].atoms[ai];
    } 
*/
  }  
  return spotted;
}  

Float
StatPostProcess::getAverageYield( FProcessMolecule fpm) const
{
  return Float(getYieldSum(fpm))/trajData.size();
}  

Float
StatPostProcess::getAverageYieldProgress( FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"average_Yield.dat");
  {
    std::ofstream fo((std::string(ofilename)+".plt").c_str());
    fo << "reset\n#set yrange [0:*]\nplot '" << ofilename << "' with impulses";
    fo.close();
  }  

  std::ofstream fo(ofilename);
  Float avg = 0.0;
  size_t i;
  for(i = 0; i < trajData.size();i++)
  {
    avg += getYield(i,fpm);
    fo << i+1 << " " << avg/(i+1) << "\n";
  }  
  fo.close();
  return avg/i;
}  






Float
StatPostProcess::getTotalEnergyOfSputtered( FProcessMolecule fpm) const
{
  Float totalSpotted = 0;
  for(size_t i = 0; i < trajData.size();i++)
  {
    totalSpotted += getEnergyOfSputtered(i,fpm);
  }  
  return totalSpotted;
}  

Float
StatPostProcess::getEnergyOfSputtered(size_t trajIndex, FProcessMolecule fpm) const
{
  Float spotted = 0;
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const Molecule& mol = td.molecules[mi];
    if (!fpm(mol)) continue;
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      spotted += SQR(atom.V.module())*atom.M/2.0;
    } 
  }  
  return spotted;
}  

Float
StatPostProcess::getAverageEnergyOfSputtered( FProcessMolecule fpm) const
{
  return Float(getTotalEnergyOfSputtered(fpm))/trajData.size();
}  



/*
void
StatPostProcess::buildSpottedMoleculesTotal()
{
  cout << "Building molecules started." << std::endl;
  cerr << "Building molecules started." << std::endl;
  for(size_t i = 0; i < trajData.size();i++)
  {
    buildSpottedMolecules(i);
  }  
  cout << "Building molecules done." << std::endl;
  cerr << "Building molecules done." << std::endl;
}  
*/
/*
#include "dosdir/dosdir.h"


#ifndef __WIN32__
#define _finddata_t dd_ffblk
#define name dd_name
#define attrib dd_attribs
#define _A_SUBDIR DD_DIREC
#define _findfirst dd_findfirst
#define _findnext dd_findnext
#define _findclose dd_findclose
#endif


void
findIntermediateStates(std::string trajDir,std::vector<std::string>& states)
{
//  TRACE("in find!");
  char stateFileWildCard[10000];
  sprintf(stateFileWildCard,"%smde0*",trajDir.c_str());
  TRACE(stateFileWildCard);
  struct _finddata_t ff;
  int done = 0;
  int findhandle;
  findhandle = _findfirst(stateFileWildCard, &ff ,DD_NORMAL);
  if (findhandle != -1)
  {
    while (!done)
    {
//      if ((ff.attrib & _A_SUBDIR) && strcmp(ff.name,".") && strcmp(ff.name,".."))
      {
        states.push_back(ff.name);
      }  
      done = _findnext(findhandle,&ff);
    }
//    _findclose(findhandle);
  }  
  sort(states.begin(), states.end());
}  
*/

void
removeDuplicates(std::vector<std::string>& states)
{
  size_t i;
/*
  for(i = 0; i < states.size(); i++)
    TRACE(states[i]);
*/
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
/*
  for(i = 0; i < states.size(); i++)
    TRACE(states[i]);

  throw;
*/
}

#include <dirent.h>

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

void
StatPostProcess::buildSpottedMolecules(mdtk::SimLoop& state,size_t trajIndex)
{
  using mdtk::Exception;    

  cout << "Building molecules for state ..." << std::endl;
  {
//TRACE(states[trajIndex]->atoms_.size());
    for(size_t atomIndex = 0; atomIndex < state.atoms_.size(); atomIndex++)
    {
      mdtk::Atom &atom = *(state.atoms_[atomIndex]);
      if (atom.coords.z < SPOTTED_DISTANCE)
      {
        {
//TRACE(atom.globalIndex);
          bool account_atom = true;
//          TRACE(trajData[trajIndex].molecules.size());
          for(size_t mi = 0; mi < trajData[trajIndex].molecules.size(); mi++)
          {
//TRACE(mi);
            if (trajData[trajIndex].molecules[mi].hasAtom(atom))
            {
              account_atom = false;
              break;
            }  
          }  
//TRACE(account_atom);
          if (account_atom) 
          {
//            TRACE("Building molecule.");
            Molecule molecule;
//            molecule.trajectory = trajIndex;//trajIndex;
//TRACE(molecule.trajectory);            
            molecule.buildFromAtom(atom,state,SPOTTED_DISTANCE);
//TRACE("After build from atom.");
//TRACE(molecule.atoms.size());
            if (molecule.atoms.size() > 0 && molecule.getVelocity().z < 0.0)
            {
              cout << "Adding molecule." << std::endl;
              trajData[trajIndex].molecules.push_back(molecule);
            }  
          }  
        }  
      }  
    }  
  }  

// InnnerCluster() ----------------------
  cout << "Building cluster dynamics for state ..." << std::endl;
  {
    trajData[trajIndex].clusterDynamics.PBC = state.getPBC();
//TRACE(states[trajIndex]->atoms_.size());
    for(size_t atomIndex = 0; atomIndex < state.atoms_.size(); atomIndex++)
    {
      mdtk::Atom &atom = *(state.atoms_[atomIndex]);
//      if (atom.ID != mdtk::Cu_EL) continue;
      if (!(atom.tag & ATOMTAG_CLUSTER)) continue;
      trajData[trajIndex].clusterDynamics.atomTrajectories.push_back(AtomTrajectory(atom));
    }
  }
// ~InnnerCluster() ----------------------

// Projectile() ----------------------
  cout << "Building projectile dynamics for state ..." << std::endl;
  {
    trajData[trajIndex].projectileDynamics.PBC = state.getPBC();
//TRACE(states[trajIndex]->atoms_.size());
    for(size_t atomIndex = 0; atomIndex < state.atoms_.size(); atomIndex++)
    {
      mdtk::Atom &atom = *(state.atoms_[atomIndex]);
      if (atom.ID != mdtk::Ar_EL) continue;
//      if (!(atom.tag & ATOMTAG_CLUSTER)) continue;
      trajData[trajIndex].projectileDynamics.atomTrajectories.push_back(AtomTrajectory(atom));
    }
  }
// ~Projectile() ----------------------

  TrajData& td = trajData[trajIndex];//trajData[traj];
//  if (td.molecules.size() > 0)
  {

    mdtk::SimLoop* mde_init = new mdtk::SimLoop();    
    mde_init->allowToFreePotentials = true;
    mde_init->allowToFreeAtoms = true;
    setupPotentials(*mde_init);
    std::string mde_init_filename = td.trajDir+"mde_init";
    TRACE(mde_init_filename);
    yaatk::text_ifstream fi(mde_init_filename.c_str()); 
mde_init->initNLafterLoading = false;
    mde_init->loadFromStream(fi/*,false*/);
UPD_NL_REF_POT(mde_init);
    setTags(mde_init);
    fi.close(); 

    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
      {
        const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
        const mdtk::Atom& atom_init = *(mde_init->atoms_[atom.globalIndex]);
        td.molecules[mi].atoms_init.push_back(atom_init); 
          td.molecules[mi].atoms_init[td.molecules[mi].atoms_init.size()-1].container = &dummy_ac;
        REQUIRE(td.molecules[mi].atoms[ai].globalIndex == td.molecules[mi].atoms_init[ai].globalIndex);
      }
      REQUIRE(td.molecules[mi].atoms.size() == td.molecules[mi].atoms_init.size());
    }  

// InnnerCluster() ----------------------
    for(size_t clusterAtomIndex = 0; clusterAtomIndex < td.clusterDynamics.atomTrajectories.size(); clusterAtomIndex++)
    {
      const mdtk::Atom &atom = td.clusterDynamics.atomTrajectories[clusterAtomIndex].endCheckPoint;
      const mdtk::Atom &atom_init = *(mde_init->atoms_[atom.globalIndex]);
      td.clusterDynamics.atomTrajectories[clusterAtomIndex].beginCheckPoint = atom_init;
    }
// ~InnnerCluster() ----------------------

// Projectile() ----------------------
    TRACE(td.projectileDynamics.atomTrajectories.size());
//    REQUIRE(td.projectileDynamics.atomTrajectories.size() == 1);
    for(size_t projectileAtomIndex = 0; projectileAtomIndex < td.projectileDynamics.atomTrajectories.size(); projectileAtomIndex++)
    {
      const mdtk::Atom &atom = td.projectileDynamics.atomTrajectories[projectileAtomIndex].endCheckPoint;
      const mdtk::Atom &atom_init = *(mde_init->atoms_[atom.globalIndex]);
      td.projectileDynamics.atomTrajectories[projectileAtomIndex].beginCheckPoint = atom_init;
    }
// ~ProjectileCluster() ----------------------


    delete mde_init;
  }  

//  if (td.molecules.size() > 0)
  {
    std::vector<std::string> interStates;
    findIntermediateStates(td.trajDir,interStates);


      mdtk::SimLoop* mde_inter = new mdtk::SimLoop();    
      mde_inter->allowToFreePotentials = true;
      mde_inter->allowToFreeAtoms = true;
      setupPotentials(*mde_inter);

{
    std::string trajFinalName = td.trajDir+"mde_init"; //was mde_final
    cout << "Loading state " << trajFinalName << std::endl;

    yaatk::text_ifstream fi(trajFinalName.c_str()); 
//    YAATK_FSTREAM_CREATE(ifstream,fi,savedStateNames[i].c_str()); 

mde_inter->initNLafterLoading = false;
    mde_inter->loadFromStream(fi);
/*
UPD_NL_REF_POT(mde_inter);
    setTags(mde_inter); NOT NEEDED
*/
//if (td.molecules.size() > 0)
//mde_inter->initNLafterLoading = true;


    fi.close();
}


    for(int stateIndex = interStates.size()-1; stateIndex >= 0; stateIndex--)
    {
//      TRACE(interStates[stateIndex]);

      std::string mde_inter_filename = td.trajDir+interStates[stateIndex];
      //to remove .GZ simply resize
      mde_inter_filename.resize(mde_inter_filename.size()-3);
      TRACE(mde_inter_filename);
/*
      YAATK_IFSTREAM_CREATE_ZIPPED_OPT(std::ifstream,fi,mde_inter_filename.c_str(),std::ios::binary); 
      mde_inter->loadFromStreamXVA_bin(fi);
*/
      yaatk::text_ifstream fi(mde_inter_filename.c_str()); 
      mde_inter->loadFromStreamXVA(fi/*,false*/);
UPD_NL_REF_POT(mde_inter);
      setTags(mde_inter);
      fi.close(); 

/*
      mde_inter->fpot.NL_init(mde_inter->atoms_);
      NeighbourList_Update(&(mde_inter->fpot),mde_inter->atoms_,SPOTTED_DISTANCE);
*/

//TRACE("***** loaded OK");

      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        bool formedNowOrEarlier = true;
        bool escapedNowOrEarlier = true;
        Molecule& molecule = td.molecules[mi];
        for(size_t ai = 0; ai < molecule.atoms.size(); ai++)
        {
          /*const */mdtk::Atom& atom_i = molecule.atoms[ai];
          /*const */mdtk::Atom& atom_i_inter = *(mde_inter->atoms_[atom_i.globalIndex]);
          for(size_t aj = 0; aj < molecule.atoms.size(); aj++)
          if (ai != aj)
          {
            /*const */mdtk::Atom& atom_j = molecule.atoms[aj];
            /*const */mdtk::Atom& atom_j_inter = *(mde_inter->atoms_[atom_j.globalIndex]);
/*
            if (!molecule_inter.hasSameAtoms(molecule))
            {
              alreadyFormed = false;
              break;
            }  
*/
/*
TRACE(atom_i.globalIndex);
TRACE(atom_j.globalIndex);
TRACE(atom_i_inter.globalIndex);
TRACE(atom_j_inter.globalIndex);

TRACE(&dummy_ac);
TRACE(atom_i.container);
TRACE(atom_i.container->getPBC());
TRACE(atom_j.container->getPBC());
*/
            Float distance_ij = 
              REF_POT_OF(mde_inter->fpot)->r_vec_module_no_touch(atom_i,atom_j);
            Float distance_ij_inter = 
              REF_POT_OF(mde_inter->fpot)->r_vec_module_no_touch(atom_i_inter,atom_j_inter);
              
            if (
                (
                molecule.Rc(atom_i      ,atom_j      ) >= distance_ij &&
                molecule.Rc(atom_i_inter,atom_j_inter) <  distance_ij_inter
                )||
                (
                molecule.Rc(atom_i      ,atom_j      ) <  distance_ij &&
                molecule.Rc(atom_i_inter,atom_j_inter) >= distance_ij_inter
                )
               )
             {
               formedNowOrEarlier = false;
               break;
             }  
          }
        }

//TRACE("***** 1st loop OK");

//        if (formedNowOrEarlier)
        for(size_t ai = 0; ai < molecule.atoms.size(); ai++)
        {
          const mdtk::Atom& atom_i = molecule.atoms[ai];
//          const mdtk::Atom& atom_i_inter = *(mde_inter->atoms_[atom_i.globalIndex]);
          Molecule molecule_inter;
          molecule_inter.buildFromAtom(*(mde_inter->atoms_[atom_i.globalIndex]),*mde_inter,SPOTTED_DISTANCE);

          if (molecule_inter.atoms.size() == 0)
          {
            escapedNowOrEarlier = false;
            break;
          }  
          if (!(molecule_inter.getVelocity().z < 0.0))
          {
            escapedNowOrEarlier = false;
            break;
          }  
          for(size_t ak_inter = 0; ak_inter < molecule_inter.atoms.size(); ak_inter++)
          {  
            const mdtk::Atom& atom_k_inter = molecule_inter.atoms[ak_inter];
            if (atom_k_inter.coords.z >= SPOTTED_DISTANCE)
            {
              escapedNowOrEarlier = false;
              break;
            }  
          }
        }
        if (formedNowOrEarlier) molecule.formationTime 
           = (mde_inter->simTime>0.01*mdtk::ps)?mde_inter->simTime:0.0*mdtk::ps;
        if (escapedNowOrEarlier) molecule.escapeTime
           = (mde_inter->simTime>0.01*mdtk::ps)?mde_inter->simTime:0.0*mdtk::ps;

      }  

//TRACE("***** 2nd loop OK");

// InnnerCluster() ----------------------
    for(size_t clusterAtomIndex = 0; clusterAtomIndex < td.clusterDynamics.atomTrajectories.size(); clusterAtomIndex++)
    {
      const mdtk::Atom &atom = td.clusterDynamics.atomTrajectories[clusterAtomIndex].endCheckPoint;
      const mdtk::Atom &atom_inter = *(mde_inter->atoms_[atom.globalIndex]);
      td.clusterDynamics.atomTrajectories[clusterAtomIndex].checkPoints.push_back(atom_inter);
    }
// ~InnnerCluster() ----------------------

// Projectile() ----------------------
//    REQUIRE(td.projectileDynamics.atomTrajectories.size() == 1);
    TRACE(td.projectileDynamics.atomTrajectories.size());
    for(size_t projectileAtomIndex = 0; projectileAtomIndex < td.projectileDynamics.atomTrajectories.size(); projectileAtomIndex++)
    {
      const mdtk::Atom &atom = td.projectileDynamics.atomTrajectories[projectileAtomIndex].endCheckPoint;
      const mdtk::Atom &atom_inter = *(mde_inter->atoms_[atom.globalIndex]);
      td.projectileDynamics.atomTrajectories[projectileAtomIndex].checkPoints.push_back(atom_inter);
    }
// ~Projectile() ----------------------


//TRACE("***** 3rd loop OK");

    }  
      delete mde_inter;
  }  

  cout << "Building molecules for state done." << std::endl;
}  

void
StatPostProcess::buildClusterFragmentsFromDyn()
{
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    TrajData& td = trajData[trajIndex];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const mdtk::Atom &atom = td.clusterDynamics.atomTrajectories[atomIndex].endCheckPoint;
      {
        {
          bool account_atom = true;
          for(size_t mi = 0; mi < td.clusterDynamics.fragments.size(); mi++)
          {
            if (td.clusterDynamics.fragments[mi].hasAtom(atom))
            {
              account_atom = false;
              break;
            }  
          }  
          if (account_atom) 
          {
            Fragment molecule;
            molecule.buildFromAtom(atom,td.clusterDynamics);
//            if (molecule.atoms.size() > 0 && molecule.getVelocity().z < 0.0)
            {
//              cout << "Adding molecule." << std::endl;
              td.clusterDynamics.fragments.push_back(molecule);
            }  
          }  
        }  
      }
    }  
  }  

}
  
void
StatPostProcess::printMoleculesTotal() const
{
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].molecules.size() > 0)
  {
    cout << "Molecules for trajectory " << traj << 
     " ("  << trajData[traj].trajDir << ") " << " :\n";
    printMolecules(traj);
    cout << std::endl;
  }  
}  

void
StatPostProcess::printMolecules(size_t trajIndex) const
{
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    cout << "Molecule #" << mi << " : ";
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      cout << "( Z = " << atom.Z/mdtk::e << " , gi = " << atom.globalIndex << " ) ";
    } 
    cout << "\n            formation time = " << td.molecules[mi].formationTime/mdtk::ps << " ps"
         << ", escape time = " << td.molecules[mi].escapeTime/mdtk::ps << " ps";
    if (td.molecules[mi].hasSubstrateAtoms() && td.molecules[mi].hasClusterAtoms())
      cout << "\nThis is HETEROGENEOUS molecule!!!";
    cout << std::endl;
  }  
}  


#define FO_COMMON "\nset terminal png giant size 2200,2000 enhanced\n"\
                  "set format y \"%3g\"\n" "set format x \"%3g\"\n"

#define STOPPING_MIN_DR_ACCOUNTED 1.0*Ao

void
StatPostProcess::printStoppingHist() const
{
  using namespace mdtk;

  yaatk::mkdir("_stoppingHist");
  yaatk::chdir("_stoppingHist");

  using mdtk::Exception;    

  const char *filename = "StoppingHist.dat";

  {
    std::ofstream fo((std::string(filename)+".plt").c_str());
    fo << "reset\nset xyplane at 0.0\n";
    fo << "splot \'" << filename << "\' with impulses\n";
    fo << "pause -1 \"Press Return\"\n";
    fo.close();
  }  
  std::ofstream fo(filename);


    const Float xmin = 0.0*mdtk::Ao;
    const Float xmax = trajData[0].clusterDynamics.PBC.x;
    const Float ymin = 0.0*mdtk::Ao;
    const Float ymax = trajData[0].clusterDynamics.PBC.y;

    const int nx = 20;
    const int ny = 15;

  fo << "# xmin = " << xmin << " Ao" << std::endl
     << "# xmax = " << xmax << " Ao" << std::endl
     << "# number of xbins = " << nx << std::endl;
  fo << "# ymin = " << ymin << " Ao" << std::endl
     << "# ymax = " << ymax << " Ao" << std::endl
     << "# number of ybins = " << ny << std::endl;
    gsl_histogram2d * h = gsl_histogram2d_alloc (nx,ny);
    gsl_histogram2d_set_ranges_uniform (h, xmin, xmax, ymin, ymax);

//  size_t atomsTotal = 0.0;

  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
//      Float x = atraj.beginCheckPoint.coords.x;
//      Float y = atraj.beginCheckPoint.coords.y;
      Float x = atraj.endCheckPoint.coords.x;
      Float y = atraj.endCheckPoint.coords.y;
      Vector3D PBC = td.clusterDynamics.PBC;
      Float de = 1e-10*Ao;
      if (x <= 0.0)   x = de;
      if (x >= PBC.x) x = PBC.x - de;
      if (y <= 0.0)   y = de;
      if (y >= PBC.y) y = PBC.y - de;

      if ((atraj.endCheckPoint.coords-atraj.beginCheckPoint.coords).module() < STOPPING_MIN_DR_ACCOUNTED) continue;

//      atomsTotal++;

      gsl_histogram2d_increment (h, x, y);
    }
  }

  for(int i = 0; i < nx; i++)
  for(int j = 0; j < ny; j++)
  {
    double lower, upper;
    gsl_histogram2d_get_xrange (h, i, &lower, &upper);
    fo << (lower+upper)/2.0 << " ";
    gsl_histogram2d_get_yrange (h, j, &lower, &upper);
    fo << (lower+upper)/2.0 << " ";
    fo << gsl_histogram2d_get(h,i,j)/trajData.size()/*Float(atomsTotal)*/ << std::endl;
  }  

  fo.close();

{
  const char *filename_bars = "StoppingHistBars.dat";
  {
    std::ofstream fo_bars((std::string(filename_bars)+".plt").c_str());
    fo_bars << "reset\n";
//    fo_bars << "set surface\n ";
//    fo_bars << "set contour surface\n ";
    fo_bars << "set xyplane at 0.0\n";
    fo_bars << "set xrange [0:"<< trajData[0].clusterDynamics.PBC.x/mdtk::Ao << "]\n";
    fo_bars << "set yrange [0:"<< trajData[0].clusterDynamics.PBC.y/mdtk::Ao << "]\n";
    fo_bars << "set xtics " << trajData[0].clusterDynamics.PBC.x/mdtk::Ao/nx << "\n";
    fo_bars << "set ytics " << trajData[0].clusterDynamics.PBC.y/mdtk::Ao/ny << "\n";
//    fo_bars << "set format xy \"%5g\"\n";
    fo_bars << "set xlabel \"x,Ao\"\n";
    fo_bars << "set ylabel \"y,Ao\"\n";
    fo_bars << "set zlabel \"Probability of stopping\"\n";
    fo_bars << "set hidden3d\n";
    fo_bars << "splot \'" << filename_bars << "\' with lines lc 1\n";
    fo_bars << "pause -1 \"Press Return\"\n";
    fo_bars.close();
  }  

  {
    std::ofstream fo_bars((std::string(filename_bars)+".png.plt").c_str());
    fo_bars << "reset\n";
//    fo_bars << "set surface\n ";
//    fo_bars << "set contour surface\n ";
    fo_bars << "set xyplane at 0.0\n";
    fo_bars << "set xrange [0:"<< trajData[0].clusterDynamics.PBC.x/mdtk::Ao << "]\n";
    fo_bars << "set yrange [0:"<< trajData[0].clusterDynamics.PBC.y/mdtk::Ao << "]\n";
    fo_bars << "set xtics " << trajData[0].clusterDynamics.PBC.x/mdtk::Ao/nx << "\n";
    fo_bars << "set ytics " << trajData[0].clusterDynamics.PBC.y/mdtk::Ao/ny << "\n";
//    fo_bars << "set format xy \"%5g\"\n";
    fo_bars << "set xlabel \"x,Ao\"\n";
    fo_bars << "set ylabel \"y,Ao\"\n";
    fo_bars << "set zlabel \"Probability of stopping\"\n";
    fo_bars << "set hidden3d\n";

  fo_bars << "set output \""<< (std::string(filename_bars)+".png.plt") << ".png\"";
  fo_bars << FO_COMMON;

    fo_bars << "splot \'" << filename_bars << "\' with lines lc 1\n";
    fo_bars << "pause -1 \"Press Return\"\n";
    fo_bars.close();
  }  

  std::ofstream fo_bars(filename_bars);
  Float de = 1e-2*mdtk::Ao;

/*
  for(int j = 0; j < ny; j++)
  {
     double y[2];
     gsl_histogram2d_get_yrange (h, j, y, y+1);
     fo_bars << -2.0*de << " " << y[0]/mdtk::Ao << " " << 0.0 << std::endl;
     fo_bars << -2.0*de << " " << y[1]/mdtk::Ao << " " << 0.0 << std::endl;
  }
  for(int j = 0; j < ny; j++)
  {
     double y[2];
     gsl_histogram2d_get_yrange (h, j, y, y+1);
     fo_bars << de << " " << y[0]/mdtk::Ao << " " << 0.0 << std::endl;
     fo_bars << de << " " << y[1]/mdtk::Ao << " " << 0.0 << std::endl;
  }
*/

  for(int i = 0; i < nx; i++)
  {
    double x[2];
    gsl_histogram2d_get_xrange (h, i, x, x+1);
      x[0] += de;
      x[1] -= de;
    for(int xi = 0; xi < 2; xi++)
    {
    for(int j = 0; j < ny; j++)
    {
      double y[2];
      gsl_histogram2d_get_yrange (h, j, y, y+1);
      double val = gsl_histogram2d_get(h,i,j);
      val /= trajData.size();//atomsTotal;

     fo_bars << x[xi]/mdtk::Ao << " " << y[0]/mdtk::Ao << " " << val << std::endl;
     fo_bars << x[xi]/mdtk::Ao << " " << y[1]/mdtk::Ao << " " << val << std::endl;

/*
      fo_bars << xlower << " " << ylower << " " << 0 << std::endl;
      fo_bars << xlower << " " << yupper << " " << 0 << std::endl;
      fo_bars << xupper << " " << yupper << " " << 0 << std::endl;
      fo_bars << xupper << " " << ylower << " " << 0 << std::endl;
      fo_bars << xlower << " " << ylower << " " << 0 << std::endl;

      fo_bars << xlower << " " << ylower << " " << val << std::endl;
      fo_bars << xlower << " " << yupper << " " << val << std::endl;
      fo_bars << xupper << " " << yupper << " " << val << std::endl;
      fo_bars << xupper << " " << ylower << " " << val << std::endl;
      fo_bars << xlower << " " << ylower << " " << val << std::endl;

      fo_bars << xlower << " " << ylower << " " << 0 << std::endl;
*/
    }
    fo_bars << std::endl;
    }
  }  

/*
  for(int j = 0; j < ny; j++)
  {
     double y[2];
     gsl_histogram2d_get_yrange (h, j, y, y+1);
     fo_bars << trajData[0].clusterDynamics.PBC.x+de << " " << y[0]/mdtk::Ao << " " << 0.0 << std::endl;
     fo_bars << trajData[0].clusterDynamics.PBC.x+de << " " << y[1]/mdtk::Ao << " " << 0.0 << std::endl;
  }
  for(int j = 0; j < ny; j++)
  {
     double y[2];
     gsl_histogram2d_get_yrange (h, j, y, y+1);
     fo_bars << trajData[0].clusterDynamics.PBC.x+2.0*de << " " << y[0]/mdtk::Ao << " " << 0.0 << std::endl;
     fo_bars << trajData[0].clusterDynamics.PBC.x+2.0*de << " " << y[1]/mdtk::Ao << " " << 0.0 << std::endl;
  }
*/

  fo_bars.close();
}

  gsl_histogram2d_free (h);

  yaatk::chdir("..");
}


void
StatPostProcess::printClusterAtomsByAzimuth(const int n) const
{
  using namespace mdtk;

  yaatk::mkdir("_byAzimuth");
  yaatk::chdir("_byAzimuth");

  using mdtk::Exception;    

  char ofilename[1024];
  sprintf(ofilename,"ClusterAtomsCount_by_azimuth_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

//  Float atomsTotal = 0.0;

  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      mdtk::Vector3D v = atraj.endCheckPoint.V;

      for(size_t i = 0; i < atraj.checkPoints.size(); i++)
        if (atraj.checkPoints[i].V.module() > v.module()) v = atraj.checkPoints[i].V;

//      TRACE(v.module());
//      if (v.module() < 300000) continue;
      if ((atraj.endCheckPoint.coords-atraj.beginCheckPoint.coords).module() < STOPPING_MIN_DR_ACCOUNTED) continue;

//      atomsTotal++;

      Float polar = atan2(v.y,v.x)/mdtk::Deg;
      gsl_histogram_accumulate (h, polar, 1.0/Float(trajData.size())   );
    }
  }

  

  saveHistogram_new(h, ofilename);
  gsl_histogram_free (h);


  yaatk::chdir("..");
}

void
StatPostProcess::printClusterAtomsByAzimuthPos(const int n, bool halfshift, Float accountedDist) const
{
  const mdtk::Vector3D& PBC = trajData[0].clusterDynamics.PBC;

  mdtk::Vector3D CenterOfCluster = mdtk::Vector3D(0.0,0.0,0.0);//PBC/2.0
  {
    const TrajData& td = trajData[0];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      mdtk::Vector3D initCoords = atraj.beginCheckPoint.coords;
      CenterOfCluster += initCoords;
    }    
    CenterOfCluster /= td.clusterDynamics.atomTrajectories.size();
    TRACE(td.clusterDynamics.atomTrajectories.size());
  }

  TRACE(PBC/2.0/mdtk::Ao);
  TRACE(CenterOfCluster/mdtk::Ao);

  using namespace mdtk;

  yaatk::mkdir("_byAzimuthPos");
  yaatk::chdir("_byAzimuthPos");

  using mdtk::Exception;    

  char ofilename[1024];
  sprintf(ofilename,"Atoms_per_traj_%05d%s-ad%.2f.dat",
    n,halfshift?"-halfshift":"",accountedDist/mdtk::Ao);
  char ofilename_count[1024];
  sprintf(ofilename_count,"Atoms_%05d%s-ad%.2f.dat",
    n,halfshift?"-halfshift":"",accountedDist/mdtk::Ao);
  char ofilename_dist[1024];
  sprintf(ofilename_dist,"Displ_%05d%s-ad%.2f.dat",
    n,halfshift?"-halfshift":"",accountedDist/mdtk::Ao);
  char ofilename_dist_per_atom[1024];
  sprintf(ofilename_dist_per_atom,"Displ_per_atom_%05d%s-ad%.2f.dat",
    n,halfshift?"-halfshift":"",accountedDist/mdtk::Ao);

//  gsl_histogram * h       = gsl_histogram_alloc (n);
  gsl_histogram * h_dist  = gsl_histogram_alloc (n);
  gsl_histogram * h_count = gsl_histogram_alloc (n);

  const Float minPolar = -180.0 +(halfshift?((360.0/n)/2.0):0.0);
  const Float maxPolar = +180.0 +(halfshift?((360.0/n)/2.0):0.0);

//  gsl_histogram_set_ranges_uniform (h,       minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (h_dist,  minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (h_count, minPolar, maxPolar);

//  Float atomsTotal = 0.0;

int coordsCut = 0;

  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
//      if (!td.clusterDynamics.isAtomInNMer(atraj.endCheckPoint,12)) continue;
      mdtk::Vector3D finCoords = atraj.endCheckPoint.coords;

      if ((atraj.endCheckPoint.coords-atraj.beginCheckPoint.coords).module() < accountedDist) continue;
      if (atraj.endCheckPoint.coords.z < -1.0*Ao) continue;

        for(size_t c = 0; c < 3; c++)
        if (PBC.X(c) < MDTK_MAX_PBC)
        {
          if (finCoords.X(c) < 0)         {finCoords.X(c) = 0.0;coordsCut++; }
          if (finCoords.X(c) >= PBC.X(c)) {finCoords.X(c) = PBC.X(c);coordsCut++; }
        }

      finCoords -= CenterOfCluster;

      Float polar = atan2(finCoords.y,finCoords.x)/mdtk::Deg;

      if (halfshift) if (polar < minPolar) polar += 360.0;

//      gsl_histogram_accumulate (h, polar, 1.0/Float(trajData.size())

// /Float(trajData[0].clusterDynamics.atomTrajectories.size())

//);
      gsl_histogram_accumulate (h_count, polar, 1.0);

      mdtk::Vector3D xdist = finCoords-atraj.beginCheckPoint.coords;
      xdist.z = 0.0;

      gsl_histogram_accumulate (h_dist, polar, xdist.module()/mdtk::Ao);
    }
  }
  
/*
  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
//    double lower, upper;
//    gsl_histogram_get_range (h, i, &lower, &upper);
 //   Float ang = (lower+upper)/2.0;
    gsl_histogram_set(gsl_histogram_get(h,i)/Float(trajData.size())/Float(trajData[0].clusterDynamics.atomTrajectories.size()));
  }  
*/
//  saveHistogram_polar(h, ofilename);
//  saveHistogram_polar(h_count, ofilename_count);


  gsl_histogram_scale(h_dist, 1.0/Float(trajData.size())/Float(trajData[0].clusterDynamics.atomTrajectories.size()));


  saveHistogram_polar(h_dist, ofilename_dist);

{
  gsl_histogram * h_count_clone = gsl_histogram_clone (h_count);

    for(int i = 0; i < n; i++)
    {
      double lower, upper;
      gsl_histogram_get_range (h_count_clone, i, &lower, &upper);
      if (gsl_histogram_get(h_count_clone, i) == 0.0)
        gsl_histogram_accumulate(h_count_clone, (lower+upper)/2.0, 1.0);
    }  

  gsl_histogram_div(h_dist, h_count_clone);
  gsl_histogram_free (h_count_clone);
  saveHistogram_polar(h_dist, ofilename_dist_per_atom);
}

  gsl_histogram_scale(h_count, 1.0/Float(trajData.size()) /* /Float(trajData[0].clusterDynamics.atomTrajectories.size())*/);
  saveHistogram_polar(h_count, ofilename);

//  gsl_histogram_free (h);
  gsl_histogram_free (h_dist);
  gsl_histogram_free (h_count);


{
  std::ofstream fo_co_cut("coordsCut.dat");
  fo_co_cut << Float(coordsCut)/trajData.size()/trajData[0].clusterDynamics.atomTrajectories.size() << "\n";
  fo_co_cut.close();
}


  yaatk::chdir("..");
}


void
StatPostProcess::printClusterDynamicsTotal() const
{
  yaatk::mkdir("_clusterdyn");
  yaatk::chdir("_clusterdyn");
  std::ofstream fo_global_plt("all.plt");
  mdtk::Float border = 5.0*mdtk::Ao;
  {
  using namespace mdtk;
  REQUIRE(trajData.size() > 0);
  }
  const mdtk::Vector3D& PBC = trajData[0].clusterDynamics.PBC;
  fo_global_plt << "reset\nset xrange [" << -border/mdtk::Ao << ":" << (PBC.x+border)/mdtk::Ao << "]\n";
  fo_global_plt <<        "set yrange [" << -border/mdtk::Ao << ":" << (PBC.y+border)/mdtk::Ao << "]\n"
  "set key off\n";
  fo_global_plt << "set output \""<< "all.plt" << ".png\"";
  fo_global_plt << FO_COMMON;
  fo_global_plt << "\nplot ";
  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    char local_plt_name[100];
    sprintf(local_plt_name,"%07lu.plt",traj);
    std::ofstream fo_local_plt(local_plt_name);
    fo_local_plt << "reset\nset xrange [" << -border/mdtk::Ao << ":" << (PBC.x+border)/mdtk::Ao << "]\n";
    fo_local_plt <<        "set yrange [" << -border/mdtk::Ao << ":" << (PBC.y+border)/mdtk::Ao << "]\n"
    "set key off\n";
    fo_local_plt << "set output \""<< local_plt_name << ".png\"";
    fo_local_plt << FO_COMMON;
    fo_local_plt << "\nplot ";
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      char traj_dat_fn[100];
      sprintf(traj_dat_fn,"-");
      fo_global_plt << "\'" << traj_dat_fn << "\' with lines";
      if (traj != trajData.size()-1 || atomIndex != td.clusterDynamics.atomTrajectories.size()-1)
        fo_global_plt << ", \\\n";
      else
        fo_global_plt << "\n";
      fo_local_plt << "\'" << traj_dat_fn << "\' with lines";
      if (atomIndex != td.clusterDynamics.atomTrajectories.size()-1)
        fo_local_plt << ", \\\n";
      else
        fo_local_plt << "\n";
    }
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      for(size_t i = 0; i < atraj.checkPoints.size(); i++)
      {
        fo_local_plt << atraj.checkPoints[i].coords.x/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.y/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.z/mdtk::Ao << "\n";
      }
      fo_local_plt << "e\n";
    }
    fo_local_plt.close();
  }  


  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      for(size_t i = 0; i < atraj.checkPoints.size(); i++)
      {
        fo_global_plt << atraj.checkPoints[i].coords.x/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.y/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.z/mdtk::Ao << "\n";
      }
      fo_global_plt << "e\n";
    }
  }  


  fo_global_plt << "\n";
  fo_global_plt.close();
  yaatk::chdir("..");
}  


void
StatPostProcess::printClusterDynamicsRSQ(bool xyOnly) const
{
if (!xyOnly)
{
  yaatk::mkdir("_clusterRSQ");
  yaatk::chdir("_clusterRSQ");
}
else
{
  yaatk::mkdir("_clusterRSQ_xy");
  yaatk::chdir("_clusterRSQ_xy");
}
  std::ofstream fo_global_plt("rsq.all.plt");
  {
  using namespace mdtk;
  REQUIRE(trajData.size() > 0);
  REQUIRE(trajData[0].clusterDynamics.atomTrajectories.size() > 0);
  }
  const mdtk::Vector3D& PBC = trajData[0].clusterDynamics.PBC;
/*
TRACE(PBC);
TRACE(PBC.X(0) < MDTK_MAX_PBC);
TRACE(PBC.X(1) < MDTK_MAX_PBC);
TRACE(PBC.X(2) < MDTK_MAX_PBC);
*/
  fo_global_plt << "set output \""<< "rsq.all.plt" << ".png\"";
  fo_global_plt << FO_COMMON;
  fo_global_plt << "\nplot \'-\' with lines\n";

//  std::vector<Float> RSQ(trajData[0].clusterDynamics.atomTrajectories.size());

  int t_size = trajData[0].clusterDynamics.atomTrajectories[0].checkPoints.size();
//  fo_global_plt << (t_size-t_size)*0.1 << " " << 0.0/trajData.size() << "\n";

//  TRACE(t_size);

  for(int t = t_size-1; t >= 0 ; t--)
  {
    Float RSQ = 0.0;
    for(size_t traj = 0; traj < trajData.size(); traj++)
    {
      const TrajData& td = trajData[traj];
//      TRACE(traj);
//      TRACE(td.trajDir);
      for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
      {
        const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];

        mdtk::Vector3D currentCoords;
        if (t < atraj.checkPoints.size())
          currentCoords = atraj.checkPoints[t].coords;
        else
        {
          currentCoords = atraj.checkPoints[atraj.checkPoints.size()-1].coords;
          TRACE("checkPoints overrun!!!");
          TRACE(t);
          TRACE(atraj.checkPoints.size());
        }

//      TRACE(atomIndex);
//        TRACE(atraj.checkPoints.size());

        for(size_t c = 0; c < 3; c++)
        if (PBC.X(c) < MDTK_MAX_PBC)
        {
          if (currentCoords.X(c) < 0)         {currentCoords.X(c) = 0.0; }
          if (currentCoords.X(c) >= PBC.X(c)) {currentCoords.X(c) = PBC.X(c); }
        }

        RSQ += SQR(currentCoords.x/mdtk::Ao-atraj.beginCheckPoint.coords.x/mdtk::Ao);
        RSQ += SQR(currentCoords.y/mdtk::Ao-atraj.beginCheckPoint.coords.y/mdtk::Ao);
        if (!xyOnly)
        RSQ += SQR(currentCoords.z/mdtk::Ao-atraj.beginCheckPoint.coords.z/mdtk::Ao);

//if (RSQ>1e10) {TRACE(t);TRACE(traj);TRACE(atomIndex);TRACE(currentCoords/mdtk::Ao);TRACE(atraj.beginCheckPoint.coords/mdtk::Ao);TRACE(atraj.checkPoints[t]);TRACE(atraj.beginCheckPoint);throw;}

      }
    }  
    fo_global_plt << (t_size-t-1)*0.1 << " " << RSQ/trajData.size() << "\n";
  }
  fo_global_plt << "e\n";

  fo_global_plt << "\n";
  fo_global_plt.close();

  yaatk::chdir("..");
}  


void
StatPostProcess::printProjectileDynamics() const
{
{
  yaatk::mkdir("_projectile");
  yaatk::chdir("_projectile");
}
  std::ofstream fo_global_plt("energy.all.plt");
  {
  using namespace mdtk;
  REQUIRE(trajData.size() > 0);
  REQUIRE(trajData[0].projectileDynamics.atomTrajectories.size() > 0);
  }
//  const mdtk::Vector3D& PBC = trajData[0].projectileDynamics.PBC;
  fo_global_plt << "reset\nset output \""<< "energy.all.plt" << ".png\"";
  fo_global_plt << FO_COMMON;
  fo_global_plt << "\nplot \'-\' with lines\n";

  const int n_overcount = 10;
  const int n = n_overcount + 10;
  gsl_histogram * h = gsl_histogram_alloc (n);
  Float c = 6.708*mdtk::Ao;
  const Float layer = c/2.0;
  gsl_histogram_set_ranges_uniform (h, -c/4.0-layer*n_overcount, -c/4.0+layer*(n-n_overcount));

  gsl_histogram * h_count = gsl_histogram_clone (h);

  int t_size = trajData[0].projectileDynamics.atomTrajectories[0].checkPoints.size();

  for(int t = t_size-1; t >= 0 ; t--)
  {
    Float energy = 0.0;
    for(size_t traj = 0; traj < trajData.size(); traj++)
    {
      const TrajData& td = trajData[traj];
      for(size_t atomIndex = 0; atomIndex < td.projectileDynamics.atomTrajectories.size(); atomIndex++)
      {
        const AtomTrajectory& atraj = td.projectileDynamics.atomTrajectories[atomIndex];

        const mdtk::Atom* currentAtom;
        if (t < atraj.checkPoints.size())
          currentAtom = &(atraj.checkPoints[t]);
        else
        {
          currentAtom = &(atraj.checkPoints[atraj.checkPoints.size()-1]);
          TRACE("checkPoints overrun!!! (Projectile)");
          TRACE(t);
          TRACE(atraj.checkPoints.size());
        }

        Float ecur = SQR(currentAtom->V.module())*currentAtom->M/2.0;
        energy += ecur;;

        if (t_size-t-1 > 2 && currentAtom->coords.z < 0)
        { 
          gsl_histogram_accumulate (h,       -30.0*mdtk::Ao, ecur);
          gsl_histogram_accumulate (h_count, -30.0*mdtk::Ao, 1.0);
        }
        else
        { 
          gsl_histogram_accumulate (h,       currentAtom->coords.z, ecur);
          gsl_histogram_accumulate (h_count, currentAtom->coords.z, 1.0);
        }
      }
    }  
    fo_global_plt << (t_size-t-1)*0.1 << " " << energy/mdtk::eV/trajData.size() << "\n";
  }
  fo_global_plt << "e\n";

  fo_global_plt << "\n";
  fo_global_plt.close();

//  gsl_histogram_scale(h, 1.0/trajData.size());
    for(int i = 0; i < n; i++)
    {
      double lower, upper;
      gsl_histogram_get_range (h_count, i, &lower, &upper);
      if (gsl_histogram_get(h_count, i) == 0.0)
        gsl_histogram_accumulate(h_count, (lower+upper)/2.0, 1.0);
    }  

  gsl_histogram_div  (h, h_count);


  {
    std::ofstream fo("energy_by_depth.all.plt");
    fo << "reset\nset output \""<< "energy_by_depth.all.plt" << ".png\"";
    fo << FO_COMMON;
    fo << "\nplot \'-\' with boxes\n";

    std::ofstream fod("energy_by_depth_delta.all.plt");
    fod << "reset\nset output \""<< "energy_by_depth_delta.all.plt" << ".png\"";
    fod << FO_COMMON;
    fod << "\nplot \'-\' with boxes\n";

    for(int i = 0; i < n; i++)
    {
      double lower, upper;
      gsl_histogram_get_range (h, i, &lower, &upper);
      fo << (lower+upper)/2.0/mdtk::Ao << " " << gsl_histogram_get(h,i)/mdtk::eV << std::endl;
    }  
    for(int i = 1; i < n; i++)
    {
      double lower, upper;
      gsl_histogram_get_range (h, i, &lower, &upper);
      fod << (lower+upper)/2.0/mdtk::Ao << " " << (gsl_histogram_get(h,i)-gsl_histogram_get(h,i-1))/mdtk::eV << std::endl;
    }  

    fo << "e\n";
    fo.close();    

    fod << "e\n";
    fod.close();    
  }
  

  gsl_histogram_free (h);
  gsl_histogram_free (h_count);

  yaatk::chdir("..");
}  




void
StatPostProcess::printProjectileStopping() const
{
{
  yaatk::mkdir("_projectileStopping");
  yaatk::chdir("_projectileStopping");
}

  std::vector<Float> depth;

  {
    for(size_t traj = 0; traj < trajData.size(); traj++)
    {
      const TrajData& td = trajData[traj];
      for(size_t atomIndex = 0; atomIndex < td.projectileDynamics.atomTrajectories.size(); atomIndex++)
      {
        const AtomTrajectory& atraj = td.projectileDynamics.atomTrajectories[atomIndex];
        depth.push_back(atraj.endCheckPoint.coords.z/mdtk::Ao);
//        Float ecur = SQR(atraj.checkPoints[t].V.module())*atraj.checkPoints[t].M/2.0;

      }
    }  
  }

REQUIRE(trajData.size()>0);
REQUIRE(trajData[0].projectileDynamics.atomTrajectories.size()>0);

  {
    depthHist2file("implants_by_depth.dat",depth, 1.0/trajData.size()/trajData[0].projectileDynamics.atomTrajectories.size());
  }  


  yaatk::chdir("..");
}  

void
StatPostProcess::printAtomTransitions() const
{

TRACE(ceil(1.1));
TRACE(ceil(0.1));
TRACE(ceil(-0.9));
TRACE(ceil(-1.9));

{
  yaatk::mkdir("_atomTransitions");
  yaatk::chdir("_atomTransitions");
}

  int shift = +4;
  int n = 25+shift;
  std::vector<std::vector<int> > tra;
  tra.resize(n);
  for(size_t i = 0; i < tra.size(); i++)
    tra[i].resize(n);

  Float c = 2.547*mdtk::Ao;//3.61*mdtk::Ao;
//  Float firstLayerBorder = -c/4.0;

int totalTransitionActual = 0;
int totalTransitionCandidates = 0;

int min1 = 10000;
int max1 = -10000;
int min2 = 10000;
int max2 = -10000;

  {
    for(size_t traj = 0; traj < trajData.size(); traj++)
    {
      const TrajData& td = trajData[traj];
//      TRACE(td.trans.translations.size());
      for(size_t tIndex = 0; tIndex < td.trans.translations.size(); tIndex++)
      {
        const Translation& t = td.trans.translations[tIndex];

        if (fabs(t.start.coords.z-t.end.coords.z) >= (c/2.0)*0.8)
          totalTransitionActual++;

        totalTransitionCandidates++;

        int layer_start;
        int layer_end;
        if (t.start.coords.z < -c/4.0) 
        {
          layer_start = /*0;*/-shift;
        }
        else
        {
          layer_start = ceil((t.start.coords.z + c/4.0)/(c/2.0));
          REQUIRE(layer_start > 0);
        }
        {
          layer_end = ceil((t.end.coords.z + c/4.0)/(c/2.0));
        }
        if (layer_end+shift < 0)
        {
          layer_end = -shift;
        }
        if (layer_end+shift >= n)
        {
          TRACE("Abnormal layer_end!!!");
          TRACE(layer_end+shift);
          TRACE(traj);
          TRACE(t.start.globalIndex);
          TRACE(t.start.coords.z/mdtk::Ao);
          TRACE(t.end.globalIndex);
          TRACE(t.end.coords.z/mdtk::Ao);
          layer_end = n-1-shift;
        }
//        TRACE(layer_start);
//        TRACE(layer_end);
        REQUIRE(layer_start+shift >= 0);
        REQUIRE(layer_end+shift >= 0);
        REQUIRE(layer_start+shift < n);
        REQUIRE(layer_end+shift < n);

        tra[layer_start+shift][layer_end+shift]++;

        if (layer_start+shift < min1) min1 = layer_start+shift;
        if (layer_start+shift > max1) max1 = layer_start+shift;
        if (layer_end+shift > max2) max2 = layer_end+shift;
        if (layer_end+shift < min2) min2 = layer_end+shift;
      }
    }  
  }

  TRACE(totalTransitionActual);
  TRACE(totalTransitionCandidates);

int totalTransitionSum = 0;
int diaTransitionSum = 0;

  std::ofstream fos("k.sh");
  fos << "konwert utf8-koi8u ./trans.plt -o ./trans-koi8u.plt\ngnuplot ./trans-koi8u.plt\n"
         "pstopnm -xborder=0 -yborder=0 -xsize=1000 -ysize=1000 -portrait trans.eps\n"
         "pnmtopng trans.eps001.ppm > trans.eps001.png\n";
  fos.close();

  system("chmod +x ./k.sh");

  std::ofstream foplt("trans.plt");
  foplt << "reset\n";
  foplt << "set xrange [" << -shift+1 << ":"<< "*" << "]\n";
  foplt << "set yrange [0:*]\n";
  foplt << "set xtics 1\n";
  foplt << "set grid xtics\n";
  foplt << "set xlabel \"Номер шару атомів\"\n";
  foplt << "set ylabel \"Частота переміщення\"\n";
  foplt << "set encoding koi8u\n";
  foplt << "set output \"trans.eps\"\n";
//  foplt << "set terminal png giant font \"c:/windows/fonts/arialuni.ttf\" 28 size 1000,1000 enhanced \\\n xffffff x000000 x000000 \\\n x000000 x000000 x000000 x000000 \\\n x000000 x000000 x000000 x000000\n";
  foplt << "set terminal postscript eps size 8cm, 10cm \"Arial,16\" enhanced\n";
  foplt << "plot \\\n";
{
  for(size_t i = /*4*/min1; i <= max1; i++)
  {
    for(size_t j = min2; j <= max2; j++)
    {
      if (tra[i][j] != 0.0)
      {
      int indi = (int(i)-shift);
      int indj = (int(j)-shift);
        foplt << "\'-\' notitle with lines lw 2 lt "; // smooth csplines lines
        if (indi > indj) foplt << "1";
        if (indi < indj) foplt << "2";
        if (indi == indj) foplt << "3";
        foplt << ", \\";
        foplt << "\n";
      }
    }
  }
  foplt << "0 notitle with lines lw 1\n";
}

  std::ofstream fo("trans_table.dat");
  fo << setw(5) << 0 << " ";
  for(size_t i = min2; i <= max2; i++)
    fo << setw(5) << int(i)-shift << " ";
  fo << "\n";
  for(size_t i = /*4*/min1; i <= max1; i++)
  {
    diaTransitionSum += tra[i][i];
    fo << setw(5) << int(i)-shift << " ";
    for(size_t j = min2; j <= max2; j++)
    {
      char s[100];
      Float val = Float(tra[i][j])/trajData.size()/13.0;
      if (val != 0.0)
        sprintf(s,"%5.2f ",val);
      else
        sprintf(s,"%5g ",val);
      fo << s;
      totalTransitionSum += tra[i][j];

      if (tra[i][j] != 0.0)
      {
      Float delta = 0;//((int(i)-shift)>(int(j)-shift))?0.1:0.2;
      int indmax = ((int(i)-shift)>(int(j)-shift))?(int(i)-shift):(int(j)-shift);
      int indmin = ((int(i)-shift)>(int(j)-shift))?(int(j)-shift):(int(i)-shift);

//      if (indmin == -shift) indmin = indmax-2*(indmax-indmin);

//      foplt << indmin-0.1          << " " << "0" << "\n";
      foplt << indmin+delta              << "     " << "0" << "\n";
//      foplt << indmin+delta              << " "     << val << "\n";
      foplt << (indmin+indmax)/2.0 << " "     << val << "\n";
//      foplt << indmax-delta              << " "     << val << "\n";
      foplt << indmax-delta              << "     " << "0" << "\n";
//      foplt << indmax+0.1          << " " << "0" << "\n";
      foplt << "e\n";
      }
    }
    fo << "\n";
  }
  fo << "\n";
  fo.close();

  foplt.close();

  TRACE(totalTransitionSum-diaTransitionSum);
  TRACE(totalTransitionSum);
  TRACE("----");
  TRACE(diaTransitionSum);

REQUIRE(trajData.size()>0);

  yaatk::chdir("..");
}  





void
StatPostProcess::printClusterFragments() const
{
  yaatk::mkdir("fragments");
  yaatk::chdir("fragments");

TRACE(trajData.size());

  {
  using namespace mdtk;
  REQUIRE(trajData.size() > 0);
  }

  std::vector<Float> freqOfFragments(trajData[0].clusterDynamics.atomTrajectories.size()+1);

TRACE(freqOfFragments.size());

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
TRACE(trajIndex);
TRACE(td.clusterDynamics.fragments.size());
    for(size_t frIndex = 0; frIndex < td.clusterDynamics.fragments.size(); frIndex++)
    {
      const Fragment& fr = td.clusterDynamics.fragments[frIndex];
TRACE(frIndex);
TRACE(fr.atoms.size());
      (freqOfFragments[fr.atoms.size()]) += 1.0;
    }
  }


{
  {
    std::ofstream fo("fragfreqs.dat.plt");
    fo << "reset\nset yrange [0:*]\nplot \'" << "fragfreqs.dat" << "\' with impulses\n";
    fo.close();
  }  

  {
    std::ofstream fo("fragfreqs.dat");
    for(size_t massIndex = 0; massIndex < freqOfFragments.size(); massIndex++)
    {
      fo << massIndex << " " 
         << freqOfFragments[massIndex]/trajData.size() << std::endl;
    } 
    fo.close();
  }
  {
    std::ofstream fo("fragcounts.dat");
    for(size_t massIndex = 0; massIndex < freqOfFragments.size(); massIndex++)
    {
      fo << massIndex << " " 
         << freqOfFragments[massIndex] << std::endl;
    } 
    fo.close();
  }
}

/*
  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      for(size_t i = 0; i < atraj.checkPoints.size(); i++)
      {
        fo_global_plt << atraj.checkPoints[i].coords.x/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.y/mdtk::Ao << " " 
                    << atraj.checkPoints[i].coords.z/mdtk::Ao << "\n";
      }
    }
  }  
*/
  yaatk::chdir("..");
}  


/*
void
StatPostProcess::printClusterDynamicsTotal() const
{
  yaatk::mkdir("_clusterdyn");
  yaatk::chdir("_clusterdyn");
  std::ofstream fo_global_plt("all.plt");
  mdtk::Float border = 5.0*mdtk::Ao;
  {
  using namespace mdtk;
  REQUIRE(trajData.size() > 0);
  }
  const mdtk::Vector3D& PBC = trajData[0].clusterDynamics.PBC;
  fo_global_plt << "reset\nset xrange [" << -border << ":" << PBC.x+border << "]\n";
  fo_global_plt <<        "set yrange [" << -border << ":" << PBC.y+border << "]\n"
  "set key off\n";
  fo_global_plt << "set output \""<< "all.plt" << ".png\"";
  fo_global_plt << FO_COMMON;
  fo_global_plt << "\nplot ";
  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = trajData[traj];
    char local_plt_name[100];
    sprintf(local_plt_name,"%07lu.plt",traj);
    std::ofstream fo_local_plt(local_plt_name);
    fo_local_plt << "reset\nset xrange [" << -border << ":" << PBC.x+border << "]\n";
    fo_local_plt <<        "set yrange [" << -border << ":" << PBC.y+border << "]\n"
    "set key off\n";
    fo_local_plt << "set output \""<< local_plt_name << ".png\"";
    fo_local_plt << FO_COMMON;
    fo_local_plt << "\nplot ";
    for(size_t atomIndex = 0; atomIndex < td.clusterDynamics.atomTrajectories.size(); atomIndex++)
    {
      char traj_dat_fn[100];
      sprintf(traj_dat_fn,"%07lu.%07lu.dat",traj,atomIndex);
      fo_global_plt << "\'" << traj_dat_fn << "\' with lines";
      if (traj != trajData.size()-1 || atomIndex != td.clusterDynamics.atomTrajectories.size()-1)
        fo_global_plt << ", \\\n";
      else
        fo_global_plt << "";
      fo_local_plt << "\'" << traj_dat_fn << "\' with lines";
      if (atomIndex != td.clusterDynamics.atomTrajectories.size()-1)
        fo_local_plt << ", \\\n";
      else
        fo_local_plt << "";
      std::ofstream fo_traj_dat(traj_dat_fn);
      const AtomTrajectory& atraj = td.clusterDynamics.atomTrajectories[atomIndex];
      for(size_t i = 0; i < atraj.checkPoints.size(); i++)
        fo_traj_dat << atraj.checkPoints[i].coords.x << " " 
                    << atraj.checkPoints[i].coords.y << " " 
                    << atraj.checkPoints[i].coords.z << "\n";
      fo_traj_dat.close();
    }
    fo_local_plt.close();
  }  
  fo_global_plt << "\n";
  fo_global_plt.close();
  yaatk::chdir("..");
}  
*/

void
StatPostProcess::printClusterDynamics(size_t /*trajIndex*/) const
{
/*
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    cout << "Molecule #" << mi << " : ";
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      cout << "( Z = " << atom.Z/mdtk::e << " , gi = " << atom.globalIndex << " ) ";
    } 
    cout << "\n            formation time = " << td.molecules[mi].formationTime/mdtk::ps << " ps"
         << ", escape time = " << td.molecules[mi].escapeTime/mdtk::ps << " ps";
    cout << std::endl;
  }  
*/
}  


void
StatPostProcess::printEmphasizedTotal() const
{
  system("mkdir emphasized");
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].molecules.size() > 0)
  {
    printEmphasized(traj);
  }  
}

void
StatPostProcess::printEmphasized(size_t trajIndex) const
{
  char fileName[1024];
  std::string stateFilename = yaatk::extractLastItem(trajData[trajIndex].trajDir);//+"_";
  replace(stateFilename.begin(),stateFilename.end(),'\\','.');
  replace(stateFilename.begin(),stateFilename.end(),'/','.');
  replace(stateFilename.begin(),stateFilename.end(),':','.');
  sprintf(fileName,"emphasized"DIR_DELIMIT_STR"%s.emphasized",stateFilename.c_str());
//  sprintf(fileName,"emphasized.txt.%03d",trajIndex);
//  TRACE(fileName);
  std::ofstream fo(fileName);
  
  std::vector<size_t> atomIndexes;

  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      atomIndexes.push_back(atom.globalIndex);
    } 
  }  

  fo << atomIndexes.size() << "\n";
  for(size_t ai = 0; ai < atomIndexes.size(); ai++)
  {
    fo << atomIndexes[ai] << "\n";
  } 
  
  fo.close();
}  

void
StatPostProcess::spottedTotalMDE() const
{
  std::ofstream fo("spots");
  fo << getYieldSum(&ProcessAll) << std::endl;
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].molecules.size() > 0)
  {
  const TrajData& td = trajData[traj];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      fo << atom.V << "\n";
    } 
  }  
  }  
  fo.close();
}  


#define MDEPP_BYDEPTH_DIR "_by_depth"

void
StatPostProcess::spottedByDepth() const
{
  system("mkdir "MDEPP_BYDEPTH_DIR);

  chdir(MDEPP_BYDEPTH_DIR);

  std::vector<Float> depth;
  std::vector<Float> depth_H;
  std::vector<Float> depth_C;
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].molecules.size() > 0)
  {
    const TrajData& td = trajData[traj];

    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
      {
//        const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
        const mdtk::Atom& atom_init = td.molecules[mi].atoms_init[ai];
        depth.push_back(atom_init.coords.z/mdtk::Ao);
        if (atom_init.ID == mdtk::H_EL)
          depth_H.push_back(atom_init.coords.z/mdtk::Ao);
        if (atom_init.ID == mdtk::C_EL)
          depth_C.push_back(atom_init.coords.z/mdtk::Ao);
        
//        fo << atom.V << "\n";
      } 
    }  
  }  
  {
    saveDepth(depth,"spotted_depths_unsorted.dat");
    sort(depth.begin(),depth.end());
    saveDepth(depth,"spotted_depths.dat");
  }  
  {
    saveDepth(depth_C,"spotted_depths_C_unsorted.dat");
    sort(depth_C.begin(),depth_C.end());
    saveDepth(depth_C,"spotted_depths_C.dat");
  }  
  {
    saveDepth(depth_H,"spotted_depths_H_unsorted.dat");
    sort(depth_H.begin(),depth_H.end());
    saveDepth(depth_H,"spotted_depths_H.dat");
  }  
  {
    depthHist2file("spots_by_depth.dat",depth);
    depthHist2file("spots_by_depth_H.dat",depth_H);
    depthHist2file("spots_by_depth_C.dat",depth_C);
  }  

  chdir("..");
}  

void saveMassSpectrum(std::vector<MassCount>& massSpectrum,const char *filename,Float scale)
{
  {
    std::ofstream fo((std::string(filename)+".plt").c_str());
    fo << "reset\nset yrange [0:*]\nplot \'" << filename << "\' with impulses\n";
    fo.close();
  }  

    std::ofstream fo(filename);
    for(size_t massIndex = 0; massIndex < massSpectrum.size(); massIndex++)
    {
      fo << massSpectrum[massIndex].amuMass << " " 
         << massSpectrum[massIndex].count*scale << std::endl;
    } 
    fo.close();
}

void saveMassSpectrumWoProjNorm(std::vector<MassCount>& massSpectrum,const char *filename, int trajnum)
{
  {
    std::ofstream fo((std::string(filename)+".plt").c_str());
    fo << "reset\nset yrange [0:*]\nplot \'" << filename << "\' with impulses\n";
    fo.close();
  }  

    std::ofstream fo(filename);
    for(size_t massIndex = 0; massIndex < massSpectrum.size(); massIndex++)
    {
      if (massSpectrum[massIndex].amuMass == 40.0) continue;
      fo << massSpectrum[massIndex].amuMass << " " 
         << massSpectrum[massIndex].count/Float(trajnum) << std::endl;
    } 
    fo.close();
}

/*
void
StatPostProcess::buildMassSpectrum_() const
{
  std::vector<MassCount>  massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    Float moleculeMass = td.molecules[mi].getAMUMass();
    bool massAlreadyAccounted = false;
    for(size_t accountIndex = 0; accountIndex < massSpectrum.size(); accountIndex++)
    {
      if (massSpectrum[accountIndex].amuMass == moleculeMass)
      {
        massSpectrum[accountIndex].count++;
        massAlreadyAccounted = true;
        break;
      }  
    }
    if (!massAlreadyAccounted)
    {
      massSpectrum.push_back(MassCount(moleculeMass));      
      massSpectrum[massSpectrum.size()-1].count = 1;
    }  
  }
  }    

  { 
    saveMassSpectrum(massSpectrum,"mass_spectrum_unsorted.dat");
  }  
  {
    sort(massSpectrum.begin(),massSpectrum.end());
  }  
  { 
    saveMassSpectrum(massSpectrum,"mass_spectrum.dat");
    saveMassSpectrumWoProjNorm(massSpectrum,"mass_spectrum_wo_proj_norm.dat",trajData.size());
  }  
}
*/
void
StatPostProcess::buildMassSpectrum(FProcessMolecule fpm) const
{
  std::vector<MassCount>  massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    if (!fpm(td.molecules[mi])) continue;
    Float moleculeMass = td.molecules[mi].getAMUMass();
    bool massAlreadyAccounted = false;
    for(size_t accountIndex = 0; accountIndex < massSpectrum.size(); accountIndex++)
    {
      if (massSpectrum[accountIndex].amuMass == moleculeMass)
      {
        massSpectrum[accountIndex].count++;
        massAlreadyAccounted = true;
        break;
      }  
    }
    if (!massAlreadyAccounted)
    {
      massSpectrum.push_back(MassCount(moleculeMass));      
      massSpectrum[massSpectrum.size()-1].count = 1;
    }  
  }
  }    

  { 
    saveMassSpectrum(massSpectrum,"mass_spectrum_unsorted.dat",1.0/trajData.size());
  }  
  {
    sort(massSpectrum.begin(),massSpectrum.end());
  }  
  { 
    saveMassSpectrum(massSpectrum,"mass_spectrum.dat",1.0/trajData.size());
//    saveMassSpectrumWoProjNorm(massSpectrum,"mass_spectrum_wo_proj_norm.dat",trajData.size());
  }  
}


void saveSpottedByMass(std::vector<MassSpotted>& massSpectrum,const char *filename)
{
    std::ofstream fo(filename);
    fo << massSpectrum.size() << std::endl;
    for(size_t massIndex = 0; massIndex < massSpectrum.size(); massIndex++)
    {
      fo << massSpectrum[massIndex].getAMUMass() << std::endl;
      MassSpotted& spotted = massSpectrum[massIndex];
      fo << spotted.species.size() << std::endl;
      for(size_t mi = 0; mi < spotted.species.size(); mi++)
      {
        fo << spotted.species[mi].getVelocity() << std::endl;
        spotted.species[mi].printGlobalIndexes(fo);
      }  
    } 
    fo.close();
}  

void
StatPostProcess::spottedTotalByMass() const
{
  std::vector<MassSpotted>  massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const Molecule& mol = td.molecules[mi];
    Float moleculeMass = mol.getAMUMass();
    bool massAlreadyAccounted = false;
    for(size_t accountIndex = 0; accountIndex < massSpectrum.size(); accountIndex++)
    {
      if (massSpectrum[accountIndex].getAMUMass() == moleculeMass)
      {
        massSpectrum[accountIndex].species.push_back(mol);
        massAlreadyAccounted = true;
        break;
      }  
    }
    if (!massAlreadyAccounted)
    {
      massSpectrum.push_back(MassSpotted(mol));      
    }  
  }
  }    
  { 
    saveSpottedByMass(massSpectrum,"spots_by_mass_unsorted.dat");
  }  
  {
    sort(massSpectrum.begin(),massSpectrum.end());
  }  
  { 
    saveSpottedByMass(massSpectrum,"spots_by_mass.dat");
  }  
}

void
StatPostProcess::setSpottedDistanceFromInit()// const
{
  using mdtk::Exception;    

  REQUIRE(trajData.size() > 0);

  std::string mde_init_filename = trajData[0].trajDir+"mde_init";

  mdtk::SimLoop* mde_init = new mdtk::SimLoop();    
  mde_init->allowToFreePotentials = true;
  mde_init->allowToFreeAtoms = true;
  setupPotentials(*mde_init);
  yaatk::text_ifstream fi(mde_init_filename.c_str()); 
mde_init->initNLafterLoading = false;
  mde_init->loadFromStream(fi/*,false*/);
  setTags(mde_init);
mde_init->initNLafterLoading = true;
  fi.close(); 
  mdtk::AtomsContainer& mde_init_atoms = mde_init->atoms_;
  Float minInitZ = 1000000.0*mdtk::Ao;
  for(size_t ai = 0; ai < mde_init_atoms.size(); ai++)
  {
    const mdtk::Atom& init_atom = *(mde_init_atoms[ai]);
    if (init_atom.coords.z < minInitZ && !isProjectileAtom(init_atom))
      minInitZ = init_atom.coords.z;
  }  
  SPOTTED_DISTANCE = minInitZ - 0.05*mdtk::Ao;
  SPOTTED_DISTANCE = -4.55*mdtk::Ao;
  SPOTTED_DISTANCE = -0.0*mdtk::Ao;
  TRACE(SPOTTED_DISTANCE/mdtk::Ao);
  delete mde_init;
}  

std::string
StatPostProcess::buildAtomByEnergy(const Float energyStep, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atom_by_energy_%.2f.dat", energyStep);

  const Float minEnergy =     0.0;
  const Float maxEnergy =  1000.0;
  const int n = (maxEnergy-minEnergy)/(energyStep);

  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, minEnergy, maxEnergy);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        gsl_histogram_accumulate (h, energy/mol.atoms.size(), mol.atoms.size());
      }
    }    

  }  

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

void
StatPostProcess::histEnergyByPolar(gsl_histogram* h, bool byAtom, FProcessMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
//        TRACE(polar);
        if (!byAtom)
          gsl_histogram_accumulate (h, polar, energy /* /mol.atoms.size()*/);
        else
          gsl_histogram_accumulate (h, polar, energy /mol.atoms.size());
      }
    }    
  }  
}  

void
StatPostProcess::histAtomsCountByPolar(gsl_histogram* h, FProcessMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
/*
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
*/
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
        gsl_histogram_accumulate (h, polar, mol.atoms.size()       /Float(trajData.size())    ); //!!!!!!!!!!
      }
    }    
  }  
}  


void
StatPostProcess::histEnergyByPolarByAtomsInRange(gsl_histogram* h, FProcessMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  int n = gsl_histogram_bins(h);

  gsl_histogram * hEnergy = gsl_histogram_alloc (n);
  gsl_histogram * hCount  = gsl_histogram_alloc (n);

  gsl_histogram_set_ranges_uniform (hEnergy, minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (hCount, minPolar, maxPolar);


  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
//        TRACE(polar);
        gsl_histogram_accumulate (hEnergy, polar, energy /* /mol.atoms.size()*/);
        gsl_histogram_accumulate (hCount, polar, mol.atoms.size());
      }
    }    
  }  

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);
  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    if (gsl_histogram_get(hCount,i) > 0.0)
    gsl_histogram_accumulate(h,
           (lower+upper)/2.0,
           gsl_histogram_get(hEnergy,i)/gsl_histogram_get(hCount,i));
  }  

  gsl_histogram_free (hEnergy);
  gsl_histogram_free (hCount);
}  

std::string
StatPostProcess::buildEnergyByPolar(const int n, bool byAtom, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_polar_%s%05d.dat", byAtom?"by_atom_":"",n);
  
  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByPolar(h, byAtom, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

std::string
StatPostProcess::buildAtomsCountByPolar(const int n, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atomsCount_by_polar_%05d.dat",
    n);
  
  gsl_histogram * h = gsl_histogram_alloc (n);

  histAtomsCountByPolar(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

std::string
StatPostProcess::buildEnergyByPolarByAtomsInRange(const int n, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_polar_by_atomsInRange_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByPolarByAtomsInRange(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

void
StatPostProcess::histEnergyByAzimuth(gsl_histogram *h, bool byAtom, FProcessMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(v.y,v.x)/mdtk::Deg;
//        TRACE(polar);
        if (!byAtom)
          gsl_histogram_accumulate (h, polar, energy /* /mol.atoms.size()*/);
        else
          gsl_histogram_accumulate (h, polar, energy  /mol.atoms.size());
      }
    }    
  }  
}  

void
StatPostProcess::histAtomsCountByAzimuth(gsl_histogram *h,  FProcessMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
/*
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
*/
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(v.y,v.x)/mdtk::Deg;
        gsl_histogram_accumulate (h, polar, mol.atoms.size()   /Float(trajData.size())   );
      }
    }    
  }  
}  


void
StatPostProcess::histEnergyByAzimuthByAtomsInRange(gsl_histogram *h, FProcessMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  int n = gsl_histogram_bins(h);

  gsl_histogram * hEnergy = gsl_histogram_alloc (n);
  gsl_histogram * hCount  = gsl_histogram_alloc (n);

  gsl_histogram_set_ranges_uniform (hEnergy, minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (hCount, minPolar, maxPolar);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if (!fpm(mol)) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        mdtk::Vector3D v = mol.getVelocity();
        Float polar = atan2(v.y,v.x)/mdtk::Deg;
//        TRACE(polar);
        gsl_histogram_accumulate (hEnergy, polar, energy /* /mol.atoms.size()*/);
        gsl_histogram_accumulate (hCount, polar, mol.atoms.size());
      }
    }    
  }  

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);
  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    if (gsl_histogram_get(hCount,i) > 0.0)
    gsl_histogram_accumulate(h,
           (lower+upper)/2.0,
           gsl_histogram_get(hEnergy,i)/gsl_histogram_get(hCount,i));
  }  

  gsl_histogram_free (hEnergy);
  gsl_histogram_free (hCount);

}  

std::string
StatPostProcess::buildEnergyByAzimuth(const int n, bool byAtom, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_azimuth_%s%05d.dat",
    byAtom?"by_atom_":"",n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByAzimuth(h, byAtom, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

std::string
StatPostProcess::buildAtomsCountByAzimuth(const int n, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atomsCount_by_azimuth_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histAtomsCountByAzimuth(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  


std::string
StatPostProcess::buildEnergyByAzimuthByAtomsInRange(const int n, FProcessMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_azimuth_by_atomsInRange_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByAzimuthByAtomsInRange(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}  

/*
#define MDEPP_ANGULAR_DIR "_angular"

void
StatPostProcess::buildAngular() const
{
  system("mkdir "MDEPP_ANGULAR_DIR);

//  system("cd "MDEPP_ANGULAR_DIR);
  chdir(MDEPP_ANGULAR_DIR);

  bool accountSputtered = false;
  do{
    bool accountBackScattered = false;
    do{
//      if (!accountSputtered || !accountBackScattered)
      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!%s%sAtomByEnergy.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo(pltFileName);
        fo << "reset\nset xrange [*:20]\n";
//        for(Float energyStep = 0.05; energyStep <= 10.0; energyStep += 0.05)  
        for(Float energyStep = 0.1; energyStep <= 10.0; energyStep += 0.1)  
        {
          datFileName = 
            buildAtomByEnergy(energyStep,accountSputtered,accountBackScattered);
          fo << "reset\nset xrange [*:20]\nplot \'" << datFileName << "\' with histeps\n";
          fo << "pause -1 \"Press Return\" \n";
        }  
        fo.close();
      }  

      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!%s%sEnergyByPolar.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo(pltFileName);
        sprintf(pltFileName,"!!!%s%sEnergyByPolarByAtom.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo_by_atom(pltFileName);
        sprintf(pltFileName,"!!!%s%sEnergyByPolarByAtomsInRange.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo_by_range(pltFileName);
        fo << "reset\n";
        fo_by_atom << "reset\n";
        fo_by_range << "reset\n";
//        for(int n = 1; n <= 180; n++)
        for(int n = 1; n <= 90; n++)
        {
          datFileName = 
            buildEnergyByPolar(n,false,accountSputtered,accountBackScattered);
          fo         << "plot \'" << datFileName << "\' with histeps\n";
          fo         << "pause -1 \"Press Return\" \n";
          datFileName = 
            buildEnergyByPolar(n,true,accountSputtered,accountBackScattered);
          fo_by_atom << "plot \'" << datFileName << "\' with histeps\n";
          fo_by_atom << "pause -1 \"Press Return\" \n";
          datFileName = 
            buildEnergyByPolarByAtomsInRange(n,accountSputtered,accountBackScattered);
          fo_by_range << "plot \'" << datFileName << "\' with histeps\n";
          fo_by_range << "pause -1 \"Press Return\" \n";
        }
        fo.close();
        fo_by_atom.close();
        fo_by_range.close();
      }    

      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!%s%sEnergyByAzimuth.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo(pltFileName);
        sprintf(pltFileName,"!!!%s%sEnergyByAzimuthByAtom.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo_by_atom(pltFileName);
        sprintf(pltFileName,"!!!%s%sEnergyByAzimuthByAtomsInRange.plt",
          accountSputtered?"S":"_", accountBackScattered?"B":"_");
        std::ofstream fo_by_range(pltFileName);
        fo << "reset\n";
        fo_by_atom << "reset\n";
        fo_by_range << "reset\n";
        for(int n = 1; n <= 180; n++)
//        for(int n = 1; n <= 360; n++)
//        for(int n = 1; n <= 360*2; n++)
        {
          datFileName = 
            buildEnergyByAzimuth(n,false,accountSputtered,accountBackScattered);
          fo         << "plot \'" << datFileName << "\' with histeps\n";
          fo         << "pause -1 \"Press Return\" \n";
          datFileName = 
            buildEnergyByAzimuth(n,true,accountSputtered,accountBackScattered);
          fo_by_atom << "plot \'" << datFileName << "\' with histeps\n";
          fo_by_atom << "pause -1 \"Press Return\" \n";
          datFileName = 
            buildEnergyByAzimuthByAtomsInRange(n,accountSputtered,accountBackScattered);
          fo_by_range << "plot \'" << datFileName << "\' with histeps\n";
          fo_by_range << "pause -1 \"Press Return\" \n";
        }  
        fo.close();
        fo_by_atom.close();
        fo_by_range.close();
      }  

      if (accountBackScattered) break;
      accountBackScattered = true;
    }while (1);
    if (accountSputtered) break;
    accountSputtered = true;
  }while (1);

//  system("cd ""..");
  chdir("..");
}  
*/

#define MDEPP_ANGULAR_DIR2 "_angular2"

void
StatPostProcess::buildAngular2(FProcessMolecule fpm) const
{
  system("mkdir "MDEPP_ANGULAR_DIR2);

//  system("cd "MDEPP_ANGULAR_DIR);
  chdir(MDEPP_ANGULAR_DIR2);
/*
  bool accountSputtered = false;
  do{
    bool accountBackScattered = false;
    do{
      if (accountSputtered || accountBackScattered)
*/
{
      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!AtomByEnergy.plt");
        std::ofstream fo(pltFileName);
        fo << "reset\nset xrange [0:10]\nset yrange [0:*]\nset format y \"%10g\"\n";
//        for(Float energyStep = 0.05; energyStep <= 10.0; energyStep += 0.05)  
        for(Float energyStep = 0.1; energyStep <= 10.0; energyStep += 0.1)  
        {
          datFileName = "-";
          fo << "reset\nset xrange [0:10]\nset yrange [0:*]\nset format y \"%10g\"\nplot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildAtomByEnergy(energyStep,fpm);
          fo << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo << tempStr << "\n";}
            fi.close();
          }  
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo << "e\n";
          fo << "pause -1 \"Press Return\" \n";
        }  
        fo.close();
      }  

      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!EnergyByPolar.plt");
        std::ofstream fo(pltFileName);
        sprintf(pltFileName,"!!!AtomsCountByPolar.plt");
        std::ofstream fo_count(pltFileName);
        sprintf(pltFileName,"!!!EnergyByPolarByAtom.plt");
        std::ofstream fo_by_atom(pltFileName);
        sprintf(pltFileName,"!!!EnergyByPolarByAtomsInRange.plt");
        std::ofstream fo_by_range(pltFileName);
        fo << "reset\nset xrange [0:90]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_count << "reset\nset xrange [0:90]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_by_atom << "reset\nset xrange [0:90]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_by_range << "reset\nset xrange [0:90]\nset yrange [0:*]\nset format y \"%10g\"\n";
//        for(int n = 1; n <= 180; n++)
        for(int n = 1; n <= 90; n++)
        {
          datFileName = "-";
          fo         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByPolar(n,false,fpm);
          fo << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo << tempStr << "\n";}
            fi.close();
          }  
if (n!=9)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo << "e\n";
          fo         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_count         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildAtomsCountByPolar(n,fpm);
          fo_count << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_count << tempStr << "\n";}
            fi.close();
          }  
if (n!=9)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_count << "e\n";
          fo_count         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_by_atom         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByPolar(n,true,fpm);
          fo_by_atom << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_by_atom << tempStr << "\n";}
            fi.close();
          }  
if (n!=9)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_by_atom << "e\n";
          fo_by_atom         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_by_range         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByPolarByAtomsInRange(n,fpm);
          fo_by_range << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_by_range << tempStr << "\n";}
            fi.close();
          }  
if (n!=9)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_by_range << "e\n";
          fo_by_range         << "pause -1 \"Press Return\" \n";
        }
        fo.close();
        fo_count.close();
        fo_by_atom.close();
        fo_by_range.close();
      }    

      {
        std::string datFileName;
        char pltFileName[1000];
        sprintf(pltFileName,"!!!EnergyByAzimuth.plt");
        std::ofstream fo(pltFileName);
        sprintf(pltFileName,"!!!AtomsCountByAzimuth.plt");
        std::ofstream fo_count(pltFileName);
        sprintf(pltFileName,"!!!EnergyByAzimuthByAtom.plt");
        std::ofstream fo_by_atom(pltFileName);
        sprintf(pltFileName,"!!!EnergyByAzimuthByAtomsInRange.plt");
        std::ofstream fo_by_range(pltFileName);
        fo << "reset\nset xrange [-180:180]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_count << "reset\nset xrange [-180:180]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_by_atom << "reset\nset xrange [-180:180]\nset yrange [0:*]\nset format y \"%10g\"\n";
        fo_by_range << "reset\nset xrange [-180:180]\nset yrange [0:*]\nset format y \"%10g\"\n";
//        for(int n = 1; n <= 180; n++)
        for(int n = 1; n <= 360; n++)
//        for(int n = 1; n <= 360*2; n++)
        {
          datFileName = "-";
          fo         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByAzimuth(n,false,fpm);
          fo << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo << tempStr << "\n";}
            fi.close();
          }  
if (n!=36)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo << "e\n";
          fo         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_count         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildAtomsCountByAzimuth(n,fpm);
          fo_count << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_count << tempStr << "\n";}
            fi.close();
          }  
if (n!=36)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_count << "e\n";
          fo_count         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_by_atom         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByAzimuth(n,true,fpm);
          fo_by_atom << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_by_atom << tempStr << "\n";}
            fi.close();
          }  
if (n!=36)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_by_atom << "e\n";
          fo_by_atom         << "pause -1 \"Press Return\" \n";

          datFileName = "-";
          fo_by_range         << "plot \'" << datFileName << "\' with histeps";
          datFileName = 
            buildEnergyByAzimuthByAtomsInRange(n,fpm);
          fo_by_range << " title \"" << datFileName << "\"\n";
          {
            char tempStr[1000];
            std::ifstream fi(datFileName.c_str());
            while(fi.getline(tempStr, 1000-1, '\n')) {fo_by_range << tempStr << "\n";}
            fi.close();
          }  
if (n!=36)
          std::remove(datFileName.c_str());
          std::remove((datFileName+".plt").c_str());
          fo_by_range << "e\n";
          fo_by_range         << "pause -1 \"Press Return\" \n";

        }  
        fo.close();
        fo_count.close();
        fo_by_atom.close();
        fo_by_range.close();
      }  
}
/*
      if (accountBackScattered) break;
      accountBackScattered = true;
    }while (1);
    if (accountSputtered) break;
    accountSputtered = true;
  }while (1);
*/
//  system("cd ""..");
  chdir("..");
}  


#define MDEPP_BYTIME_DIR "_by_time"

void
StatPostProcess::buildByTime(FProcessMolecule fpm) const
{
  system("mkdir "MDEPP_BYTIME_DIR);

//  system("cd "MDEPP_ANGULAR_DIR);
  chdir(MDEPP_BYTIME_DIR);

  char ofilename[1024];

  sprintf(ofilename,"mass_by_formation_time.dat");

  std::string byFormationTimeDat(ofilename);
  std::ofstream foByFormation(byFormationTimeDat.c_str());
  std::ofstream foByFormationPlt((byFormationTimeDat+".plt").c_str());
  foByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << byFormationTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_formation_time.dat");

  std::string EnergyByFormationTimeDat(ofilename);
  std::ofstream foEnergyByFormation(EnergyByFormationTimeDat.c_str());
  std::ofstream foEnergyByFormationPlt((EnergyByFormationTimeDat+".plt").c_str());
  foEnergyByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByFormationTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_atom_by_formation_time.dat");

  std::string EnergyByAtomByFormationTimeDat(ofilename);
  std::ofstream foEnergyByAtomByFormation(EnergyByAtomByFormationTimeDat.c_str());
  std::ofstream foEnergyByAtomByFormationPlt((EnergyByAtomByFormationTimeDat+".plt").c_str());
  foEnergyByAtomByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByAtomByFormationTimeDat << "\' with points\n";

  const Float minTime =    -0.05 /* *mdtk::ps*/;
  const Float maxTime =     6.05 /* *mdtk::ps*/;
  const int n = ceil((maxTime-minTime)/(0.1 /* *mdtk::ps*/));

  gsl_histogram * h_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_by_formation, minTime, maxTime);
  gsl_histogram * h_count_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_count_by_formation, minTime, maxTime);
  gsl_histogram * h_energy_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_energy_by_formation, minTime, maxTime);

  sprintf(ofilename,"mass_by_escape_time.dat");

  std::string byEscapeTimeDat(ofilename);
  std::ofstream foByEscape(byEscapeTimeDat.c_str());
  std::ofstream foByEscapePlt((byEscapeTimeDat+".plt").c_str());
  foByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << byEscapeTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_escape_time.dat");

  std::string EnergyByEscapeTimeDat(ofilename);
  std::ofstream foEnergyByEscape(EnergyByEscapeTimeDat.c_str());
  std::ofstream foEnergyByEscapePlt((EnergyByEscapeTimeDat+".plt").c_str());
  foEnergyByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByEscapeTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_atom_by_escape_time.dat");

  std::string EnergyByAtomByEscapeTimeDat(ofilename);
  std::ofstream foEnergyByAtomByEscape(EnergyByAtomByEscapeTimeDat.c_str());
  std::ofstream foEnergyByAtomByEscapePlt((EnergyByAtomByEscapeTimeDat+".plt").c_str());
  foEnergyByAtomByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByAtomByEscapeTimeDat << "\' with points\n";


  gsl_histogram * h_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_by_escape, minTime, maxTime);
  gsl_histogram * h_count_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_count_by_escape, minTime, maxTime);
  gsl_histogram * h_energy_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_energy_by_escape, minTime, maxTime);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const Molecule& m = td.molecules[mi];
//      if (m.hasProjectileAtoms()) continue;
      if (!fpm(m)) continue;

      Float atomsCount = m.atoms.size();
      Float moleculeMass = m.getAMUMass();
      Float formationTime = m.formationTime /mdtk::ps;
      Float escapeTime = m.escapeTime /mdtk::ps;
      mdtk::Vector3D velocity = m.getVelocity();
      Float kineticEnergy = 0.5*moleculeMass*mdtk::amu*SQR(velocity.module())/mdtk::eV;
      foByFormation << formationTime << " " << moleculeMass << "\n";
      foByEscape    << escapeTime    << " " << moleculeMass << "\n";
      foEnergyByFormation << formationTime << " " << kineticEnergy << "\n";
      foEnergyByEscape    << escapeTime    << " " << kineticEnergy << "\n";
      foEnergyByAtomByFormation << formationTime << " " << kineticEnergy/atomsCount << "\n";
      foEnergyByAtomByEscape    << escapeTime    << " " << kineticEnergy/atomsCount << "\n";
      gsl_histogram_accumulate (h_by_formation, formationTime, moleculeMass);
      gsl_histogram_accumulate (h_by_escape,       escapeTime, moleculeMass);
      gsl_histogram_accumulate (h_count_by_formation, formationTime, atomsCount);
      gsl_histogram_accumulate (h_count_by_escape,       escapeTime, atomsCount);
      gsl_histogram_accumulate (h_energy_by_formation, formationTime, kineticEnergy);
      gsl_histogram_accumulate (h_energy_by_escape,       escapeTime, kineticEnergy);
    }
  }    

  foByEscape.close();
  foByEscapePlt.close();
  foEnergyByEscape.close();
  foEnergyByEscapePlt.close();
  foEnergyByAtomByEscape.close();
  foEnergyByAtomByEscapePlt.close();

  foByFormation.close();
  foByFormationPlt.close();
  foEnergyByFormation.close();
  foEnergyByFormationPlt.close();
  foEnergyByAtomByFormation.close();
  foEnergyByAtomByFormationPlt.close();

  sprintf(ofilename,"mass_by_escape_time_HIST.dat");
  saveHistogram(h_by_escape,ofilename);
  gsl_histogram_free(h_by_escape);

  sprintf(ofilename,"mass_by_formation_time_HIST.dat");
  saveHistogram(h_by_formation,ofilename);
  gsl_histogram_free(h_by_formation);

  sprintf(ofilename,"atomscount_by_escape_time_HIST.dat");
  saveHistogram(h_count_by_escape,ofilename);
  gsl_histogram_free(h_count_by_escape);

  sprintf(ofilename,"energy_by_escape_time_HIST.dat");
  saveHistogram(h_energy_by_escape,ofilename);
  gsl_histogram_free(h_energy_by_escape);

  sprintf(ofilename,"atomscount_by_formation_time_HIST.dat");
  saveHistogram(h_count_by_formation,ofilename);
  gsl_histogram_free(h_count_by_formation);

  sprintf(ofilename,"energy_by_formation_time_HIST.dat");
  saveHistogram(h_energy_by_formation,ofilename);
  gsl_histogram_free(h_energy_by_formation);
  
  chdir("..");
}

/*
void
StatPostProcess::buildByTime() const
{
  bool accountSputtered = false;
  do{
    bool accountBackScattered = false;
    do{
//      if (!accountSputtered || !accountBackScattered)

      buildByTime(fbm);

      if (accountBackScattered) break;
      accountBackScattered = true;
    }while (1);
    if (accountSputtered) break;
    accountSputtered = true;
  }while (1);
}  
*/
//#define SPOTTED_DISTANCE (-10.0*Ao)
//double SPOTTED_DISTANCE = (-5.0*mdtk::Ao);

}  

/*
  std::ofstream fo(ofilename);
  {
    std::ofstream fo((std::string(ofilename)+".plt").c_str());
    fo << "reset\nset xrange [*:20]\nplot \'" << ofilename << "\' with histeps\n";
    fo.close();
  }  

  const Float minEnergy =     0.0;
  const Float maxEnergy =  1000.0;
  const int n = (maxEnergy-minEnergy)/(energyStep);

  fo << "# min energy = " << minEnergy << " eV" << std::endl
     << "# max energy = " << maxEnergy << " eV" << std::endl
     << "# number of bins = " << n << std::endl;
  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, minEnergy, maxEnergy);

  {
    for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
    {
      const TrajData& td = trajData[trajIndex];
      for(size_t mi = 0; mi < td.molecules.size(); mi++)
      {
        const Molecule& mol = td.molecules[mi];
        if ( mol.hasProjectileAtoms() && !accountBackScattered) continue;
        if (!mol.hasProjectileAtoms() && !accountSputtered    ) continue;
        Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
        gsl_histogram_accumulate (h, energy/mol.atoms.size(), mol.atoms.size());
      }
    }    

  }  

  for(size_t i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    fo << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << std::endl;
  }  
  gsl_histogram_free (h);

  fo.close();
*/
