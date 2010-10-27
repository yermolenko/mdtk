/* 
   Molecular dynamics postprocessor, main classes, header

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

#ifndef mdtk_StatPostProcess_hpp
#define mdtk_StatPostProcess_hpp

#include <iostream>
#include <stdlib.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>

using namespace std;

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>

#include <mdtk/Atom.hpp>
#include <mdtk/SimLoop.hpp>

#include "tools.hpp"

#include "Molecule.hpp"
#include "Fragment.hpp"

namespace mdepp
{

inline
void
updateNeighborLists(mdtk::SimLoop* ml)
{
  ml->updateGlobalIndexes();
  REF_POT_OF(ml->fpot)->NL_init(ml->atoms_);
  REF_POT_OF(ml->fpot)->NL_UpdateIfNeeded(ml->atoms_); /// !!!!!!!!!!!!!
}

using namespace mdtk;

class Translation
{
public:
  mdtk::Atom start;
  mdtk::Atom end;
  Translation(mdtk::Atom& as, mdtk::Atom& ae)
  :start(as),end(ae)
  {
  }
  Translation()
  :start(),end()
  {
  }
  void saveToStream(std::ostream& os) const
  {
    os << start << "\n";
    os << end << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
    is >> start;
    is >> end;
  }  
};

class Translations
{
public:
  std::vector<Translation> translations;
  Translations()
  :translations()
  {
  }
  void saveToStream(std::ostream& os) const
  {
    os << translations.size() << "\n";
    for(size_t i = 0; i < translations.size(); i++)
      translations[i].saveToStream(os);
  }  
  void loadFromStream(std::istream& is)
  {
    size_t trajCount;
    is >> trajCount;
    translations.resize(trajCount);
    for(size_t i = 0; i < translations.size(); i++)
      translations[i].loadFromStream(is);
  }  
};

class StatPostProcess
{
public:
  typedef bool (*FProcessMolecule)(const Molecule&);
  FProcessMolecule testProcessMolecule;
  static
  bool ProcessAll(const Molecule&)
  {
    return true;
  }
  static
  bool ProcessProjectile(const Molecule& mol)
  {
    return mol.hasOnlyProjectileAtoms();
  }
  static
  bool ProcessCluster(const Molecule& mol)
  {
    return mol.hasOnlyClusterAtoms();
  }
  static
  bool ProcessSubstrate(const Molecule& mol)
  {
    return mol.hasOnlySubstrateAtoms();
  }
  static
  bool ProcessClusterAndSubstrate(const Molecule& mol)
  {
//    return mol.hasSubstrateAtoms() || mol.hasClusterAtoms();
    return mol.hasOnlySubstrateOrClusterAtoms();
  }
  struct TrajData
  {
    std::string trajDir;
    std::vector<Molecule> molecules;
    ClusterDynamics clusterDynamics;
    ProjectileDynamics projectileDynamics;
    Translations trans;
    TrajData() : 
      trajDir(),
      molecules(),
      clusterDynamics(), projectileDynamics(), trans() {;}
    void saveToStream(std::ostream& os) const
    {
      os << trajDir << "\n";
      os << molecules.size() << "\n";
      for(size_t i = 0; i < molecules.size(); i++)
        molecules[i].saveToStream(os);
        
//      os << SPOTTED_DISTANCE << "\n";
      clusterDynamics.saveToStream(os);
      projectileDynamics.saveToStream(os);
      trans.saveToStream(os);
    }  
    void loadFromStream(std::istream& is)
    {
      is >> trajDir;
      size_t sz;
      is >> sz;
      molecules.resize(sz);
      for(size_t i = 0; i < molecules.size(); i++)
        molecules[i].loadFromStream(is);
        
//      is >> SPOTTED_DISTANCE;
      clusterDynamics.loadFromStream(is);
      projectileDynamics.loadFromStream(is);
      trans.loadFromStream(is);
    }  
  };
  std::vector<TrajData> trajData;
  double SPOTTED_DISTANCE;
  StatPostProcess(std::vector<std::string>& savedStateNames)
   :testProcessMolecule(&ProcessAll),
    trajData(),
    SPOTTED_DISTANCE(-5.0*mdtk::Ao)
  {
    using mdtk::Exception;
    
    std::vector<_SavedStateSortStruct> sorted;
    for(size_t i = 0; i < savedStateNames.size(); i++)
      sorted.push_back(savedStateNames[i]);
  
    sort(sorted.begin(), sorted.end());

    for(size_t i = 0; i < sorted.size(); i++)
    {
      trajData.push_back(TrajData());
      trajData.back().trajDir = sorted[i].fullTrajDirName;
    }  
  
    REQUIRE(savedStateNames.size() == trajData.size());

    setSpottedDistanceFromInit();

    for(size_t i = 0; i < trajData.size(); i++)
      TRACE(trajData[i].trajDir);
  }  
  StatPostProcess()
   :trajData(),
    SPOTTED_DISTANCE(-5.0*mdtk::Ao)
  {
  }  
  
  int   getAboveSpottedHeight(mdtk::SimLoop&) const; 

  int   getYield(size_t trajIndex, FProcessMolecule fpm) const;
  int   getYieldSum( FProcessMolecule fpm) const;

  Float   getEnergyOfSputtered(size_t trajIndex, FProcessMolecule fpm) const;
  Float   getTotalEnergyOfSputtered( FProcessMolecule fpm) const;
  Float   getAverageEnergyOfSputtered( FProcessMolecule fpm) const;

  Float getAverageYield( FProcessMolecule fpm) const;
  Float getAverageYieldProgress( FProcessMolecule fpm) const;

  enum StateType{STATE_FINAL,STATE_INIT,STATE_INTER};

  void  buildSputteredMolecules(mdtk::SimLoop&,size_t trajIndex, StateType s);
  void  buildClusterDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);
  void  buildProjectileDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);
  void  buildDummyDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);

//  void  buildTransitions(mdtk::SimLoop&,size_t trajIndex, StateType s);

  void  printMolecules(size_t trajIndex) const;
  void  printMoleculesTotal() const;

  void  printClusterDynamics(size_t trajIndex) const;
  void  printClusterDynamicsTotal() const;

  void  printClusterDynamicsRSQ(bool xyOnly = false) const;

  void  printProjectileDynamics() const;
  void  printProjectileStopping() const;

  void  printAtomTransitions() const;

  void  buildClusterFragmentsFromDyn();
  void  printClusterFragments() const;

  void  printStoppingHist() const;
  void  printClusterAtomsByAzimuth(const int n = 36) const;
  void  printClusterAtomsByAzimuthPos(const int n = 36, bool halfshift = false, Float accountedDist = 0.75*2.46*mdtk::Ao) const;
//0.75*(2.46*mdtk::Ao*cos(30.0*mdtk::Deg)*2.0)

  void  printEmphasized(size_t trajIndex) const;
  void  printEmphasizedTotal() const;
  
  void  spottedTotalMDE() const;
  void  spottedTotalByMass() const;

  void  spottedByDepth() const;
  
  void  buildMassSpectrum(FProcessMolecule fpm = &ProcessAll) const;
  
  std::string  buildAtomByEnergy(const Float energyStep, FProcessMolecule fpm) const;
  
  std::string  buildEnergyByPolar(const int n, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  std::string  buildAtomsCountByPolar(const int n, FProcessMolecule fpm) const;
  std::string  buildEnergyByAzimuth(const int n, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  std::string  buildAtomsCountByAzimuth(const int n, FProcessMolecule fpm) const;

  std::string  buildEnergyByPolarByAtomsInRange(const int n, FProcessMolecule fpm) const;
  std::string  buildEnergyByAzimuthByAtomsInRange(const int n, FProcessMolecule fpm) const;

  void  histEnergyByPolar(gsl_histogram* h, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByPolar(gsl_histogram* h, FProcessMolecule fpm) const;
  void  histEnergyByAzimuth(gsl_histogram* h, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByAzimuth(gsl_histogram* h, FProcessMolecule fpm) const;

  void  histEnergyByPolarByAtomsInRange(gsl_histogram* h, FProcessMolecule fpm) const;
  void  histEnergyByAzimuthByAtomsInRange(gsl_histogram* h, FProcessMolecule fpm) const;

  void  buildAngular2(FProcessMolecule fpm) const;

  void  buildByTime( FProcessMolecule fpm) const;

  void  setSpottedDistanceFromInit();// const;

  void saveToStream(std::ostream& os) const
  {
    os << trajData.size() << "\n";
    for(size_t i = 0; i < trajData.size(); i++)
      trajData[i].saveToStream(os);

    os << SPOTTED_DISTANCE << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
    size_t sz;
    is >> sz;
    trajData.resize(sz);
    for(size_t i = 0; i < trajData.size(); i++)
      trajData[i].loadFromStream(is);

    is >> SPOTTED_DISTANCE;
  }  
  void loadFromStreamADD(std::istream& is)
  {
    size_t sz;
    is >> sz;
    size_t oldSize = trajData.size();
    TRACE(oldSize);
    trajData.resize(oldSize+sz);
    TRACE(trajData.size());
    for(size_t i = oldSize+0; i < trajData.size(); i++)
      trajData[i].loadFromStream(is);

    is >> SPOTTED_DISTANCE;
  }  
  void execute();
};

struct MassCount
{
  Float amuMass;
  int count;
  MassCount(Float amuMassVal = 0) : amuMass(amuMassVal), count(0) {;}
  friend int operator<(const MassCount& left, const MassCount& right);
  friend int operator<(MassCount& left, MassCount& right);
};  

inline
int operator<(const MassCount& left, const MassCount& right)
{
  return left.amuMass < right.amuMass;
}

inline
int operator<(MassCount& left, MassCount& right)
{
  return left.amuMass < right.amuMass;
}


struct MassSpotted
{
  std::vector<Molecule> species;
  Float getAMUMass()const{return species[0].getAMUMass();};
  MassSpotted(Molecule& m):species(){species.push_back(m);}
  MassSpotted(const Molecule& m):species(){species.push_back(m);}
  friend int operator<(const MassSpotted& left, const MassSpotted& right);
  friend int operator<(MassSpotted& left, MassSpotted& right);
};  

inline
int operator<(const MassSpotted& left, const MassSpotted& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
int operator<(MassSpotted& left, MassSpotted& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

void saveSpottedByMass(std::vector<MassSpotted>& massSpectrum);


}  


#endif
