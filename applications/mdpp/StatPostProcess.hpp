/* 
   Molecular dynamics postprocessor, main classes, header

   Copyright (C) 2007, 2008, 2009, 2010, 2011 Oleksandr Yermolenko
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
#include <map>
#include <algorithm>
#include <fstream>

#include <mdtk/Atom.hpp>
#include <mdtk/SimLoop.hpp>

#include "tools.hpp"

#include "NeighbourList.hpp"

#include "ClassicMolecule.hpp"
#include "Fragment.hpp"

#include "Fullerene.hpp"

namespace mdepp
{

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
  typedef bool (*FProcessClassicMolecule)(const ClassicMolecule&);
  FProcessClassicMolecule testProcessClassicMolecule;
  static
  bool ProcessAll(const ClassicMolecule&)
  {
    return true;
  }
  static
  bool ProcessProjectile(const ClassicMolecule& mol)
  {
    return mol.hasProjectileAtoms();
  }
  static
  bool ProcessCluster(const ClassicMolecule& mol)
  {
    return mol.hasClusterAtoms();
  }
  static
  bool ProcessFullerene(const ClassicMolecule& mol)
  {
    return mol.hasFullereneAtoms();
  }
  static
  bool ProcessSubstrate(const ClassicMolecule& mol)
  {
    return mol.hasSubstrateAtoms();
  }
  static
  bool ProcessClusterAndSubstrate(const ClassicMolecule& mol)
  {
//    return mol.hasSubstrateAtoms() || mol.hasClusterAtoms();
    return mol.hasOnlySubstrateOrClusterAtoms();
  }
  struct TrajData
  {
    std::string trajDir;
    std::vector<ClassicMolecule> molecules;
    ClusterDynamics clusterDynamics;
    ProjectileDynamics projectileDynamics;
    Translations trans;
    std::map <Float,Fullerene> trajFullerene;
    Vector3D PBC;
    TrajData() : 
      trajDir(),
      molecules(),
      clusterDynamics(), projectileDynamics(), trans(),
      trajFullerene(), PBC() {;}
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

      os << trajFullerene.size() << "\n";
      std::map< Float, Fullerene >::const_iterator i;
      for( i = trajFullerene.begin(); i != trajFullerene.end() ; ++i )
	os << i->first << "\t " << i->second << "\n";

      os << PBC << "\n";
    }  
    void loadFromStream(std::istream& is)
    {
      is >> trajDir;
      size_t sz, i;
      is >> sz;
      molecules.resize(sz);
      for(i = 0; i < molecules.size(); i++)
        molecules[i].loadFromStream(is);
        
//      is >> SPOTTED_DISTANCE;
      clusterDynamics.loadFromStream(is);
      projectileDynamics.loadFromStream(is);
      trans.loadFromStream(is);

      is >> sz;
      for(i = 0; i < sz; ++i)
      {
	Float t;
	Fullerene f;
	is >> t >> f;
	trajFullerene[t] = f;
      }

      is >> PBC;
    }  
  };
  std::vector<TrajData> trajData;
  double SPOTTED_DISTANCE;
  StatPostProcess(std::vector<std::string>& savedStateNames)
   :testProcessClassicMolecule(&ProcessAll),
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

    TRACE(trajData.size());

    for(size_t i = 0; i < trajData.size(); i++)
      TRACE(trajData[i].trajDir);
  }  
  StatPostProcess()
   :trajData(),
    SPOTTED_DISTANCE(-5.0*mdtk::Ao)
  {
  }  
  
  int   getAboveSpottedHeight(mdtk::SimLoop&) const; 

  int   getYield(size_t trajIndex, FProcessClassicMolecule fpm) const;
  int   getYieldSum( FProcessClassicMolecule fpm) const;

  Float   getEnergyOfSputtered(size_t trajIndex, FProcessClassicMolecule fpm) const;
  Float   getTotalEnergyOfSputtered( FProcessClassicMolecule fpm) const;
  Float   getAverageEnergyOfSputtered( FProcessClassicMolecule fpm) const;

  Float getAverageYield( FProcessClassicMolecule fpm) const;
  Float getAverageYieldProgress( FProcessClassicMolecule fpm) const;

  enum StateType{STATE_FINAL,STATE_INIT,STATE_INTER};

  void  buildSputteredClassicMolecules(mdtk::SimLoop&,size_t trajIndex, StateType s, NeighbourList& nl);
  void  buildClusterDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s, NeighbourList& nl);
  void  buildProjectileDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);
  void  buildDummyDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);

//  void  buildTransitions(mdtk::SimLoop&,size_t trajIndex, StateType s);

  void  printClassicMolecules(size_t trajIndex) const;
  void  printClassicMoleculesTotal() const;

  void  printFullereneInfo(size_t trajIndex) const;
  void  printFullereneInfo() const;

  void  plotFullereneLandings(bool endo, const std::string rotDir, Float integralThreshold = 3.0*Ao) const;
  void  plotFullereneImplantDepth(bool endo, const std::string rotDir, Float integralThreshold = 3.0*Ao) const;

  bool  isThereAnythingToPlot(bool endo, const std::string rotDir) const;

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
  
  void  buildMassSpectrum(FProcessClassicMolecule fpm = &ProcessAll) const;
  
  std::string  buildAtomByEnergy(const Float energyStep, FProcessClassicMolecule fpm) const;
  
  std::string  buildEnergyByPolar(const int n, bool byAtom = false, FProcessClassicMolecule fpm = &ProcessAll) const;
  std::string  buildAtomsCountByPolar(const int n, FProcessClassicMolecule fpm) const;
  std::string  buildEnergyByAzimuth(const int n, bool byAtom = false, FProcessClassicMolecule fpm = &ProcessAll) const;
  std::string  buildAtomsCountByAzimuth(const int n, FProcessClassicMolecule fpm) const;

  std::string  buildEnergyByPolarByAtomsInRange(const int n, FProcessClassicMolecule fpm) const;
  std::string  buildEnergyByAzimuthByAtomsInRange(const int n, FProcessClassicMolecule fpm) const;

  void  histEnergyByPolar(gsl_histogram* h, bool byAtom = false, FProcessClassicMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByPolar(gsl_histogram* h, FProcessClassicMolecule fpm) const;
  void  histEnergyByAzimuth(gsl_histogram* h, bool byAtom = false, FProcessClassicMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByAzimuth(gsl_histogram* h, FProcessClassicMolecule fpm) const;

  void  histEnergyByPolarByAtomsInRange(gsl_histogram* h, FProcessClassicMolecule fpm) const;
  void  histEnergyByAzimuthByAtomsInRange(gsl_histogram* h, FProcessClassicMolecule fpm) const;

  void  buildAngular2(FProcessClassicMolecule fpm) const;

  void  buildByTime( FProcessClassicMolecule fpm) const;

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
  std::vector<ClassicMolecule> species;
  Float getAMUMass()const{return species[0].getAMUMass();};
  MassSpotted(ClassicMolecule& m):species(){species.push_back(m);}
  MassSpotted(const ClassicMolecule& m):species(){species.push_back(m);}
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
