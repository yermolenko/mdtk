/* 
   Molecular dynamics postprocessor, main classes, header

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Oleksandr Yermolenko
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
    std::map <Float,AtomGroup> trajProjectile;
    std::map <Float,AtomGroup> trajCluster;
    Vector3D PBC;
    TrajData() :
      trajDir(),
      molecules(),
      trajProjectile(),
      trajCluster(),
      PBC() {;}
    void saveToStream(std::ostream& os) const
    {
      os << trajDir << "\n";
      os << molecules.size() << "\n";
      for(size_t i = 0; i < molecules.size(); i++)
        molecules[i].saveToStream(os);

      {
        os << trajProjectile.size() << "\n";
        std::map< Float, AtomGroup >::const_iterator i;
        for( i = trajProjectile.begin(); i != trajProjectile.end() ; ++i )
          os << i->first << "\t " << i->second << "\n";
      }

      {
        os << trajCluster.size() << "\n";
        std::map< Float, AtomGroup >::const_iterator i;
        for( i = trajCluster.begin(); i != trajCluster.end() ; ++i )
          os << i->first << "\t " << i->second << "\n";
      }

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

      is >> sz;
      for(i = 0; i < sz; ++i)
      {
	Float t;
	AtomGroup f;
	is >> t >> f;
	trajProjectile[t] = f;
      }

      is >> sz;
      for(i = 0; i < sz; ++i)
      {
	Float t;
	AtomGroup f;
	is >> t >> f;
	trajCluster[t] = f;
      }

      is >> PBC;
    }
  };
  std::vector<TrajData> trajData;
  double SPOTTED_DISTANCE;
  struct Id
  {
    std::string str;
    mdtk::ElementID clusterElement;
    size_t clusterSize;
    mdtk::ElementID ionElement;
    Float ionEnergy;
    Id(std::string s);
    Id():str(),
         clusterElement(),
         clusterSize(),
         ionElement(),
         ionEnergy(){};
    void saveToStream(std::ostream& os) const
      {
        os << str << "\n";
        os << clusterElement << "\n";
        os << clusterSize << "\n";
        os << ionElement << "\n";
        os << ionEnergy << "\n";
      }
    void loadFromStream(std::istream& is)
      {
        is >> str;
        int ID;
        is >> ID; clusterElement = ElementID(ID);
        is >> clusterSize;
        is >> ID; ionElement = ElementID(ID);
        is >> ionEnergy;
      }
  }id;
  StatPostProcess(std::string trajsetDir)
   :testProcessClassicMolecule(&ProcessAll),
    trajData(),
    SPOTTED_DISTANCE(-5.0*mdtk::Ao),
    id(yaatk::extractItemFromEnd(trajsetDir,1))
  {
    using mdtk::Exception;

    std::vector<std::string> savedStateNames;
    mdepp::FProcessTrajectory fpt = mdepp::trajProcess_Custom2;
    mdepp::addTrajDirNames(savedStateNames,trajsetDir.c_str(),fpt);
    std::sort(savedStateNames.begin(),savedStateNames.end());

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
    SPOTTED_DISTANCE(-5.0*mdtk::Ao),
    id()
  {
  }

  void  setSpottedDistanceFromInit();

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
  void  buildDummyDynamics(mdtk::SimLoop&,size_t trajIndex, StateType s);

  void  printClassicMolecules(size_t trajIndex) const;
  void  printClassicMoleculesTotal() const;

  void  printFullereneInfo(size_t trajIndex) const;
  void  printFullereneInfo() const;

  void  printCoefficients() const;

  void  buildMassSpectrum(FProcessClassicMolecule fpm = &ProcessAll) const;

  void  spottedByDepth() const;

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

  void  buildAngular(FProcessClassicMolecule fpm) const;

  void  buildByTime( FProcessClassicMolecule fpm) const;

  void  saveToStream(std::ostream& os) const
  {
    os << trajData.size() << "\n";
    for(size_t i = 0; i < trajData.size(); i++)
      trajData[i].saveToStream(os);

    os << SPOTTED_DISTANCE << "\n";

    id.saveToStream(os);
  }
  void loadFromStream(std::istream& is)
  {
    size_t sz;
    is >> sz;
    trajData.resize(sz);
    for(size_t i = 0; i < trajData.size(); i++)
      trajData[i].loadFromStream(is);

    is >> SPOTTED_DISTANCE;

    id.loadFromStream(is);
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
  void removeBadTrajectories();
};

}

#endif
