/* 
   Molecular dynamics postprocessor, main classes, header

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015,
   2016 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

extern std::string initialWorkingDirectory;

class StatPostProcess
{
public:
  typedef bool (*FProcessClassicMolecule)(const ClassicMolecule&);
  FProcessClassicMolecule testProcessClassicMolecule;
  static bool ProcessAll(const ClassicMolecule&);
  static bool ProcessProjectile(const ClassicMolecule& mol);
  static bool ProcessCluster(const ClassicMolecule& mol);
  static bool ProcessFullerene(const ClassicMolecule& mol);
  static bool ProcessSubstrate(const ClassicMolecule& mol);
  static bool ProcessClusterAndSubstrate(const ClassicMolecule& mol);
  static std::string FProcessClassicMoleculeToString(FProcessClassicMolecule fpm);

  struct TrajData
  {
    std::string trajDir;
    std::vector<ClassicMolecule> molecules;
    std::map <Float,AtomGroup> trajProjectile;
    std::map <Float,AtomGroup> trajCluster;

    Vector3D PBC;

    TrajData();
    void saveToStream(std::ostream& os) const;
    void loadFromStream(std::istream& is);
    void execute(mdtk::SimLoop& state);

    enum StateType{STATE_FINAL,STATE_INIT,STATE_INTER};
    void  buildSputteredClassicMolecules(
      mdtk::SimLoop& state,
      StateType s,
      NeighbourList& nl);
    void  buildDummyDynamics(
      mdtk::SimLoop& state,
      StateType s);

    std::map <Float,AtomGroup> getTrajWithPartialSnapshots(
      const std::map<Float,AtomGroup>& trajectorySeed) const;

    double SPOTTED_DISTANCE;
    void  setSpottedDistanceFromInit();

    int   getAboveSpottedHeight(mdtk::SimLoop&) const;

    bool hasIntactClusterSputtering() const;
    static std::map <Float,Float> plot_Ekin_t(const std::map <Float,AtomGroup>&);

    int yield(FProcessClassicMolecule fpm) const;
  };
  static std::string getCacheFilename(std::string);
  static std::string cacheDir;
  static std::map<std::string,TrajData> ramCache;
  std::vector<TrajData*> trajData;

  typedef bool (*TrajFilter)(const TrajData&);
  static TrajFilter trajFilter;
  static bool TrajFilterProcessAll(const TrajData&);
  static bool TrajFilterProcessIntactClusterOnly(const TrajData&);

  struct Id
  {
    std::string str;
    mdtk::ElementID clusterElement;
    size_t clusterSize;
    mdtk::ElementID ionElement;
    Float ionEnergy;
    Id(std::string s);
    Id();
    void saveToStream(std::ostream& os) const;
    void loadFromStream(std::istream& is);
  }id;

  static mdtk::SimLoop stateTemplate;

  StatPostProcess(const std::vector<std::string> trajDirs);
  // StatPostProcess();
  static size_t instanceCounter;
  virtual ~StatPostProcess();

  Float yield(FProcessClassicMolecule fpm);

  int   getYield(size_t trajIndex, FProcessClassicMolecule fpm) const;
  Float getYieldNormalizedByClusterSize(size_t trajIndex, FProcessClassicMolecule fpm) const;
  int   getYieldSum( FProcessClassicMolecule fpm) const;
  // int   getYieldMax( FProcessClassicMolecule fpm) const;
  Float getYieldAverage( FProcessClassicMolecule fpm) const;
  Float getYieldAverageProgress( FProcessClassicMolecule fpm) const;

  // Float getYieldMass(size_t trajIndex, FProcessClassicMolecule fpm) const;
  // Float getYieldMassSum( FProcessClassicMolecule fpm) const;
  // Float getYieldMassMax( FProcessClassicMolecule fpm) const;
  // Float getYieldMassAverage( FProcessClassicMolecule fpm) const;

  int   getYieldFragmentsCount(size_t trajIndex, FProcessClassicMolecule fpm) const;
  Float getYieldFragmentsCountNormalizedByClusterSize(size_t trajIndex, FProcessClassicMolecule fpm) const;
  // int   getYieldFragmentsCountSum( FProcessClassicMolecule fpm) const;
  // int   getYieldFragmentsCountMax( FProcessClassicMolecule fpm) const;
  // Float getYieldFragmentsCountAverage( FProcessClassicMolecule fpm) const;

  // int   getFragmentsCount(size_t trajIndex, FProcessClassicMolecule fpm) const;
  // int   getFragmentsCountSum( FProcessClassicMolecule fpm) const;
  // int   getFragmentsCountMax( FProcessClassicMolecule fpm) const;
  // Float getFragmentsCountAverage( FProcessClassicMolecule fpm) const;

  Float   getEnergyOfSputtered(size_t trajIndex, FProcessClassicMolecule fpm) const;
  Float   getTotalEnergyOfSputtered( FProcessClassicMolecule fpm) const;
  Float   getAverageEnergyOfSputtered( FProcessClassicMolecule fpm) const;

  void  printClassicMolecules(size_t trajIndex) const;
  void  printClassicMoleculesTotal() const;

  void  printFullereneInfo(size_t trajIndex) const;
  void  printFullereneInfo() const;

  void  printCoefficients() const;

  std::map<ClassicMolecule, size_t>  buildMassSpectrum(FProcessClassicMolecule fpm = &ProcessAll) const;

  void  spottedByDepth() const;

  typedef double (*FMoleculeAttribute)(const ClassicMolecule&);
  static double moleculeEnergy(const ClassicMolecule&);
  static double moleculeEnergyInEV(const ClassicMolecule&);
  static double moleculeMass(const ClassicMolecule&);
  static double moleculeMassInAMU(const ClassicMolecule&);
  static double moleculeCount(const ClassicMolecule&);
  static double moleculeAtomsCount(const ClassicMolecule&);
  static double moleculeEnergyByAtom(const ClassicMolecule&);
  static std::string FMoleculeAttributeToString(FMoleculeAttribute fma);

  enum AngleType{ANGLE_POLAR,ANGLE_AZIMUTH};
  std::map<Float, Float> distByAngle(
    AngleType angleType,
    const int n,
    FMoleculeAttribute fma,
    FProcessClassicMolecule fpm) const;

  Float maxMoleculeAttribute(
    FMoleculeAttribute fma, FProcessClassicMolecule moleculeFilter) const;
  static Float suggestedBinWidth(
    FMoleculeAttribute fma, FProcessClassicMolecule moleculeFilter,
    ElementID ionElement, size_t clusterSize, ElementID clusterElement);

  std::map<Float, Float> distBy(
    FMoleculeAttribute histFunc,
    Float binWidth,
    Float histMin,
    Float histMax,
    FMoleculeAttribute fma,
    FProcessClassicMolecule fpm
    ) const;

  static std::map<Float, Float>
  divideHistograms(
    std::map<Float, Float>& h1,
    std::map<Float, Float>& h2);

  void  buildByTime( FProcessClassicMolecule fpm) const;
};

bool
isAmongSputtered(const AtomGroup& atoms, const std::vector<ClassicMolecule>& molecules);

class BadTrajectoryException
{
  std::string _msg;
public:
  BadTrajectoryException() : _msg("Trajectory is bad.") { }
  BadTrajectoryException(const char* msg) : _msg(msg) { }
  BadTrajectoryException(std::string msg) : _msg(msg) { }
  const char* what() const
  {
    return _msg.c_str();
  }
};

}

#endif
