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


namespace mdepp
{


#define ATOMTAG_FIXED 1<<0
#define ATOMTAG_SUBSTRATE 1<<1
#define ATOMTAG_CLUSTER   1<<2
#define ATOMTAG_NOTAG 0



#define REF_POT_OF(ML_FPOT) ML_FPOT.potentials[0]

inline
void
updateNeighborLists(mdtk::SimLoop* ml)
{
    ml->updateGlobalIndexes();
    REF_POT_OF(ml->fpot)->NL_init(ml->atoms_);
    REF_POT_OF(ml->fpot)->NL_UpdateIfNeeded(ml->atoms_); /// !!!!!!!!!!!!!
}

inline
void
setTags(mdtk::SimLoop* ml)
{
  for(size_t i = 0; i < ml->atoms_.size(); i++)
  {
    mdtk::Atom& atom = *(ml->atoms_[i]);
    atom.tag = 0;
    if (atom.M > 1000.0*mdtk::amu) atom.tag |= ATOMTAG_FIXED;
    if (atom.ID == mdtk::Cu_EL) atom.tag |= ATOMTAG_CLUSTER;
    if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) atom.tag |= ATOMTAG_SUBSTRATE;
//    if (atom.ID == mdtk::Ar_EL) atom.tag |= ATOMTAG_PROJECTILE;
  }
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


  using namespace mdtk;

//extern double SPOTTED_DISTANCE;

extern mdtk::AtomsContainer dummy_ac;

inline
bool  isProjectileAtom(const mdtk::Atom& atom)// const
{
  if (atom.ID == mdtk::Ar_EL) return true;
//  if (atom.tag & ATOMTAG_CLUSTER) return true;
  return false;
}  

inline
bool  isClusterAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::Cu_EL) return true;
  if (atom.tag & ATOMTAG_CLUSTER) return true;
  return false;
}  

inline
bool  isSubstrateAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) return true;
  if (atom.tag & ATOMTAG_SUBSTRATE) return true;
  return false;
}  

class Molecule
{
private:
  enum {ECOUNT = 4};
  enum {C = 0};
  enum {H = 1};
  enum {Ar = 2};
  enum {Cu = 3};
  Float Rc_[ECOUNT][ECOUNT];
  size_t e2i(const mdtk::Atom &atom) const
  {
    using namespace mdtk;
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      case Ar_EL : return Ar; break;
      case Cu_EL : return Cu; break;
      default : throw Exception("Molecule::e2i() : unknown element");
    };  
  }  
public:
  Float Rc(const mdtk::Atom &atom1,const mdtk::Atom &atom2) const
  {
    return Rc_[e2i(atom1)][e2i(atom2)];
  }  
  std::set<mdtk::ElementID> handledElements;
  bool
  isHandled(mdtk::Atom& atom) const
  {
    if (handledElements.find(atom.ID) != handledElements.end())
      return true;
    else
      return false;
  }  
public:
  Float formationTime;
  Float escapeTime;
  std::vector<mdtk::Atom> atoms;
  std::vector<mdtk::Atom> atoms_init;
//  AtomsContainer atoms;
//  AtomsContainer atoms_init;
  int  trajectory_dummy;
/*
Molecule( const Molecule &C )
 :atoms(),atoms_init()
{
  initParams();    

  formationTime = C.formationTime;
  escapeTime = C.escapeTime;
  trajectory = C.trajectory;
  atoms = C.atoms;
  atoms_init = C.atoms_init;
  TRACE(atoms.size());
  for(size_t i = 0; i < atoms.size(); i++)
  {
    TRACE(i);
    atoms[i].container = &dummy_ac;
    TRACE("1");
    atoms_init[i].container = &dummy_ac;
    TRACE("2");
  }
}

inline
Molecule&
operator =(const Molecule &C) 
{
  if (this == &C) return *this;

  formationTime = C.formationTime;
  escapeTime = C.escapeTime;
  trajectory = C.trajectory;
  atoms = C.atoms;
  atoms_init = C.atoms_init;
//  TRACE(atoms.size());
  for(size_t i = 0; i < atoms.size(); i++)
  {
    atoms[i].container = &dummy_ac;
    atoms_init[i].container = &dummy_ac;
  }

  return *this;
}
*/
  void initParams()
  {
    using namespace mdtk;
    handledElements.insert(H_EL);
    handledElements.insert(C_EL);
    handledElements.insert(Ar_EL);
    handledElements.insert(Cu_EL);

    for(size_t e1 = 0; e1 < ECOUNT; e1++)
      for(size_t e2 = 0; e2 < ECOUNT; e2++)
        Rc_[e1][e2] = 0.001*Ao;

    Rc_[C][C] = 2.0*Ao;
    Rc_[H][H] = 1.7*Ao;
    Rc_[C][H] = 1.8*Ao;
      Rc_[H][C] = Rc_[C][H];

    Rc_[Cu][Cu] = 3.0*Ao;


    Rc_[Cu][H] = 2.75*Ao;
      Rc_[H][Cu] = Rc_[Cu][H];

    Rc_[Cu][C] = 2.75*Ao;
      Rc_[C][Cu] = Rc_[Cu][C];
      
//    Rc_[Ar][Ar] = 0.001*Ao;
//      Rc_[Ar][ H] = Rc_[ H][Ar] = Rc_[ C][Ar] = Rc_[Ar][ C] = Rc_[Ar][Ar];
  }

  void addAtom(mdtk::Atom& a) {atoms.push_back(a);atoms[atoms.size()-1].container = &dummy_ac;}
  void buildFromAtom(mdtk::Atom&, mdtk::SimLoop& ml,double SPOTTED_DISTANCE);
  bool hasAtom(mdtk::Atom&) const;
  Molecule():handledElements(),formationTime(-1),escapeTime(-1),
    atoms(/*0*/),atoms_init(),trajectory_dummy(-1)
  {
    initParams();    
  }
  ~Molecule(){;}
/*  
  bool  hasBackScatteredAtoms() const
  {
    return hasProjectileAtoms();
  }  
*/  

  bool  hasProjectileAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isProjectileAtom(atoms[ai])) return true;
    return false;
  }  
  bool  hasClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isClusterAtom(atoms[ai])) return true;
    return false;
  }  
  bool  hasSubstrateAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isSubstrateAtom(atoms[ai])) return true;
    return false;
  }

  bool  hasOnlyProjectileAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isProjectileAtom(atoms[ai])) return false;
    return true;
  }  
  bool  hasOnlyClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isClusterAtom(atoms[ai])) return false;
    return true;
  }  
  bool  hasOnlySubstrateAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isSubstrateAtom(atoms[ai])) return false;
    return true;
  }
  bool  hasOnlySubstrateOrClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!(isSubstrateAtom(atoms[ai]) || isClusterAtom(atoms[ai]))) return false;
    return true;
  }
  Float getAMUMass() const
  {
    Float moleculeMass = 0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      moleculeMass += atom.M;
    } 
    return mdtk::academic_round(moleculeMass/mdtk::amu);
  }
  mdtk::Vector3D getVelocity() const
  {
    using mdtk::Exception;    
    
    REQUIRE(atoms.size() > 0);
    mdtk::Vector3D sumOfP = 0.0;
    Float sumOfM = 0.0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      sumOfM += atom.M;
      sumOfP += atom.V*atom.M;
    };
    return sumOfP/sumOfM;    
  }  
  void printGlobalIndexes(std::ostream& fo) const
  {
    fo << atoms.size() << std::endl;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      fo << atom.globalIndex << std::endl;
    } 
  } 
  void saveToStream(std::ostream& os) const
  {
    os << formationTime << "\n";
    os << escapeTime << "\n";
    os << atoms.size() << "\n";
    for(size_t i = 0; i < atoms.size(); i++)
      os << atoms[i] << "\n";
    os << atoms_init.size() << "\n";
    for(size_t i = 0; i < atoms_init.size(); i++)
      os << atoms_init[i] << "\n";
    os << trajectory_dummy << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
//    ERRTRACE("Loading Molecule....");
    is >> formationTime;
    is >> escapeTime;
    size_t sz;
    is >> sz;
    atoms.resize(sz);
    for(size_t i = 0; i < atoms.size(); i++)
    {
//      ERRTRACE(i);
      is >> atoms[i];
    }
    is >> sz;
    atoms_init.resize(sz);
    for(size_t i = 0; i < atoms_init.size(); i++)
      is >> atoms_init[i];
    is >> trajectory_dummy;
  }  

  friend int operator<(const Molecule& left, const Molecule& right);
  friend int operator<(Molecule& left, Molecule& right);
};  


inline
int operator<(const Molecule& left, const Molecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
int operator<(Molecule& left, Molecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}


class AtomTrajectory
{
public:
  mdtk::Atom beginCheckPoint;
  mdtk::Atom endCheckPoint;
  std::vector<mdtk::Atom> checkPoints;
  AtomTrajectory(mdtk::Atom& finalCheckPoint)
  :beginCheckPoint(),endCheckPoint(finalCheckPoint),checkPoints()
  {
//    checkPoints.push_back(finalCheckPoint);
  }
  AtomTrajectory()
  :beginCheckPoint(),endCheckPoint(),checkPoints()
  {
  }
  void saveToStream(std::ostream& os) const
  {
    os << beginCheckPoint << "\n";
    os << endCheckPoint << "\n";
    os << checkPoints.size() << "\n";
    for(size_t i = 0; i < checkPoints.size(); i++)
      os << checkPoints[i] << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
//    ERRTRACE("Loading AtomTrajectory....");

    is >> beginCheckPoint;
    is >> endCheckPoint;
    size_t sz;
    is >> sz;
    checkPoints.resize(sz);
    for(size_t i = 0; i < checkPoints.size(); i++)
      is >> checkPoints[i];
  }  
};

class ClusterDynamics;

class Fragment
{
private:
public:
  std::set<mdtk::ElementID> handledElements;
  bool
  isHandled(const mdtk::Atom& atom) const
  {
    if (handledElements.find(atom.ID) != handledElements.end() && isClusterAtom(atom))
      return true;
    else
      return false;
  }  
public:
  std::vector<mdtk::Atom> atoms;
//  std::vector<mdtk::Atom> atoms_init;
  int  trajectory_dummy;

  void initParams()
  {
    using namespace mdtk;
    handledElements.insert(Cu_EL);
  }

  void addAtom(mdtk::Atom& a) {atoms.push_back(a);atoms[atoms.size()-1].container = &dummy_ac;}
  void buildFromAtom(const mdtk::Atom&, const ClusterDynamics& cd);
  bool hasAtom(const mdtk::Atom&) const;
  Fragment():handledElements(),
    atoms(/*0*/),/*atoms_init(),*/trajectory_dummy(-1)
  {
    initParams();    
  }
  ~Fragment(){;}
  Float getAMUMass() const
  {
    Float moleculeMass = 0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      moleculeMass += atom.M;
    } 
    return mdtk::academic_round(moleculeMass/mdtk::amu);
  }
  mdtk::Vector3D getVelocity() const
  {
    using mdtk::Exception;    
    
    REQUIRE(atoms.size() > 0);
    mdtk::Vector3D sumOfP = 0.0;
    Float sumOfM = 0.0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      sumOfM += atom.M;
      sumOfP += atom.V*atom.M;
    };
    return sumOfP/sumOfM;    
  }  
  void printGlobalIndexes(std::ostream& fo) const
  {
    fo << atoms.size() << std::endl;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      fo << atom.globalIndex << std::endl;
    } 
  } 
  void saveToStream(std::ostream& os) const
  {
    os << atoms.size() << "\n";
    for(size_t i = 0; i < atoms.size(); i++)
      os << atoms[i] << "\n";
/*
    os << atoms_init.size() << "\n";
    for(size_t i = 0; i < atoms_init.size(); i++)
      os << atoms_init[i] << "\n";
*/
    os << trajectory_dummy << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
    size_t sz;
    is >> sz;
    atoms.resize(sz);
    for(size_t i = 0; i < atoms.size(); i++)
    {
//      ERRTRACE(i);
      is >> atoms[i];
    }
/*
    is >> sz;
    atoms_init.resize(sz);
    for(size_t i = 0; i < atoms_init.size(); i++)
      is >> atoms_init[i];
*/
    is >> trajectory_dummy;
  }  

  friend int operator<(const Fragment& left, const Fragment& right);
  friend int operator<(Fragment& left, Fragment& right);
};  


inline
int operator<(const Fragment& left, const Fragment& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
int operator<(Fragment& left, Fragment& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
bool
Fragment::hasAtom(const mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
      found = true;
  return found;
}  


class ClusterDynamics
{
public:
  mdtk::Vector3D PBC;
  std::vector<AtomTrajectory> atomTrajectories;
  std::vector<Fragment> fragments;
  ClusterDynamics()
  :PBC(),atomTrajectories(),fragments()
  {
  }
  void saveToStream(std::ostream& os) const
  {
    os << PBC << "\n";
    os << atomTrajectories.size() << "\n";
    for(size_t i = 0; i < atomTrajectories.size(); i++)
      atomTrajectories[i].saveToStream(os);
  }  
  void loadFromStream(std::istream& is)
  {
//    ERRTRACE("Loading ClusterDynamics....");
    is >> PBC;
    size_t trajCount;
    is >> trajCount;
    atomTrajectories.resize(trajCount);
    for(size_t i = 0; i < atomTrajectories.size(); i++)
      atomTrajectories[i].loadFromStream(is);
  }  
  bool isAtomInNMer(const mdtk::Atom& atom, const unsigned int N) const
  {
//    TRACE(atom.globalIndex);
    for(size_t mi = 0; mi < fragments.size(); mi++)
    {
//      TRACE(fragments[mi].atoms_init.size());
      if (fragments[mi].hasAtom(atom) && (fragments[mi].atoms.size() == N))
      {
        return true;
      }  
    }  

    return false;
  }
};


class ProjectileDynamics
{
public:
  mdtk::Vector3D PBC;
  std::vector<AtomTrajectory> atomTrajectories;
  ProjectileDynamics()
  :PBC(),atomTrajectories()
  {
  }
  void saveToStream(std::ostream& os) const
  {
    os << PBC << "\n";
    os << atomTrajectories.size() << "\n";
    for(size_t i = 0; i < atomTrajectories.size(); i++)
      atomTrajectories[i].saveToStream(os);
  }  
  void loadFromStream(std::istream& is)
  {
//    ERRTRACE("Loading ClusterDynamics....");
    is >> PBC;
    size_t trajCount;
    is >> trajCount;
    atomTrajectories.resize(trajCount);
    for(size_t i = 0; i < atomTrajectories.size(); i++)
      atomTrajectories[i].loadFromStream(is);
  }  
};


inline
void
Fragment::buildFromAtom(const mdtk::Atom& a, const ClusterDynamics& cd)
{
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);atoms[atoms.size()-1].container = &dummy_ac;

  for(size_t atomIndex = 0; atomIndex < cd.atomTrajectories.size(); atomIndex++)
  {
    const mdtk::Atom& nb_a = cd.atomTrajectories[atomIndex].endCheckPoint;
    if (!isHandled(nb_a)) continue;
    Float distance = (a.coords-nb_a.coords).module();//REF_POT_OF(ml.fpot)->r_vec_module_no_touch(a,nb_a);//sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
//    if (/*Rc(a,nb_a)*/1.50*(2.46*mdtk::Ao*cos(30.0*mdtk::Deg)*2.0) >= distance)
//    if (/*Rc(a,nb_a)*/1.50*(2.46*mdtk::Ao) >= distance)
//    if (/*Rc(a,nb_a)*/1.50*(3.61*mdtk::Ao/2.0) >= distance)
    if (/*Rc(a,nb_a)*/1.00*(3.61*mdtk::Ao) >= distance)
    {
//      if (nb_a.coords.z < SPOTTED_DISTANCE)
//      {
        buildFromAtom(nb_a,cd);
//      }
//      else
//      {
//        atoms.clear();
//      }  
//      if (atoms.size() == 0) break;
    }  
  }  
}  


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
    int aboveSpottedHeight;
    std::vector<Molecule> molecules;
    ClusterDynamics clusterDynamics;
    ProjectileDynamics projectileDynamics;
    Translations trans;
    TrajData() : 
      trajDir(),
      aboveSpottedHeight(0),
      molecules(),
      clusterDynamics(), projectileDynamics(), trans() {;}
    void saveToStream(std::ostream& os) const
    {
      os << trajDir << "\n";
      os << aboveSpottedHeight << "\n";
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
      is >> aboveSpottedHeight;
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

/*  
  int   getAboveSpottedHeight(mdtk::SimLoop&) const; // depends on states !!!!!
  int   getAboveSpottedHeightTotal() const; 
*/
  int   getYield(size_t trajIndex, FProcessMolecule fpm) const;
  int   getYieldSum( FProcessMolecule fpm) const;

  Float   getEnergyOfSputtered(size_t trajIndex, FProcessMolecule fpm) const;
  Float   getTotalEnergyOfSputtered( FProcessMolecule fpm) const;
  Float   getAverageEnergyOfSputtered( FProcessMolecule fpm) const;

  Float getAverageYield( FProcessMolecule fpm) const;
  Float getAverageYieldProgress( FProcessMolecule fpm) const;

  void  buildSpottedMolecules(mdtk::SimLoop&,size_t trajIndex);

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
  
  std::string
  buildAtomByEnergy(const Float energyStep, FProcessMolecule fpm) const;
  
  std::string
  buildEnergyByPolar(const int n, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  std::string
  buildAtomsCountByPolar(const int n, FProcessMolecule fpm) const;
  std::string
  buildEnergyByAzimuth(const int n, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  std::string
  buildAtomsCountByAzimuth(const int n, FProcessMolecule fpm) const;

  std::string
  buildEnergyByPolarByAtomsInRange(const int n, FProcessMolecule fpm) const;
  std::string
  buildEnergyByAzimuthByAtomsInRange(const int n, FProcessMolecule fpm) const;

  void  histEnergyByPolar(gsl_histogram* h, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByPolar(gsl_histogram* h, FProcessMolecule fpm) const;
  void  histEnergyByAzimuth(gsl_histogram* h, bool byAtom = false, FProcessMolecule fpm = &ProcessAll) const;
  void  histAtomsCountByAzimuth(gsl_histogram* h, FProcessMolecule fpm) const;

  void  histEnergyByPolarByAtomsInRange(gsl_histogram* h, FProcessMolecule fpm) const;
  void  histEnergyByAzimuthByAtomsInRange(gsl_histogram* h, FProcessMolecule fpm) const;
/*
  void  buildEnergyByPolarByAtom(const int n,
    bool accountSputtered = true, bool accountBackScattered = true) const;
  void  buildEnergyByAzimuthByAtom(const int n,
    bool accountSputtered = true, bool accountBackScattered = true) const;
*/  

//  void  buildAngular() const;
  void  buildAngular2(FProcessMolecule fpm) const;

//  void  buildByTime() const;
  void  buildByTime( FProcessMolecule fpm) const;

//  void  setSpottedDistanceFromInit(std::string mde_init_filename) const;
  void  setSpottedDistanceFromInit();// const;

  void saveToStream(std::ostream& os) const
  {
//    PostProcess::saveToStream(os);
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
