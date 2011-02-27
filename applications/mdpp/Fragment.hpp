#ifndef mdpp_Fragment_hpp
#define mdpp_Fragment_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>

#include "ClassicMolecule.hpp"

namespace mdepp
{

using namespace mdtk;

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

  void addAtom(mdtk::Atom& a) {atoms.push_back(a);}
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

}

#endif
