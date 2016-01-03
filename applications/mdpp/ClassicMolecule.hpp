/*
   The Simple Molecule class (header file).

   Copyright (C) 2010, 2011, 2012, 2014, 2015 Oleksandr Yermolenko
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

#ifndef mdpp_ClassicMolecule_hpp
#define mdpp_ClassicMolecule_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <map>
#include <algorithm>

#include "NeighbourList.hpp"

namespace mdepp
{

using namespace mdtk;

//#define REF_POT_OF(ML_FPOT) ML_FPOT.potentials[0]

inline
void
setTags(mdtk::SimLoop& ml)
{
  return; // disable retagging
  for(size_t i = 0; i < ml.atoms.size(); i++)
  {
    mdtk::Atom& atom = ml.atoms[i];
    atom.clearTags();
//    if (atom.M > 1000.0*mdtk::amu) atom.tag |= ATOMTAG_FIXED;
    if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL)
    {
      atom.tag(ATOMTAG_SUBSTRATE);
    }
    if (atom.ID == mdtk::Ar_EL || atom.ID == mdtk::Xe_EL)
    {
      atom.tag(ATOMTAG_PROJECTILE);
    }
    if (atom.ID == mdtk::Cu_EL || atom.ID == mdtk::Au_EL)
    {
      atom.tag(ATOMTAG_CLUSTER);
    }
  }
}

inline
bool  isProjectileAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::Ar_EL) return true;
  if (atom.hasTag(ATOMTAG_PROJECTILE)) return true;
  return false;
}  

inline
bool  isClusterAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::Cu_EL) return true;
  if (atom.hasTag(ATOMTAG_CLUSTER)) return true;
  return false;
}  

inline
bool  isSubstrateAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) return true;
  if (atom.hasTag(ATOMTAG_SUBSTRATE)) return true;
  return false;
}  

inline
bool  isFullereneAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) return true;
  if (atom.hasTag(ATOMTAG_FULLERENE)) return true;
  return false;
}  

class ClassicMolecule
{
private:
  enum {ECOUNT = 6};
  enum {C = 0};
  enum {H = 1};
  enum {Ar = 2};
  enum {Xe = 3};
  enum {Cu = 4};
  enum {Au = 5};
  Float Rc_[ECOUNT][ECOUNT];
  size_t e2i(const mdtk::Atom &atom) const
  {
    using namespace mdtk;
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      case Ar_EL : return Ar; break;
      case Xe_EL : return Xe; break;
      case Cu_EL : return Cu; break;
      case Au_EL : return Au; break;
      default : throw Exception("ClassicMolecule::e2i() : unknown element");
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
ClassicMolecule( const ClassicMolecule &C )
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
ClassicMolecule&
operator =(const ClassicMolecule &C) 
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
    handledElements.insert(Xe_EL);
    handledElements.insert(Cu_EL);
    handledElements.insert(Au_EL);

    for(size_t e1 = 0; e1 < ECOUNT; e1++)
      for(size_t e2 = 0; e2 < ECOUNT; e2++)
        Rc_[e1][e2] = 0.001*Ao;

    Rc_[C][C] = 3.40*Ao*1.5;
    Rc_[H][H] = 2.65*Ao*1.5;
    Rc_[C][H] = 0.5*(Rc_[C][C] + Rc_[H][H]);
      Rc_[H][C] = Rc_[C][H];

    Rc_[Cu][Cu] = 3.0*Ao*1.5;


    Rc_[Cu][H] = 0.5*(Rc_[Cu][Cu] + Rc_[H][H]);
      Rc_[H][Cu] = Rc_[Cu][H];

    Rc_[Cu][C] = 0.5*(Rc_[Cu][Cu] + Rc_[C][C]);
      Rc_[C][Cu] = Rc_[Cu][C];


      Rc_[Au][Au] = 3.0*Ao*1.5;// 2.6*Ao*1.5


    Rc_[Au][H] = 0.5*(Rc_[Au][Au] + Rc_[H][H]);
      Rc_[H][Au] = Rc_[Au][H];

    Rc_[Au][C] = 0.5*(Rc_[Au][Au] + Rc_[C][C]);
      Rc_[C][Au] = Rc_[Au][C];

//    Rc_[Ar][Ar] = 0.001*Ao;
//      Rc_[Ar][ H] = Rc_[ H][Ar] = Rc_[ C][Ar] = Rc_[Ar][ C] = Rc_[Ar][Ar];
  }

  void addAtom(mdtk::Atom& a) {atoms.push_back(a);}
  void buildFromAtom(mdtk::Atom&, NeighbourList& nl,double SPOTTED_DISTANCE);
  bool hasAtom(mdtk::Atom&) const;
  ClassicMolecule():handledElements(),formationTime(-1),escapeTime(-1),
    atoms(/*0*/),atoms_init(),trajectory_dummy(-1)
  {
    initParams();    
  }
  ~ClassicMolecule(){;}
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
  bool  hasFullereneAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isFullereneAtom(atoms[ai])) return true;
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
  bool  hasOnlyFullereneAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isFullereneAtom(atoms[ai])) return false;
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
    return mdtk::academic_round(moleculeMass/mdtk::amu/10)*10;
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
  std::string formula() const
  {
    std::ostringstream os;
    std::map<ElementID, size_t> formula;
    for(size_t ai = 0; ai < atoms.size(); ai++)
      formula[atoms[ai].ID]++;
    std::map<ElementID, size_t>::iterator i = formula.begin();
    while (i != formula.end())
    {
      os << ElementIDtoString(i->first);
      if (i->second > 1)
        os << i->second;
      i++;
    }
    return os.str();
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
//    ERRTRACE("Loading ClassicMolecule....");
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

  friend int operator<(const ClassicMolecule& left, const ClassicMolecule& right);
  friend int operator<(ClassicMolecule& left, ClassicMolecule& right);
};  


inline
int operator<(const ClassicMolecule& left, const ClassicMolecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
int operator<(ClassicMolecule& left, ClassicMolecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}


}

#endif
