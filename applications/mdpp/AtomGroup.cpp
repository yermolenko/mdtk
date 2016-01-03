/*
   The AtomGroup class.

   Copyright (C) 2010, 2012, 2015 Oleksandr Yermolenko
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

#include "AtomGroup.hpp"
#include "Molecule.hpp"

namespace mdepp
{

AtomGroup::AtomGroup()
 :atoms()
{
}

AtomGroup::~AtomGroup()
{
}

AtomGroup::AtomGroup(const AtomGroup &c)
 :atoms()
{
  atoms = c.atoms;
}

AtomGroup::AtomGroup(const ClassicMolecule &c)
 :atoms()
{
  atoms = c.atoms;
}

AtomGroup& 
AtomGroup::operator =(const AtomGroup &c) 
{
  if (this == &c) return *this;
  atoms = c.atoms;
  return *this;
}

void
AtomGroup::put(std::ostream& os) const
{
  os << atoms.size() << "\n";
  for(size_t i = 0; i < atoms.size(); i++)
    os << atoms[i] << "\n";
}

void
AtomGroup::get(std::istream& is)
{
  size_t sz;
  is >> sz;
  atoms.resize(sz);
  for(size_t i = 0; i < atoms.size(); i++)
    is >> atoms[i];
}

std::ostream&
operator<<(std::ostream& os, const AtomGroup& v)
{
  v.put(os);
  return os;
}

std::istream&
operator>>(std::istream& is, AtomGroup& v)
{
  v.get(is);
  return is;
}

bool
AtomGroup::hasAtom(const mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
      found = true;
  return found;
}  

void
AtomGroup::addAtom(const mdtk::Atom& a)
{
  atoms.push_back(a);
//  atoms.back().container = NULL;
}

void
AtomGroup::update(const mdtk::SimLoop& ml)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    mdtk::Atom& a = atoms[i];
    REQUIRE(ac.size() > a.globalIndex);
    const mdtk::Atom& aUpd = ac[a.globalIndex];
    a.V = aUpd.V;
    a.coords = aUpd.coords;
  }
}

void
AtomGroup::update(
  const mdtk::SnapshotList::SelectedAtomSnapshotList& atomSnapshotList,
  const mdtk::SnapshotList& sn)
{
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    bool atomFoundInSnapshot = false;
    mdtk::Atom& a = atoms[ai];
    for(size_t i = 0; i < atomSnapshotList.size(); i++)
    {
      if (a.globalIndex == sn.atomsSelectedForSaving[i])
      {
        const mdtk::SnapshotList::AtomSnapshot& as = atomSnapshotList[i];
        as.restoreToAtom(a);
        atomFoundInSnapshot = true;
        break;
      }
    }
    REQUIRE(atomFoundInSnapshot);
  }
}

void
AtomGroup::build(const mdtk::SimLoop& ml, 
		 const double SPOTTED_DISTANCE)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = ac[i];
    if (a.coords.z < SPOTTED_DISTANCE)
      addAtom(a);
  }
}  

void
AtomGroup::buildByTag(const mdtk::SimLoop& ml,
                      const unsigned int tag)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = ac[i];
    if (a.hasTag(tag))
      addAtom(a);
  }
}

Molecule
AtomGroup::molecule(size_t atomIndex) const
{
  REQUIRE(atomIndex < atoms.size());
  Molecule m;
  m.buildFromAtom(atoms[atomIndex],*this);
  return m;
}

Molecule
AtomGroup::molecule(const mdtk::Atom& a) const
{
  REQUIRE(hasAtom(a));
  return molecule(a.globalIndex);
}

bool
AtomGroup::isMolecule() const
{
  return maxMolecule().size() == atoms.size();
}

bool
AtomGroup::isMetalCluster() const
{
  if (!isMolecule()) return false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].ID != Cu_EL &&
        atoms[i].ID != Ag_EL &&
        atoms[i].ID != Au_EL) return false;

  return true;
}

bool
AtomGroup::isFullerene() const
{
  if (!isMolecule()) return false;
  if (atoms.size() != 60) return false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].ID != C_EL) return false;

  return true;
}

Molecule
AtomGroup::maxMolecule() const
{
  Molecule mMax;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Molecule m = molecule(i);
    if (m.size() > mMax.size())
      mMax = m;
  }
  return mMax;
}

Float
AtomGroup::mass() const
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    moleculeMass += atom.M;
  } 
  return moleculeMass;
}

mdtk::Vector3D
AtomGroup::velocity() const
{
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

mdtk::Vector3D
AtomGroup::massCenter() const
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.coords*atom.M;
  };
  return sumOfP/sumOfM;    
}  

Float
AtomGroup::kineticEnergy() const
{
  Float e = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    e += atom.M*SQR(atom.V.module())/2.0;
  }
  return e;
}

}
