/*
   The AtomGroup class (header file).

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

#ifndef mdpp_AtomGroup_hpp
#define mdpp_AtomGroup_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <mdtk/SnapshotList.hpp>
#include <algorithm>
#include <vector>
#include "ClassicMolecule.hpp"

namespace mdepp
{
  using namespace mdtk;

  class Molecule;

class AtomGroup
{
public:
  std::vector<mdtk::Atom> atoms;

  AtomGroup();
  virtual ~AtomGroup();
  AtomGroup(const AtomGroup &c);
  AtomGroup(const ClassicMolecule &c);

  AtomGroup& operator =(const AtomGroup &c);

  virtual void get(std::istream& is);
  virtual void put(std::ostream& os) const;

  void addAtom(const mdtk::Atom& a);
  bool hasAtom(const mdtk::Atom&) const;

  size_t size() const {return atoms.size();};

  void build(const mdtk::SimLoop& ml, 
	     const double SPOTTED_DISTANCE = +1000.0*mdtk::Ao);
  void buildByTag(const mdtk::SimLoop& ml,
                  const unsigned int tag);
  void update(const mdtk::SimLoop& ml);
  void update(
    const mdtk::SnapshotList::SelectedAtomSnapshotList& atomSnapshotList,
    const mdtk::SnapshotList& sn);

  Molecule molecule(size_t atomIndex) const;
  Molecule molecule(const mdtk::Atom& a) const;
  Molecule maxMolecule() const;
  bool isMolecule() const;
  bool isMonomer() const { return atoms.size() == 1;}
  bool isMetalCluster() const;
  bool isFullerene() const;
  mdtk::Vector3D massCenter() const;

  Float          mass() const;
  mdtk::Vector3D velocity() const;
  Float          kineticEnergy() const;

  friend std::istream&  operator>> (std::istream& is, AtomGroup& vec);
  friend std::ostream&  operator<< (std::ostream& os, const AtomGroup& vec);
};

std::istream& operator>> (std::istream& is, AtomGroup& vec);
std::ostream& operator<< (std::ostream& os, const AtomGroup& vec);

}

#endif
