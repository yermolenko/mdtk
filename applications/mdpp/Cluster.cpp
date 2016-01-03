/*
   The Cluster class.

   Copyright (C) 2010, 2012 Oleksandr Yermolenko
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

#include "Cluster.hpp"

namespace mdepp
{

Cluster::Cluster()
 :AtomGroup()
{
}

Cluster::~Cluster()
{
}

Cluster::Cluster(const Cluster &c)
 :AtomGroup(c)
{
}

Cluster& 
Cluster::operator =(const Cluster &c) 
{
  if (this == &c) return *this;
  AtomGroup::operator=(c);
  return *this;
}

void
Cluster::put(std::ostream& os) const
{
  AtomGroup::put(os);
}

void
Cluster::get(std::istream& is)
{
  AtomGroup::get(is);
}

void
Cluster::build(const mdtk::SimLoop& ml)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = ac[i];
    if (a.ID == Cu_EL && a.coords.z < -3.615*Ao*1.5)
      addAtom(a);
    if (a.ID == Au_EL && a.coords.z < -3.615*Ao*1.5)
      addAtom(a);
  }
//  REQUIRE(atoms.size() == 13);
}  

}
