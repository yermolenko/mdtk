/*
   The MDPP NeighbourList class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2012 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include "NeighbourList.hpp"

#include "ClassicMolecule.hpp"

namespace mdepp
{

  void
  NeighbourList::update(mdtk::AtomsArray& atoms)
  {
    REQUIRE(Rc > 0.0);
    size_t i,j,N;
    N = atoms.size();
    for(i = 0; i < N; ++i)
    {
      mdtk::Atom& atom_i = atoms[i];
      if (
	atom_i.coords.z > 3.615*mdtk::Ao
	&& !(atom_i.hasTag(ATOMTAG_FULLERENE))
	&& !(atom_i.hasTag(ATOMTAG_CLUSTER))
	) continue;
      mdtk::AtomRefsContainer& nl_ = nl[i];
    
      nl_.clear();
      nl_.reserve(50);

      for(j = 0; j < N; ++j)
      {
	mdtk::Atom& atom_j = atoms[j];
	if (j != i && depos(atom_i,atom_j).module() < Rc) 
          nl_.push_back(&atom_j);
      }
    }
  }

}

