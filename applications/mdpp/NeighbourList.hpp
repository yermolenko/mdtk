/*
   The MDPP NeighbourList class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Oleksandr
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

#ifndef mdpp_NeighbourList_hpp
#define mdpp_NeighbourList_hpp

#include <vector>
#include <mdtk/Atom.hpp>
#include <mdtk/Vector3D.hpp>
#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>
#include <mdtk/SimLoop.hpp>

namespace mdepp
{
class NeighbourList
{
public:
  mdtk::Float Rc;
  std::vector<mdtk::AtomsContainer> nl;
  NeighbourList(mdtk::AtomsContainer& atoms, const mdtk::Float Rc_ = 5.0*mdtk::Ao)
    : Rc(Rc_), nl(atoms.size())
  {
    update(atoms);
  }  

  void update(mdtk::AtomsContainer&);
};

}

#endif

