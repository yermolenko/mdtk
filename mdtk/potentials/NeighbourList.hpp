/*
   The NeighbourList class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012
   Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#ifndef mdtk_NeighbourList_hpp
#define mdtk_NeighbourList_hpp

#include <vector>
#include <mdtk/Atom.hpp>
#include <mdtk/AtomsContainer.hpp>
#include <mdtk/Vector3D.hpp>
#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>


#define MDTK_NB_RESERVE_ADD 5

namespace mdtk
{

class FGeneral;

class NeighbourList
{
  const FGeneral* fpot;
  bool ListUpdateRequested;
public:
  Float Rcutoff;
  std::vector<AtomRefsContainer> nl;
  std::vector<Vector3D> displacements;
  NeighbourList(const FGeneral* pot)
   : fpot(pot), ListUpdateRequested(true),
     Rcutoff(0.0),
     nl(), displacements()
  {
  }  

  void init(AtomsArray& atoms)
  {
    nl.clear();
    displacements.clear();
    nl.resize(atoms.size());

    displacements.resize(atoms.size());
    ListUpdateRequested = true;
  }  

  void requestUpdate() {ListUpdateRequested = true;};
  void checkRequestUpdate(AtomsArray& atoms)
  {
    ListUpdateRequested = MovedTooMuch(atoms);
  };
  bool MovedTooMuch(AtomsArray&);
  void Update(AtomsArray&);
  void UpdateIfNeeded(AtomsArray& atoms)
  {
    if (ListUpdateRequested)
    {
      Update(atoms);
      ListUpdateRequested = false;
    }
  }  
  

  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    YAATK_FSTREAM_WRITE(os,ListUpdateRequested,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    YAATK_FSTREAM_READ(is,ListUpdateRequested,smode);
    ListUpdateRequested = true;
  }  
};

}

#endif

