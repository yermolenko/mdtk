/*
   The NeighbourList class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013
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

#include <mdtk/potentials/NeighbourList.hpp>
#include <mdtk/potentials/FGeneral.hpp>

namespace mdtk
{

#define NLSKIN_FACTOR 0.44

void
NeighbourList::Update(AtomsArray& atoms_, std::vector<NeighbourList*>& nlObjectsToUpdate)
{
  if (nlObjectsToUpdate.size() == 0)
    return;

  {
    int NL2UPD = nlObjectsToUpdate.size();
    TRACE(NL2UPD);
  }

  size_t N = atoms_.size();

  for(size_t i = 0; i < N; i++)
  {
    for(size_t nloi = 0; nloi < nlObjectsToUpdate.size(); ++nloi)
    {
      NeighbourList& nlObject = *(nlObjectsToUpdate[nloi]);

      nlObject.displacements[i] = Vector3D(0,0,0);

      AtomRefsContainer& nl_ = nlObject.nl[i];

      size_t nl_size_prev = nl_.size();
      nl_.clear();
      nl_.reserve(nl_size_prev+MDTK_NB_RESERVE_ADD);
    }
  }

  Float range_squared_max = 0.0;
  for(size_t nloi = 0; nloi < nlObjectsToUpdate.size(); ++nloi)
  {
    NeighbourList& nlObject = *(nlObjectsToUpdate[nloi]);

    Float range_squared = SQR((1.0+NLSKIN_FACTOR)*nlObject.Rcutoff);

    if (range_squared_max < range_squared)
      range_squared_max = range_squared;
  }
  REQUIRE(range_squared_max > 0.0);
//  TRACE(sqrt(range_squared_max)/Ao);

  for(size_t i = 0; i < N; i++)
  {
    Atom& atom_i = atoms_[i];

    for(size_t j = i+1; j < N; j++)
    {
      Atom& atom_j = atoms_[j];

      Float dij_squared = depos(atom_i,atom_j).module_squared();

      if (dij_squared > range_squared_max)
        continue;

      for(size_t nloi = 0; nloi < nlObjectsToUpdate.size(); ++nloi)
      {
        NeighbourList& nlObject = *(nlObjectsToUpdate[nloi]);

        if (!nlObject.fpot->isHandled(atom_i)) continue;
        if (!nlObject.fpot->isHandled(atom_j)) continue;
//      if (!nlObject.fpot->isHandledPair(atom_i,atom_j)) continue;

//        REQUIRE(nlObject.Rcutoff > 0.0);
//        PRINT("NL update\n");
        Float range_squared = SQR((1.0+NLSKIN_FACTOR)*nlObject.Rcutoff);

        if (dij_squared < range_squared)
        {
          nlObject.nl[i].push_back(&atom_j);
          nlObject.nl[j].push_back(&atom_i);
        }
      }
    }
  }
}

bool
NeighbourList::MovedTooMuch(AtomsArray& atoms_)
{
  REQUIRE(Rcutoff > 0.0);

  Float disp1, disp2 ,disp;
  disp1 = 0.0;
  disp2 = 0.0;
  int i, atoms_count;
  atoms_count = atoms_.size();
  for (i = 0; i < atoms_count; i++)
  {
   disp = displacements[i].module();
   if (disp >= disp1)
   {        
      disp2 = disp1;
      disp1 = disp;                
   }
   else
   if (disp >= disp2)
   {   
      disp2 = disp;
   }
  }
  return (disp1 + disp2 > NLSKIN_FACTOR*Rcutoff);
}


}

