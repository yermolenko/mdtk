/*
   The NeighbourList class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr
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

#include <mdtk/potentials/NeighbourList.hpp>
#include <mdtk/potentials/FGeneral.hpp>

namespace mdtk
{

#define NLSKIN_FACTOR 0.44

void
NeighbourList::Update(AtomsContainer& atoms_)
{
  REQUIRE(Rcutoff > 0.0);
  PVLOG("NL update\n");
  Float range_squared,dij_squared;
  range_squared = SQR((1.0+NLSKIN_FACTOR)*Rcutoff);
  int i,j,N;
  N = atoms_.size();
  for(i = 0; i < N; i++)
  {
    AtomsContainer& nl_ = nl[i];
    AtomsContainer& nl_with_self_ = nl_with_self[i];
    
    Atom& atom_i = *(atoms_[i]);
    size_t nl_size_prev = nl_.size();
    nl_.clear();
    nl_.reserve(nl_size_prev+MDTK_NB_RESERVE_ADD /*50+1*/);
    size_t nl_with_self_size_prev = nl_with_self_.size();
    nl_with_self_.clear();
    nl_with_self_.reserve(nl_with_self_size_prev+MDTK_NB_RESERVE_ADD /*50+1*/);

    for(j = 0; j < N; j++)
    {
      Atom& atom_j = *(atoms_[j]);
      if (j != i)
      {
        dij_squared = fpot->r_vec_module_squared_no_touch(atom_i,atom_j);
        if (dij_squared < range_squared) 
        {
          nl_.push_back(&atom_j);
          nl_with_self_.push_back(&atom_j);
        }
      }
      else
      {
        nl_with_self_.push_back(&atom_j);
      }  
    }
    
    displacements[i] = Vector3D(0,0,0);
  }
}

bool
NeighbourList::MovedTooMuch(AtomsContainer& atoms_)
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

