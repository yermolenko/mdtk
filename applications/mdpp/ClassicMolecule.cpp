/*
   The Simple Molecule class.

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

#include "ClassicMolecule.hpp"

namespace mdepp
{

bool
ClassicMolecule::hasAtom(mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
    {
      found = true;
      break;
    }
  return found;
}  

void
ClassicMolecule::buildFromAtom(mdtk::Atom& a, NeighbourList& nl,double SPOTTED_DISTANCE)
{
//TRACE("adding atom");
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);

//  TRACE(a.fixed);
//  TRACE(atoms[atoms.size()-1].fixed);
  

//TRACE(a.globalIndex);  
//TRACE(a.ID);  
//TRACE("1");  
//TRACE(a.nl(&(ml.fpot)).size());
  for(size_t i = 0; i < nl.nl[a.globalIndex].size(); i++)
  {
//TRACE("2");
    mdtk::Atom& nb_a = *nl.nl[a.globalIndex][i];
//    nb_a.container = &dummy_ac;
//       a.container = &dummy_ac;
    if (!isHandled(nb_a)) continue;
    Float distance = depos(a,nb_a).module();//sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
    if (Rc(a,nb_a) >= distance)
    {
      if (nb_a.coords.z < SPOTTED_DISTANCE)
      {
//TRACE("Next dive");        
        buildFromAtom(nb_a,nl,SPOTTED_DISTANCE);
      }
      else
      {
//TRACE("atoms.clear");        
        atoms.clear();
      }  
      if (atoms.size() == 0) break;
    }  
  }  
}  

}
