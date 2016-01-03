/*
   The Fragment class.

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

#include "Fragment.hpp"

namespace mdepp
{

void
Fragment::buildFromAtom(const mdtk::Atom& a, const ClusterDynamics& cd)
{
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);

  for(size_t atomIndex = 0; atomIndex < cd.atomTrajectories.size(); atomIndex++)
  {
    const mdtk::Atom& nb_a = cd.atomTrajectories[atomIndex].endCheckPoint;
    if (!isHandled(nb_a)) continue;
    Float distance = (a.coords-nb_a.coords).module();//REF_POT_OF(ml.fpot)->r_vec_module_no_touch(a,nb_a);//sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
//    if (/*Rc(a,nb_a)*/1.50*(2.46*mdtk::Ao*cos(30.0*mdtk::Deg)*2.0) >= distance)
//    if (/*Rc(a,nb_a)*/1.50*(2.46*mdtk::Ao) >= distance)
//    if (/*Rc(a,nb_a)*/1.50*(3.61*mdtk::Ao/2.0) >= distance)
    if (/*Rc(a,nb_a)*/1.00*(3.61*mdtk::Ao) >= distance)
    {
//      if (nb_a.coords.z < SPOTTED_DISTANCE)
//      {
        buildFromAtom(nb_a,cd);
//      }
//      else
//      {
//        atoms.clear();
//      }  
//      if (atoms.size() == 0) break;
    }  
  }  
}  



}
