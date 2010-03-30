/*
   The Ziegler-Biersack-Littmark interatomic potential implementation.

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

#include "FBZL.hpp"
#include <iostream>

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt; 
using std::pow; 

FBZL::FBZL(Rcutoff rcutoff):
  FPairwise(rcutoff),
  AB_(0.53e-8)
{
  handledElementPairs.insert(std::make_pair(Ar_EL,C_EL));
  handledElementPairs.insert(std::make_pair(Ar_EL,H_EL));
  handledElementPairs.insert(std::make_pair(C_EL,Ar_EL));
  handledElementPairs.insert(std::make_pair(H_EL,Ar_EL));

  handledElementPairs.insert(std::make_pair(Ar_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,Ar_EL));
  handledElementPairs.insert(std::make_pair(Ar_EL,Ar_EL));
}

Float
FBZL::F11(Atom& a1, Atom& a2)
{
  Vector3D r = r_vec(a1,a2);
  Float R = r.module();

  if (R > getRcutoff()) return 0.0;

  Float ZA = a1.Z; Float ZB = a2.Z;

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));

  Float Y=R/AS;

  return  ZA*ZB/R*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y));
}

Vector3D
FBZL::dF11(Atom& a1, Atom& a2, Atom& da)
{
  Vector3D r = r_vec(a1,a2);
  Float R = r.module();
  Float ZA = a1.Z; Float ZB = a2.Z;

  if (R > getRcutoff()) return Vector3D(0.0,0.0,0.0);

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));
  Float Y=R/AS;
  Float Der =
          -ZA*ZB/(R*R)*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y))-
          ZA*ZB/(R*AS)*(0.18175*3.1998*exp(-3.1998*Y)+
          0.50986*0.94229*exp(-0.94229*Y)+0.28022*0.4029*exp(-0.4029*Y)+
          0.02817*0.20162*exp(-0.20162*Y));

  return  Der*dr_vec_module(a1,a2,da);
}

Float
FBZL::operator()(AtomsContainer& gl)
{
  Float Ei = 0;
for(size_t i = 0; i < gl.size(); i++)
{
  Atom& atom = *(gl[i]);
  AtomsContainer& nl = NL(atom);
  Index j;
  for(j = 0; j < nl.size(); j++)
  {
    Atom &atom_j = *(nl[j]);
    if (atom.globalIndex > atom_j.globalIndex) continue;
    if (isHandledPair(atom,atom_j))
    if (&atom != &atom_j)
    {
      Ei += F11(atom,atom_j)*f(atom,atom_j);
    }  
  }  
}  
  return Ei;
}

Vector3D
FBZL::grad(Atom &atom,AtomsContainer&)
{
  Index j;
  Vector3D dEi(0.0,0.0,0.0);
  for(j = 0; j < NL(atom).size(); j++)
  {
    Atom &atom_j = *(NL(atom)[j]);
    if (isHandledPair(atom,atom_j))
    if (&atom != &atom_j)
    {
      dEi += dF11(atom,atom_j,atom)*f(atom,atom_j)
             +F11(atom,atom_j)    *df(atom,atom_j,atom);
    }  
  }  
  return dEi;
}


} // namespace mdtk

