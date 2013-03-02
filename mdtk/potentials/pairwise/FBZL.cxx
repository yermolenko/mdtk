/*
   The Ziegler-Biersack-Littmark interatomic potential implementation.

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
  std::vector<ElementID> elements;
  elements.push_back(H_EL);
  elements.push_back(C_EL);
  elements.push_back(Cu_EL);
  elements.push_back(Au_EL);
  elements.push_back(Ag_EL);

  std::vector<ElementID> ions;
  ions.push_back(Ar_EL);
  ions.push_back(Xe_EL);

  for(size_t i = 0; i < elements.size(); ++i)
    handledElements.insert(elements[i]);

  for(size_t i = 0; i < ions.size(); ++i)
    handledElements.insert(ions[i]);

  for(size_t i = 0; i < ions.size(); ++i)
    for(size_t j = 0; j < ions.size(); ++j)
    {
      handledElementPairs.insert(std::make_pair(ions[i],ions[j]));
      handledElementPairs.insert(std::make_pair(ions[j],ions[i]));
    }

  for(size_t i = 0; i < ions.size(); ++i)
    for(size_t j = 0; j < elements.size(); ++j)
    {
      handledElementPairs.insert(std::make_pair(ions[i],elements[j]));
      handledElementPairs.insert(std::make_pair(elements[j],ions[i]));
    }
}

Float
FBZL::F11(AtomsPair& ij, const Float V)
{
  Float R = ij.r();

  if (R > getRcutoff()) return 0.0;

  Float ZA = ij.atom1.Z; Float ZB = ij.atom2.Z;

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));

  Float Y=R/AS;

  if (V != 0.0)
  {
    Float Der =
          -ZA*ZB/(R*R)*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y))-
          ZA*ZB/(R*AS)*(0.18175*3.1998*exp(-3.1998*Y)+
          0.50986*0.94229*exp(-0.94229*Y)+0.28022*0.4029*exp(-0.4029*Y)+
          0.02817*0.20162*exp(-0.20162*Y));

    ij.r(Der*V);
  }

  return  ZA*ZB/R*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y));
}

Float
FBZL::operator()(AtomsArray& gl)
{
  Float Ei = 0;
for(size_t i = 0; i < gl.size(); i++)
{
  Atom& atom = gl[i];
  AtomRefsContainer& nl = NL(atom);
  for(size_t j = 0; j < nl.size(); j++)
  {
    Atom &atom_j = *(nl[j]);
    if (atom.globalIndex > atom_j.globalIndex) continue;
    if (isHandledPair(atom,atom_j))
    if (&atom != &atom_j)
    {
      if (!probablyAreNeighbours(atom,atom_j)) continue;
      AtomsPair ij(atom,atom_j,R(0),R(1));
      Float f = ij.f();
      Float F11Val = F11(ij,f);
//      if (V != 0)
      ij.f(F11Val);
      Ei += F11Val*f;
    }
  }
}
  return Ei;
}

} // namespace mdtk
