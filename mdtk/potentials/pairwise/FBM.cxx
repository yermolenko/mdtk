/*
   The Born-Mayer interatomic potential implementation.

   Copyright (C) 2004, 2005, 2009, 2012, 2013 Oleksandr Yermolenko
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

#include "FBM.hpp"
#include <iostream>

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt;
using std::pow;

FBM::FBM(Rcutoff rcutoff):
  FPairwise(rcutoff)
{
  handledElements.insert(Cu_EL);
  handledElements.insert(Ar_EL);

  handledElementPairs.insert(std::make_pair(Ar_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,Ar_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,Cu_EL));

  A3[Cu_EL][Cu_EL] = 22.565*1000.0*eV;
    A3[Cu_EL][Cu_EL] = A3[Cu_EL][Cu_EL];
  A4[Cu_EL][Cu_EL] = 50.88/(10.0*Ao);
    A4[Cu_EL][Cu_EL] = A4[Cu_EL][Cu_EL];

  A3[Cu_EL][Ar_EL] = 59.874*1000.0*eV;
    A3[Ar_EL][Cu_EL] = A3[Cu_EL][Ar_EL];
  A4[Cu_EL][Ar_EL] = 72.0/(10.0*Ao);
    A4[Ar_EL][Cu_EL] = A4[Cu_EL][Ar_EL];
}

Float
FBM::F11(AtomsPair& ij, const Float V)
{
  Float R = ij.r();

  if (V != 0.0)
  {
    Float Der = -A3[ij.atom1.ID][ij.atom2.ID]*A4[ij.atom1.ID][ij.atom2.ID]*exp(-A4[ij.atom1.ID][ij.atom2.ID]*R);

    ij.r(Der*V);
  }

  return  A3[ij.atom1.ID][ij.atom2.ID]*exp(-A4[ij.atom1.ID][ij.atom2.ID]*R);
}

Float
FBM::operator()(AtomsArray& gl)
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
