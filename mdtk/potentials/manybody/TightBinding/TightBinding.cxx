/*
   Implementation of the many-body interatomic potential for copper.
   See [G. Betz, W. Husinsky, Nucl. Instr. and Meth. B 102, 281 (1995)]

   Copyright (C) 2006, 2007, 2008, 2009, 2012, 2013 Oleksandr
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

#include "TightBinding.hpp"
#include <algorithm>

namespace mdtk
{

Float
TightBinding::operator()(AtomsArray& gl)
{
  Float Ei = 0;
  for(size_t ii = 0; ii < gl.size(); ii++)
  {
    Atom &atom_i = gl[ii];
    if (isHandled(atom_i))
    {
      Ei += F(atom_i);
    }

    for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
    {
      Atom &atom_j = *(NL(atom_i)[jj]);
      if (atom_i.globalIndex > atom_j.globalIndex) continue;
      if (&atom_i != &atom_j)
      {
        if (!probablyAreNeighbours(atom_i,atom_j)) continue;
        AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));
        Ei += Phi(ij);
      }
    }
  }

  return Ei;
}

inline
Float
TightBinding::Phi(AtomsPair& ij)
{
  Float fvar = ij.f();

#ifdef TightBinding_OPTIMIZED
  if (fvar == 0.0) return 0.0;
#endif

  Float r = ij.r();

#ifndef  EAM_HANDLE_SHORTRANGE
  Spline& spline = *(this->spline);
  if (r < spline.x1())
  {
// if (V != 0)
    {
      Float Der = -BM_B*BM_A*exp(-BM_B*r);
      ij.r(Der);
    }
    return  BM_A*exp(-BM_B*r);
  }
  else
  {
    if (r < spline.x2())
    {
//    if (V != 0)
      ij.r(spline.der(r));
      return spline(r);
    }
  }
#endif

  Float Val = Phi0_*exp(-alpha_*r);

// if (V != 0)
  {
    Float Der = Val*(-alpha_);
    ij.r(Der*fvar);
    ij.f(Val);
  }

  return fvar*Val;
}

inline
Float
TightBinding::g(AtomsPair& ij, const Float V)
{
  Float fvar = ij.f();

#ifdef TightBinding_OPTIMIZED
  if (fvar == 0.0) return 0.0;
#endif

  Float r = ij.r();

#ifndef  EAM_HANDLE_SHORTRANGE
//  if (r < 1.029*Ao) return 0.0;
#endif

  Float Val = exp(-beta_*r);

// if (V != 0)
  {
    Float Der = Val*(-beta_);
    ij.r(Der*fvar*V);
    ij.f(Val*V);
  }

  return fvar*Val;
}

inline
Float
TightBinding::F(Atom &atom1)
{
  Float rhovar = rho(atom1);
  REQUIRE(rhovar >= 0.0);
  if (rhovar != 0)
  {
    rho(atom1,-c_/(2.0*sqrt(rhovar)));
  }
  return -c_*sqrt(rhovar);
}

Float
TightBinding::rho(Atom &atom_i, const Float V)
{
  Float rhoij = 0.0;
  for(size_t j = 0; j < NL(atom_i).size(); j++)
  {
    Atom& atom_j = *(NL(atom_i)[j]);
    if (!probablyAreNeighbours(atom_i,atom_j)) continue;
    AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));
    {
      rhoij += g(ij,V);
    }
  }
  return rhoij;
}

TightBinding::TightBinding():
  FManybody()
{
  handledElements.insert(Cu_EL);
  handledElementPairs.insert(std::make_pair(Cu_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,DUMMY_EL));
  handledElementPairs.insert(std::make_pair(DUMMY_EL,Cu_EL));

  setupPotential();

  nl.Rcutoff = getRcutoff();
}

void
TightBinding::setupPotential()
{
  alpha_ = 42.87/(10.0*Ao);
  beta_  = 18.00/(10.0*Ao);
  c_     = 12.17*eV;
  Phi0_  = 9.892 *1000.0 *eV;
  R_[0]  = 5.0*Ao;
  R_[1]  = 5.5*Ao;

  BM_A = 22.565*1000.0*eV;
  BM_B = 50.88/(10.0*Ao);
  fillR_concat_();

  PRINT("TightBinding interatomic potential configured.\n");
}

void
TightBinding::fillR_concat_()
{
  Float r;

  Atom atom1; atom1.ID = Cu_EL; atom1.setAttributesByElementID();
  Atom atom2; atom2.ID = Cu_EL; atom2.setAttributesByElementID();

  Float       x[2]; x[0] = 1.0*Ao; x[1] = 1.2*Ao;
  Float       v[2];
  Float    dvdx[2];
//  Float    d2vdxdx[2];
//  d2vdxdx[0] = 0;
//  d2vdxdx[1] = 0;

  REQUIRE(x[1] < R_[0]);

  {
    r = x[0];

    Float VShortRange = 0.0;
    Float DerVShortRange = 0.0;
    {
      VShortRange=   BM_A*exp(-BM_B*r);
      DerVShortRange = -BM_B*BM_A*exp(-BM_B*r);
    }
    v[0] = VShortRange;
    dvdx[0] = DerVShortRange;

    r = x[1];

    Float VLongRange = 0.0;
    Float DerVLongRange = 0.0;
    {
      VLongRange = Phi0_*exp(-alpha_*r);
      DerVLongRange = -alpha_*Phi0_*exp(-alpha_*r);
    }

    v[1] = VLongRange;
    dvdx[1] = DerVLongRange;
  }

  spline = new Spline(x,v,dvdx/*,d2vdxdx*/);
}

}
