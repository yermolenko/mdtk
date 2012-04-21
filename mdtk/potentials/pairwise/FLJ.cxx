/*
   The Lennard-Jones interatomic potential implementation.

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

#include "FLJ.hpp"
#include <iostream>

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt;
using std::pow;

FLJ::FLJ(Rcutoff rcutoff):
  FPairwise(rcutoff)
{
  handledElements.insert(Cu_EL);
  handledElements.insert(C_EL);
  handledElements.insert(H_EL);

  handledElementPairs.insert(std::make_pair(Cu_EL,C_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,H_EL));
  handledElementPairs.insert(std::make_pair(C_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(H_EL,Cu_EL));

//  See [D.E. Ellis et al., Materials Science in Semiconductor
//  Processing 3, 123 (2000)]

  Float C6_CuC = 41.548*eV*pow(Ao,6);
  Float C12_CuC = 2989.105*eV*pow(Ao,12);

  sigma_[Cu][Cu] = 0.0*Ao;
  sigma_[C][C]   = 0.0*Ao;
  sigma_[Cu][C]  = pow(C12_CuC/C6_CuC,1.0/6.0);
    sigma_[C][Cu] = sigma_[Cu][C];
  sigma_[Cu][Cu] = 0.0*Ao;
  sigma_[H][H]   = 0.0*Ao;
  sigma_[Cu][H]  = 2.049067052*Ao;
    sigma_[H][Cu] = sigma_[Cu][H];

  zeta_[Cu][Cu] = 0.0*eV;
  zeta_[C][C]   = 0.0*eV;
  zeta_[Cu][C]  = pow(C6_CuC,2)/(4.0*C12_CuC);
    zeta_[C][Cu] = zeta_[Cu][C];
  zeta_[Cu][Cu] = 0.0*eV;
  zeta_[H][H]   = 0.0*eV;
  zeta_[Cu][H]  = 0.01*eV;
    zeta_[H][Cu] = zeta_[Cu][H];

  fillR_concat_();
}

void
FLJ::fillR_concat_()
{
  Float r;

  AtomsArray atoms(3);
  atoms[0].ID = Cu_EL;
  atoms[1].ID = C_EL;
  atoms[2].ID = H_EL;
  atoms.setAttributesByElementID();

for(size_t i = 0; i < atoms.size(); i++)
for(size_t j = 0; j < atoms.size(); j++)
{
Atom& atom1 = atoms[i];
Atom& atom2 = atoms[j];

         Float       x[2]; x[0] = 1.2*Ao; x[1] = 1.9*Ao;
         Float       v[2];
         Float    dvdx[2];
//         Float    d2vdxdx[2];
//d2vdxdx[0] = 0;
//d2vdxdx[1] = 0;

{
  AtomsPair ij(atom1,atom2,false);

  r = x[0];

  Float VZBL = 0.0;
  Float DerVZBL = 0.0;
  {
  Float AB_ = 0.53e-8;

  Float ZA = atom1.Z; Float ZB = atom2.Z;

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));

  Float Y=r/AS;

  VZBL=   ZA*ZB/r*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y));

  DerVZBL =
          -ZA*ZB/(r*r)*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y))-
          ZA*ZB/(r*AS)*(0.18175*3.1998*exp(-3.1998*Y)+
          0.50986*0.94229*exp(-0.94229*Y)+0.28022*0.4029*exp(-0.4029*Y)+
          0.02817*0.20162*exp(-0.20162*Y));

  }
  v[0] = VZBL;
  dvdx[0] = DerVZBL/* = -1.5*/;

  r = x[1];

  Float VLJ = 0.0;
  Float DerVLJ = 0.0;
  {
    Float sigma_ij=     this->sigma(ij);
    Float zeta_ij =     this->zeta(ij);

    Float s_div_r    = sigma_ij/r;
    Float s_div_r_6  = s_div_r;
    int i;
    for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
    Float s_div_r_12 = s_div_r_6*s_div_r_6;
    DerVLJ = 4.0*zeta_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
    VLJ = 4.0*zeta_ij*(s_div_r_12-s_div_r_6);
  }

  v[1] = VLJ;
  dvdx[1] = DerVLJ/* = 0*/;
}

splines[e2i(atom1)][e2i(atom2)] = new Spline(x,v,dvdx/*,d2vdxdx*/);

}

}

Float
FLJ::VLJ(AtomsPair& ij)
{
  Float r = ij.r();

#ifndef  LJ_HANDLE_SHORTRANGE
  Spline& spline = *(splines[e2i(ij.atom1)][e2i(ij.atom2)]);
  if (r < spline.x1())
  {
  Float AB_ = 0.53e-8;

  Float ZA = ij.atom1.Z; Float ZB = ij.atom2.Z;

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));

  Float Y=r/AS;

//  if (V != 0.0)
  {
    Float Der =
          -ZA*ZB/(r*r)*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y))-
          ZA*ZB/(r*AS)*(0.18175*3.1998*exp(-3.1998*Y)+
          0.50986*0.94229*exp(-0.94229*Y)+0.28022*0.4029*exp(-0.4029*Y)+
          0.02817*0.20162*exp(-0.20162*Y));

    ij.r(Der);
  }

  return  ZA*ZB/r*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y));
  }
  else
  {
    if (r < spline.x2())
    {
//  if (V != 0.0)
      ij.r(spline.der(r));

      return spline(r);
    }
  }
#endif

  Float sigma_ij=     this->sigma(ij);
  Float zeta_ij =     this->zeta(ij);

  Float s_div_r    = sigma_ij/r;

  Float s_div_r_6  = s_div_r;
  int i;
  for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
  Float s_div_r_12 = s_div_r_6*s_div_r_6;

  Float Val = 4.0*zeta_ij*(s_div_r_12-s_div_r_6);


  Float f = ij.f();

// if (V != 0)
  {
    Float Der = 4.0*zeta_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
    ij.r(Der*f)+ij.f(Val);
  }

  return f*Val;
}

Float
FLJ::operator()(AtomsArray& gl)
{
  Float Ei = 0;
for(size_t i = 0; i < gl.size(); i++)
{
  Atom& atom = gl[i];
  AtomRefsContainer& nl = NL(atom);
  Index j;
  for(j = 0; j < nl.size(); j++)
  {
    Atom &atom_j = *(nl[j]);
    if (atom.globalIndex > atom_j.globalIndex) continue;
    if (isHandledPair(atom,atom_j))
    if (&atom != &atom_j)
    {
      if (!probablyAreNeighbours(atom,atom_j)) continue;
      AtomsPair ij(atom,atom_j,R(0),R(1));
      Ei += VLJ(ij);
    }
  }
}
  return Ei;
}

} // namespace mdtk
