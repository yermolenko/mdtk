/*
   The Lennard-Jones interatomic potential implementation.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013,
   2015 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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
  FPairwise(rcutoff),
  splines()
{
  std::vector<ElementID> elements1;
  elements1.push_back(H_EL);
  elements1.push_back(C_EL);

  std::vector<ElementID> elements2;
  elements2.push_back(Cu_EL);
//  elements2.push_back(Ag_EL);
  elements2.push_back(Au_EL);

  for(size_t i = 0; i < elements1.size(); ++i)
    handledElements.insert(elements1[i]);

  for(size_t i = 0; i < elements2.size(); ++i)
    handledElements.insert(elements2[i]);

  for(size_t i = 0; i < elements1.size(); ++i)
    for(size_t j = 0; j < elements2.size(); ++j)
    {
      handledElementPairs.insert(std::make_pair(elements1[i],elements2[j]));
      handledElementPairs.insert(std::make_pair(elements2[j],elements1[i]));
    }

  if (0) // disable this set of parameters
  {
//  See [D.E. Ellis et al., Materials Science in Semiconductor
//  Processing 3, 123 (2000)]

    Float C6_CuC = 41.548*eV*pow(Ao,6);
    Float C12_CuC = 2989.105*eV*pow(Ao,12);

    sigma_[Cu][Cu] = 0.0*Ao;
    sigma_[C][C]   = 0.0*Ao;
    sigma_[Cu][C]  = pow(C12_CuC/C6_CuC,1.0/6.0); // ~ 2.039*Ao
      sigma_[C][Cu] = sigma_[Cu][C];
    sigma_[H][H]   = 0.0*Ao;
    sigma_[Cu][H]  = 2.049067052*Ao;
      sigma_[H][Cu] = sigma_[Cu][H];

    epsilon_[Cu][Cu] = 0.0*eV;
    epsilon_[C][C]   = 0.0*eV;
    epsilon_[Cu][C]  = pow(C6_CuC,2)/(4.0*C12_CuC); // ~ 0.14438*eV
      epsilon_[C][Cu] = epsilon_[Cu][C];
    epsilon_[H][H]   = 0.0*eV;
    epsilon_[Cu][H]  = 0.01*eV;
      epsilon_[H][Cu] = epsilon_[Cu][H];

//  See [Arnaud Delcorte, Barbara J. Garrison, Nuclear Instruments and
//  Methods in Physics Research B 269, 1572 (2011)]

    sigma_[Au][Au] = 0.0*Ao;
    sigma_[C][C]   = 0.0*Ao;
    sigma_[Au][C]  = 3.172*Ao;
      sigma_[C][Au] = sigma_[Au][C];
    sigma_[H][H]   = 0.0*Ao;
    sigma_[Au][H]  = 2.746*Ao;
      sigma_[H][Au] = sigma_[Au][H];

    epsilon_[Au][Au] = 0.0*eV;
    epsilon_[C][C]   = 0.0*eV;
    epsilon_[Au][C]  = 0.00277*eV;
      epsilon_[C][Au] = epsilon_[Au][C];
    epsilon_[H][H]   = 0.0*eV;
    epsilon_[Au][H]  = 0.00179*eV;
      epsilon_[H][Au] = epsilon_[Au][H];
  }

  {
//  Lorentz–Berthelot mixing rules 1

//  S.-P. Huang et al. / Surface Science 545 (2003) 163
    sigma_[Cu][Cu] = 3.05*Ao;
    epsilon_[Cu][Cu] = 0.165565*eV;
    sigma_[C][C]   = 3.4*Ao;
    epsilon_[C][C]   = 0.002413*eV;
    sigma_[Au][Au] = 2.569*Ao;
    epsilon_[Au][Au] = 0.408*eV;

//  M.P. Allen, D.J. Tildesley, 1987
    sigma_[H][H]   = 2.81*Ao;
    epsilon_[H][H]   = 8.6*K*kb;

    sigma_[Cu][C]  = (sigma_[Cu][Cu] + sigma_[C][C])/2; // should be ~ 3.225*Ao
    sigma_[C][Cu] = sigma_[Cu][C];
    sigma_[Cu][H]  = (sigma_[Cu][Cu] + sigma_[H][H])/2;
    sigma_[H][Cu] = sigma_[Cu][H];

    epsilon_[Cu][C]  = sqrt(epsilon_[Cu][Cu]*epsilon_[C][C]); // should be ~ 0.019996*eV
    epsilon_[C][Cu] = epsilon_[Cu][C];
    epsilon_[Cu][H]  = sqrt(epsilon_[Cu][Cu]*epsilon_[H][H]);
    epsilon_[H][Cu] = epsilon_[Cu][H];

    sigma_[Au][C]  = (sigma_[Au][Au] + sigma_[C][C])/2; // should be ~ 2.985*Ao
    sigma_[C][Au] = sigma_[Au][C];
    sigma_[Au][H]  = (sigma_[Au][Au] + sigma_[H][H])/2;
    sigma_[H][Au] = sigma_[Au][H];

    epsilon_[Au][C]  = sqrt(epsilon_[Au][Au]*epsilon_[C][C]); // should be ~ 0.033244*eV
    epsilon_[C][Au] = epsilon_[Au][C];
    epsilon_[Au][H]  = sqrt(epsilon_[Au][Au]*epsilon_[H][H]);
    epsilon_[H][Au] = epsilon_[Au][H];

    if (0)
    {
      TRACE(sigma_[C][Cu]/Ao);
      TRACE(epsilon_[C][Cu]/eV);
      TRACE(sigma_[C][Au]/Ao);
      TRACE(epsilon_[C][Au]/eV);
      TRACE(sigma_[H][Cu]/Ao);
      TRACE(epsilon_[H][Cu]/eV);
      TRACE(sigma_[H][Au]/Ao);
      TRACE(epsilon_[H][Au]/eV);
      exit(1);
    }
  }

  fillR_concat_();
}

FLJ::~FLJ()
{
  for(size_t i = 0; i < ECOUNT; i++)
    for(size_t j = 0; j < ECOUNT; j++)
      if (splines[i][j])
        delete splines[i][j];
}

void
FLJ::fillR_concat_()
{
  Float r;

  AtomsArray atoms(ECOUNT);
  atoms[Cu].ID = Cu_EL;
  atoms[H].ID = C_EL;
  atoms[C].ID = H_EL;
  atoms[Au].ID = Au_EL;
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
    Float epsilon_ij =     this->epsilon(ij);

    Float s_div_r    = sigma_ij/r;
    Float s_div_r_6  = s_div_r;
    int i;
    for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
    Float s_div_r_12 = s_div_r_6*s_div_r_6;
    DerVLJ = 4.0*epsilon_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
    VLJ = 4.0*epsilon_ij*(s_div_r_12-s_div_r_6);
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
  Float epsilon_ij =     this->epsilon(ij);

  Float s_div_r    = sigma_ij/r;

  Float s_div_r_6  = s_div_r;
  int i;
  for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
  Float s_div_r_12 = s_div_r_6*s_div_r_6;

  Float Val = 4.0*epsilon_ij*(s_div_r_12-s_div_r_6);


  Float f = ij.f();

// if (V != 0)
  {
    Float Der = 4.0*epsilon_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
    ij.r(Der*f);
    ij.f(Val);
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
  for(size_t j = 0; j < nl.size(); j++)
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
