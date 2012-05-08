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
  std::vector<ElementID> elements1;
  elements1.push_back(H_EL);
  elements1.push_back(C_EL);

  std::vector<ElementID> elements2;
  elements2.push_back(Cu_EL);
//  elements2.push_back(Ag_EL);
//  elements2.push_back(Au_EL);

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
    Float sigma_ij=     this->sigma(atom1,atom2);
    Float zeta_ij =     this->zeta(atom1,atom2);

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
FLJ::VLJ(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef FLJ_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

#ifndef  LJ_HANDLE_SHORTRANGE
  Spline& spline = *(splines[e2i(atom1)][e2i(atom2)]);
  if (r < spline.x1())
  {
  Float AB_ = 0.53e-8;

  Float ZA = atom1.Z; Float ZB = atom2.Z;

  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));

  Float Y=r/AS;

  if (ontouch_enabled) r_vec_touch_only(atom1,atom2);

  return  ZA*ZB/r*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y));
  }
  else
  {
    if (r < spline.x2())
    {
      if (ontouch_enabled) r_vec_touch_only(atom1,atom2);

      return spline(r);
    }
  }
#endif

  Float sigma_ij=     this->sigma(atom1,atom2);
  Float zeta_ij =     this->zeta(atom1,atom2);

  Float s_div_r    = sigma_ij/r;
  
  Float s_div_r_6  = s_div_r;
  int i;
  for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
  Float s_div_r_12 = s_div_r_6*s_div_r_6;

  return fvar*4.0*zeta_ij*(s_div_r_12-s_div_r_6);
}

Vector3D
FLJ::dVLJ(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef FLJ_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

#ifndef  LJ_HANDLE_SHORTRANGE
  Spline& spline = *(splines[e2i(atom1)][e2i(atom2)]);
  if (r < spline.x1())
  {
  Float AB_ = 0.53e-8;

  Float ZA = atom1.Z; Float ZB = atom2.Z;
  Float AS=8.8534e-1*AB_/(pow(ZA/e,Float(0.23))+pow(ZB/e,Float(0.23)));
  Float Y=r/AS;
  Float Der =
          -ZA*ZB/(r*r)*(0.18175*exp(-3.1998*Y)+
          0.50986*exp(-0.94229*Y)+
          0.28022*exp(-0.4029*Y)+0.02817*exp(-0.20162*Y))-
          ZA*ZB/(r*AS)*(0.18175*3.1998*exp(-3.1998*Y)+
          0.50986*0.94229*exp(-0.94229*Y)+0.28022*0.4029*exp(-0.4029*Y)+
          0.02817*0.20162*exp(-0.20162*Y));

  return Der*drmodvar;
  }
  else
  {
    if (r < spline.x2()) return spline.der(r)*drmodvar;
  }
#endif

  Float sigma_ij=     this->sigma(atom1,atom2);
  Float zeta_ij =     this->zeta(atom1,atom2);

  Float s_div_r    = sigma_ij/r;
  Float s_div_r_6  = s_div_r;
  int i;
  for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
  Float s_div_r_12 = s_div_r_6*s_div_r_6;
  Float Der = 4.0*zeta_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
#ifdef FLJ_OPTIMIZED  
  if (dfvar == 0.0) return Der*f(atom1,atom2)*drmodvar;
#endif
  Float Val = 4.0*zeta_ij*(s_div_r_12-s_div_r_6);
  
  return  Der*f(atom1,atom2)*drmodvar+Val*dfvar;
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
      Ei += VLJ(atom,atom_j);
    }  
  }  
}  
  return Ei;
}

Vector3D
FLJ::grad(Atom &atom,AtomsArray&)
{
  Index j;
  Vector3D dEi(0.0,0.0,0.0);
  for(j = 0; j < NL(atom).size(); j++)
  {
    Atom &atom_j = *(NL(atom)[j]);
    if (isHandledPair(atom,atom_j))
    if (&atom != &atom_j)
    {
      dEi += dVLJ(atom,atom_j,atom);
    }  
  }  
  return dEi;
}


} // namespace mdtk


