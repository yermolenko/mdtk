/*
   The Lennard-Jones interatomic potential implementation.

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
  handledElementPairs.insert(std::make_pair(Cu_EL,C_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,H_EL));
  handledElementPairs.insert(std::make_pair(C_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(H_EL,Cu_EL));

  sigma_[Cu][Cu] = 0.0*Ao;
  sigma_[C][C]   = 0.0*Ao;
  sigma_[Cu][C]  = 2.049067052*Ao;
    sigma_[C][Cu] = sigma_[Cu][C];
  sigma_[Cu][Cu] = 0.0*Ao;
  sigma_[H][H]   = 0.0*Ao;
  sigma_[Cu][H]  = 2.049067052*Ao;
    sigma_[H][Cu] = sigma_[Cu][H];
  
  zeta_[Cu][Cu] = 0.0*eV;
  zeta_[C][C]   = 0.0*eV;
  zeta_[Cu][C]  = 0.05*eV;
    zeta_[C][Cu] = zeta_[Cu][C];
  zeta_[Cu][Cu] = 0.0*eV;
  zeta_[H][H]   = 0.0*eV;
  zeta_[Cu][H]  = 0.01*eV;
    zeta_[H][Cu] = zeta_[Cu][H];
}

Float
FLJ::VLJ(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef FLJ_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

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
FLJ::operator()(AtomsContainer& gl)
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
      Ei += VLJ(atom,atom_j);
    }  
  }  
}  
  return Ei;
}

Vector3D
FLJ::grad(Atom &atom,AtomsContainer&)
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


