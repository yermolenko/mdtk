/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The Lennard-Jones part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013, 2015
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

#include "AIREBO_LJ.hpp"
#include <algorithm>

#include <fstream>

namespace mdtk
{

Float
AIREBO::operator()(AtomsArray& gl)
{
  cleanup_Cij();
  fill_Cij(gl);

  Float Ei = 0;
  for(size_t ii = 0; ii < gl.size(); ii++)
  {
    Atom &atom_i = gl[ii];
    if (isHandled(atom_i))
    {
      AtomRefsContainer& nli = NL(atom_i);
      for(size_t jj = 0; jj < nli.size(); jj++)
      {
        Atom &atom_j = *(nli[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        if (&atom_i != &atom_j)
        {
          if (!probablyAreNeighbours(atom_i,atom_j)) continue;
          AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));

          Float V = 1.0;

          Float C = Cij(ij);

          if (C > 0.0) // C !=0.0
          {
            Float SrSb = StrStb(ij);
            Float LJ = VLJ(ij,SrSb*C*V);
            if (V != 0)
            {
              StrStb(ij,C*LJ*V);
              Cij(ij,SrSb*LJ*V);
            }
            Ei += SrSb*C*LJ;
          }
        }
      }
    }
  }

  return Ei;
}

inline
Float
AIREBO::StrStb(AtomsPair& ij, const Float V)
{
  Float Val = 1.0;

  Float rij = ij.r();
  Float rij_min = RLJ(0,ij);
  Float rij_max = RLJ(1,ij);
  Float Str = S(rij,rij_min,rij_max);
  if (Str != 0.0)
  {
    Float bij = BijAsterix(ij);
    Float bij_min = b(0,ij);
    Float bij_max = b(1,ij);
    Float Stb = S(bij,bij_min,bij_max);

    if (V != 0.0)
    {
      Float dStr = dS(rij,rij_min,rij_max);
      if (dStr != 0.0)
        ij.r(-dStr*V);

      Float dStb = dS(bij,bij_min,bij_max);
      if (dStr != 0.0 && Stb != 0.0)
        ij.r(dStr*Stb*V);
      if (dStb != 0.0)
        BijAsterix(ij,Str*dStb*V);
    }

    Val += Str*Stb-Str;
  }

//  Val = 1.0+Str*(Stb-1.0);
  return Val;
}

inline
Float
AIREBO::S(Float arg, Float arg_min, Float arg_max) const
{
  REQUIRE(arg_min < arg_max);

  if (arg<arg_min)
  {
    return 1.0;
  }
  else if (arg>arg_max)
  {
    return 0.0;
  }
  else
  {
    Float t = (arg-arg_min)/(arg_max-arg_min);
    return 1.0-t*t*(3.0-2.0*t);
  }
}

inline
Float
AIREBO::dS(Float arg, Float arg_min, Float arg_max) const
{
  REQUIRE(arg_min < arg_max);

  if (arg<arg_min)
  {
    return 0.0;
  }
  else if (arg>arg_max)
  {
    return 0.0;
  }
  else
  {
    Float t = (arg-arg_min)/(arg_max-arg_min);
    return (-2.0*t*(3.0-2.0*t)+2.0*t*t)/(arg_max-arg_min);
  }
}

inline
Float
AIREBO::VLJ(AtomsPair& ij, const Float V)
{
  Float fvar = ij.f();

#ifdef AIREBO_OPTIMIZED
  if (fvar == 0.0) return 0.0;
#endif

  Float r = ij.r();

  Float sigma_ij=     this->sigma(ij);
  Float zeta_ij =     this->zeta(ij);

  Float s_div_r    = sigma_ij/r;
  Float s_div_r_6  = s_div_r;
  int i;
  for(i = 2; i <= 6; i++) s_div_r_6 *= s_div_r;
  Float s_div_r_12 = s_div_r_6*s_div_r_6;

  Float Val = 4.0*zeta_ij*(s_div_r_12-s_div_r_6);

  if (V != 0.0)
  {
    Float Der = 4.0*zeta_ij*(-12.0*s_div_r_12/r+6.0*s_div_r_6/r);
    ij.r(Der*fvar*V);
    ij.f(Val*V);
  }

  return fvar*Val;
}

Float
AIREBO::BijAsterix(AtomsPair& ij, const Float V)
{
  Float bij = 0.0;

  AtomsPair ijAsterix(ij.atom1,ij.atom2,
                      rebo.R(0,ij.atom1,ij.atom2),rebo.R(1,ij.atom1,ij.atom2),
                      rebo.R(0,ij.atom1,ij.atom2));
  bij = rebo.Baver(ijAsterix, V);

  return bij;
}

AIREBO::AIREBO(CREBO* crebo):
  FManybody()
  ,rebo(*crebo)
  ,CA()
{
  setupPotential();

  nl.Rcutoff = getRcutoff();
}

void
AIREBO::setupPotential()
{
  handledElements.insert(H_EL);
  handledElements.insert(C_EL);
  handledElementPairs.insert(std::make_pair(H_EL,C_EL));
  handledElementPairs.insert(std::make_pair(C_EL,H_EL));
  handledElementPairs.insert(std::make_pair(H_EL,H_EL));
  handledElementPairs.insert(std::make_pair(C_EL,C_EL));

  sigma_[C][C] = 3.40*Ao;
  sigma_[H][H] = 2.65*Ao;
  sigma_[C][H] = 0.5*(sigma_[C][C]+sigma_[H][H]);
    sigma_[H][C] = sigma_[C][H];
  
  zeta_[C][C] = 0.00284*eV;
  zeta_[H][H] = 0.00150*eV;
  zeta_[C][H] = sqrt(zeta_[C][C]*zeta_[H][H]);
    zeta_[H][C] = zeta_[C][H];
  
  RLJ_[C][C][0] = sigma_[C][C];
  RLJ_[H][H][0] = sigma_[H][H];
  RLJ_[C][H][0] = sigma_[C][H];
    RLJ_[H][C][0] = RLJ_[C][H][0];

  RLJ_[C][C][1] = sigma_[C][C]*pow(2.0,1.0/6.0);
  RLJ_[H][H][1] = sigma_[H][H]*pow(2.0,1.0/6.0);
  RLJ_[C][H][1] = sigma_[C][H]*pow(2.0,1.0/6.0);
    RLJ_[H][C][1] = RLJ_[C][H][1];

  R_[C][C][0] = /*7.0*Ao;*/ 5.0*Ao;
  R_[H][H][0] = /*6.0*Ao;*/ 4.0*Ao;
  R_[C][H][0] = /*6.5*Ao;*/ 4.5*Ao;
    R_[H][C][0] = R_[C][H][0];

  REQUIRE(R_[C][C][0] > RLJ_[C][C][1]);
  REQUIRE(R_[H][H][0] > RLJ_[H][H][1]);
  REQUIRE(R_[C][H][0] > RLJ_[C][H][1]);
  REQUIRE(R_[H][C][0] > RLJ_[H][C][1]);

  R_[C][C][1] = R_[C][C][0]+0.5*Ao;
  R_[H][H][1] = R_[H][H][0]+0.5*Ao;
  R_[C][H][1] = R_[C][H][0]+0.5*Ao;
    R_[H][C][1] = R_[C][H][1];

  b_[C][C][0] = 0.77;
  b_[H][H][0] = 0.32;
  b_[C][H][0] = 0.75;
    b_[H][C][0] = b_[C][H][0];

  b_[C][C][1] = 0.81;
  b_[H][H][1] = 0.42;
  b_[C][H][1] = 0.90;
    b_[H][C][1] = b_[C][H][1];

#ifdef AIREBO_USING_BRENNER
  b_[C][C][0] = 0.77;
  b_[H][H][0] = 0.32;
  b_[C][H][0] = 0.75;//0.80
    b_[H][C][0] = b_[C][H][0];

  b_[C][C][1] = 0.89;//modified
  b_[H][H][1] = 0.42;//0.56(0.60)
  b_[C][H][1] = 0.90;//0.85
    b_[H][C][1] = b_[C][H][1];
#endif

  PRINT("AIREBO::LJ interatomic potential configured.\n");
}

void
AIREBO::fill_Cij(AtomsArray& gl)
{
  CA.resize(gl.size());

  for(size_t i1 = 0; i1 < gl.size(); i1++)
  {
    Atom &atom1 = gl[i1];
    if (isHandled(atom1))
    {
      AtomRefsContainer& nl1 = rebo.NL(atom1);
      for(size_t i2 = 0; i2 < nl1.size(); i2++)
      {
        Atom &atom2 = *(nl1[i2]);
        if (!rebo.probablyAreNeighbours(atom1,atom2)) continue;
        CEelement& pa
          = CA[atom1.globalIndex][atom2.globalIndex];
        Float w = pa.first;
        AtomsPair p1(atom1,atom2,
                     rebo.R(0,atom1,atom2),
                     rebo.R(1,atom1,atom2));
        Float new_w = p1.f();
        if (new_w > w || (w > 0.0 && new_w == w && pa.second.size() > 1))
        {
          for(size_t pj = 0; pj < pa.second.size(); ++pj)
            delete pa.second[pj];
          pa.second.clear();
          pa.first = new_w;
          pa.second.push_back(new AtomsPair(p1));
        }
        AtomRefsContainer& nl2 = rebo.NL(atom2);
        for(size_t i3 = 0; i3 < nl2.size(); i3++)
        {
          Atom &atom3 = *(nl2[i3]);
          if (!rebo.probablyAreNeighbours(atom2,atom3)) continue;
          CEelement& pa
            = CA[atom1.globalIndex][atom3.globalIndex];
          Float w = pa.first;
          AtomsPair p2(atom2,atom3,
                       rebo.R(0,atom2,atom3),
                       rebo.R(1,atom2,atom3));
          Float new_w = p1.f()*p2.f();
          if (new_w > w || (w > 0.0 && new_w == w && pa.second.size() > 2))
          {
            for(size_t pj = 0; pj < pa.second.size(); ++pj)
              delete pa.second[pj];
            pa.second.clear();
            pa.first = new_w;
            pa.second.push_back(new AtomsPair(p1));
            pa.second.push_back(new AtomsPair(p2));
          }
          AtomRefsContainer& nl3 = rebo.NL(atom3);
          for(size_t i4 = 0; i4 < nl3.size(); i4++)
          {
            Atom &atom4 = *(nl3[i4]);
            if (!rebo.probablyAreNeighbours(atom3,atom4)) continue;
            CEelement& pa
              = CA[atom1.globalIndex][atom4.globalIndex];
            Float w = pa.first;
            AtomsPair p3(atom3,atom4,
                         rebo.R(0,atom3,atom4),
                         rebo.R(1,atom3,atom4));
            Float new_w = p1.f()*p2.f()*p3.f();
            if (new_w > w/* || (w > 0.0 && new_w == w && pa.second.size() > 3)*/)
            {
              for(size_t pj = 0; pj < pa.second.size(); ++pj)
                delete pa.second[pj];
              pa.second.clear();
              pa.first = new_w;
              pa.second.push_back(new AtomsPair(p1));
              pa.second.push_back(new AtomsPair(p2));
              pa.second.push_back(new AtomsPair(p3));
            }
          }
        }
      }
    }
  }
}

void
AIREBO::cleanup_Cij()
{
  for(size_t i = 0; i < CA.size(); i++)
  {
    for(CArray::iterator ca = CA[i].begin(); ca != CA[i].end(); ++ca)
      for(size_t pj = 0; pj < ca->second.second.size(); ++pj)
        delete ca->second.second[pj];
    CA[i].clear();
  }
}

inline
Float
AIREBO::Cij(AtomsPair& ij, const Float V)
{
  CArray::iterator cel
    = CA[ij.atom1.globalIndex].find(ij.atom2.globalIndex);

  Float w = 0.0;

  if (cel != CA[ij.atom1.globalIndex].end())
  {
    w = cel->second.first;

    if (V != 0.0)
    {
      for(size_t pi = 0; pi < cel->second.second.size(); ++pi)
      {
        Float m = 1.0;
        for(size_t pj = 0; pj < cel->second.second.size(); ++pj)
          if (pi != pj)
            m *= (cel->second.second)[pj]->f();
        (cel->second.second)[pi]->f(-m*V);
      }
    }
  }

  return 1 - w;
}

}
