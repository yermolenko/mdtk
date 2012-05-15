/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The Lennard-Jones part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011, 2012 Oleksandr
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

#include "AIREBO_LJ.hpp"
#include <algorithm>

#include <fstream>

namespace mdtk
{

Float
AIREBO::operator()(AtomsArray& gl)
{
  std::vector<std::vector<AtomPair> > backup = rebo.pairs;
  if (gl.size() != rebo.pairs.size()) rebo.pairs.resize(gl.size());
  size_t iir;
  for(iir = 0; iir < gl.size(); iir++)
  {
    size_t prevSize = rebo.pairs[iir].size();
    rebo.pairs[iir].clear();
    rebo.pairs[iir].reserve(prevSize+FMANYBODY_PAIRS_RESERVE_ADD);
  }

  Float Ei = 0;
  if (gl.size() != pairs.size()) pairs.resize(gl.size());
  size_t ii;
  for(ii = 0; ii < gl.size(); ii++)
  {
    size_t prevSize = pairs[ii].size();
    pairs[ii].clear();
    pairs[ii].reserve(prevSize+FMANYBODY_PAIRS_RESERVE_ADD);
  }
  for(ii = 0; ii < gl.size(); ii++)
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
{
  setupPotential();

  nl.Rcutoff = getRcutoff();
}

void
AIREBO::setupPotential()
{
  PTRACE("Setup AIREBO");

  handledElements.insert(H_EL);
  handledElements.insert(C_EL);
  handledElementPairs.insert(std::make_pair(H_EL,C_EL));
  handledElementPairs.insert(std::make_pair(C_EL,H_EL));

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
}

inline
Float
AIREBO::Cij(AtomsPair& ij, const Float V)
{
  std::vector<std::pair<Atom*,Vector3D> > grads_bak;

#define RESTORE_GRADIENTS                                       \
  {                                                             \
    for(size_t x = 0; x < grads_bak.size(); ++x)                \
      grads_bak[x].first->grad = grads_bak[x].second;           \
    grads_bak.clear();                                          \
  }

#define BACKUP_GRADIENT(ATOM)                                           \
  {                                                                     \
    grads_bak.push_back(std::pair<Atom*,Vector3D>(&ATOM,ATOM.grad));    \
  }

  bool found = false;
  Float wmax = 0.0;
  {
    AtomsPair ij_rebo(ij.atom1,ij.atom2,rebo.R(0,ij),rebo.R(1,ij));
    Float new_wmax = ij_rebo.f();
    if (new_wmax > wmax)
    {
      wmax  = new_wmax;
      if (V != 0.0)
      {
        RESTORE_GRADIENTS;
        BACKUP_GRADIENT(ij_rebo.atom1);
        BACKUP_GRADIENT(ij_rebo.atom2);
        REQUIRE(grads_bak.size()==2);

        ij_rebo.f(-V);
      }
      if (wmax >= 1.0)
        found = true;
    }
    AtomRefsContainer& nli = rebo.NL(ij_rebo.atom1);
    for(size_t k = 0; k < nli.size() && !found; k++)
    {
      Atom &atom_k = *(nli[k]);
      if (!rebo.probablyAreNeighbours(ij_rebo.atom1,atom_k)) continue;
      if (&atom_k != &ij_rebo.atom1 && &atom_k != &ij_rebo.atom2)
      {
        AtomsPair ik_rebo(ij_rebo.atom1,atom_k,
                          rebo.R(0,ij_rebo.atom1,atom_k),
                          rebo.R(1,ij_rebo.atom1,atom_k));
        if (rebo.probablyAreNeighbours(ij_rebo.atom2,atom_k))
        {
          AtomsPair kj_rebo(atom_k,ij_rebo.atom2,
                            rebo.R(0,ij_rebo.atom2,atom_k),
                            rebo.R(1,ij_rebo.atom2,atom_k));
          Float new_wmax = ik_rebo.f()*kj_rebo.f();
          if (new_wmax > wmax)
          {
            wmax  = new_wmax;
            if (V != 0.0)
            {
              RESTORE_GRADIENTS;
              BACKUP_GRADIENT(ik_rebo.atom1);
              BACKUP_GRADIENT(ik_rebo.atom2);
              BACKUP_GRADIENT(kj_rebo.atom2);
              REQUIRE(grads_bak.size()==3);

              ik_rebo.f(-kj_rebo.f()*V);
              kj_rebo.f(-ik_rebo.f()*V);
            }
            if (wmax >= 1.0)
            {
              found = true;
              break;
            }
          }
        }
        AtomRefsContainer& nlj = rebo.NL(ij_rebo.atom2);
        for(Index l = 0; l < nlj.size() && !found; l++)
        {
          Atom &atom_l = *(nlj[l]);
          if (!rebo.probablyAreNeighbours(ij_rebo.atom2,atom_l)) continue;
          if (&atom_l != &atom_k && &atom_l != &ij_rebo.atom1 && &atom_l != &ij_rebo.atom2)
          {
            if (rebo.probablyAreNeighbours(atom_k,atom_l))
            {
              AtomsPair lj_rebo(atom_l,ij_rebo.atom2,
                                rebo.R(0,ij_rebo.atom2,atom_l),
                                rebo.R(1,ij_rebo.atom2,atom_l));
              AtomsPair kl_rebo(atom_k,atom_l,
                                rebo.R(0,atom_l,atom_k),
                                rebo.R(1,atom_l,atom_k));
              Float new_wmax = ik_rebo.f()*kl_rebo.f()*lj_rebo.f();
              if (new_wmax > wmax)
              {
                wmax  = new_wmax;
                if (V != 0.0)
                {
                  RESTORE_GRADIENTS;
                  BACKUP_GRADIENT(ik_rebo.atom1);
                  BACKUP_GRADIENT(ik_rebo.atom2);
                  BACKUP_GRADIENT(kl_rebo.atom2);
                  BACKUP_GRADIENT(lj_rebo.atom2);
                  REQUIRE(grads_bak.size()==4);

                  ik_rebo.f(-kl_rebo.f()*lj_rebo.f()*V);
                  kl_rebo.f(-ik_rebo.f()*lj_rebo.f()*V);
                  lj_rebo.f(-ik_rebo.f()*kl_rebo.f()*V);
                }
                if (wmax >= 1.0)
                {
                  found = true;
                  break;
                }
              }
            }
          }
        }
      }
    }
  }
  REQUIRE(1.0-wmax >= 0);

  return 1.0-wmax;
}

}
