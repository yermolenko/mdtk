/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The Lennard-Jones part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr
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

  Float  AIREBO::buildPairs(AtomsContainer& gl)
  {
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
      Atom &atom_i = *(gl[ii]);
      if (isHandled(atom_i))
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        std::pair<int,int> sample_pair(atom_i.globalIndex,atom_j.globalIndex);
        if (isHandled(atom_j))
        if (&atom_i != &atom_j)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
#ifdef AIREBO_OPTIMIZED_EVEN_BETTER
        if (Cij(atom_i, atom_j) > 0.0 /* C!=0 */)
#endif
        {
    currentPairPtr = &sample_pair;
    ontouch_enabled = true;

      {
        Float C = Cij(atom_i, atom_j);
#ifdef AIREBO_OPTIMIZED
          if (C > 0.0) // C !=0.0
          {
#endif
        Float rij = r_vec_module(atom_i, atom_j);
        Float rij_min = RLJ(0,atom_i, atom_j);
        Float rij_max = RLJ(1,atom_i, atom_j);
        Float Str = S(rij,rij_min,rij_max);
#ifdef AIREBO_OPTIMIZED
          if (Str != 0.0)
          {
#endif
        Float bij = BijAsterix(atom_i, atom_j);
        Float bij_min = b(0,atom_i, atom_j);
        Float bij_max = b(1,atom_i, atom_j);
        Float Stb = S(bij,bij_min,bij_max);
        Ei += (1.0+Str*(Stb-1.0))*C*VLJ(atom_i,atom_j);
#ifdef AIREBO_OPTIMIZED
          } // Str!=0
          else
          {
        Ei += (1.0              )*C*VLJ(atom_i,atom_j);
          }  
#endif 
#ifdef AIREBO_OPTIMIZED
          } // C!=0
#endif 

      }  
    currentPairPtr = NULL;
    ontouch_enabled = false;
        }  
      }  
    }  
    return Ei;
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
AIREBO::VLJ(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef AIREBO_OPTIMIZED  
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

inline
Vector3D
AIREBO::dVLJ(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef AIREBO_OPTIMIZED  
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
#ifdef AIREBO_OPTIMIZED  
  if (dfvar == 0.0) return Der*f(atom1,atom2)*drmodvar;
#endif
  Float Val = 4.0*zeta_ij*(s_div_r_12-s_div_r_6);
  
  return  Der*f(atom1,atom2)*drmodvar+Val*dfvar;
}  

Float
AIREBO::operator()(AtomsContainer& gl)
{
  Float Ei = 0.0;
  Ei += ELJ(gl);
  return Ei;
}  



Vector3D
AIREBO::grad(Atom &atom,AtomsContainer &gl)
{
  Vector3D dEi = 0.0;
  dEi += dELJ(atom, gl);
  return  dEi;
}  

Float
AIREBO::BijAsterix(Atom &atom1,Atom &atom2)
{
  Float bij = 0.0;
  
  set_r_vec_exception(atom1,atom2,CREBO::R(0,atom1,atom2));
  bij = Baver(atom1, atom2);
  cleat_r_vec_exception();
  
  return bij;
}

Vector3D
AIREBO::dBijAsterix(Atom &/*atom1*/,Atom &/*atom2*/, Atom &/*datom*/)
{
  return 0.0; // TODO
}   

Float
AIREBO::ELJ(AtomsContainer& gl)
{
  return buildPairs(gl);
}  

Vector3D
AIREBO::dELJ(Atom &atom,AtomsContainer &gl)
{
  Vector3D dEi = 0.0;

  Index i;

  if (isHandled(atom))
  {
    std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];
    
    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = *(gl[acnt[i].first]);
      if (isHandled(atom_i))
      {
        Atom &atom_j = *(gl[acnt[i].second]);

        REQUIREM(&atom_j != &atom_i,"must be (&atom_j != &atom_i)");
        if (isHandled(atom_j))
        {
          Vector3D dEij = 0.0; dEij = 0.0;

          Float C = Cij(atom_i, atom_j);
#ifdef AIREBO_OPTIMIZED
          if (C > 0.0) // (C != 0.0)
          {
#endif
          Vector3D dC = 0.0;
#ifdef AIREBO_OPTIMIZED
            if (C < 1.0)
            {
#endif
          dC = dCij(atom_i, atom_j, atom);
#ifdef AIREBO_OPTIMIZED
            }   
            else
            {
          dC = 0.0;
            };
#endif
          Float rij = r_vec_module(atom_i, atom_j);
          Float rij_min = RLJ(0,atom_i, atom_j);
          Float rij_max = RLJ(1,atom_i, atom_j);
          Float Str = S(rij,rij_min,rij_max);
          Vector3D dStr = dS(rij,rij_min,rij_max)*dr_vec_module(atom_i, atom_j,atom);
#ifdef AIREBO_OPTIMIZED
          if (Str != 0.0)
          {
#endif
          Float bij = BijAsterix(atom_i, atom_j);
          Float bij_min = b(0,atom_i, atom_j);
          Float bij_max = b(1,atom_i, atom_j);
          Float Stb = S(bij,bij_min,bij_max);
#ifdef AIREBO_OPTIMIZED
          Float dSdBij  = dS(bij,bij_min,bij_max);
          Vector3D dStb = 0.0;
          if (dSdBij != 0.0)
          {
            dStb = dSdBij*dBijAsterix(atom_i, atom_j,atom);
          }
#else
          Vector3D dStb = dS(bij,bij_min,bij_max)*dBijAsterix(atom_i, atom_j,atom);
#endif

          Float VLJ_val = VLJ(atom_i,atom_j);

          dEij = (dStr*(Stb-1.0)+Str*dStb)  * C* VLJ_val+
                 (1.0+Str*(Stb-1.0))        *dC* VLJ_val+
                 (1.0+Str*(Stb-1.0))        * C*dVLJ(atom_i,atom_j,atom);

#ifdef AIREBO_OPTIMIZED
          } // Str != 0.0
          else
          {
          dEij =  0.0+
                 (1.0              )        *dC* VLJ(atom_i,atom_j)+
                 (1.0              )        * C*dVLJ(atom_i,atom_j,atom);
          }  
#endif

#ifdef AIREBO_OPTIMIZED
          } // C!=0
#endif
          dEi += dEij;
        }  
      }
    }    
    
  }  

  return  dEi;
}  



AIREBO::AIREBO(CREBO* crebo):
#ifdef AIREBO_USING_BRENNER
  CREBO(CREBO::POTENTIAL2)
#endif
#ifdef AIREBO_USING_REBO
  CREBO(CREBO::POTENTIAL1)
#endif
  ,rebo(*crebo)
  ,nl(this)
{
  setupPotential();

  nl.Rcutoff = getRcutoff();
}

void
AIREBO::setupPotential()
{
  PTRACE("Setup AIREBO");

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
AIREBO::Cij(Atom &atom1,Atom &atom2)
{
  Float wmax = 0.0;
  Float new_wmax = rebo.f(atom1,atom2);
  if (new_wmax > wmax)
  {
    wmax  = new_wmax;
if (wmax >= 1.0)  return 0.0/*1.0-wmax*/;
  }
  for(Index k = 0; k < rebo.NL(atom1).size(); k++)
  {
    Atom &atom_k = *(rebo.NL(atom1)[k]);
    if (isHandled(atom_k))
    if (&atom_k != &atom1 && &atom_k != &atom2)
    {
      {
        Float new_wmax = rebo.f(atom1,atom_k)*rebo.f(atom_k,atom2);
        if (new_wmax > wmax)
        {
          wmax  = new_wmax;
if (wmax >= 1.0)  return 0.0/*1.0-wmax*/;
        }  
      }  
      for(Index l = 0; l < rebo.NL(atom2).size(); l++)
      {
        Atom &atom_l = *(rebo.NL(atom2)[l]); 
        if (isHandled(atom_l))
        if (&atom_l != &atom_k && &atom_l != &atom1 && &atom_l != &atom2)
        {
          {
            Float new_wmax = rebo.f(atom1,atom_k)*rebo.f(atom_k,atom_l)*rebo.f(atom_l,atom2);
            if (new_wmax > wmax)
            {
               wmax  = new_wmax;
if (wmax >= 1.0)  return 0.0/*1.0-wmax*/;
            }  
          }  
        }  
      }
    }
  }    
  return 1.0-wmax;
}
   
inline
Vector3D
AIREBO::dCij(Atom &atom1,Atom &atom2, Atom &datom)
{
  Float wmax = 0.0;
  Vector3D dwmax = 0.0;
  Float new_wmax = rebo.f(atom1,atom2);
  if (new_wmax > wmax)
  {
    dwmax  = rebo.df(atom1,atom2,datom);
     wmax  = new_wmax;
if (wmax >= 1.0)  {REQUIRE(dwmax==0);return -dwmax;}
  }  
  for(Index k = 0; k < rebo.NL(atom1).size(); k++)
  {
    Atom &atom_k = *(rebo.NL(atom1)[k]);
    if (isHandled(atom_k))
    if (&atom_k != &atom1 && &atom_k != &atom2)
    {
      {
        Float new_wmax = rebo.f(atom1,atom_k)*rebo.f(atom_k,atom2);
        if (new_wmax > wmax)
        {
          dwmax  = rebo.df(atom1,atom_k,datom)*rebo. f(atom_k,atom2)+
                   rebo. f(atom1,atom_k)      *rebo.df(atom_k,atom2,datom);
           wmax  = new_wmax;
if (wmax >= 1.0)  {REQUIRE(dwmax==0);return -dwmax;}
        }
      }
      for(Index l = 0; l < rebo.NL(atom2).size(); l++)
      {
        Atom &atom_l = *(rebo.NL(atom2)[l]); 
        if (isHandled(atom_l))
        if (&atom_l != &atom_k && &atom_l != &atom1 && &atom_l != &atom2)
        {
          {
            Float new_wmax = rebo.f(atom1,atom_k)*rebo.f(atom_k,atom_l)*rebo.f(atom_l,atom2);
            if (new_wmax > wmax)
            {
              dwmax  = rebo.df(atom1,atom_k,datom)*rebo. f(atom_k,atom_l)      *rebo. f(atom_l,atom2)+
                       rebo. f(atom1,atom_k)      *rebo.df(atom_k,atom_l,datom)*rebo. f(atom_l,atom2)+
                       rebo. f(atom1,atom_k)      *rebo. f(atom_k,atom_l)      *rebo.df(atom_l,atom2,datom);
               wmax  = new_wmax;
if (wmax >= 1.0)  {REQUIRE(dwmax==0);return -dwmax;}
            }  
          }  
        }  
      }
    }
  }    
  return -dwmax;
}   

}


