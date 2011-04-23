/*
   Implementation of the many-body interatomic potential for copper.
   See [G. Betz, W. Husinsky, Nucl. Instr. and Meth. B 102, 281 (1995)]

   Copyright (C) 2006, 2007, 2008, 2009 Oleksandr Yermolenko
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

#include "TightBinding.hpp"
#include <algorithm>

namespace mdtk
{

  Float  TightBinding::buildPairs(AtomsContainer& gl)
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
      {
      std::pair<int,int> sample_pair(atom_i.globalIndex,DUMMY_EL);
    currentPairPtr = &sample_pair;
    ontouch_enabled = true;
        Ei += F(atom_i);
/*
TRACE(atom_i.globalIndex);
TRACE(F(atom_i)/eV);
*/
    currentPairPtr = NULL;
    ontouch_enabled = false;
      };

      if (isHandled(atom_i))
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        std::pair<int,int> sample_pair(atom_i.globalIndex,atom_j.globalIndex);
        if (isHandled(atom_j))
        if (&atom_i != &atom_j)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
        {

    currentPairPtr = &sample_pair;
    ontouch_enabled = true;
        Ei += Phi(atom_i,atom_j);
/*
TRACE(atom_j.globalIndex);
TRACE(Phi(atom_i,atom_j)/eV);
*/
    currentPairPtr = NULL;
    ontouch_enabled = false;
        }  
      }  
    }  

return Ei;
  }  

Vector3D
TightBinding::grad(Atom &atom,AtomsContainer&gl)
{
  Index i;
  
  Vector3D dEi(0.0,0.0,0.0);

  if (isHandled(atom))
  {

    std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];

    TRACE(atom.globalIndex);
    TRACE(acnt.size());

    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = *(gl[acnt[i].first]);
      
      if (isHandled(atom_i))
      {
      if (acnt[i].second == DUMMY_EL)
      {
        dEi += dF(atom_i,atom);
/*
TRACE(atom.globalIndex);
TRACE(atom_i.globalIndex);
TRACE(dF(atom_i,atom));
*/
        continue;
      }  
      };
      
      if (isHandled(atom_i))
      {
        Atom &atom_j = *(gl[acnt[i].second]);

        REQUIREM(&atom_j != &atom_i,"must be (&atom_j != &atom_i)");
        if (isHandled(atom_j))
        {
          dEi += dPhi(atom_i,atom_j,atom);
/*
TRACE(atom_j.globalIndex);
TRACE(dPhi(atom_i,atom_j,atom));
*/
        }  
      }
    }    
  }  

  return  dEi;
}  


inline
Float
TightBinding::Phi(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef TightBinding_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r;
  if (ontouch_enabled)
    r  = r_vec_module_no_touch(atom1,atom2);
  else
    r  = r_vec_module(atom1,atom2);

#ifndef  EAM_HANDLE_SHORTRANGE
  if (r < 1.029*Ao) return 0.0;
#endif

  if (ontouch_enabled) r_vec_touch_only(atom1,atom2);

  return fvar*Phi0_*exp(-alpha_*r);
}  

inline
Vector3D
TightBinding::dPhi(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef TightBinding_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

#ifndef  EAM_HANDLE_SHORTRANGE
  if (r < 1.029*Ao) return 0.0;
#endif

  return Phi0_*exp(-alpha_*r)*(dfvar-alpha_*drmodvar*f(atom1,atom2));
}

inline
Float
TightBinding::g(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef TightBinding_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r;
  if (ontouch_enabled)
    r  = r_vec_module_no_touch(atom1,atom2);
  else
    r  = r_vec_module(atom1,atom2);

#ifndef  EAM_HANDLE_SHORTRANGE
//  if (r < 1.029*Ao) return 0.0;
#endif

  if (ontouch_enabled) r_vec_touch_only(atom1,atom2);

  return fvar*exp(-beta_*r);
}  

inline
Vector3D
TightBinding::dg(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef TightBinding_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

#ifndef  EAM_HANDLE_SHORTRANGE
//  if (r < 1.029*Ao) return 0.0;
#endif

  return exp(-beta_*r)*(dfvar-beta_*drmodvar*f(atom1,atom2));
}

  
inline
Float
TightBinding::F(Atom &atom1)
{
  Float rhovar = rho(atom1);
  REQUIRE(rhovar >= 0.0);
  return -c_*sqrt(rhovar);
}

inline
Vector3D
TightBinding::dF(Atom &atom1, Atom &datom)
{
  Float rhovar = rho(atom1);
  REQUIRE(rhovar >= 0.0);
  Vector3D drhovar = drho(atom1,datom);
  if (rhovar == 0.0 || drhovar == 0.0) return 0.0;
  return -c_/(2.0*sqrt(rhovar))*drhovar;
}

Float
TightBinding::rho(Atom &atom_i)
{
  Index j;
  Float rhoij = 0.0;
  for(j = 0; j < NL(atom_i).size(); j++)
  {
    Atom& atom_j = *(NL(atom_i)[j]);
    if (/*atom_j.globalIndex > atom_i.globalIndex &&*/ isHandled(atom_j))
    {
      rhoij += g(atom_i,atom_j);
    }  
  }  
  return rhoij;
}

Vector3D
TightBinding::drho(Atom &atom_i, Atom &datom)
{
  Vector3D Derrho = 0.0;
  for(Index j = 0; j < NL(atom_i).size(); j++)
  {
    Atom& atom_j = *(NL(atom_i)[j]);
    if (/*atom_j.globalIndex > atom_i.globalIndex &&*/ isHandled(atom_j))
    {
      Derrho += dg(atom_i,atom_j,datom);
    }  
  }  
  return Derrho;
}

Float
TightBinding::operator()(AtomsContainer& gl)
{
  return buildPairs(gl);
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
  PTRACE("Setup TightBinding");

  alpha_ = 42.87/(10.0*Ao);
  beta_  = 18.00/(10.0*Ao);
  c_     = 12.17*eV;
  Phi0_  = 9.892 *1000.0 *eV;
  R_[0]  = 5.0*Ao;
  R_[1]  = 5.5*Ao;
}  

}


