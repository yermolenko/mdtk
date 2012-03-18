/*
   Implementation of the many-body interatomic potential for copper,
   gold, silver and their alloys.
   See [G.J. Ackland and V. Vitek, Phys. Rev. B 41, 10324 (1990)]

   Copyright (C) 2007, 2008, 2009, 2012 Oleksandr Yermolenko
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

#include "Ackland.hpp"
#include <algorithm>

namespace mdtk
{

//#define ACKLAND_ZBL_CORRETION (-22.67812*eV)
#define ACKLAND_ZBL_CORRETION (-40.0*eV)

  Float  Ackland::buildPairs(AtomsArray& gl)
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
      Atom &atom_i = gl[ii];
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
//        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
        {

    currentPairPtr = &sample_pair;
    ontouch_enabled = true;
        Ei += Phi(atom_i,atom_j);
/*
TRACE(atom_j.globalIndex);
TRACE(r_vec_module_no_touch(atom_i,atom_j));
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
Ackland::grad(Atom &atom,AtomsArray& gl)
{
  Index i;
  
  Vector3D dEi(0.0,0.0,0.0);

  if (isHandled(atom))
  {

    std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];

    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = gl[acnt[i].first];
      
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
        Atom &atom_j = gl[acnt[i].second];

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
Ackland::Phi(Atom &atom1,Atom &atom2) // V
{
  Float r;
  if (ontouch_enabled)
    r  = r_vec_module_no_touch(atom1,atom2);
  else
    r  = r_vec_module(atom1,atom2);

#ifndef  Ackland_HANDLE_SHORTRANGE
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

#ifdef Ackland_OPTIMIZED  
  if (r >= rk(1,atom1,atom2))  return 0.0;
#endif

  if (ontouch_enabled) r_vec_touch_only(atom1,atom2);

  Float V = 0.0;

  size_t nk = (atom1.ID == atom2.ID)?6:3;

  for(size_t k = 1; k <= nk; k++)
  {
    Float rt = rk(k,atom1,atom2)-r;
    if (rt <= 0) continue;
    V += ak(k,atom1,atom2)*rt*rt*rt;
//TRACE(ak(k,atom1,atom2)*rt*rt*rt/eV);
  }

  return V;
}  

inline
Vector3D
Ackland::dPhi(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef Ackland_OPTIMIZED  
  if (drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

#ifndef  Ackland_HANDLE_SHORTRANGE
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

#ifdef Ackland_OPTIMIZED  
  if (r >= rk(1,atom1,atom2))  return 0.0;
#endif

  Float dV = 0.0;
  Float V = 0.0;

  size_t nk = (atom1.ID == atom2.ID)?6:3;

  for(size_t k = 1; k <= nk; k++)
  {
    Float rt = rk(k,atom1,atom2)-r;
    if (rt <= 0) continue;
    Float akval = ak(k,atom1,atom2);
    V  += akval*rt*rt*rt;
    dV += akval*3.0*rt*rt*(-1.0);
  }

  return dV*drmodvar;
}

inline
Float
Ackland::PhiCap(size_t a1_id, size_t a2_id, Float r) const
{
#ifdef Ackland_OPTIMIZED  
  if (r >= Rk_[a1_id][a2_id][1])  return 0.0;
#endif
/*
  TRACE(a1_id);
  TRACE(a2_id);
  TRACE(r);
*/
  Float val = 0.0;
  for(size_t k = 1; k <= 2; k++)
  {
//    TRACE(Rk_[a1_id][a2_id][k]);
//    TRACE(Ak_[a1_id][a2_id][k]);

    Float rt = Rk_[a1_id][a2_id][k]-r;
    if (rt <= 0) continue;
    val += Ak_[a1_id][a2_id][k]*rt*rt*rt;
  }

//  TRACE(val);

  return val;
}

inline
Float
Ackland::dPhiCap(size_t a1_id, size_t a2_id, Float r) const
{
#ifdef Ackland_OPTIMIZED  
  if (r >= Rk_[a1_id][a2_id][1])  return 0.0;
#endif

  Float val = 0.0;
  for(size_t k = 1; k <= 2; k++)
  {
    Float rt = Rk_[a1_id][a2_id][k]-r;
    if (rt <= 0) continue;
    val += Ak_[a1_id][a2_id][k]*3.0*rt*rt*(-1.0);
  }
  return val;
}

inline
Float
Ackland::g(Atom &atom1,Atom &atom2)  //PhiBig
{
  Float r;
  if (ontouch_enabled)
    r  = r_vec_module_no_touch(atom1,atom2);
  else
    r  = r_vec_module(atom1,atom2);

  Float PhiCapVal = 0.0;

  size_t a1_id = e2i(atom1);
  size_t a2_id = e2i(atom2);

  if (a1_id != a2_id)
  {
    bool touch = true;
    if (r >= Rk_[a1_id][a1_id][1]) touch = false;
    if (r >= Rk_[a2_id][a2_id][1]) touch = false;
    if (ontouch_enabled && touch) r_vec_touch_only(atom1,atom2);
    PhiCapVal = sqrt(PhiCap(a1_id, a1_id, r)*PhiCap(a2_id, a2_id, r));
  }
  else
  {
    bool touch = true;
    if (r >= Rk_[a1_id][a2_id][1]) touch = false;
    if (ontouch_enabled && touch) r_vec_touch_only(atom1,atom2);
    PhiCapVal = PhiCap(a1_id, a2_id, r);
  }

  return PhiCapVal;
}

inline
Vector3D
Ackland::dg(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef Ackland_OPTIMIZED  
  if (drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  size_t a1_id = e2i(atom1);
  size_t a2_id = e2i(atom2);

  if (a1_id != a2_id)
  {
    Float PhiCapVal1 = PhiCap(a1_id, a1_id, r);
    Float PhiCapVal2 = PhiCap(a2_id, a2_id, r);
//    TRACE(PhiCapVal1);
//    TRACE(PhiCapVal2);
    REQUIRE(PhiCapVal1*PhiCapVal2>=0);
    if (PhiCapVal1*PhiCapVal2<=0) return 0.0;
    Float dPhiCapVal1 = dPhiCap(a1_id, a1_id, r);
    Float dPhiCapVal2 = dPhiCap(a2_id, a2_id, r);
    Float dsqrtPhiCapVal = 0.5/sqrt(PhiCapVal1*PhiCapVal2)*(dPhiCapVal1*PhiCapVal2+PhiCapVal1*dPhiCapVal2);
    return dsqrtPhiCapVal*drmodvar;
  }
  else
  {
    Float dPhiCapVal = dPhiCap(a1_id, a2_id, r);
    return dPhiCapVal*drmodvar;
  }
}
  
inline
Float
Ackland::F(Atom &atom1)
{
//  return 0;

  Float rhovar = rho(atom1);
  REQUIRE(rhovar >= 0.0);
  return -c_*sqrt(rhovar);
}

inline
Vector3D
Ackland::dF(Atom &atom1, Atom &datom)
{
//  return 0;

  Float rhovar = rho(atom1);
  REQUIRE(rhovar >= 0.0);
  Vector3D drhovar = drho(atom1,datom);
  if (rhovar == 0.0 || drhovar == 0.0) return 0.0;
  return -c_/(2.0*sqrt(rhovar))*drhovar;
}

Float
Ackland::rho(Atom &atom_i)
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
Ackland::drho(Atom &atom_i, Atom &datom)
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
Ackland::operator()(AtomsArray& gl)
{
  return buildPairs(gl);
}  

Ackland::Ackland():
  FManybody()
{
  handledElements.insert(Cu_EL);
  handledElementPairs.insert(std::make_pair(Cu_EL,Cu_EL));
  handledElementPairs.insert(std::make_pair(Cu_EL,DUMMY_EL));
  handledElementPairs.insert(std::make_pair(DUMMY_EL,Cu_EL));

  handledElements.insert(Ag_EL);
  handledElementPairs.insert(std::make_pair(Ag_EL,Ag_EL));
  handledElementPairs.insert(std::make_pair(Ag_EL,DUMMY_EL));
  handledElementPairs.insert(std::make_pair(DUMMY_EL,Ag_EL));

  handledElements.insert(Au_EL);
  handledElementPairs.insert(std::make_pair(Au_EL,Au_EL));
  handledElementPairs.insert(std::make_pair(Au_EL,DUMMY_EL));
  handledElementPairs.insert(std::make_pair(DUMMY_EL,Au_EL));

  handledElementPairs.insert(std::make_pair(Cu_EL,Ag_EL));
  handledElementPairs.insert(std::make_pair(Ag_EL,Cu_EL));

  handledElementPairs.insert(std::make_pair(Cu_EL,Au_EL));
  handledElementPairs.insert(std::make_pair(Au_EL,Cu_EL));

  handledElementPairs.insert(std::make_pair(Ag_EL,Au_EL));
  handledElementPairs.insert(std::make_pair(Au_EL,Ag_EL));

  setupPotential();

  nl.Rcutoff = getRcutoff();
}

void
Ackland::setupPotential()
{
  PTRACE("Setup Ackland");

//  Float koe = 1.3067977/61.73525861;

/*
  alpha_ = 42.87/(10.0*Ao);
  beta_  = 18.00/(10.0*Ao);
*/
  c_     = +1.0/* *eV */;//12.17*eV;
/*
  Phi0_  = 9.892 *1000.0 *eV;
*/
  Float a;

  rk_[Ag][Au][1] = rk_[Au][Ag][1] = 4.40856*Ao;
  rk_[Ag][Au][2] = rk_[Au][Ag][2] = 3.46970*Ao;
  rk_[Ag][Au][3] = rk_[Au][Ag][3] = 3.00000*Ao;
  ak_[Ag][Au][1] = ak_[Au][Ag][1] = 0.003566503525*eV/(Ao*Ao*Ao);
  ak_[Ag][Au][2] = ak_[Au][Ag][2] = 1.015075326*eV/(Ao*Ao*Ao);
  ak_[Ag][Au][3] = ak_[Au][Ag][3] = 0.00000000*eV/(Ao*Ao*Ao);

  rk_[Cu][Au][1] = rk_[Au][Cu][1] = 4.309816*Ao;
  rk_[Cu][Au][2] = rk_[Au][Cu][2] = 4.047479*Ao;
  rk_[Cu][Au][3] = rk_[Au][Cu][3] = 3.297946*Ao;
  ak_[Cu][Au][1] = ak_[Au][Cu][1] = -0.08554551660*eV/(Ao*Ao*Ao);
  ak_[Cu][Au][2] = ak_[Au][Cu][2] = 0.192835880*eV/(Ao*Ao*Ao);
  ak_[Cu][Au][3] = ak_[Au][Cu][3] = 0.75932286*eV/(Ao*Ao*Ao);

  rk_[Ag][Cu][1] = rk_[Cu][Ag][1] = 4.15800*Ao;
  rk_[Ag][Cu][2] = rk_[Cu][Ag][2] = 3.15700*Ao;
  rk_[Ag][Cu][3] = rk_[Cu][Ag][3] = 3.00000*Ao;
  ak_[Ag][Cu][1] = ak_[Cu][Ag][1] = 0.03533037752*eV/(Ao*Ao*Ao);
  ak_[Ag][Cu][2] = ak_[Cu][Ag][2] = 0.8518466353*eV/(Ao*Ao*Ao);
  ak_[Ag][Cu][3] = ak_[Cu][Ag][3] = 0.0000000000*eV/(Ao*Ao*Ao);

//-----------------
  a = 3.615*Ao;

  ak_[Cu][Cu][1] = ak_[Cu][Cu][1] = 61.73525861*eV/(a*a*a);
  ak_[Cu][Cu][2] = ak_[Cu][Cu][2] = -108.18467800*eV/(a*a*a);
  ak_[Cu][Cu][3] = ak_[Cu][Cu][3] = 57.00053948*eV/(a*a*a);
  ak_[Cu][Cu][4] = ak_[Cu][Cu][4] = -12.88796578*eV/(a*a*a);
  ak_[Cu][Cu][5] = ak_[Cu][Cu][5] = 39.16381901*eV/(a*a*a);
  ak_[Cu][Cu][6] = ak_[Cu][Cu][6] = 0.00000000*eV/(a*a*a);

  Ak_[Cu][Cu][1] = Ak_[Cu][Cu][1] = 10.03718305*(eV*eV) /(a*a*a);
  Ak_[Cu][Cu][2] = Ak_[Cu][Cu][2] = 17.06363299*(eV*eV) /(a*a*a);

  rk_[Cu][Cu][1] = rk_[Cu][Cu][1] = 1.225*a;
  rk_[Cu][Cu][2] = rk_[Cu][Cu][2] = 1.202*a;
  rk_[Cu][Cu][3] = rk_[Cu][Cu][3] = 1.154*a;
  rk_[Cu][Cu][4] = rk_[Cu][Cu][4] = 1.050*a;
  rk_[Cu][Cu][5] = rk_[Cu][Cu][5] = 0.866*a;
  rk_[Cu][Cu][6] = rk_[Cu][Cu][6] = 0.707*a;

  Rk_[Cu][Cu][1] = Rk_[Cu][Cu][1] = 1.225*a;
  Rk_[Cu][Cu][2] = Rk_[Cu][Cu][2] = 0.990*a;

//-----------------
  a = 4.086*Ao;

  ak_[Ag][Ag][1] = ak_[Ag][Ag][1] = 20.368404*eV/(a*a*a);
  ak_[Ag][Ag][2] = ak_[Ag][Ag][2] = -102.36075*eV/(a*a*a);
  ak_[Ag][Ag][3] = ak_[Ag][Ag][3] = 94.31277*eV/(a*a*a);
  ak_[Ag][Ag][4] = ak_[Ag][Ag][4] = -6.220051*eV/(a*a*a);
  ak_[Ag][Ag][5] = ak_[Ag][Ag][5] = 31.080889*eV/(a*a*a);
  ak_[Ag][Ag][6] = ak_[Ag][Ag][6] = 175.56047*eV/(a*a*a);

  Ak_[Ag][Ag][1] = Ak_[Ag][Ag][1] = 1.458761*(eV*eV) /(a*a*a);
  Ak_[Ag][Ag][2] = Ak_[Ag][Ag][2] = 42.946555*(eV*eV) /(a*a*a);

  rk_[Ag][Ag][1] = rk_[Ag][Ag][1] = 1.2247449*a;
  rk_[Ag][Ag][2] = rk_[Ag][Ag][2] = 1.1547054*a;
  rk_[Ag][Ag][3] = rk_[Ag][Ag][3] = 1.1180065*a;
  rk_[Ag][Ag][4] = rk_[Ag][Ag][4] = 1.0000000*a;
  rk_[Ag][Ag][5] = rk_[Ag][Ag][5] = 0.8660254*a;
  rk_[Ag][Ag][6] = rk_[Ag][Ag][6] = 0.7071068*a;

  Rk_[Ag][Ag][1] = Rk_[Ag][Ag][1] = 1.2247449*a;
  Rk_[Ag][Ag][2] = Rk_[Ag][Ag][2] = 1.0000000*a;

//-----------------
  a = 4.078*Ao;

  ak_[Au][Au][1] = ak_[Au][Au][1] = 29.059066*eV/(a*a*a);
  ak_[Au][Au][2] = ak_[Au][Au][2] = -153.14779*eV/(a*a*a);
  ak_[Au][Au][3] = ak_[Au][Au][3] = 148.17881*eV/(a*a*a);
  ak_[Au][Au][4] = ak_[Au][Au][4] = -22.20508*eV/(a*a*a);
  ak_[Au][Au][5] = ak_[Au][Au][5] = 72.71465*eV/(a*a*a);
  ak_[Au][Au][6] = ak_[Au][Au][6] = 199.26269*eV/(a*a*a);

  Ak_[Au][Au][1] = Ak_[Au][Au][1] = 21.930125*(eV*eV) /(a*a*a);
  Ak_[Au][Au][2] = Ak_[Au][Au][2] = 284.99631*(eV*eV) /(a*a*a);

  rk_[Au][Au][1] = rk_[Au][Au][1] = 1.2247449*a;
  rk_[Au][Au][2] = rk_[Au][Au][2] = 1.1547054*a;
  rk_[Au][Au][3] = rk_[Au][Au][3] = 1.1180065*a;
  rk_[Au][Au][4] = rk_[Au][Au][4] = 1.0000000*a;
  rk_[Au][Au][5] = rk_[Au][Au][5] = 0.8660254*a;
  rk_[Au][Au][6] = rk_[Au][Au][6] = 0.7071068*a;

  Rk_[Au][Au][1] = Rk_[Au][Au][1] = 1.1180065*a;
  Rk_[Au][Au][2] = Rk_[Au][Au][2] = 0.8660254*a;
/*
  R_[0]  = 5.0*Ao;
  R_[1]  = 5.5*Ao;
*/

/*
  for(int i1 = 0; i1 < 3; i1++)
    for(int i2 = 0; i2 < 3; i2++)
      for(int i3 = 0; i3 < 7; i3++)
        ak_[i1][i2][i3] *= koe;

  for(int i1 = 0; i1 < 3; i1++)
    for(int i2 = 0; i2 < 3; i2++)
      for(int i3 = 0; i3 < 3; i3++)
        Ak_[i1][i2][i3] *= koe;
*/
fillR_concat_();
}  

void
Ackland::fillR_concat_()
{
  Float r;

  AtomsArray atoms(3);
  atoms[0].ID = Cu_EL;
  atoms[1].ID = Ag_EL;
  atoms[2].ID = Au_EL;
  atoms.setAttributesByElementID();

for(size_t i = 0; i < atoms.size(); i++)
for(size_t j = 0; j < atoms.size(); j++)
{
Atom& atom1 = atoms[i];
Atom& atom2 = atoms[j];

         Float       x[2]; x[0] = 1.0*Ao; x[1] = 1.4*Ao;
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

  Float VAckland = 0.0;
  Float DerVAckland = 0.0;
  {

  size_t nk = (atom1.ID == atom2.ID)?6:3;

  for(size_t k = 1; k <= nk; k++)
  {
    Float rt = rk(k,atom1,atom2)-r;
    if (rt <= 0) continue;
    VAckland += ak(k,atom1,atom2)*rt*rt*rt;
    DerVAckland += ak(k,atom1,atom2)*3.0*rt*rt*(-1.0);
  }
  }

  v[1] = VAckland;
  dvdx[1] = DerVAckland/* = 0*/;
}

splines[e2i(atom1)][e2i(atom2)] = new Spline(x,v,dvdx/*,d2vdxdx*/);

}

for(size_t i = 0; i < ECOUNT; i++)
for(size_t j = 0; j < ECOUNT; j++)
{
//R_concat_[i][j] = 1.45*Ao;
//TRACE(R_concat_[i][j]/Ao);
}

}

}


