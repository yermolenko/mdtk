/*
   The generalized interatomic potential class (header file).

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

#ifndef mdtk_FGeneral_hpp
#define mdtk_FGeneral_hpp

#include <mdtk/config.hpp>
#include <set>
#include <utility>
#include <mdtk/Atom.hpp>
#include <mdtk/potentials/NeighbourList.hpp>
#include <mdtk/potentials/AtomsPair.hpp>
#include <mdtk/consts.hpp>

#include <limits>
#include <cmath>

#define FGENERAL_CHECKS

#define FGENERAL_OPTIMIZED

namespace mdtk
{

using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::sqrt;
using std::acos;
using std::FILE;

class FGeneral
{
public:
  virtual Float getRcutoff() const = 0;

  std::set<ElementID> handledElements;
  std::set<std::pair<ElementID,ElementID> > handledElementPairs;
  bool isHandled(Atom& atom) const;
  bool isHandledPair(Atom& atom1, Atom& atom2) const;
protected:
  NeighbourList nl;
public:
  AtomRefsContainer& NL(const Atom& atom)
  {
    return nl.nl[atom.globalIndex];
  }

  virtual
  void NL_checkRequestUpdate(AtomsArray& atoms) {nl.checkRequestUpdate(atoms);}
  virtual
  void NL_UpdateIfNeeded(AtomsArray& atoms) {nl.UpdateIfNeeded(atoms);}
  virtual
  void NL_init(AtomsArray& atoms) {nl.init(atoms);}

  virtual
  void incDisplacement(Atom& atom, Vector3D inc)
  {
    nl.displacements[atom.globalIndex] += inc;
  }

  FGeneral();
  virtual ~FGeneral(){;}

  virtual Float operator()(AtomsArray&) = 0;
private:
public:
  Float SinTheta(AtomsPair& ij, AtomsPair& ik);
  Float CosTheta(AtomsPair& ij, AtomsPair& ik, const Float V = 0.0);
  Vector3D dCosTheta(const AtomsPair& ij, const AtomsPair& ik, const Atom &datom) const;

  Float CosDihedral(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V = 0.0);
  Vector3D dCosDihedral(const AtomsPair& ij, const AtomsPair& ik, const AtomsPair& jl, const Atom &datom) const;

  virtual
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    nl.SaveToStream(os,smode);
  }
  virtual
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    nl.LoadFromStream(is,smode);
  }
};

inline
bool
FGeneral::isHandled(Atom& atom) const
{
  if (handledElements.find(atom.ID) != handledElements.end())
    return true;
  else
    return false;
}

inline
bool
FGeneral::isHandledPair(Atom& atom1, Atom& atom2) const
{
  if (handledElementPairs.find(std::make_pair(atom1.ID,atom2.ID)) != handledElementPairs.end() )
    return true;
  else
    return false;
}

#define CosT_epsilon 1e-7

inline
Float
FGeneral::CosTheta(AtomsPair& ij, AtomsPair& ik, const Float V)
{
  Float CosT = scalarmul(ij.rv,ik.rv)/(ij.r_*ik.r_);

  REQUIRE(CosT>=-1.0-CosT_epsilon && CosT<=+1.0+CosT_epsilon);

  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  if (V != 0.0)
  {
    ij.atom1.grad += dCosTheta(ij,ik,ij.atom1)*V;
    ij.atom2.grad += dCosTheta(ij,ik,ij.atom2)*V;
    ik.atom2.grad += dCosTheta(ij,ik,ik.atom2)*V;
  }

  return CosT;
}

inline
Float
FGeneral::SinTheta(AtomsPair& ij, AtomsPair& ik)
{
  Float SinT = vectormul(ij.rv,ik.rv).module()/(ij.r_*ik.r_);

  REQUIRE(SinT>=-1.0-CosT_epsilon && SinT<=+1.0+CosT_epsilon);

  if (SinT<-1.0) SinT = -1.0;
  if (SinT>+1.0) SinT = +1.0;

  return SinT;
}

inline
Vector3D
FGeneral::dCosTheta(const AtomsPair& ij, const AtomsPair& ik, const Atom &datom) const
{
  Float e1 = scalarmul(ij.rv,ik.rv);
  Float e2 = ij.r_*ik.r_;

  Vector3D dscalarmul;
  {
    if (&ij.atom1 == &datom)
      dscalarmul = -ij.rv-ik.rv;
    else if (&ij.atom2 == &datom)
      dscalarmul = ik.rv;
    else if (&ik.atom2 == &datom)
      dscalarmul = ij.rv;
    else
      dscalarmul = 0.0;
  };

  Vector3D de1 = -dscalarmul;
  Vector3D de2 = ij.dr(datom)*ik.r_+ij.r_*ik.dr(datom);

  return (de1*e2-e1*de2)/SQR(e2);
}

inline
Float
FGeneral::CosDihedral(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V)
{
  AtomsPair ji = -ij;

  Vector3D ejik = vectormul(ji.rv,ik.rv);
  Vector3D eijl = vectormul(ij.rv,jl.rv);

  Float scmul = scalarmul(ejik,eijl);
  Float modprod = (ejik.module()*eijl.module());

  Float CosT;

  if (/*scmul == 0.0 &&*/ modprod == 0.0)
    CosT = 0.0;
  else

  CosT = scmul/modprod;

  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  if (V != 0.0)
  {
    ij.atom1.grad += dCosDihedral(ij,ik,jl,ij.atom1)*V;
    ij.atom2.grad += dCosDihedral(ij,ik,jl,ij.atom2)*V;
    ik.atom2.grad += dCosDihedral(ij,ik,jl,ik.atom2)*V;
    jl.atom2.grad += dCosDihedral(ij,ik,jl,jl.atom2)*V;
  }

  return CosT;
}


#define dCosDihedral_COMMON_INC \
\
Float de1_x = dscalar_dx_;\
Float de2_1_x = 0.5*pow(ejik_module_squared,-0.5)*dejik_module_squared_dx_;\
Float de2_2_x = 0.5*pow(eijl_module_squared,-0.5)*deijl_module_squared_dx_;\
Float de2_x = de2_1_x*e2_2 + e2_1*de2_2_x;\
\
Float de1_y = dscalar_dy_;\
Float de2_1_y = 0.5*pow(ejik_module_squared,-0.5)*dejik_module_squared_dy_;\
Float de2_2_y = 0.5*pow(eijl_module_squared,-0.5)*deijl_module_squared_dy_;\
Float de2_y = de2_1_y*e2_2 + e2_1*de2_2_y;\
\
Float de1_z = dscalar_dz_;\
Float de2_1_z = 0.5*pow(ejik_module_squared,-0.5)*dejik_module_squared_dz_;\
Float de2_2_z = 0.5*pow(eijl_module_squared,-0.5)*deijl_module_squared_dz_;\
Float de2_z = de2_1_z*e2_2 + e2_1*de2_2_z;\
\
  if (e2 == 0.0)\
    return 0.0;\
  else \
\
    return Vector3D\
    (\
      (de1_x*e2-e1*de2_x)/SQR(e2),\
      (de1_y*e2-e1*de2_y)/SQR(e2),\
      (de1_z*e2-e1*de2_z)/SQR(e2)\
    );\


inline
Vector3D
FGeneral::dCosDihedral(const AtomsPair& ij, const AtomsPair& ik, const AtomsPair& jl, const Atom &datom) const
{
  Vector3D ij_rv = ij.rv;
  Float  ij_rv_x = ij_rv.x;
  Float  ij_rv_y = ij_rv.y;
  Float  ij_rv_z = ij_rv.z;
  Vector3D ik_rv = ik.rv;
  Float  ik_rv_x = ik_rv.x;
  Float  ik_rv_y = ik_rv.y;
  Float  ik_rv_z = ik_rv.z;

  Vector3D ji_rv = (-ij).rv;
  Float  ji_rv_x = ji_rv.x;
  Float  ji_rv_y = ji_rv.y;
  Float  ji_rv_z = ji_rv.z;
  Vector3D jl_rv = jl.rv;
  Float  jl_rv_x = jl_rv.x;
  Float  jl_rv_y = jl_rv.y;
  Float  jl_rv_z = jl_rv.z;

  Vector3D ki_rv = (-ik).rv;
  Float  ki_rv_x = ki_rv.x;
  Float  ki_rv_y = ki_rv.y;
  Float  ki_rv_z = ki_rv.z;

  Vector3D lj_rv = (-jl).rv;
  Float  lj_rv_x = lj_rv.x;
  Float  lj_rv_y = lj_rv.y;
  Float  lj_rv_z = lj_rv.z;

//  Vector3D il_rv = r_vec(atom_i,atom_l);
  AtomsPair il(ij.atom1,jl.atom2,true);
  Vector3D il_rv = il.rv;//ij_rv+jl_rv;
  Float  il_rv_x = il_rv.x;
  Float  il_rv_y = il_rv.y;
  Float  il_rv_z = il_rv.z;

//  Vector3D jk_rv = r_vec(atom_j,atom_k);
  AtomsPair jk(ij.atom2,ik.atom2,true);
  Vector3D jk_rv = jk.rv;//ji_rv+ik_rv;
  Float  jk_rv_x = jk_rv.x;
  Float  jk_rv_y = jk_rv.y;
  Float  jk_rv_z = jk_rv.z;

//  Vector3D kj_rv = r_vec(atom_k,atom_j);
  AtomsPair kj = -jk;
  Vector3D kj_rv = kj.rv;//ki_rv+ij_rv;
  Float  kj_rv_x = kj_rv.x;
  Float  kj_rv_y = kj_rv.y;
  Float  kj_rv_z = kj_rv.z;

//  Vector3D li_rv = r_vec(atom_l,atom_i);
  AtomsPair li = -il;
  Vector3D li_rv = li.rv;//lj_rv+ji_rv;
  Float  li_rv_x = li_rv.x;
  Float  li_rv_y = li_rv.y;
  Float  li_rv_z = li_rv.z;

Float e1  = (ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x);
Float ejik_module_squared = pow(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y,2.0)+pow(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z,2.0)+pow(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x,2.0);
Float eijl_module_squared = pow(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y,2.0)+pow(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z,2.0)+pow(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x,2.0);
Float e2_1 = sqrt(ejik_module_squared);
Float e2_2 = sqrt(eijl_module_squared);
Float e2 = e2_1*e2_2;

Float dscalar_dx_,dscalar_dy_,dscalar_dz_;
Float dejik_module_squared_dx_,
      dejik_module_squared_dy_,
      dejik_module_squared_dz_,
      deijl_module_squared_dx_,
      deijl_module_squared_dy_,
      deijl_module_squared_dz_;

  if (&ij.atom1 == &datom)
  {
dscalar_dx_ = jk_rv_z*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*lj_rv_z+kj_rv_y*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x)+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*jl_rv_y;
      dscalar_dy_ = kj_rv_z*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*
jl_rv_z+jk_rv_x*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x)+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*lj_rv_x;
      dscalar_dz_ = jk_rv_y*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*
lj_rv_y+kj_rv_x*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*jl_rv_x;
      dejik_module_squared_dx_ = 2.0*(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*jk_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*kj_rv_y;
      dejik_module_squared_dy_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*kj_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*jk_rv_x;
      dejik_module_squared_dz_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*jk_rv_y+2.0*(ji_rv_z
*ik_rv_x-ji_rv_x*ik_rv_z)*kj_rv_x;
      deijl_module_squared_dx_ = 2.0*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)*lj_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*jl_rv_y;
      deijl_module_squared_dy_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*jl_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*lj_rv_x;
      deijl_module_squared_dz_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*lj_rv_y+2.0*(ij_rv_z
*jl_rv_x-ij_rv_x*jl_rv_z)*jl_rv_x;

dCosDihedral_COMMON_INC;
  }
  else if (&ij.atom2 == &datom)
  {
dscalar_dx_ = ki_rv_z*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*il_rv_z+ik_rv_y*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x)+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*li_rv_y;
      dscalar_dy_ = ik_rv_z*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*
li_rv_z+ki_rv_x*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x)+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*il_rv_x;
      dscalar_dz_ = ki_rv_y*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*
il_rv_y+ik_rv_x*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*li_rv_x;
      dejik_module_squared_dx_ = 2.0*(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*ki_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*ik_rv_y;
      dejik_module_squared_dy_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ik_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*ki_rv_x;
      dejik_module_squared_dz_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ki_rv_y+2.0*(ji_rv_z
*ik_rv_x-ji_rv_x*ik_rv_z)*ik_rv_x;
      deijl_module_squared_dx_ = 2.0*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)*il_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*li_rv_y;
      deijl_module_squared_dy_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*li_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*il_rv_x;
      deijl_module_squared_dz_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*il_rv_y+2.0*(ij_rv_z
*jl_rv_x-ij_rv_x*jl_rv_z)*li_rv_x;

dCosDihedral_COMMON_INC;
  }
  else if (&ik.atom2 == &datom)
  {
dscalar_dx_ = ij_rv_z*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)+ji_rv_y*(ij_rv_x*jl_rv_y-ij_rv_y*jl_rv_x);
      dscalar_dy_ = ji_rv_z*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+ij_rv_x*(ij_rv_x*jl_rv_y-ij_rv_y*
jl_rv_x);
      dscalar_dz_ = ij_rv_y*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)+ji_rv_x*(ij_rv_z*jl_rv_x-ij_rv_x*
jl_rv_z);
      dejik_module_squared_dx_ = 2.0*(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*ij_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*ji_rv_y;
      dejik_module_squared_dy_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ji_rv_z+2.0*(ji_rv_x
*ik_rv_y-ji_rv_y*ik_rv_x)*ij_rv_x;
      dejik_module_squared_dz_ = 2.0*(ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ij_rv_y+2.0*(ji_rv_z
*ik_rv_x-ji_rv_x*ik_rv_z)*ji_rv_x;
      deijl_module_squared_dx_ = 0.0;
      deijl_module_squared_dy_ = 0.0;
      deijl_module_squared_dz_ = 0.0;

dCosDihedral_COMMON_INC;
  }
  else if (&jl.atom2 == &datom)
  {
dscalar_dx_ = (ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*ji_rv_z+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*ij_rv_y;
      dscalar_dy_ = (ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ij_rv_z+(ji_rv_x*ik_rv_y-ji_rv_y*ik_rv_x)*
ji_rv_x;
      dscalar_dz_ = (ji_rv_y*ik_rv_z-ji_rv_z*ik_rv_y)*ji_rv_y+(ji_rv_z*ik_rv_x-ji_rv_x*ik_rv_z)*
ij_rv_x;
      dejik_module_squared_dx_ = 0.0;
      dejik_module_squared_dy_ = 0.0;
      dejik_module_squared_dz_ = 0.0;
      deijl_module_squared_dx_ = 2.0*(ij_rv_z*jl_rv_x-ij_rv_x*jl_rv_z)*ji_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*ij_rv_y;
      deijl_module_squared_dy_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*ij_rv_z+2.0*(ij_rv_x
*jl_rv_y-ij_rv_y*jl_rv_x)*ji_rv_x;
      deijl_module_squared_dz_ = 2.0*(ij_rv_y*jl_rv_z-ij_rv_z*jl_rv_y)*ji_rv_y+2.0*(ij_rv_z
*jl_rv_x-ij_rv_x*jl_rv_z)*ij_rv_x;

dCosDihedral_COMMON_INC;
  }
  else
    return 0.0;
}


} // namespace mdtk

#endif
