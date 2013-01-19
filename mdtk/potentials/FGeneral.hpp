/*
   The generalized interatomic potential class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013
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
  bool isHandled(const Atom& atom) const;
  bool isHandledPair(const Atom& atom1, const Atom& atom2) const;
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

  Float CosDihedral(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V = 0.0);

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
FGeneral::isHandled(const Atom& atom) const
{
  if (handledElements.find(atom.ID) != handledElements.end())
    return true;
  else
    return false;
}

inline
bool
FGeneral::isHandledPair(const Atom& atom1, const Atom& atom2) const
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
  Float e1 = scalarmul(ij.rv,ik.rv);
  Float e2 = ij.r_*ik.r_;

  Float CosT = e1/e2;

  REQUIRE(CosT>=-1.0-CosT_epsilon && CosT<=+1.0+CosT_epsilon);

  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  if (V != 0.0)
  {
    {
      Vector3D de1 = ij.rv+ik.rv;
      Vector3D de2 = ij.dr(ij.atom1)*ik.r_+ij.r_*ik.dr(ij.atom1);

      Vector3D dCosTheta = (de1*e2-e1*de2)/SQR(e2);

      ij.atom1.grad += dCosTheta*V;
    }
    {
      Vector3D de1 = -ik.rv;
      Vector3D de2 = ij.dr(ij.atom2)*ik.r_;

      Vector3D dCosTheta = (de1*e2-e1*de2)/SQR(e2);

      ij.atom2.grad += dCosTheta*V;
    }
    {
      Vector3D de1 = -ij.rv;
      Vector3D de2 = ij.r_*ik.dr(ik.atom2);

      Vector3D dCosTheta = (de1*e2-e1*de2)/SQR(e2);

      ik.atom2.grad += dCosTheta*V;
    }
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
Float
FGeneral::CosDihedral(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V)
{
  Vector3D ejik = vectormul(-(ij.rv),ik.rv);
  Vector3D eijl = vectormul(ij.rv,jl.rv);

  Float e1  = scalarmul(ejik,eijl);
  Float ejik_module_squared = ejik.module_squared();
  Float eijl_module_squared = eijl.module_squared();
  Float e2_1 = sqrt(ejik_module_squared);
  Float e2_2 = sqrt(eijl_module_squared);
  Float e2 = e2_1*e2_2;

  Float CosT = (e2 == 0.0)?0.0:e1/e2;

  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  if (V != 0.0 && e2 != 0.0)
  {
    AtomsPair il(ij.atom1,jl.atom2,true);
    AtomsPair jk(ij.atom2,ik.atom2,true);

    Float dejik_module_squared_dx_,
          dejik_module_squared_dy_,
          dejik_module_squared_dz_,
          deijl_module_squared_dx_,
          deijl_module_squared_dy_,
          deijl_module_squared_dz_;

    Vector3D de1;
    Vector3D de2;

    Vector3D dCosDihedral = 0.0;

    Float e2_2_e2_1 = e2_2/e2_1;
    Float e2_1_e2_2 = e2_1/e2_2;

    {
      de1.x =  jk.rv.z*eijl.y-ejik.y*jl.rv.z-jk.rv.y*eijl.z+ejik.z*jl.rv.y;
      de1.y = -jk.rv.z*eijl.x+ejik.x*jl.rv.z+jk.rv.x*eijl.z-ejik.z*jl.rv.x;
      de1.z =  jk.rv.y*eijl.x-ejik.x*jl.rv.y-jk.rv.x*eijl.y+ejik.y*jl.rv.x;
      dejik_module_squared_dx_ =  ejik.y*jk.rv.z-ejik.z*jk.rv.y;
      dejik_module_squared_dy_ = -ejik.x*jk.rv.z+ejik.z*jk.rv.x;
      dejik_module_squared_dz_ =  ejik.x*jk.rv.y-ejik.y*jk.rv.x;
      deijl_module_squared_dx_ = -eijl.y*jl.rv.z+eijl.z*jl.rv.y;
      deijl_module_squared_dy_ =  eijl.x*jl.rv.z-eijl.z*jl.rv.x;
      deijl_module_squared_dz_ = -eijl.x*jl.rv.y+eijl.y*jl.rv.x;

      de2.x = e2_2_e2_1*dejik_module_squared_dx_ + e2_1_e2_2*deijl_module_squared_dx_;
      de2.y = e2_2_e2_1*dejik_module_squared_dy_ + e2_1_e2_2*deijl_module_squared_dy_;
      de2.z = e2_2_e2_1*dejik_module_squared_dz_ + e2_1_e2_2*deijl_module_squared_dz_;
      dCosDihedral = (de1*e2-e1*de2)/SQR(e2);

      ij.atom1.grad += dCosDihedral*V;
    }
    {
      de1.x = -ik.rv.z*eijl.y+ejik.y*il.rv.z+ik.rv.y*eijl.z-ejik.z*il.rv.y;
      de1.y =  ik.rv.z*eijl.x-ejik.x*il.rv.z-ik.rv.x*eijl.z+ejik.z*il.rv.x;
      de1.z = -ik.rv.y*eijl.x+ejik.x*il.rv.y+ik.rv.x*eijl.y-ejik.y*il.rv.x;
      dejik_module_squared_dx_ = -ejik.y*ik.rv.z+ejik.z*ik.rv.y;
      dejik_module_squared_dy_ =  ejik.x*ik.rv.z-ejik.z*ik.rv.x;
      dejik_module_squared_dz_ = -ejik.x*ik.rv.y+ejik.y*ik.rv.x;
      deijl_module_squared_dx_ =  eijl.y*il.rv.z-eijl.z*il.rv.y;
      deijl_module_squared_dy_ = -eijl.x*il.rv.z+eijl.z*il.rv.x;
      deijl_module_squared_dz_ =  eijl.x*il.rv.y-eijl.y*il.rv.x;

      de2.x = e2_2_e2_1*dejik_module_squared_dx_ + e2_1_e2_2*deijl_module_squared_dx_;
      de2.y = e2_2_e2_1*dejik_module_squared_dy_ + e2_1_e2_2*deijl_module_squared_dy_;
      de2.z = e2_2_e2_1*dejik_module_squared_dz_ + e2_1_e2_2*deijl_module_squared_dz_;
      dCosDihedral = (de1*e2-e1*de2)/SQR(e2);

      ij.atom2.grad += dCosDihedral*V;
    }
    {
      de1.x =  ij.rv.z*eijl.y-ij.rv.y*eijl.z;
      de1.y = -ij.rv.z*eijl.x+ij.rv.x*eijl.z;
      de1.z =  ij.rv.y*eijl.x-ij.rv.x*eijl.y;
      dejik_module_squared_dx_ =  ejik.y*ij.rv.z-ejik.z*ij.rv.y;
      dejik_module_squared_dy_ = -ejik.x*ij.rv.z+ejik.z*ij.rv.x;
      dejik_module_squared_dz_ =  ejik.x*ij.rv.y-ejik.y*ij.rv.x;

      de2.x = e2_2_e2_1*dejik_module_squared_dx_;
      de2.y = e2_2_e2_1*dejik_module_squared_dy_;
      de2.z = e2_2_e2_1*dejik_module_squared_dz_;
      dCosDihedral = (de1*e2-e1*de2)/SQR(e2);

      ik.atom2.grad += dCosDihedral*V;
    }
    {
      de1.x = -ejik.y*ij.rv.z+ejik.z*ij.rv.y;
      de1.y =  ejik.x*ij.rv.z-ejik.z*ij.rv.x;
      de1.z = -ejik.x*ij.rv.y+ejik.y*ij.rv.x;
      deijl_module_squared_dx_ = -eijl.y*ij.rv.z+eijl.z*ij.rv.y;
      deijl_module_squared_dy_ =  eijl.x*ij.rv.z-eijl.z*ij.rv.x;
      deijl_module_squared_dz_ = -eijl.x*ij.rv.y+eijl.y*ij.rv.x;

      de2.x = e2_1_e2_2*deijl_module_squared_dx_;
      de2.y = e2_1_e2_2*deijl_module_squared_dy_;
      de2.z = e2_1_e2_2*deijl_module_squared_dz_;
      dCosDihedral = (de1*e2-e1*de2)/SQR(e2);

      jl.atom2.grad += dCosDihedral*V;
    }
  }

  return CosT;
}

} // namespace mdtk

#endif
