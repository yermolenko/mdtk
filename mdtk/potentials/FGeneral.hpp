/*
   The generalized interatomic potential class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr
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

#ifndef mdtk_FGeneral_hpp
#define mdtk_FGeneral_hpp

#include <mdtk/config.hpp>
#include <set>
#include <utility>
#include <mdtk/Atom.hpp>
#include <mdtk/potentials/NeighbourList.hpp>
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
  virtual bool probablyAreNeighbours(Atom& atom1, Atom& atom2) {return true;}
  AtomsContainer& NL(Atom& atom)
  {
    return nl.nl[atom.globalIndex];
  }  
  AtomsContainer& NL_with_self(Atom& atom)
  {
    return nl.nl_with_self[atom.globalIndex];
  }  

  virtual
  void NL_checkRequestUpdate(AtomsContainer& atoms) {nl.checkRequestUpdate(atoms);}
  virtual
  void NL_UpdateIfNeeded(AtomsContainer& atoms) {nl.UpdateIfNeeded(atoms);}
  virtual
  void NL_init(AtomsContainer& atoms) {nl.init(atoms);}

  virtual
  void incDisplacement(Atom& atom, Vector3D inc)
  {
    nl.displacements[atom.globalIndex] += inc;
  }  

  FGeneral();
  virtual ~FGeneral(){;}

  virtual Float operator()(AtomsContainer&) = 0;
  virtual Vector3D grad(Atom &a1,AtomsContainer&) = 0;
private:
public:
  bool ontouch_enabled;
  virtual void onTouch(Atom& a)  = 0;

  Vector3D  r_vec(Atom &atom1, Atom &atom2) ;  
  Vector3D  r_vec_no_touch(Atom &atom1, Atom &atom2) const ;  
  void      r_vec_touch_only(Atom &atom1, Atom &atom2) ;  

  Float r_vec_module(Atom &atom1,Atom &atom2) ; 
  Float r_vec_module_squared(Atom &atom1,Atom &atom2) ; 
  Float r_vec_module_no_touch(Atom &atom1,Atom &atom2) ; 
  Float r_vec_module_squared_no_touch(Atom &atom1,Atom &atom2) ; 
    Vector3D dr_vec_module(Atom &atom1,Atom &atom2, Atom &datom) ; 
  Float rsqr_vec_module(Atom &atom1,Atom &atom2) ; 
    Vector3D drsqr_vec_module(Atom &atom1,Atom &atom2, Atom &datom) ; 
  Float CosTheta(Atom &atom_i,Atom &atom_j,Atom &atom_k) ;
    Vector3D dCosTheta(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) ;
  Float CosThetaJIK(Atom &atom_i,Atom &atom_j,Atom &atom_k) ;
    Vector3D dCosThetaJIK(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) ;

  Float CosDihedral(Atom &atom_i,Atom &atom_j,Atom &atom_k,Atom &atom_l) ;
    Vector3D dCosDihedral(Atom &atom_i,Atom &atom_j,Atom &atom_k,Atom &atom_l, Atom &datom) ;

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
public:  
  struct R_VEC_EX{
    Atom* atom1;
    Atom* atom2;
    Float new_r_vec;
    Atom* atom1_bak;
    Atom* atom2_bak;
    R_VEC_EX():atom1(NULL),atom2(NULL),new_r_vec(0.0),
      atom1_bak(NULL),atom2_bak(NULL){;}
  }r_vec_ex;
  void disable_r_vec_ex()
  {
    r_vec_ex.atom1_bak = r_vec_ex.atom1;
    r_vec_ex.atom2_bak = r_vec_ex.atom2;
    r_vec_ex.atom1 = NULL;
    r_vec_ex.atom2 = NULL;
  }
  void enable_r_vec_ex()
  {
    r_vec_ex.atom1 = r_vec_ex.atom1_bak;
    r_vec_ex.atom2 = r_vec_ex.atom2_bak;
    r_vec_ex.atom1_bak = NULL;
    r_vec_ex.atom2_bak = NULL;
  }
  void set_r_vec_exception(Atom& a1,Atom& a2,Float new_r_vec)
  {
    r_vec_ex.atom1 = &a1;
    r_vec_ex.atom2 = &a2;
    r_vec_ex.new_r_vec = new_r_vec;
    REQUIRE(r_vec_ex.atom1 != r_vec_ex.atom2);
    REQUIRE(new_r_vec != 0.0);
  }  
  void cleat_r_vec_exception()
  {
    r_vec_ex.atom1 = NULL;
    r_vec_ex.atom2 = NULL;
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


inline
Vector3D
FGeneral::r_vec(Atom &atom1, Atom &atom2) 
{
  if (ontouch_enabled)
  {
    onTouch(atom1);
    onTouch(atom2);
  }  
  return r_vec_no_touch(atom1, atom2);
}

inline
void
FGeneral::r_vec_touch_only(Atom &atom1, Atom &atom2) 
{
  if (ontouch_enabled)
  {
    onTouch(atom1);
    onTouch(atom2);
  }  
}

inline
Float
FGeneral::r_vec_module_no_touch(Atom &atom1,Atom &atom2) 
{
#ifdef FGENERAL_CHECKS  
  REQUIREM(&atom1 != &atom2,"r_vec_module_no_touch: &atom1 == &atom2");
#endif  
  return r_vec_no_touch(atom1,atom2).module();
}   

inline
Float
FGeneral::r_vec_module_squared_no_touch(Atom &atom1,Atom &atom2) 
{
#ifdef FGENERAL_CHECKS  
  REQUIREM(&atom1 != &atom2,"r_vec_module_no_touch: &atom1 == &atom2");
#endif  
  return r_vec_no_touch(atom1,atom2).module_squared();
}   

inline
Vector3D
FGeneral::r_vec_no_touch(Atom &atom1, Atom &atom2) const
{
  Vector3D val;

  val = depos(atom1,atom2);

  if ( (r_vec_ex.atom1 != NULL) &&
       (
       (r_vec_ex.atom1 == &atom1 && r_vec_ex.atom2 == &atom2) ||
       (r_vec_ex.atom1 == &atom2 && r_vec_ex.atom2 == &atom1)
       )
     )
    return val.normalized()*r_vec_ex.new_r_vec;

  return val;
}  

inline
Float
FGeneral::r_vec_module(Atom &atom1,Atom &atom2) 
{
#ifdef FGENERAL_CHECKS  
  REQUIREM(&atom1 != &atom2,"r_vec_module: &atom1 == &atom2");
#endif  
  return r_vec(atom1,atom2).module();
}   

inline
Float
FGeneral::r_vec_module_squared(Atom &atom1,Atom &atom2) 
{
#ifdef FGENERAL_CHECKS  
  REQUIREM(&atom1 != &atom2,"r_vec_module: &atom1 == &atom2");
#endif  
  return r_vec(atom1,atom2).module_squared();
}   

inline
Vector3D
FGeneral::dr_vec_module(Atom &atom1,Atom &atom2, Atom &datom) 
{
#ifdef FGENERAL_OPTIMIZED  
  Vector3D dvar = drsqr_vec_module(atom1,atom2,datom);
  if (dvar != 0.0) 
    return 0.5/r_vec_module(atom1,atom2)*dvar;
  else
    return 0.0;
#else    
  return 0.5/r_vec_module(atom1,atom2)*drsqr_vec_module(atom1,atom2,datom);
#endif  
}   

inline
Float
FGeneral::rsqr_vec_module(Atom &atom1,Atom &atom2) 
{
  return SQR(r_vec(atom1,atom2).module());
}   

inline
Vector3D
FGeneral::drsqr_vec_module(Atom &atom1,Atom &atom2, Atom &datom) 
{
#ifdef FGENERAL_CHECKS  
  REQUIREM(&atom1 != &atom2,"drsqr_vec_module: &atom1 == &atom2");
#endif  
  if (&atom1 == &datom)
    return 2.0*r_vec(atom1,atom2);
  else if (&atom2 == &datom) 
    return -2.0*r_vec(atom1,atom2);
  else
    return 0.0;
}   

#define CosT_epsilon 1e-7

inline
Float
FGeneral::CosTheta(Atom &atom_i,Atom &atom_j,Atom &atom_k) 
{
  Vector3D rij = r_vec(atom_i,atom_j);
  Vector3D rik = r_vec(atom_i,atom_k);

  Float CosT = scalarmul(rij,rik)/(rij.module()*rik.module());

  REQUIRE(CosT>=-1.0-CosT_epsilon && CosT<=+1.0+CosT_epsilon);
    
  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  return CosT;
}

inline
Vector3D
FGeneral::dCosTheta(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) 
{
  Vector3D rij = r_vec(atom_i,atom_j);
  Vector3D rik = r_vec(atom_i,atom_k);

  Float e1 = scalarmul(rij,rik);
  Float e2 = rij.module()*rik.module();

  Vector3D drij = dr_vec_module(atom_i,atom_j,datom);
  Vector3D drik = dr_vec_module(atom_i,atom_k,datom);

  Vector3D dscalarmul(0.0,0.0,0.0);
  {
    if (&atom_i == &datom)
      dscalarmul = -rij-rik;
    else if (&atom_j == &datom)
      dscalarmul = rik;
    else if (&atom_k == &datom)
      dscalarmul = rij;
    else
      dscalarmul = 0.0;
  };

  Vector3D de1 = -dscalarmul;
  Vector3D de2 = drij*rik.module()+rij.module()*drik;

  return (de1*e2-e1*de2)/SQR(e2);
}

inline
Float
FGeneral::CosThetaJIK(Atom &atom_i,Atom &atom_j,Atom &atom_k) 
{
  Vector3D rji = r_vec(atom_j,atom_i);
  Vector3D rki = r_vec(atom_k,atom_i);

  Float CosT = scalarmul(rji,rki)/(rji.module()*rki.module());

  REQUIRE(CosT>=-1.0-CosT_epsilon && CosT<=+1.0+CosT_epsilon);
    
  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

  return CosT;
}

inline
Vector3D
FGeneral::dCosThetaJIK(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) 
{
  Vector3D rji = r_vec(atom_j,atom_i);
  Vector3D rki = r_vec(atom_k,atom_i);

  Float e1 = scalarmul(rji,rki);
  Float e2 = rji.module()*rki.module();

  Vector3D drji = dr_vec_module(atom_j,atom_i,datom);
  Vector3D drki = dr_vec_module(atom_k,atom_i,datom);

  Vector3D dscalarmul(0.0,0.0,0.0);
  {
    if (&atom_i == &datom)
      dscalarmul = -rji-rki;
    else if (&atom_j == &datom)
      dscalarmul = rki;
    else if (&atom_k == &datom)
      dscalarmul = rji;
    else
      dscalarmul = 0.0;
  };

  Vector3D de1 = dscalarmul;
  Vector3D de2 = drji*rki.module()+rji.module()*drki;

  return (de1*e2-e1*de2)/SQR(e2);
}

inline
Float
FGeneral::CosDihedral(Atom &atom_i,Atom &atom_j,Atom &atom_k,Atom &atom_l) 
{
  Vector3D rij = r_vec(atom_i,atom_j);
  Vector3D rji = r_vec(atom_j,atom_i);
  Vector3D rik = r_vec(atom_i,atom_k);
  Vector3D rjl = r_vec(atom_j,atom_l);

  Vector3D ejik = vectormul(rji,rik);
  Vector3D eijl = vectormul(rij,rjl);

  Float scmul = scalarmul(ejik,eijl);
  Float modprod = (ejik.module()*eijl.module());

  Float CosT;

  if (/*scmul == 0.0 &&*/ modprod == 0.0)
    CosT = 0.0;
  else 

  CosT = scmul/modprod;
    
  if (CosT<-1.0) CosT = -1.0;
  if (CosT>+1.0) CosT = +1.0;

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
FGeneral::dCosDihedral(Atom &atom_i,Atom &atom_j,Atom &atom_k,Atom &atom_l, Atom &datom) 
{
  Vector3D rij = r_vec(atom_i,atom_j);
    Float  rij_x = rij.x;
    Float  rij_y = rij.y;
    Float  rij_z = rij.z;
  Vector3D rik = r_vec(atom_i,atom_k);
    Float  rik_x = rik.x;
    Float  rik_y = rik.y;
    Float  rik_z = rik.z;

  Vector3D rji = r_vec(atom_j,atom_i);
    Float  rji_x = rji.x;
    Float  rji_y = rji.y;
    Float  rji_z = rji.z;
  Vector3D rjl = r_vec(atom_j,atom_l);
    Float  rjl_x = rjl.x;
    Float  rjl_y = rjl.y;
    Float  rjl_z = rjl.z;

  Vector3D rki = r_vec(atom_k,atom_i);
    Float  rki_x = rki.x;
    Float  rki_y = rki.y;
    Float  rki_z = rki.z;
/*
  Vector3D rkl = r_vec(atom_k,atom_l);
    Float  rkl_x = rkl.x;
    Float  rkl_y = rkl.y;
    Float  rkl_z = rkl.z;
*/

  Vector3D rlj = r_vec(atom_l,atom_j);
    Float  rlj_x = rlj.x;
    Float  rlj_y = rlj.y;
    Float  rlj_z = rlj.z;
/*
  Vector3D rlk = r_vec(atom_l,atom_k);
    Float  rlk_x = rlk.x;
    Float  rlk_y = rlk.y;
    Float  rlk_z = rlk.z;
*/

//  Vector3D ril = r_vec(atom_i,atom_l); 
  Vector3D ril = rij+rjl;
    Float  ril_x = ril.x;
    Float  ril_y = ril.y;
    Float  ril_z = ril.z;

//  Vector3D rjk = r_vec(atom_j,atom_k); 
  Vector3D rjk = rji+rik; 
    Float  rjk_x = rjk.x;
    Float  rjk_y = rjk.y;
    Float  rjk_z = rjk.z;

//  Vector3D rkj = r_vec(atom_k,atom_j); 
  Vector3D rkj = rki+rij; 
    Float  rkj_x = rkj.x;
    Float  rkj_y = rkj.y;
    Float  rkj_z = rkj.z;

//  Vector3D rli = r_vec(atom_l,atom_i); 
  Vector3D rli = rlj+rji; 
    Float  rli_x = rli.x;
    Float  rli_y = rli.y;
    Float  rli_z = rli.z;

Float e1  = (rji_y*rik_z-rji_z*rik_y)*(rij_y*rjl_z-rij_z*rjl_y)+(rji_z*rik_x-rji_x*rik_z)*(rij_z*rjl_x-rij_x*rjl_z)+(rji_x*rik_y-rji_y*rik_x)*(rij_x*rjl_y-rij_y*rjl_x);
Float ejik_module_squared = pow(rji_y*rik_z-rji_z*rik_y,2.0)+pow(rji_z*rik_x-rji_x*rik_z,2.0)+pow(rji_x*rik_y-rji_y*rik_x,2.0);
Float eijl_module_squared = pow(rij_y*rjl_z-rij_z*rjl_y,2.0)+pow(rij_z*rjl_x-rij_x*rjl_z,2.0)+pow(rij_x*rjl_y-rij_y*rjl_x,2.0);
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

  if (&atom_i == &datom)
  {
dscalar_dx_ = rjk_z*(rij_z*rjl_x-rij_x*rjl_z)+(rji_z*rik_x-rji_x*rik_z)*rlj_z+rkj_y*(rij_x*rjl_y-rij_y*rjl_x)+(rji_x*rik_y-rji_y*rik_x)*rjl_y;
      dscalar_dy_ = rkj_z*(rij_y*rjl_z-rij_z*rjl_y)+(rji_y*rik_z-rji_z*rik_y)*
rjl_z+rjk_x*(rij_x*rjl_y-rij_y*rjl_x)+(rji_x*rik_y-rji_y*rik_x)*rlj_x;
      dscalar_dz_ = rjk_y*(rij_y*rjl_z-rij_z*rjl_y)+(rji_y*rik_z-rji_z*rik_y)*
rlj_y+rkj_x*(rij_z*rjl_x-rij_x*rjl_z)+(rji_z*rik_x-rji_x*rik_z)*rjl_x;
      dejik_module_squared_dx_ = 2.0*(rji_z*rik_x-rji_x*rik_z)*rjk_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rkj_y;
      dejik_module_squared_dy_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rkj_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rjk_x;
      dejik_module_squared_dz_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rjk_y+2.0*(rji_z
*rik_x-rji_x*rik_z)*rkj_x;
      deijl_module_squared_dx_ = 2.0*(rij_z*rjl_x-rij_x*rjl_z)*rlj_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*rjl_y;
      deijl_module_squared_dy_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*rjl_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*rlj_x;
      deijl_module_squared_dz_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*rlj_y+2.0*(rij_z
*rjl_x-rij_x*rjl_z)*rjl_x;

dCosDihedral_COMMON_INC;
  }  
  else if (&atom_j == &datom)
  {
dscalar_dx_ = rki_z*(rij_z*rjl_x-rij_x*rjl_z)+(rji_z*rik_x-rji_x*rik_z)*ril_z+rik_y*(rij_x*rjl_y-rij_y*rjl_x)+(rji_x*rik_y-rji_y*rik_x)*rli_y;
      dscalar_dy_ = rik_z*(rij_y*rjl_z-rij_z*rjl_y)+(rji_y*rik_z-rji_z*rik_y)*
rli_z+rki_x*(rij_x*rjl_y-rij_y*rjl_x)+(rji_x*rik_y-rji_y*rik_x)*ril_x;
      dscalar_dz_ = rki_y*(rij_y*rjl_z-rij_z*rjl_y)+(rji_y*rik_z-rji_z*rik_y)*
ril_y+rik_x*(rij_z*rjl_x-rij_x*rjl_z)+(rji_z*rik_x-rji_x*rik_z)*rli_x;
      dejik_module_squared_dx_ = 2.0*(rji_z*rik_x-rji_x*rik_z)*rki_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rik_y;
      dejik_module_squared_dy_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rik_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rki_x;
      dejik_module_squared_dz_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rki_y+2.0*(rji_z
*rik_x-rji_x*rik_z)*rik_x;
      deijl_module_squared_dx_ = 2.0*(rij_z*rjl_x-rij_x*rjl_z)*ril_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*rli_y;
      deijl_module_squared_dy_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*rli_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*ril_x;
      deijl_module_squared_dz_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*ril_y+2.0*(rij_z
*rjl_x-rij_x*rjl_z)*rli_x;

dCosDihedral_COMMON_INC;
  }  
  else if (&atom_k == &datom)
  {
dscalar_dx_ = rij_z*(rij_z*rjl_x-rij_x*rjl_z)+rji_y*(rij_x*rjl_y-rij_y*rjl_x);
      dscalar_dy_ = rji_z*(rij_y*rjl_z-rij_z*rjl_y)+rij_x*(rij_x*rjl_y-rij_y*
rjl_x);
      dscalar_dz_ = rij_y*(rij_y*rjl_z-rij_z*rjl_y)+rji_x*(rij_z*rjl_x-rij_x*
rjl_z);
      dejik_module_squared_dx_ = 2.0*(rji_z*rik_x-rji_x*rik_z)*rij_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rji_y;
      dejik_module_squared_dy_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rji_z+2.0*(rji_x
*rik_y-rji_y*rik_x)*rij_x;
      dejik_module_squared_dz_ = 2.0*(rji_y*rik_z-rji_z*rik_y)*rij_y+2.0*(rji_z
*rik_x-rji_x*rik_z)*rji_x;
      deijl_module_squared_dx_ = 0.0;
      deijl_module_squared_dy_ = 0.0;
      deijl_module_squared_dz_ = 0.0;

dCosDihedral_COMMON_INC;
  }  
  else if (&atom_l == &datom)
  {
dscalar_dx_ = (rji_z*rik_x-rji_x*rik_z)*rji_z+(rji_x*rik_y-rji_y*rik_x)*rij_y;
      dscalar_dy_ = (rji_y*rik_z-rji_z*rik_y)*rij_z+(rji_x*rik_y-rji_y*rik_x)*
rji_x;
      dscalar_dz_ = (rji_y*rik_z-rji_z*rik_y)*rji_y+(rji_z*rik_x-rji_x*rik_z)*
rij_x;
      dejik_module_squared_dx_ = 0.0;
      dejik_module_squared_dy_ = 0.0;
      dejik_module_squared_dz_ = 0.0;
      deijl_module_squared_dx_ = 2.0*(rij_z*rjl_x-rij_x*rjl_z)*rji_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*rij_y;
      deijl_module_squared_dy_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*rij_z+2.0*(rij_x
*rjl_y-rij_y*rjl_x)*rji_x;
      deijl_module_squared_dz_ = 2.0*(rij_y*rjl_z-rij_z*rjl_y)*rji_y+2.0*(rij_z
*rjl_x-rij_x*rjl_z)*rij_x;

dCosDihedral_COMMON_INC;
  }  
  else
    return 0.0;
}


} // namespace mdtk

#endif

