/*
   Implementation of the many-body interatomic potential for
   hydrocarbons (header file).
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009 Oleksandr Yermolenko
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

#ifndef mdtk_brenner_hpp
#define mdtk_brenner_hpp

#include <cstdlib>
#include <cctype>

#include <mdtk/Atom.hpp>
#include <mdtk/consts.hpp>
#include <mdtk/potentials/manybody/FManybody.hpp>

#include "cubici.hpp"

//#define RETURN_BRENNER_0 {return 0;};
#define RETURN_BRENNER_0 ;

#define BRENNER_OPTIMIZED
#define BRENNER_OPTIMIZED_EVEN_BETTER

namespace mdtk
{


class Brenner : public FManybody
{
  friend class AIREBO;
public:
  Float VR(Atom &atom1,Atom &atom2); 
    Vector3D dVR(Atom &atom1,Atom &atom2, Atom &datom); 
  Float VR_Exp(Atom &atom1,Atom &atom2); 
    Vector3D dVR_Exp(Atom &atom1,Atom &atom2, Atom &datom); 

  Float VA(Atom &atom1,Atom &atom2); 
    Vector3D dVA(Atom &atom1,Atom &atom2, Atom &datom); 
  Float VA_Exp(Atom &atom1,Atom &atom2); 
    Vector3D dVA_Exp(Atom &atom1,Atom &atom2, Atom &datom); 
private:

  Float Baver(Atom &atom1,Atom &atom2); 
    Vector3D dBaver(Atom &atom1,Atom &atom2, Atom &datom); 
  Float B(Atom &atom1,Atom &atom2); 
    Vector3D dB(Atom &atom1,Atom &atom2, Atom &datom); 
  Float D(Atom &atom1,Atom &atom2); 
    Vector3D dD(Atom &atom1,Atom &atom2, Atom &datom); 
  Float ExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k); 
    Vector3D dExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom); 
  Float G(Atom &atom_i,Atom &atom_j,Atom &atom_k); 
    Vector3D dG(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom); 
      Float dGdCT(Float CosT) const;

  Float Nt(Atom &atom_i, Atom &atom_j); 
    Vector3D dNt(Atom &atom_i, Atom &atom_j, Atom &datom); 
  Float NH(Atom &atom_i, Atom &atom_j); 
    Vector3D dNH(Atom &atom_i, Atom &atom_j, Atom &datom); 
  Float NC(Atom &atom_i, Atom &atom_j); 
    Vector3D dNC(Atom &atom_i, Atom &atom_j, Atom &datom); 

  Float Nconj(Atom &atom1,Atom &atom2); 
    Vector3D dNconj(Atom &atom1,Atom &atom2, Atom &datom); 
  Float F(Float x) const; 
    Float dF(Float x) const; 

  Float HN_CC_num(Float a1, Float a2) const;
    Float HN_CC_dH_num(Float a1, Float a2) const;
    Float HN_CC_dC_num(Float a1, Float a2) const;
  Float HN_CH_num(Float a1, Float a2) const;
    Float HN_CH_dH_num(Float a1, Float a2) const;
    Float HN_CH_dC_num(Float a1, Float a2) const;

  Float HN(Atom &atom_i, Atom &atom_j);
    Vector3D dHN(Atom &atom_i, Atom &atom_j, Atom &datom);
  
  Float FN_num(Float a1, Float a2, Float a3) const; 
    Float FN_dNconj_num(Float a1, Float a2, Float a3) const;
    Float FN_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float FN_dNt_j_num(Float a1, Float a2, Float a3) const;

  Float FN(Atom &atom_i, Atom &atom_j);
    Vector3D dFN(Atom &atom_i, Atom &atom_j, Atom &datom); 
  
public:
  enum ParamSet{POTENTIAL1,POTENTIAL2} paramSet;  

  virtual Float operator()(AtomsContainer&);
  virtual Vector3D grad(Atom &,AtomsContainer&);

  Brenner(ParamSet parSet = POTENTIAL1);
//  virtual 
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}

private:
  enum {ECOUNT = 2};
  enum {C = 0};
  enum {H = 1};

  Float Re_[ECOUNT][ECOUNT];
  Float De_[ECOUNT][ECOUNT];
  Float beta_[ECOUNT][ECOUNT];
  Float S_[ECOUNT][ECOUNT];
  Float delta_[ECOUNT][ECOUNT];
  Float R_[ECOUNT][ECOUNT][2];

  Float alpha_[ECOUNT][ECOUNT][ECOUNT];
  
  Float ao_;
  Float co_;
  Float do_;

  Float G_HH_;

  void setupPotential1();
  void setupPotential2();
  
  size_t e2i(Atom &atom) const
  {
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      default : throw Exception("e2i() : unknown element");
    };  
  }  

  Float Re(Atom &atom1,Atom &atom2) const
  {
    return Re_[e2i(atom1)][e2i(atom2)];
  }  
  Float De(Atom &atom1,Atom &atom2) const
  {
    return De_[e2i(atom1)][e2i(atom2)];
  }  
  Float beta(Atom &atom1,Atom &atom2) const
  {
    return beta_[e2i(atom1)][e2i(atom2)];
  }  
  Float S(Atom &atom1,Atom &atom2) const
  {
    return S_[e2i(atom1)][e2i(atom2)];
  }  
  Float delta(Atom &atom1,Atom &atom2) const
  {
    return delta_[e2i(atom1)][e2i(atom2)];
  }  
    Float delta(Atom &atom) const
    {
      return delta_[e2i(atom)][e2i(atom)];
    }  
//  virtual
  Float R(int i,Atom &atom1,Atom &atom2) const
  {
    return R_[e2i(atom1)][e2i(atom2)][i];
  }  
  Float alpha(Atom &atom1,Atom &atom2,Atom &atom3) const
  {
    return alpha_[e2i(atom1)][e2i(atom2)][e2i(atom3)];
  }  

  FuncH_CC funcH_CC;
  FuncH_CH funcH_CH;

  FuncF funcF;
public:
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FManybody::SaveToStream(os,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FManybody::LoadFromStream(is,smode);
  }  

  Float  buildPairs(AtomsContainer& gl);

Float
f(Atom &atom1,Atom &atom2)
{
  Float r;
  if (ontouch_enabled)
  {
    r  = r_vec_module_no_touch(atom1,atom2);
  }
  else
  {
    r  = r_vec_module(atom1,atom2);
  }

  Float R1;
  Float R2;

  if (r<(R1=R(0,atom1,atom2)))
  {
    return 1.0;
  }
  else if (r>(R2=R(1,atom1,atom2)))
  {
    return 0.0;
  }
  else
  {
    if (ontouch_enabled) r_vec_touch_only(atom1,atom2);
    return (1.0+cos(M_PI*(r-R1)
            /(R2-R1)))/2.0;
  }  
}

Vector3D
df(Atom &atom1,Atom &atom2, Atom &datom)
{
  if (&datom != &atom1 && &datom != &atom2) return 0.0;

  Float r = r_vec_module(atom1,atom2);

  Float R1;
  Float R2;
  
  if (r<(R1=R(0,atom1,atom2)))
  {
    return 0.0;
  }
  else if (r>(R2=R(1,atom1,atom2)))
  {
    return 0.0;
  }
  else
  {
#ifdef FGENERAL_OPTIMIZED  
    Vector3D dvar = dr_vec_module(atom1,atom2,datom);
    if (dvar != 0.0)
      return (-(M_PI/(R2-R1))*sin(M_PI*(r-R1)/(R2-R1)))/2.0
             *dvar;
    else
      return 0.0;
#else
    return (-(M_PI/(R2-R1))*sin(M_PI*(r-R1)/(R2-R1)))/2.0
           *dr_vec_module(atom1,atom2,datom);
#endif
  }     
}


};

  


}

#endif


