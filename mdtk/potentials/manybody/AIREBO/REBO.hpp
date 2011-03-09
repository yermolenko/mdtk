/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The REBO potential part (header file).
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

#ifndef mdtk_REBO_h
#define mdtk_REBO_h

#include <cstdlib>
#include <cctype>

#include <mdtk/Atom.hpp>
#include <mdtk/consts.hpp>
#include <mdtk/potentials/manybody/FManybody.hpp>

#include "cubici.hpp"

#define REBO_DIHEDRAL

//#define RETURN_REBO_0 {return 0;};
#define RETURN_REBO_0 ;

#define REBO_OPTIMIZED
#define REBO_OPTIMIZED_EVEN_BETTER

namespace mdtk
{


class REBO : public FManybody
{
  friend class AIREBO;
public:
  Float Sprime(Float arg, Float arg_min, Float arg_max) const;
  Float dSprime(Float arg, Float arg_min, Float arg_max) const;

  Float fprime(Atom &atom1,Atom &atom2);
  Vector3D dfprime(Atom &atom1,Atom &atom2, Atom &datom);

  Float VR(Atom &atom1,Atom &atom2); 
    Vector3D dVR(Atom &atom1,Atom &atom2, Atom &datom); 
  Float VR_Exp(Atom &atom1,Atom &atom2); 
    Vector3D dVR_Exp(Atom &atom1,Atom &atom2, Atom &datom); 

  Float VA(Atom &atom1,Atom &atom2); 
    Vector3D dVA(Atom &atom1,Atom &atom2, Atom &datom); 
protected:
  Float Baver(Atom &atom1,Atom &atom2); 
    Vector3D dBaver(Atom &atom1,Atom &atom2, Atom &datom); 
  Float p_sigma_pi(Atom &atom1,Atom &atom2); 
    Vector3D dp_sigma_pi(Atom &atom1,Atom &atom2, Atom &datom); 
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

  Float P_CC_num(Float a1, Float a2) const;
    Float P_CC_dH_num(Float a1, Float a2) const;
    Float P_CC_dC_num(Float a1, Float a2) const;
  Float P_CH_num(Float a1, Float a2) const;
    Float P_CH_dH_num(Float a1, Float a2) const;
    Float P_CH_dC_num(Float a1, Float a2) const;

  Float P(Atom &atom_i, Atom &atom_j);
    Vector3D dP(Atom &atom_i, Atom &atom_j, Atom &datom);
  
  Float pi_rc_CC_num(Float a1, Float a2, Float a3) const; 
    Float pi_rc_CC_dNconj_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_CC_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_CC_dNt_j_num(Float a1, Float a2, Float a3) const;

  Float pi_rc_CH_num(Float a1, Float a2, Float a3) const; 
    Float pi_rc_CH_dNconj_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_CH_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_CH_dNt_j_num(Float a1, Float a2, Float a3) const;

  Float pi_rc_HH_num(Float a1, Float a2, Float a3) const; 
    Float pi_rc_HH_dNconj_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_HH_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float pi_rc_HH_dNt_j_num(Float a1, Float a2, Float a3) const;

  Float pi_rc(Atom &atom_i, Atom &atom_j);
    Vector3D dpi_rc(Atom &atom_i, Atom &atom_j, Atom &datom); 
public:  
  Float pi_dh(Atom &atom_i, Atom &atom_j);
    Vector3D dpi_dh(Atom &atom_i, Atom &atom_j, Atom &datom); 

  Float Tij_num(Float a1, Float a2, Float a3) const; 
    Float Tij_dNconj_num(Float a1, Float a2, Float a3) const;
    Float Tij_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float Tij_dNt_j_num(Float a1, Float a2, Float a3) const;

public:
  enum ParamSet{POTENTIAL1,POTENTIAL2} /*paramSet*/;  

  virtual Float operator()(AtomsContainer&);
  virtual Vector3D grad(Atom &,AtomsContainer&);

  REBO(ParamSet /*parSet*/ = POTENTIAL1);
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}
  bool probablyAreNeighbours(Atom& atom1, Atom& atom2)
    {
      if (r_vec_no_touch(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }

protected:
  enum {ECOUNT = 2};
  enum {C = 0};
  enum {H = 1};

  Float Q_[ECOUNT][ECOUNT];
  Float alpha_[ECOUNT][ECOUNT];
  Float A_[ECOUNT][ECOUNT];
  Float B_[3][ECOUNT][ECOUNT];
  Float beta_[3][ECOUNT][ECOUNT];
  Float rho_[ECOUNT][ECOUNT];

  Float R_[ECOUNT][ECOUNT][2];
  Float Rprime_[ECOUNT][ECOUNT][2];
  Float Ntot_[2];

  Float ao_;
  Float co_;
  Float do_;

  Float G_HH_;

  void setupPotential1();
  
  size_t e2i(Atom &atom) const
  {
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      default : throw Exception("e2i() : unknown element");
    };  
  }  

  Float Q(Atom &atom1,Atom &atom2) const
  {
    return Q_[e2i(atom1)][e2i(atom2)];
  }  
  Float alpha(Atom &atom1,Atom &atom2) const
  {
    return alpha_[e2i(atom1)][e2i(atom2)];
  }  
  Float A(Atom &atom1,Atom &atom2) const
  {
    return A_[e2i(atom1)][e2i(atom2)];
  }  
  Float B(int i, Atom &atom1,Atom &atom2) const
  {
    return B_[i-1][e2i(atom1)][e2i(atom2)];
  }  
  Float beta(int i, Atom &atom1,Atom &atom2) const
  {
    return beta_[i-1][e2i(atom1)][e2i(atom2)];
  }  
  Float rho(Atom &atom1,Atom &atom2) const
  {
    if (rho_[e2i(atom1)][e2i(atom2)] == 0.0) throw Exception("REBO rho. CC");
    return rho_[e2i(atom1)][e2i(atom2)];
  }  
//  virtual
  Float R(int i,Atom &atom1,Atom &atom2) const
  {
    return R_[e2i(atom1)][e2i(atom2)][i];
  }  
  Float Rprime(int i,Atom &atom1,Atom &atom2) const
  {
    return Rprime_[e2i(atom1)][e2i(atom2)][i];
  }  

  FuncP_CC funcP_CC;
  FuncP_CH funcP_CH;

  Func_pi_rc_CC func_pi_rc_CC;
  Func_pi_rc_CH func_pi_rc_CH;
  Func_pi_rc_HH func_pi_rc_HH;
  
  FuncG_C1 funcG_C1; 
  FuncG_C2 funcG_C2; 
  FuncG_H  funcG_H; 

  Func_Tij func_Tij;

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


