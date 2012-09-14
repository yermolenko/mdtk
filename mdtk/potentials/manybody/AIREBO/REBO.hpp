/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The REBO potential part (header file).
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

  Float fprime(AtomsPair& ij, const Float V = 0.0);

  Float VR(AtomsPair& ij, const Float V = 0.0);
  Float VR_Exp(AtomsPair& ij, const Float V = 0.0);

  Float VA(AtomsPair& ij, const Float V = 0.0);
protected:
  Float Baver(AtomsPair& ij, const Float V = 0.0);
  Float p_sigma_pi(AtomsPair& ij, const Float V = 0.0);
  Float D(AtomsPair& ij, const Float V = 0.0);
  Float ExpTerm(AtomsPair& ij, AtomsPair& ik, const Float V = 0.0);
  Float G(AtomsPair& ij, AtomsPair& ik, const Float V = 0.0);
  Float dGdCT(Float CosT) const;

  Float Nt(AtomsPair& ij, const Float V = 0.0);
  void Nt_donly(AtomsPair& ij, const Float V = 0.0);
  Float NH(AtomsPair& ij, const Float V = 0.0);
  void NH_donly(AtomsPair& ij, const Float V = 0.0);
  Float NC(AtomsPair& ij, const Float V = 0.0);
  void NC_donly(AtomsPair& ij, const Float V = 0.0);

  Float NconjSum1(AtomsPair& ij, const Float V = 0.0);
  Float NconjSum2(AtomsPair& ij, const Float V = 0.0);
  Float Nconj(AtomsPair& ij, const Float V = 0.0);
  Float F(Float x) const;
  Float dF(Float x) const;

  Float P_CC_num(Float a1, Float a2) const;
    Float P_CC_dH_num(Float a1, Float a2) const;
    Float P_CC_dC_num(Float a1, Float a2) const;
  Float P_CH_num(Float a1, Float a2) const;
    Float P_CH_dH_num(Float a1, Float a2) const;
    Float P_CH_dC_num(Float a1, Float a2) const;

  Float P(AtomsPair& ij, const Float V = 0.0);

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

  Float pi_rc(AtomsPair& ij, const Float V = 0.0);
public:
  Float pi_dh_Tij(AtomsPair& ij, const Float V = 0.0);
  Float pi_dh_sum(AtomsPair& ij, const Float V = 0.0);
  Float pi_dh(AtomsPair& ij, const Float V = 0.0);

  Float Tij_num(Float a1, Float a2, Float a3) const; 
    Float Tij_dNconj_num(Float a1, Float a2, Float a3) const;
    Float Tij_dNt_i_num(Float a1, Float a2, Float a3) const;
    Float Tij_dNt_j_num(Float a1, Float a2, Float a3) const;

public:
  enum ParamSet{POTENTIAL1,POTENTIAL2} /*paramSet*/;  

  virtual Float operator()(AtomsArray&);

  REBO(ParamSet /*parSet*/ = POTENTIAL1);
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}
  bool probablyAreNeighbours(const Atom& atom1, const Atom& atom2) const
    {
      if (depos(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
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

  size_t e2i(const Atom &atom) const
  {
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      default : throw Exception("e2i() : unknown element");
    };
  }

  Float Q(const AtomsPair& ij) const
  {
    return Q_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float alpha(const AtomsPair& ij) const
  {
    return alpha_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float A(const AtomsPair& ij) const
  {
    return A_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float B(int i, const AtomsPair& ij) const
  {
    return B_[i-1][e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float beta(int i, const AtomsPair& ij) const
  {
    return beta_[i-1][e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float rho(const AtomsPair& ij) const
  {
    if (rho_[e2i(ij.atom1)][e2i(ij.atom2)] == 0.0) throw Exception("REBO rho. CC");
    return rho_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float R(int i, const Atom &atom1, const Atom &atom2) const
  {
    return R_[e2i(atom1)][e2i(atom2)][i];
  }
  Float R(int i, const AtomsPair& ij) const
  {
    return R_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }
  Float Rprime(int i, const AtomsPair& ij) const
  {
    return Rprime_[e2i(ij.atom1)][e2i(ij.atom2)][i];
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

  std::vector<Float> Nts;
  std::vector<std::vector<Atom*> > Nts2touch;
  std::vector<Float> NCs;
  std::vector<std::vector<Atom*> > NCs2touch;
  std::vector<Float> NHs;
  std::vector<std::vector<Atom*> > NHs2touch;
  void dfThem(const std::vector<std::vector<Atom*> >& patoms, AtomsPair& ij, const Float V);
  bool areNeighbours(Atom &atom_i, Atom &atom_j);
  void countNeighbours(AtomsArray& gl);
};

}

#endif
