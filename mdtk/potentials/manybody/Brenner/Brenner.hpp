/*
   Implementation of the many-body interatomic potential for
   hydrocarbons (header file).
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009, 2012, 2013 Oleksandr
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
  Float VR(AtomsPair& ij, const Float V = 0.0);
  Float VR_Exp(AtomsPair& ij, const Float V = 0.0);

  Float VA(AtomsPair& ij, const Float V = 0.0);
  Float VA_Exp(AtomsPair& ij, const Float V = 0.0);
private:

  Float Baver(AtomsPair& ij, const Float V = 0.0);
  Float B(AtomsPair& ij, const Float V = 0.0);
  Float D(AtomsPair& ij, const Float V = 0.0);
  Float ExpTerm(AtomsPair& ij, AtomsPair& ik, const Float V = 0.0);
  Float G(AtomsPair& ij, AtomsPair& ik, const Float V = 0.0);
  Float dGdCT(Float CosT) const;

  Float Nt(AtomsPair& ij, const Float V = 0.0);
  Float NH(AtomsPair& ij, const Float V = 0.0);
  Float NC(AtomsPair& ij, const Float V = 0.0);

  Float Nconj(AtomsPair& ij, const Float V = 0.0);
  Float F(Float x) const;
  Float dF(Float x) const;

  Float HN_CC_num(Float a1, Float a2) const
    { return funcH_CC(a1,a2); }
  Float HN_CC_dH_num(Float a1, Float a2) const
    { return funcH_CC.dH(a1,a2); }
  Float HN_CC_dC_num(Float a1, Float a2) const
    { return funcH_CC.dC(a1,a2); }

  Float HN_CH_num(Float a1, Float a2) const
    { return funcH_CH(a1,a2); }
  Float HN_CH_dH_num(Float a1, Float a2) const
    { return funcH_CH.dH(a1,a2); }
  Float HN_CH_dC_num(Float a1, Float a2) const
    { return funcH_CH.dC(a1,a2); }

  Float HN(AtomsPair& ij, const Float V = 0.0);

  Float FN_num(Float a1, Float a2, Float a3) const
    { return funcF(a1,a2,a3); }
  Float FN_dNconj_num(Float a1, Float a2, Float a3) const
    { return funcF.dk(a1,a2,a3); }
  Float FN_dNt_i_num(Float a1, Float a2, Float a3) const
    { return funcF.di(a1,a2,a3); }
  Float FN_dNt_j_num(Float a1, Float a2, Float a3) const
    { return funcF.dj(a1,a2,a3); }

  Float FN(AtomsPair& ij, const Float V = 0.0);

public:
  enum ParamSet{POTENTIAL1,POTENTIAL2} paramSet;

  virtual Float operator()(AtomsArray&);

  Brenner(ParamSet parSet = POTENTIAL1);
//  virtual
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}
  bool probablyAreNeighbours(const Atom& atom1, const Atom& atom2) const
    {
      if (depos(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }
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

  size_t e2i(const Atom &atom) const
  {
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      default : throw Exception("e2i() : unknown element");
    };
  }

  Float Re(const AtomsPair& ij) const
  {
    return Re_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float De(const AtomsPair& ij) const
  {
    return De_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float beta(const AtomsPair& ij) const
  {
    return beta_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float S(const AtomsPair& ij) const
  {
    return S_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float delta(const AtomsPair& ij) const
  {
    return delta_[e2i(ij.atom1)][e2i(ij.atom2)];
  }
  Float delta(const Atom &atom) const
    {
      return delta_[e2i(atom)][e2i(atom)];
    }
//  virtual
  Float R(int i, const AtomsPair& ij) const
  {
    return R_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }
  Float R(int i, const Atom& atom1, const Atom& atom2) const
  {
    return R_[e2i(atom1)][e2i(atom2)][i];
  }
  Float alpha(const AtomsPair& ij, const AtomsPair& ik) const
  {
    return alpha_[e2i(ij.atom1)][e2i(ij.atom2)][e2i(ik.atom2)];
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
};

}

#endif
