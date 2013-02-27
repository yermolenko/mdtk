/*
   Implementation of the many-body interatomic potential for copper,
   gold, silver and their alloys (header file).
   See [G.J. Ackland and V. Vitek, Phys. Rev. B 41, 10324 (1990)]

   Copyright (C) 2007, 2008, 2009, 2012, 2013 Oleksandr Yermolenko
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

#ifndef mdtk_ACKLAND_h
#define mdtk_ACKLAND_h

#include <cstdlib>
#include <cctype>

#include <mdtk/Atom.hpp>
#include <mdtk/consts.hpp>
#include <mdtk/potentials/manybody/FManybody.hpp>

#include <mdtk/Spline.hpp>
#include <mdtk/Spline5n.hpp>
//#include <mdtk/potentials/pairwise/FBM.hpp>

//#define  Ackland_HANDLE_SHORTRANGE
#define Ackland_OPTIMIZED

namespace mdtk
{

class Ackland : public FManybody
{
private:
  Float Phi(AtomsPair& ij);
  Float F(Atom &atom1);
  Float rho(Atom &atom1, const Float V = 0.0);
  Float g(AtomsPair& ij, const Float V = 0.0);
public:
  virtual Float operator()(AtomsArray&);
  void setupPotential();

  Ackland();
  Float getRcutoff() const {return max3(rk_[Au][Au][1],rk_[Ag][Ag][1],rk_[Cu][Cu][1]);} // needs unification!!!
private:
  enum {ECOUNT = 3};
  enum {Cu = 0};
  enum {Ag = 1};
  enum {Au = 2};

  Float rk_[ECOUNT][ECOUNT][7];
  Float ak_[ECOUNT][ECOUNT][7];
  Float Rk_[ECOUNT][ECOUNT][3];
  Float Ak_[ECOUNT][ECOUNT][3];

  Spline* splines[ECOUNT][ECOUNT];

  Float PhiCap(size_t a1_id, size_t a2_id, Float r) const;
  Float dPhiCap(size_t a1_id, size_t a2_id, Float r) const;

  size_t e2i(const Atom &atom) const
  {
    switch (atom.ID)
    {
      case Cu_EL : return Cu; break;
      case Ag_EL : return Ag; break;
      case Au_EL : return Au; break;
      default : throw Exception("e2i() : unknown element");
    };
  }
  Float rk(int i, const AtomsPair& ij) const
  {
    return rk_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }
  Float ak(int i, const AtomsPair& ij) const
  {
    return ak_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }
  Float Rk(int i, const AtomsPair& ij) const
  {
    return Rk_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }
  Float Ak(int i, const AtomsPair& ij) const
  {
    return Ak_[e2i(ij.atom1)][e2i(ij.atom2)][i];
  }

  Float Heaviside(Float x)
  {
    if (x <= 0.0) return 0.0;
    return 1.0;
  }
public:
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FManybody::SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FManybody::LoadFromStream(is,smode);
  }

  void fillR_concat_();
};

}

#endif


