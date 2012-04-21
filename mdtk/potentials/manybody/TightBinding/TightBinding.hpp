/*
   Implementation of the many-body interatomic potential for copper
   (header file).
   See [G. Betz, W. Husinsky, Nucl. Instr. and Meth. B 102, 281 (1995)]

   Copyright (C) 2006, 2007, 2008, 2009, 2012 Oleksandr Yermolenko
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

#ifndef mdtk_TIGHTBINDING_h
#define mdtk_TIGHTBINDING_h

#include <cstdlib>
#include <cctype>

#include <mdtk/Atom.hpp>
#include <mdtk/consts.hpp>
#include <mdtk/potentials/manybody/FManybody.hpp>
#include <mdtk/Spline.hpp>
#include <mdtk/potentials/pairwise/FBM.hpp>

//#define  EAM_HANDLE_SHORTRANGE
#define TightBinding_OPTIMIZED

namespace mdtk
{

class TightBinding : public FManybody
{
private:
  Float Phi(AtomsPair& ij);
  Float F(Atom &atom1);
  Float rho(Atom &atom1, const Float V = 0.0);
  Float g(AtomsPair& ij, const Float V = 0.0);
public:
  virtual Float operator()(AtomsArray&);
  void setupPotential();

  TightBinding();
  Float getRcutoff() const {return R_[1];}
private:
  Spline* spline;
  void fillR_concat_();
  Float BM_A;
  Float BM_B;

  Float alpha_;
  Float beta_;
  Float c_;
  Float Phi0_;
  Float R_[2];
  Float R(int i, const Atom &, const Atom &) const
  {
    return R_[i];
  }
  Float R(int i) const
  {
    return R_[i];
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

  bool probablyAreNeighbours(const Atom& atom1, const Atom& atom2) const
    {
      if (depos(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }
};

}

#endif


