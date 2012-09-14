/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The torsional potential part (header file).
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

#ifndef mdtk_AIREBO_ETors_hpp
#define mdtk_AIREBO_ETors_hpp

#include <cstdlib>
#include <cctype>

#include <mdtk/potentials/manybody/FManybody.hpp>
#include <mdtk/potentials/manybody/AIREBO/REBO.hpp>

#define ETORS_OPTIMIZED
#define ETORS_OPTIMIZED_EVEN_BETTER

#define EREBO REBO

namespace mdtk
{

class ETors : public EREBO
{
public:
  virtual Float operator()(AtomsArray& nl);

  Float Vtors(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V);

  ETors();
//  bool probablyAreNeighbours(Atom& atom1, Atom& atom2) const; not needed because of the same R[][][] values
private:
  void setupPotential();

  enum {ECOUNT = EREBO::ECOUNT};
  enum {C = EREBO::C};
  enum {H = EREBO::H};

  Float zetaCC(const Atom& atom1, const Atom& atom2) const
  {
    return zetaCC_[EREBO::e2i(atom1)][EREBO::e2i(atom2)];
  }

  Float zetaCC_[ECOUNT][ECOUNT];
public:
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    EREBO::SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    EREBO::LoadFromStream(is,smode);
  }
};

}

#endif
