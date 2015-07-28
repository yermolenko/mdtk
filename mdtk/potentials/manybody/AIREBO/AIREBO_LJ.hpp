/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The Lennard-Jones part (header file).
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2015
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

#ifndef mdtk_AIREBO_LJ_hpp
#define mdtk_AIREBO_LJ_hpp

#include <cstdlib>
#include <cctype>
#include <map>

#include <mdtk/potentials/manybody/FManybody.hpp>
#include <mdtk/potentials/manybody/AIREBO/REBO.hpp>
#include <mdtk/potentials/manybody/Brenner/Brenner.hpp>

#define AIREBO_OPTIMIZED
#define AIREBO_OPTIMIZED_EVEN_BETTER

//#define AIREBO_USING_BRENNER
#define AIREBO_USING_REBO

#if (!defined(AIREBO_USING_BRENNER) && !defined(AIREBO_USING_REBO))
#error "Please, define AIREBO_USING_"
#endif
#if (defined(AIREBO_USING_BRENNER) && defined(AIREBO_USING_REBO))
#error "Please, define only one AIREBO_USING_"
#endif

#ifdef  AIREBO_USING_BRENNER
#define CREBO Brenner
#endif

#ifdef  AIREBO_USING_REBO
#define CREBO REBO
#endif

namespace mdtk
{

class AIREBO : public FManybody
{
  CREBO& rebo;

  Float S(Float arg, Float arg_min, Float arg_max) const;
    Float dS(Float arg, Float arg_min, Float arg_max) const;

  Float StrStb(AtomsPair& ij, const Float V = 0.0);

  Float VLJ(AtomsPair& ij, const Float V = 0.0);

  Float Vtors(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V = 0.0);

  typedef std::pair<Float,std::vector<AtomsPair*> > CEelement;
  typedef std::map<size_t, CEelement> CArray;
  std::vector<CArray> CA;
  void  fill_Cij(AtomsArray& gl);
  void  cleanup_Cij();
  Float Cij(AtomsPair& ij, const Float V = 0.0);

  Float BijAsterix(AtomsPair& ij, const Float V = 0.0);
public:
  virtual Float operator()(AtomsArray& nl);

  AIREBO(CREBO* crebo);
//  virtual
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}
  bool probablyAreNeighbours(const Atom& atom1, const Atom& atom2) const
    {
      if (depos(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }
private:
  void setupPotential();

  enum {ECOUNT = CREBO::ECOUNT};
  enum {C = CREBO::C};
  enum {H = CREBO::H};


  Float sigma(const AtomsPair& ij) const
  {
    return sigma_[rebo.e2i(ij.atom1)][rebo.e2i(ij.atom2)];
  }

  Float R(int i, const Atom &atom1, const Atom &atom2) const
  {
    return R_[rebo.e2i(atom1)][rebo.e2i(atom2)][i];
  }

  Float RLJ(int i, const AtomsPair& ij) const
  {
    return RLJ_[rebo.e2i(ij.atom1)][rebo.e2i(ij.atom2)][i];
  }

  Float b(int i, const AtomsPair& ij) const
  {
    return b_[rebo.e2i(ij.atom1)][rebo.e2i(ij.atom2)][i];
  }

  Float zeta(const AtomsPair& ij) const
  {
    return zeta_[rebo.e2i(ij.atom1)][rebo.e2i(ij.atom2)];
  }
  Float R_[ECOUNT][ECOUNT][2];
  Float RLJ_[ECOUNT][ECOUNT][2];
  Float b_[ECOUNT][ECOUNT][2];
  Float sigma_[ECOUNT][ECOUNT];
  Float zeta_[ECOUNT][ECOUNT];
public:
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FManybody::SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FManybody::LoadFromStream(is,smode);
  }
  Float  buildPairs(AtomsArray& gl);
};

}

#endif

