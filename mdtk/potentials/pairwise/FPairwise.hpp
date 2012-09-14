/*
   The generalized pairwise interatomic potential class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2012 Oleksandr
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

#ifndef mdtk_FPairwise_hpp
#define mdtk_FPairwise_hpp

#include <mdtk/potentials/FGeneral.hpp>

namespace mdtk
{

struct Rcutoff
{
  Float R[2];
  explicit Rcutoff(Float R0 = 3.1*Ao, Float R1 = 3.5*Ao)
  {
    REQUIRE(!(R0>R1));
    R[0] = R0;
    R[1] = R1;
  }
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    YAATK_FSTREAM_WRITE(os,R[0],smode);
    YAATK_FSTREAM_WRITE(os,R[1],smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    YAATK_FSTREAM_READ(is,R[0],smode);
    YAATK_FSTREAM_READ(is,R[1],smode);
  }
};

class FPairwise : public FGeneral
{
protected:
  Rcutoff rc;
  Float R(int i, const Atom &, const Atom &) const { return rc.R[i]; }
  Float R(int i) const { return rc.R[i]; }
public:
  Float getRcutoff() const {return rc.R[1];};
public:
  FPairwise(Rcutoff = Rcutoff());
  virtual void onTouch(const Atom&) {}
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::SaveToStream(os,smode);
    rc.SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::LoadFromStream(is,smode);
    rc.LoadFromStream(is,smode);
  }
  bool probablyAreNeighbours(const Atom& atom1, const Atom& atom2) const
    {
      if (depos(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }
};

} // namespace mdtk

#endif
