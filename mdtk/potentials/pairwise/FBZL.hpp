/*
   The Ziegler-Biersack-Littmark interatomic potential implementation
   (header file).

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

#ifndef mdtk_FBZL_hpp
#define mdtk_FBZL_hpp

#include "FPairwise.hpp"

namespace mdtk
{

class FBZL : public FPairwise
{
private:
  Float AB_;
public:
  Float F11(AtomsPair& ij, const Float V = 0.0);
public:
  virtual Float operator()(AtomsArray&);
  FBZL(Rcutoff = Rcutoff());

  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FPairwise::SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FPairwise::LoadFromStream(is,smode);
  }
};

}

#endif
