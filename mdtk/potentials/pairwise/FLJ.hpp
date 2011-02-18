/*
   The Lennard-Jones interatomic potential implementation
   (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr
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

#ifndef mdtk_FLJ_hpp
#define mdtk_FLJ_hpp

#include "FPairwise.hpp"
#include <mdtk/Spline.hpp>

//#define  LJ_HANDLE_SHORTRANGE

namespace mdtk
{

class FLJ : public FPairwise
{
private:
public:  
  Float VLJ(Atom &a1,Atom &a2); 
  Vector3D dVLJ(Atom &a1,Atom &a2, Atom &da); 

  enum {ECOUNT = 3};
  enum {Cu = 0};
  enum {H = 1};
  enum {C = 2};

  Float sigma_[ECOUNT][ECOUNT];
  Float zeta_[ECOUNT][ECOUNT];
  Float sigma(Atom &atom1,Atom &atom2) const
  {
    return sigma_[e2i(atom1)][e2i(atom2)];
  }  
  Float zeta(Atom &atom1,Atom &atom2) const
  {
    return zeta_[e2i(atom1)][e2i(atom2)];
  }  
  size_t e2i(Atom &atom) const
  {
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      case Cu_EL : return Cu; break;
      default : throw Exception("e2i() : unknown element");
    };  
  }  

  Spline* splines[ECOUNT][ECOUNT];
  void fillR_concat_();

public:
  virtual Float operator()(AtomsContainer&);
  virtual Vector3D grad(Atom &a1,AtomsContainer&);
	FLJ(Rcutoff = Rcutoff());

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

