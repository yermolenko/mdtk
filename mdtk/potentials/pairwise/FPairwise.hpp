/*
   The generalized pairwise interatomic potential class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Oleksandr
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
  Rcutoff rcutoff_;
  Float R(int i,Atom &,Atom &) const { return rcutoff_.R[i]; }  
public:
  Float getRcutoff() const {return rcutoff_.R[1];};
public:
  FPairwise(Rcutoff = Rcutoff());
  virtual void onTouch(Atom&) {}
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::SaveToStream(os,smode);
    rcutoff_.SaveToStream(os,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::LoadFromStream(is,smode);
    rcutoff_.LoadFromStream(is,smode);
  }  

Float
f(Atom &atom1,Atom &atom2)
{
  Float r = r_vec_module(atom1,atom2);

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
    if (R1 == R2) return 0.0;
    return (1.0+cos(M_PI*(r-R1)
            /(R2-R1)))/2.0;
  }  
}

Vector3D
df(Atom &atom1,Atom &atom2, Atom &datom)
{
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
    if (R1 == R2) return 0.0;
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

} // namespace mdtk

#endif

