/*
   Implementation of the many-body interatomic potential for copper
   (header file).
   See [G. Betz, W. Husinsky, Nucl. Instr. and Meth. B 102, 281 (1995)]

   Copyright (C) 2006, 2007, 2008, 2009 Oleksandr Yermolenko
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
  Float Phi(Atom &atom1,Atom &atom2); 
    Vector3D dPhi(Atom &atom1,Atom &atom2, Atom &datom); 
  Float F(Atom &atom1); 
    Vector3D dF(Atom &atom1, Atom &datom); 
  Float rho(Atom &atom1); 
    Vector3D drho(Atom &atom1, Atom &datom); 
  Float g(Atom &atom1,Atom &atom2); 
    Vector3D dg(Atom &atom1,Atom &atom2, Atom &datom); 
public:
  virtual Float operator()(AtomsContainer&);
  virtual Vector3D grad(Atom &,AtomsContainer&);
  void setupPotential();

  TightBinding();
  Float getRcutoff() const {return R_[1];}
private:
  Spline* spline;
  void fillR_concat_();

  Float alpha_;
  Float beta_;
  Float c_;
  Float Phi0_;
  Float R_[2];
  Float R(int i,Atom &,Atom &) const
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

  Float  buildPairs(AtomsContainer& gl);

Float
f(Atom &atom1,Atom &atom2)
{
  Float r;
  if (ontouch_enabled)
  {
    r  = r_vec_module_no_touch(atom1,atom2);
  }
  else
  {
    r  = r_vec_module(atom1,atom2);
  }

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
    if (ontouch_enabled) r_vec_touch_only(atom1,atom2);
    return (1.0+cos(M_PI*(r-R1)
            /(R2-R1)))/2.0;
  }  
}

Vector3D
df(Atom &atom1,Atom &atom2, Atom &datom)
{
  if (&datom != &atom1 && &datom != &atom2) return 0.0;

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

}

#endif


