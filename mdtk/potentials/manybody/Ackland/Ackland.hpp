/*
   Implementation of the many-body interatomic potential for copper,
   gold, silver and their alloys (header file).
   See [G.J. Ackland and V. Vitek, Phys. Rev. B 41, 10324 (1990)]

   Copyright (C) 2007, 2008, 2009, 2012 Oleksandr Yermolenko
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
  Float Phi(Atom &atom1,Atom &atom2); 
    Vector3D dPhi(Atom &atom1,Atom &atom2, Atom &datom); 
  Float F(Atom &atom1); 
    Vector3D dF(Atom &atom1, Atom &datom); 
  Float rho(Atom &atom1); 
    Vector3D drho(Atom &atom1, Atom &datom); 
  Float g(Atom &atom1,Atom &atom2); 
    Vector3D dg(Atom &atom1,Atom &atom2, Atom &datom); 
public:
  virtual Float operator()(AtomsArray&);
  virtual Vector3D grad(Atom &,AtomsArray&);
  void setupPotential();

  Ackland();
  Float getRcutoff() const {return max3(rk_[Au][Au][1],rk_[Ag][Ag][1],rk_[Cu][Cu][1]);} // needs unification!!!
private:
/*
  Float alpha_;
  Float beta_;
  Float Phi0_;
*/
  Float c_;

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

  size_t e2i(Atom &atom) const
  {
    switch (atom.ID)
    {
      case Cu_EL : return Cu; break;
      case Ag_EL : return Ag; break;
      case Au_EL : return Au; break;
      default : throw Exception("e2i() : unknown element");
    };  
  }
  Float rk(int i,Atom &atom1,Atom &atom2) const
  {
    return rk_[e2i(atom1)][e2i(atom2)][i];
  }  
  Float ak(int i,Atom &atom1,Atom &atom2) const
  {
    return ak_[e2i(atom1)][e2i(atom2)][i];
  }  
  Float Rk(int i,Atom &atom1,Atom &atom2) const
  {
    return Rk_[e2i(atom1)][e2i(atom2)][i];
  }  
  Float Ak(int i,Atom &atom1,Atom &atom2) const
  {
    return Ak_[e2i(atom1)][e2i(atom2)][i];
  }  

  Float Heaviside(Float x)
  {
    if (x <= 0.0) return 0.0;
    return 1.0;
  }
/*
  Float R_[2];
  Float R(int i,Atom &,Atom &) const
  {
    return R_[i];
  }  
*/
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
/*
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
*/
  void fillR_concat_();
};

}

#endif


