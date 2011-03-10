/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The Lennard-Jones part (header file).
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr
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

#ifndef mdtk_AIREBO_LJ_hpp
#define mdtk_AIREBO_LJ_hpp

#include <cstdlib>
#include <cctype>

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

class AIREBO : public CREBO
{
  CREBO& rebo;

  Float S(Float arg, Float arg_min, Float arg_max) const;
    Float dS(Float arg, Float arg_min, Float arg_max) const;

  Float VLJ(Atom &a1,Atom &a2); 
    Vector3D dVLJ(Atom &a1,Atom &a2, Atom &da); 

  Float Vtors(Atom &ai,Atom &aj,Atom &ak,Atom &al); 
    Vector3D dVtors(Atom &ai,Atom &aj,Atom &ak,Atom &al, Atom &da); 

  Float Cij(Atom &atom1,Atom &atom2); 
    Vector3D dCij(Atom &atom1,Atom &atom2, Atom &datom); 

  Float BijAsterix(Atom &atom1,Atom &atom2); 
    Vector3D dBijAsterix(Atom &atom1,Atom &atom2, Atom &datom); 

public:
  virtual Float operator()(AtomsContainer& nl);
  virtual Vector3D grad(Atom &atom,AtomsContainer&);

  Float ELJ(AtomsContainer&);
  Vector3D dELJ(Atom &,AtomsContainer&);

  AIREBO(CREBO* crebo);
//  virtual
  Float getRcutoff() const {return      max3(R_[C][C][1],R_[C][H][1],R_[H][H][1]);}
  bool probablyAreNeighbours(Atom& atom1, Atom& atom2)
    {
      if (r_vec_no_touch(atom1,atom2).module_squared() > SQR(R(1,atom1,atom2)))
        return false;

      return true;
    }

private:
  void setupPotential();

  enum {ECOUNT = CREBO::ECOUNT};
  enum {C = CREBO::C};
  enum {H = CREBO::H};


  Float sigma(Atom &atom1,Atom &atom2) const
  {
    return sigma_[CREBO::e2i(atom1)][CREBO::e2i(atom2)];
  }  

//  virtual
  Float R(int i,Atom &atom1,Atom &atom2) const
  {
    return R_[CREBO::e2i(atom1)][CREBO::e2i(atom2)][i];
  }  

  Float RLJ(int i,Atom &atom1,Atom &atom2) const
  {
    return RLJ_[CREBO::e2i(atom1)][CREBO::e2i(atom2)][i];
  }  

  Float b(int i,Atom &atom1,Atom &atom2) const
  {
    return b_[CREBO::e2i(atom1)][CREBO::e2i(atom2)][i];
  }  

  Float zeta(Atom &atom1,Atom &atom2) const
  {
    return zeta_[CREBO::e2i(atom1)][CREBO::e2i(atom2)];
  }  
  Float R_[ECOUNT][ECOUNT][2];
  Float RLJ_[ECOUNT][ECOUNT][2];
  Float b_[ECOUNT][ECOUNT][2];
  Float sigma_[ECOUNT][ECOUNT];
  Float zeta_[ECOUNT][ECOUNT];
public:
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    CREBO::SaveToStream(os,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    CREBO::LoadFromStream(is,smode);
  }  
  Float  buildPairs(AtomsContainer& gl);

  NeighbourList nl;

  AtomsContainer& NL(Atom& atom)
  {
    return nl.nl[atom.globalIndex];
  }  
  AtomsContainer& NL_with_self(Atom& atom)
  {
    return nl.nl_with_self[atom.globalIndex];
  }  

  void NL_checkRequestUpdate(AtomsContainer& atoms)
  {
    CREBO::NL_checkRequestUpdate(atoms);
    nl.checkRequestUpdate(atoms);
  }
  void NL_UpdateIfNeeded(AtomsContainer& atoms)
  {
    CREBO::NL_UpdateIfNeeded(atoms);
    nl.UpdateIfNeeded(atoms);
  }
  void NL_init(AtomsContainer& atoms)
  {
    CREBO::NL_init(atoms);
    nl.init(atoms);
  }
  void incDisplacement(Atom& atom, Vector3D inc)
  {
    CREBO::incDisplacement(atom,inc);
    nl.displacements[atom.globalIndex] += inc;
  }  


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

