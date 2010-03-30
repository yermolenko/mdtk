/*
   The proxy class for interatomic potentials (header file).

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

#ifndef mdtk_FProxy_hpp
#define mdtk_FProxy_hpp

#include <mdtk/potentials/FGeneral.hpp>

namespace mdtk
{

class FProxy
{
private:
public:
  std::vector<FGeneral*> potentials;
public:
  void addPotential(FGeneral* p)
  {
    potentials.push_back(p);
  }
public:
  Float operator()(AtomsContainer&);
  Vector3D grad(Atom &,AtomsContainer&);
  Float getRcutoff() const;
  FProxy();
  virtual
  ~FProxy()
  {
  }
  void freePotentials()
  {
    for(size_t i = 0; i < potentials.size(); i++)
      delete potentials[i];
  }

  void incDisplacement(Atom& atom, Vector3D inc)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      potentials[i]->incDisplacement(atom,inc);
  }  

  void NL_init(AtomsContainer& atoms)
  {
    for(size_t i = 0; i < potentials.size(); i++)
    {
      potentials[i]->NL_init(atoms);
    }
  };  

  void NL_checkRequestUpdate(AtomsContainer& atoms)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      potentials[i]->NL_checkRequestUpdate(atoms);
  }
  void NL_UpdateIfNeeded(AtomsContainer& atoms)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      potentials[i]->NL_UpdateIfNeeded(atoms);
  }

  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      potentials[i]->SaveToStream(os,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      potentials[i]->LoadFromStream(is,smode);
  }  
  bool hasNB(Atom& atom)
  {
    for(size_t i = 0; i < potentials.size(); i++)
      if (potentials[i]->NL(atom).size() > 0) return true;
    return false;
  }  
  
  Vector3D
  r_vec(Atom &atom1, Atom &atom2)
  {
    REQUIRE(potentials.size() > 0);
    return potentials[0]->r_vec_no_touch(atom1,atom2);
  }
  Float
  r_vec_module(Atom &atom1,Atom &atom2)
  {
    return r_vec(atom1,atom2).module();
  }   
  struct fakeNL
  {
    Float nlSkin() {TRACE("***CAUTION*** Called fake nlSkin()");return 0.0;}
  }nl;  

  void diagnose();
};

} 

#endif

