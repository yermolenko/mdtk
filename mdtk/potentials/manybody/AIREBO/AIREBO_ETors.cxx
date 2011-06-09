/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The torsional potential part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011 Oleksandr Yermolenko
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

#include "AIREBO_ETors.hpp"
#include <algorithm>

namespace mdtk
{

  Float  ETors::buildPairs(AtomsContainer& gl)
  {

#ifdef VERBOSE_LOG
    std::cout << "Building Tors ... " << std::flush;
#endif
    Float Ei = 0;
    if (gl.size() != pairs.size()) pairs.resize(gl.size());    
    size_t ii;
    for(ii = 0; ii < gl.size(); ii++)
    {
      size_t prevSize = pairs[ii].size();
      pairs[ii].clear();
      pairs[ii].reserve(prevSize+FMANYBODY_PAIRS_RESERVE_ADD);
    }  
    for(ii = 0; ii < gl.size(); ii++)
    {
      Atom &atom_i = *(gl[ii]);
      if (isHandled(atom_i))
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        std::pair<int,int> sample_pair(atom_i.globalIndex,atom_j.globalIndex);
        if (&atom_i != &atom_j)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
        {
    currentPairPtr = &sample_pair;
    ontouch_enabled = true;

      {
  for(Index k = 0; k < EREBO::NL(atom_i).size(); k++)
  {
    Atom &atom_k = *(EREBO::NL(atom_i)[k]);
    if (&atom_k != &atom_i && &atom_k != &atom_j)
    if (r_vec_module_no_touch(atom_k,atom_i) < EREBO::R(1,atom_k,atom_i))
    for(Index l = 0; l < EREBO::NL(atom_j).size(); l++)
    {
      Atom &atom_l = *(EREBO::NL(atom_j)[l]);
      if (&atom_l != &atom_i && &atom_l != &atom_j &&  &atom_l != &atom_k )
      {
        Float wki  = EREBO::f(atom_k,atom_i);
        Float wij  = EREBO::f(atom_i,atom_j);
        if (r_vec_module_no_touch(atom_j,atom_l) < EREBO::R(1,atom_j,atom_l))
        {
          Float wjl  = EREBO::f(atom_j,atom_l);
          Ei += wki*wij*wjl*Vtors(atom_i,atom_j,atom_k,atom_l);
        }  
      }  
    }
  }    

      }  

    currentPairPtr = NULL;
    ontouch_enabled = false;
        }  
      }  
    }  

#ifdef VERBOSE_LOG
    std::cout << "done." << std::endl;
#endif
    return Ei;
  }  


Float
ETors::Vtors(Atom &ai,Atom &aj,Atom &ak,Atom &al)
{
  if (ai.ID != C_EL || aj.ID != C_EL) return 0.0;
  Float CosDh = - CosDihedral(ai,aj,ak,al);
  return (256.0/405.0)*zetaCC(ak,al)*pow(0.5*(1.0+CosDh),5.0)-0.1*zetaCC(ak,al);
}

Vector3D
ETors::dVtors(Atom &ai,Atom &aj,Atom &ak,Atom &al, Atom &da)
{
  if (ai.ID != C_EL || aj.ID != C_EL) return 0.0;

  Float CosDh = - CosDihedral(ai,aj,ak,al);
  Vector3D dCosDh = - dCosDihedral(ai,aj,ak,al,da);

  return (256.0/405.0)*zetaCC(ak,al)*pow(0.5*(1.0+CosDh),5.0-1.0)*0.5*dCosDh;
}  
 

Float
ETors::operator()(AtomsContainer& gl)
{
  Float Ei = 0.0;
  Ei += ETor(gl);
  return Ei;
}  

Vector3D
ETors::grad(Atom &atom,AtomsContainer &gl)
{
  Vector3D dEi = 0.0;
  dEi += dETor(atom, gl);
  return  dEi;
}  

Float
ETors::ETor(AtomsContainer& gl)
{
  return buildPairs(gl);
}

Vector3D
ETors::dETor(Atom &atom,AtomsContainer &gl)
{
  Vector3D dEi = 0.0;

  Index i;

  if (isHandled(atom))
  {
    std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];
    
    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = *(gl[acnt[i].first]);
      {
        Atom &atom_j = *(gl[acnt[i].second]);

        REQUIREM(&atom_j != &atom_i,"must be (&atom_j != &atom_i)");
        {


  for(Index k = 0; k < EREBO::NL(atom_i).size(); k++)
  {
    Atom &atom_k = *(EREBO::NL(atom_i)[k]);
    if (&atom_k != &atom_i && &atom_k != &atom_j)
    if (r_vec_module_no_touch(atom_k,atom_i) < EREBO::R(1,atom_k,atom_i))
    for(Index l = 0; l < EREBO::NL(atom_j).size(); l++)
    {
      Atom &atom_l = *(EREBO::NL(atom_j)[l]);
      if (&atom_l != &atom_i && &atom_l != &atom_j &&  &atom_l != &atom_k )
      {
        Float wki     = EREBO::f(atom_k,atom_i);
        Vector3D dwki = EREBO::df(atom_k,atom_i,atom);
        Float wij     = EREBO::f(atom_i,atom_j);
        Vector3D dwij = EREBO::df(atom_i,atom_j,atom);

        if (r_vec_module_no_touch(atom_j,atom_l) < EREBO::R(1,atom_j,atom_l))
        {
          Float wjl     = EREBO::f(atom_j,atom_l);
          Vector3D dwjl = EREBO::df(atom_j,atom_l,atom);

          Float Vtors_val = Vtors(atom_i,atom_j,atom_k,atom_l);
          Vector3D dVtors_val = dVtors(atom_i,atom_j,atom_k,atom_l,atom); 

            Vector3D dEij;
            dEij = dwki* wij* wjl* Vtors_val+
                    wki*dwij* wjl* Vtors_val+
                    wki* wij*dwjl* Vtors_val+
                    wki* wij* wjl*dVtors_val;
            dEi += dEij;
        }  
      }  
    }
  }    
 
        }  
          
      }
    }    
    
  }  

  return  dEi;
}  

ETors::ETors():
  EREBO(REBO::POTENTIAL1)
{
  setupPotential();
  nl.Rcutoff = getRcutoff();
}

void
ETors::setupPotential()
{
  PTRACE("Setup ETORS");

  zetaCC_[C][C] = 0.3079*eV;
  zetaCC_[H][H] = 0.1250*eV;
  zetaCC_[C][H] = 0.1787*eV;
    zetaCC_[H][C] = zetaCC_[C][H];
}

}


