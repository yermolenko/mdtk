/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The torsional potential part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011, 2012, 2013
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

#include "AIREBO_ETors.hpp"
#include <algorithm>

namespace mdtk
{

Float
ETors::operator()(AtomsArray& gl)
{
  Float Ei = 0;
  for(size_t ii = 0; ii < gl.size(); ii++)
  {
    Atom &atom_i = gl[ii];
    if (isHandled(atom_i))
    {
      AtomRefsContainer& nli = NL(atom_i);
      for(size_t jj = 0; jj < nli.size(); jj++)
      {
        Atom &atom_j = *(nli[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        if (&atom_i != &atom_j)
        {
          if (!probablyAreNeighbours(atom_i,atom_j)) continue;

          AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));

          for(size_t k = 0; k < nli.size(); k++)
          {
            Atom &atom_k = *(nli[k]);
            if (&atom_k != &atom_i && &atom_k != &atom_j)
            {
              if (!probablyAreNeighbours(atom_i,atom_k)) continue;
              AtomsPair ki(atom_k,atom_i,R(0,atom_k,atom_i),R(1,atom_k,atom_i));
              AtomRefsContainer& nlj = NL(atom_j);
              for(size_t l = 0; l < nlj.size(); l++)
              {
                Atom &atom_l = *(nlj[l]);
                if (&atom_l != &atom_i && &atom_l != &atom_j &&  &atom_l != &atom_k )
                {
                  if (!probablyAreNeighbours(atom_j,atom_l)) continue;
                  AtomsPair jl(atom_j,atom_l,R(0,atom_j,atom_l),R(1,atom_j,atom_l));

                  if (fabs(SinTheta(ij,ki))<0.1) continue;
                  if (fabs(SinTheta(ij,jl))<0.1) continue;

                  Float V = 1.0;

                  AtomsPair ik(-ki);

                  Float wki  = ki.f();
                  Float wij  = ij.f();
                  Float wjl  = jl.f();
                  Float VtorsVal = Vtors(ij,ik,jl,wki*wij*wjl*V);

                  if (V != 0)
                  {
                    ki.f(wij*wjl*VtorsVal*V);
                    ij.f(wki*wjl*VtorsVal*V);
                    jl.f(wki*wij*VtorsVal*V);
                  }

                  Ei += wki*wij*wjl*VtorsVal;
                }
              }
            }
          }
        }
      }
    }
  }

  return Ei;
}

Float
ETors::Vtors(AtomsPair& ij, AtomsPair& ik, AtomsPair& jl, const Float V)
{
  if (ij.atom1.ID != C_EL || ij.atom2.ID != C_EL) return 0.0;
  Float CosDh = - CosDihedral(ij,ik,jl);

  Float Val = (256.0/405.0)*zetaCC(ik.atom2,jl.atom2)*pow(0.5*(1.0+CosDh),5.0)-0.1*zetaCC(ik.atom2,jl.atom2);

  if (Val != 0)
  {
    Float Der = (256.0/405.0)*zetaCC(ik.atom2,jl.atom2)*pow(0.5*(1.0+CosDh),5.0-1.0)*0.5;
    CosDihedral(ij,ik,jl,-1.0*Der*V);
  }

  return Val;
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
  zetaCC_[C][C] = 0.3079*eV;
  zetaCC_[H][H] = 0.1250*eV;
  zetaCC_[C][H] = 0.1787*eV;
    zetaCC_[H][C] = zetaCC_[C][H];

  PRINT("AIREBO::ETors interatomic potential configured.\n");
}

}


