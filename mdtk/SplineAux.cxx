/*
   Auxillary N-cubic spline classes and functions.

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

#include "SplineAux.hpp"

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt;
using std::pow;


#include <gsl/gsl_linalg.h>

void
solveLinSys(gsl_matrix* a, gsl_vector* b, gsl_vector* x)
{
  gsl_permutation* p = gsl_permutation_alloc(x->size);
  int signum;

  gsl_linalg_LU_decomp(a, p, &signum);    
  gsl_linalg_LU_solve(a, p, b, x);

  gsl_permutation_free(p);
}

void
addEquation(SplineMultiD& s, std::vector<Float>& xv, Float& fv,
   std::vector<std::vector<Float> >& xvs, std::vector<Float>& fvs)
{
  xvs.push_back(std::vector<Float>(s.terms.size()));
  fvs.push_back(fv);

  size_t j;
  for(j = 0; j < s.terms.size(); j++)
  {
    xvs[xvs.size()-1][j] = s.terms[j].woC(xv);
  }
  fvs[fvs.size()-1] = fv;
}

void
SplineMultiD::addEquations(std::vector<std::vector<Float> >& xvs, std::vector<Float>& fvs)
{
  REQUIRE(xvs.size() == fvs.size());
  REQUIRE(xvs[0].size() == terms.size());
  REQUIRE(xvs.size() == terms.size());

  gsl_matrix* va = gsl_matrix_alloc(terms.size(), terms.size());
  gsl_vector* vb = gsl_vector_alloc(terms.size());
  gsl_vector* vx = gsl_vector_alloc(terms.size());

  size_t j, eqIndex;
  for(eqIndex = 0; eqIndex < terms.size(); eqIndex++)
  {
    for(j = 0; j < terms.size(); j++)
    {
      gsl_matrix_set(va, eqIndex, j, xvs[eqIndex][j]);
    }
    gsl_vector_set(vb, eqIndex, fvs[eqIndex]);
  }

  solveLinSys(va,vb,vx);

  size_t i;
  for(i = 0; i < terms.size() ; i++)
    terms[i].c = gsl_vector_get(vx, i);

  gsl_matrix_free(va);
  gsl_vector_free(vb);
  gsl_vector_free(vx);
}


}

