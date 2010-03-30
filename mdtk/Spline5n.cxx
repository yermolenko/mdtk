/*
   The fifth-order spline class.

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

#include "Spline5n.hpp"

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt;
using std::pow;


void
Spline5n::init(
         Float       x[2],
         Float       v[2],
         Float    dvdx[2],
         Float    d2vdxdx[2]
)
{
  size_t N = 1;

  SplineMultiD Fv(N,5);

  SplineMultiD Fdvdx(N,5);

  Fdvdx.differentiate(0);

  SplineMultiD Fd2vdxdx(N,5);

  Fd2vdxdx.differentiate(0); Fd2vdxdx.differentiate(0);

/*
  f.bounds[0][0] = x[0];
  f.bounds[0][1] = x[1];
  f.bounds[1][0] = y[0];
  f.bounds[1][1] = y[1];
  f.bounds[2][0] = z[0];
  f.bounds[2][1] = z[1];
*/

std::vector<std::vector<Float> > xvs; std::vector<Float> fvs;

  std::vector<Float> xv;
  xv.resize(N);
 
  Float fv;

#define FILL(VFUN,i1) \
{ \
  xv[0] = x[i1]; \
  fv = VFUN[i1]; \
  addEquation(F##VFUN,xv,fv,xvs,fvs); \
}

for(size_t i1 = 0; i1 < 2; i1++)
  FILL(v,i1);

for(size_t i1 = 0; i1 < 2; i1++)
  FILL(dvdx,i1);

for(size_t i1 = 0; i1 < 2; i1++)
  FILL(d2vdxdx,i1);

#define FILL0(VFUN,i1) \
{ \
  xv[0] = x[i1]; \
  fv = 0.0; \
  addEquation(F##VFUN,xv,fv,xvs,fvs); \
}

  f.addEquations(xvs,fvs);


  dfdx = f;

  dfdx.differentiate(0);

  for(int i = 0;i < 2;i++)
  {
    x_[i] = x[i];
  }  

}


}

