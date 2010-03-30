/*
   The tricubic spline class.

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

#include "Spline3D.hpp"

namespace mdtk
{

using std::fabs;
using std::exp;
using std::sqrt;
using std::pow;


void
Spline3D::init(
         Float          x[2],
         Float          y[2],
         Float          z[2],
         Float    v[2][2][2],
         Float dvdx[2][2][2],
         Float dvdy[2][2][2],
         Float dvdz[2][2][2]
         )
{

{
  size_t N = 3;

  SplineMultiD Fv(N);

  SplineMultiD Fdvdx(N);
  SplineMultiD Fdvdy(N);
  SplineMultiD Fdvdz(N);

  Fdvdx.differentiate(0);
  Fdvdy.differentiate(1);
  Fdvdz.differentiate(2);

  SplineMultiD Fd2vdxdy(N);
  SplineMultiD Fd2vdxdz(N);
  SplineMultiD Fd2vdydz(N);

  Fd2vdxdy.differentiate(0);Fd2vdxdy.differentiate(1);
  Fd2vdxdz.differentiate(0);Fd2vdxdz.differentiate(2);
  Fd2vdydz.differentiate(1);Fd2vdydz.differentiate(2);

  SplineMultiD Fd3vdxdydz(N);
  Fd3vdxdydz.differentiate(0);Fd3vdxdydz.differentiate(1);Fd3vdxdydz.differentiate(2);

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

#define FILL(VFUN,i1,i2,i3) \
{ \
  xv[0] = x[i1]; \
  xv[1] = y[i2]; \
  xv[2] = z[i3]; \
  fv = VFUN[i1][i2][i3]; \
  addEquation(F##VFUN,xv,fv,xvs,fvs); \
}

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL(v,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL(dvdx,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL(dvdy,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL(dvdz,i1,i2,i3);

#define FILL0(VFUN,i1,i2,i3) \
{ \
  xv[0] = x[i1]; \
  xv[1] = y[i2]; \
  xv[2] = z[i3]; \
  fv = 0.0; \
  addEquation(F##VFUN,xv,fv,xvs,fvs); \
}

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL0(d2vdxdy,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL0(d2vdxdz,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL0(d2vdydz,i1,i2,i3);

for(size_t i1 = 0; i1 < 2; i1++)
for(size_t i2 = 0; i2 < 2; i2++)
for(size_t i3 = 0; i3 < 2; i3++)
  FILL0(d3vdxdydz,i1,i2,i3);

  f.addEquations(xvs,fvs);

}

  dfdx = f;
  dfdy = f;
  dfdz = f;

  dfdx.differentiate(0);
  dfdy.differentiate(1);
  dfdz.differentiate(2);

  for(int i = 0;i < 2;i++)
  {
    x_[i] = x[i];
    y_[i] = y[i];
    z_[i] = z[i];
  }  

}


}

