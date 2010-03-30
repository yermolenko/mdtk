/*
   The bicubic spline class header file.

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

#ifndef mdtk_Spline2D_hpp
#define mdtk_Spline2D_hpp

#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>

#include <mdtk/SplineAux.hpp>

namespace mdtk
{

class Spline2D
{
private:
  SplineMultiD f;
  SplineMultiD dfdx;
  SplineMultiD dfdy;
  Float  x_[2];
  Float  y_[2];
public:
  Float operator()(Float x, Float y) const;
  Float dx(Float x, Float y) const;
  Float dy(Float x, Float y) const;
  void init(
         Float       x[2],
         Float       y[2],
         Float    v[2][2],
         Float dvdx[2][2],
         Float dvdy[2][2]
         );
  Spline2D(
         Float       x[2],
         Float       y[2],
         Float    v[2][2],
         Float dvdx[2][2],
         Float dvdy[2][2]
         )
  : f(2), dfdx(2), dfdy(2) {init(x,y,v,dvdx,dvdy);}
  Spline2D() : f(2), dfdx(2), dfdy(2){};
  virtual ~Spline2D() {}
  
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    int i;
    f.SaveToStream(os,smode);
    dfdx.SaveToStream(os,smode);
    dfdy.SaveToStream(os,smode);

    for(i = 0;i < 2;i++)
      YAATK_FSTREAM_WRITE(os,x_[i],smode);
    for(i = 0;i < 2;i++)
      YAATK_FSTREAM_WRITE(os,y_[i],smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    int i;
    f.LoadFromStream(is,smode);
    dfdx.LoadFromStream(is,smode);
    dfdy.LoadFromStream(is,smode);

    for(i = 0;i < 2;i++)
      YAATK_FSTREAM_READ(is,x_[i],smode);
    for(i = 0;i < 2;i++)
      YAATK_FSTREAM_READ(is,y_[i],smode);
  }  
};


inline
Float
Spline2D::operator()(Float x, Float y) const
{
  std::vector<Float> xv(2);
  xv[0] = x;
  xv[1] = y;
  return f(xv);
}

inline
Float
Spline2D::dx(Float x, Float y) const
{
  std::vector<Float> xv(2);
  xv[0] = x;
  xv[1] = y;
  return dfdx(xv);
}

inline
Float
Spline2D::dy(Float x, Float y) const
{
  std::vector<Float> xv(2);
  xv[0] = x;
  xv[1] = y;
  return dfdy(xv);
}


} 

#endif

