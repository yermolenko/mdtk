/*
   The fifth-order spline class header file.

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

#ifndef mdtk_Spline5n_hpp
#define mdtk_Spline5n_hpp

#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>

#include <mdtk/SplineAux.hpp>

namespace mdtk
{

class Spline5n
{
private:
  SplineMultiD f;
  SplineMultiD dfdx;
  Float  x_[2];
public:
  Float operator()(Float x) const;
  void init(
         Float       x[2],
         Float       v[2],
         Float    dvdx[2],
         Float    d2vdxdx[2]
  );
  Float der(Float x) const;
  Spline5n(
         Float       x[2],
         Float       v[2],
         Float    dvdx[2],
         Float    d2vdxdx[2]
  ):f(1,5), dfdx(1,5)
  {
    init(x,v,dvdx,d2vdxdx);
  }
  Spline5n(Float x1 = 0.0, Float y1 = 1.0, Float dy1 = 0.0, Float d2y1 = 0.0,           Float x2 = 1.0, Float y2 = 1.0, Float dy2 = 0.0, Float d2y2 = 0.0)
  :f(1,5), dfdx(1,5)
  {
    Float       x[2];
    Float       v[2];
    Float    dvdx[2];
    Float    d2vdxdx[2];
 
    x[0] = x1;
    x[1] = x2;
    v[0] = y1;
    v[1] = y2;
    dvdx[0] = dy1;
    dvdx[1] = dy2;
    d2vdxdx[0] = d2y1;
    d2vdxdx[1] = d2y2;
    init(x,v,dvdx,d2vdxdx);
  }
 
  virtual ~Spline5n() {}

  Float x1() const {return x_[0];};
  Float x2() const {return x_[1];};

  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    f.SaveToStream(os,smode);
    dfdx.SaveToStream(os,smode);

    YAATK_FSTREAM_WRITE(os,x_[0],smode);
    YAATK_FSTREAM_WRITE(os,x_[1],smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    f.LoadFromStream(is,smode);
    dfdx.LoadFromStream(is,smode);

    YAATK_FSTREAM_READ(is,x_[0],smode);
    YAATK_FSTREAM_READ(is,x_[1],smode);
  }  
};

inline
Float
Spline5n::operator()(Float x) const
{
  std::vector<Float> xv(1);
  xv[0] = x;
  return f(xv);
}

inline
Float
Spline5n::der(Float x) const
{
  std::vector<Float> xv(1);
  xv[0] = x;
  return dfdx(xv);
}


} 

#endif

