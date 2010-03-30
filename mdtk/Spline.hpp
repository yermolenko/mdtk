/*
   The cubic spline class header file.

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

#ifndef mdtk_Spline_hpp
#define mdtk_Spline_hpp

#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>

#include <mdtk/SplineAux.hpp>

namespace mdtk
{

class Spline
{
private:
  SplineMultiD f;
  SplineMultiD dfdx;
  Float  x_[2];
public:
  Float operator()(Float x) const;
  Float der(Float x) const;
  void init(
         Float       x[2],
         Float       v[2],
         Float    dvdx[2]
  );
  Spline(
         Float       x[2],
         Float       v[2],
         Float    dvdx[2]
  ): f(1), dfdx(1) {init(x,v,dvdx);};
  virtual ~Spline() {}

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
Spline::operator()(Float x) const
{
  std::vector<Float> xv(1);
  xv[0] = x;
  return f(xv);
}

inline
Float
Spline::der(Float x) const
{
  std::vector<Float> xv(1);
  xv[0] = x;
  return dfdx(xv);
}

} 

#endif

