/*
   Implementation of the many-body interatomic potential for
   hydrocarbons. Splines (header file).
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009 Oleksandr Yermolenko
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

#ifndef mdtk_cubici_hpp
#define mdtk_cubici_hpp

#include <cstdlib>
#include <cctype>

#include "mdtk/config.hpp"

#include "mdtk/Spline2D.hpp"
#include "mdtk/Spline3D.hpp"

namespace mdtk
{

class FuncH_CC
{
  Spline2D spline[4][4];
public:
  FuncH_CC(int paramSet = 0);
  virtual ~FuncH_CC() {};
  Float operator()(Float h, Float c) const;  
  Float dH(Float h, Float c) const;  
  Float dC(Float h, Float c) const;  
};

class FuncH_CH
{
  Spline2D spline[4][4];
public:
  FuncH_CH(int paramSet = 0);
  virtual ~FuncH_CH(){};
  Float operator()(Float h, Float c) const;  
  Float dH(Float h, Float c) const;  
  Float dC(Float h, Float c) const;  
};

class FuncF
{
  Spline3D spline[4][4][5];
public:
  FuncF(int paramSet = 0);
  virtual ~FuncF(){};
  Float operator()(Float i, Float j, Float k) const;  
  Float di(Float i, Float j, Float k) const;  
  Float dj(Float i, Float j, Float k) const;  
  Float dk(Float i, Float j, Float k) const;  

  Float di_num(Float i, Float j, Float k) const;  
  Float dj_num(Float i, Float j, Float k) const;  
  Float dk_num(Float i, Float j, Float k) const;  
};

inline
Float
FuncH_CC::operator()(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].operator()(h,c);
}
    
inline
Float
FuncH_CC::dH(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dx(h,c);
}
    
inline
Float
FuncH_CC::dC(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dy(h,c);
}


inline
Float
FuncH_CH::operator()(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].operator()(h,c);
}

    
inline
Float
FuncH_CH::dH(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dx(h,c);
}

    
inline
Float
FuncH_CH::dC(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dy(h,c);
}


inline
Float
FuncF::operator()(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=4 || j>=4) return 0.0;
  if (k > (5-1)) k = (5-1);
  return spline[int(i)][int(j)][int(k)].operator()(i,j,k);
}

inline
Float
FuncF::di(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=4 || j>=4) return 0.0;
  if (k > (5-1)) k = (5-1);
  return spline[int(i)][int(j)][int(k)].dx(i,j,k);
}


inline
Float
FuncF::dj(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=4 || j>=4) return 0.0;
  if (k > (5-1)) k = (5-1);
  return spline[int(i)][int(j)][int(k)].dy(i,j,k);
}


inline
Float
FuncF::dk(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=4 || j>=4) return 0.0;
  if (k > (5-1)) {k = (5-1);/*return 0.0;*/};
  return spline[int(i)][int(j)][int(k)].dz(i,j,k);
}

}

#endif

