/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. Splines (header file).
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009 Oleksandr Yermolenko
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

#ifndef mdtk_AIREBO_cubici_hpp
#define mdtk_AIREBO_cubici_hpp

#include <cstdlib>
#include <cctype>

#include "mdtk/config.hpp"
#include "mdtk/tools.hpp"

#include "mdtk/Spline2D.hpp"
#include "mdtk/Spline3D.hpp"
#include "mdtk/Spline5n.hpp"

namespace mdtk
{

class FuncG
{
protected:
  Spline5n spline[100];
  int      splineCount;
public:
  FuncG(){;};
  virtual ~FuncG(){;};
  Float operator()(Float CosT) const;  
  Float dCosT(Float CosT) const;  

  virtual void init() = 0;
};

class FuncG_C1 : public FuncG
{
public:
  void init();
  FuncG_C1()
   :FuncG()
   {
   }  
};

class FuncG_C2 : public FuncG
{
public:
  void init();
  FuncG_C2()
   :FuncG()
   {
   }  
};

class FuncG_H : public FuncG
{
public:
  void init();
  FuncG_H()
   :FuncG()
   {
   }  
};

inline
Float
FuncG::operator()(Float CosT) const
{
  REQUIREM(CosT>=-1.0 && CosT<=+1.0,"FuncG: CosT>=-1.0 && CosT<=+1.0");
  int splineIndex;
  for(splineIndex = 0; splineIndex < splineCount; splineIndex++)
    if (CosT <= spline[splineIndex].x2()) break;
  return spline[splineIndex].operator()(CosT);
}
    
inline
Float
FuncG::dCosT(Float CosT) const
{
  REQUIREM(CosT>=-1.0 && CosT<=+1.0,"FuncG: CosT>=-1.0 && CosT<=+1.0");
  int splineIndex;
  for(splineIndex = 0; splineIndex < splineCount; splineIndex++)
    if (CosT <= spline[splineIndex].x2()) break;

  return spline[splineIndex].der(CosT);
}




class FuncP_CC
{
  Spline2D spline[4][4];
public:
  FuncP_CC(int paramSet = 0);
  virtual ~FuncP_CC(){};
  Float operator()(Float h, Float c) const;  
  Float dH(Float h, Float c) const;  
  Float dC(Float h, Float c) const;  
};

class FuncP_CH
{
  Spline2D spline[4][4];
public:
  FuncP_CH(int paramSet = 0);
  virtual ~FuncP_CH(){};
  Float operator()(Float h, Float c) const;  
  Float dH(Float h, Float c) const;  
  Float dC(Float h, Float c) const;  
};

class Func_pi_rc
{
protected:
#define i_size_pi_rc 5
#define j_size_pi_rc 5
#define k_size_pi_rc 16
  Spline3D spline[i_size_pi_rc][j_size_pi_rc][k_size_pi_rc];
public:
  Func_pi_rc()
  {
  }  
  virtual ~Func_pi_rc(){};
  virtual void init(int paramSet = 0) = 0;
  Float operator()(Float i, Float j, Float k) const;  
  Float di(Float i, Float j, Float k) const;  
  Float dj(Float i, Float j, Float k) const;  
  Float dk(Float i, Float j, Float k) const;  

  Float di_num(Float i, Float j, Float k) const;  
  Float dj_num(Float i, Float j, Float k) const;  
  Float dk_num(Float i, Float j, Float k) const;  
};

class Func_pi_rc_CC : public Func_pi_rc
{
public:
  void init(int paramSet = 0);
  Func_pi_rc_CC()
   :Func_pi_rc()
   {
   }  
};

class Func_pi_rc_CH : public Func_pi_rc
{
public:
  void init(int paramSet = 0);
  Func_pi_rc_CH()
   :Func_pi_rc()
   {
   }  
};

class Func_pi_rc_HH : public Func_pi_rc
{
public:
  void init(int paramSet = 0);
  Func_pi_rc_HH()
   :Func_pi_rc()
   {
   }  
};

class Func_Tij : public Func_pi_rc
{
public:
  void init(int paramSet = 0);
  Func_Tij()
   :Func_pi_rc()
   {
   }  
};


inline
Float
FuncP_CC::operator()(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].operator()(h,c);
}
    
inline
Float
FuncP_CC::dH(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dx(h,c);
}
    
inline
Float
FuncP_CC::dC(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dy(h,c);
}


inline
Float
FuncP_CH::operator()(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].operator()(h,c);
}

    
inline
Float
FuncP_CH::dH(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dx(h,c);
}

    
inline
Float
FuncP_CH::dC(Float h, Float c) const
{
  REQUIREM(h>=0 && c>=0,"H: h>=0 && c>=0");
  if (h>=4 || c>=4) return 0.0;
  return spline[int(h)][int(c)].dy(h,c);
}


inline
Float
Func_pi_rc::operator()(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=i_size_pi_rc-1 || j>=j_size_pi_rc-1) return 0.0;
  if (k > (k_size_pi_rc-1)) k = (k_size_pi_rc-1);
  return spline[int(i)][int(j)][int(k)].operator()(i,j,k);
}

inline
Float
Func_pi_rc::di(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=i_size_pi_rc-1 || j>=j_size_pi_rc-1) return 0.0;
  if (k > (k_size_pi_rc-1)) k = (k_size_pi_rc-1);
  return spline[int(i)][int(j)][int(k)].dx(i,j,k);
}


inline
Float
Func_pi_rc::dj(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=i_size_pi_rc-1 || j>=j_size_pi_rc-1) return 0.0;
  if (k > (k_size_pi_rc-1)) k = (k_size_pi_rc-1);
  return spline[int(i)][int(j)][int(k)].dy(i,j,k);
}


inline
Float
Func_pi_rc::dk(Float i, Float j, Float k) const
{
  REQUIREM(i>=0 && j>=0 && k>=1,"F: i>=0 && j>=0 && k>=1");
  if (i>=i_size_pi_rc-1 || j>=j_size_pi_rc-1) return 0.0;
  if (k > (k_size_pi_rc-1)) k = (k_size_pi_rc-1);
  return spline[int(i)][int(j)][int(k)].dz(i,j,k);
}

}

#endif

