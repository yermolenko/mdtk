/*
   Some useful macros and routines (header file).

   Copyright (C) 2004, 2009 Oleksandr Yermolenko
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

#ifndef mdtk_tools_hpp
#define mdtk_tools_hpp


#include <cstdio>
#include <cstdlib>

#include <string>
#include <vector>

#include <mdtk/config.hpp>

#include <mdtk/Exception.hpp>


#include "tools.h"

#include <mdtk/yaatk.hpp>


#include <string>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>

#include <vector>

namespace mdtk
{

template <class T>
Float SIGN2FLOAT(T x)
{
  return ((x<0)?(-1.0):(1.0));
}

template <class T>
Float BOOL2FLOAT(T x)
{
  return ((x)?(1.0):(0.0));
}

template <class T>
T SQR(T x)
{
  return x*x;
}

template <class T>
T SQR3(T x)
{
  return x*x*x;
}


inline 
Float 
academic_round(Float enr)
{
  int   dig;
  dig = int(fabs(enr)*100)%10;
  if (enr < 0)
  {
    if (dig >= 5)
      enr = floor(enr*10.0)/10.0;
    else
      enr = ceil(enr*10.0)/10.0;
  }
  else
  {
    if (dig >= 5)
      enr = ceil(enr*10.0)/10.0;
    else
      enr = floor(enr*10.0)/10.0;
  }
  return enr;
}

inline 
Float 
intpow(Float a, int n)
{
  switch (n)
  {
    case 0 : return 1.0;
    case 1 : return a;
    case -1 : return 1.0/(a);
    case 2 : return a*a;
    case -2 : return 1.0/(a*a);
    case 3 : return a*a*a;
    case -3 : return 1.0/(a*a*a);
    case 4 : return a*a*a*a;
    case -4 : return 1.0/(a*a*a*a);
    case 5 : return a*a*a*a*a;
    case -5 : return 1.0/(a*a*a*a*a);
    case 6 : return a*a*a*a*a*a;
    case -6 : return 1.0/(a*a*a*a*a*a);
    default: {
  Float v = 1.0;
  int abs_n = abs(n);
  for(int i = 1; i <= abs_n; i++)
    v *= a;
  if (n < 0) v = 1.0/v;
  return v;
             }
  }
}

}

#endif

