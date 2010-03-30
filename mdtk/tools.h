/*
   Some useful macros.

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

#include <string>           
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>

#include <cassert>

#include <vector>

#include "Exception.hpp"

#define max2(x,y) ((x)>(y)?(x):(y))
#define max3(x,y,z) ((max2((x),(y)))>(z)?(max2((x),(y))):(z))

#define min2(x,y) ((x)>(y)?(y):(x))
#define min3(x,y,z) ((min2((x),(y)))>(z)?(z):(min2((x),(y))))


#define powTo6(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define powTo2(x) ((x)*(x))

//#define MDTK_TRY_UNUSE_CPU

#include <mdtk/procmon.hpp>

#ifdef MDTK_TRY_UNUSE_CPU
#define MDTK_RELEASE_CPU {procmon::unUseCPU();};
#else
#define MDTK_RELEASE_CPU {;};
#endif



