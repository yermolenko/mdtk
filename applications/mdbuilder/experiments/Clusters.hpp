/*
   Building of various clusters

   Copyright (C) 2010, 2011 Oleksandr Yermolenko
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

#ifndef MDBUILDER_Clusters_HPP
#define MDBUILDER_Clusters_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

inline
void
place_C60(mdtk::SimLoop& sl)
{
  std::ifstream fcoords("C60.coords");
  REQUIRE(fcoords != 0);
  size_t atomsCount;
  fcoords >> atomsCount;
  
  for(size_t i = 0; i < atomsCount; i++)
  {
    int index, tmp1, tmp2;
    Float x,y,z;
    fcoords >> index >> tmp1 >> tmp2 >> x >> y >> z;
    place(C_EL,sl,Vector3D(x*Ao,y*Ao,z*Ao));
  }
  
  fcoords.close();
}

inline
void
build_C60_optimized(mdtk::SimLoop& sl)
{
  glLoadIdentity();
  mdbuilder::place_C60(sl);
  
  yaatk::mkdir("_tmp-C60-optimized");
  yaatk::chdir("_tmp-C60-optimized");
  setupPotentials(sl);
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.2*ps;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.thermalBath.zMin = -100000.0*Ao;
  sl.execute();
  yaatk::chdir("..");
}

}

#endif
