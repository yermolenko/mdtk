/*
   Building of fullerite stuctures

   Copyright (C) 2007, 2008, 2009, 2011, 2012 Oleksandr Yermolenko
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

#ifndef MDBUILDER_Fullerite_HPP
#define MDBUILDER_Fullerite_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

SimLoop
build_Fullerite_C60(
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  double a = 14.17*Ao
  );

}

#endif
