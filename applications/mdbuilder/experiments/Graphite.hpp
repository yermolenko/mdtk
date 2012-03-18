/*
   Building of graphite stuctures (header file)

   Copyright (C) 2008, 2009, 2011, 2012 Oleksandr Yermolenko
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

#ifndef MDBUILDER_Graphite_HPP
#define MDBUILDER_Graphite_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_Graphite_cell(
  SimLoop& sl,
  double a = 2.46*Ao,
  double b = 2.46*Ao,
  double c = 6.708*Ao,
  double gamma = 60.0*Deg
  );

void
place_Graphite_lattice(
  AtomsArray& sl,
  int a_num = 12,
  int b_num = 14,
  int c_num = 3,
  bool fixBottomCellLayer = true,
  double a = 2.46*Ao,
  double b = 2.46*Ao,
  double c = 6.708*Ao,
  double gamma = 60.0*Deg
  );

SimLoop
build_Graphite_lattice(
  int a_num = 12,
  int b_num = 14,
  int c_num = 3,
  bool fixBottomCellLayer = true,
  double a = 2.46*Ao,
  double b = 2.46*Ao,
  double c = 6.708*Ao,
  double gamma = 60.0*Deg
  );

}

#endif
