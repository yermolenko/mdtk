/*
   Building of FCC stuctures (header file)

   Copyright (C) 2008, 2011, 2012 Oleksandr Yermolenko
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

#ifndef MDBUILDER_FCC_HPP
#define MDBUILDER_FCC_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_FCC_cell(
  SimLoop& sl,
  ElementID el = Cu_EL,
  double a = 3.615*Ao,
  double b = 3.615*Ao,
  double c = 3.615*Ao
  );

void
place_FCC_lattice(
  SimLoop& sl,
  int a_num = 14,
  int b_num = 14,
  int c_num = 7,
  ElementID el = Cu_EL,
  bool fixBottomLayer = true,
  double a = 3.615*Ao,
  double b = 3.615*Ao,
  double c = 3.615*Ao
  );

SimLoop
build_FCC_lattice(
  int a_num = 14,
  int b_num = 14,
  int c_num = 7,
  ElementID el = Cu_EL,
  bool fixBottomLayer = true,
  double a = 3.615*Ao,
  double b = 3.615*Ao,
  double c = 3.615*Ao
  );

void
place_Generic_FCC_cell(
  SimLoop& sl,
  const SimLoop sl_element,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  );

void
place_Generic_NegFCC_cell(
  SimLoop& sl,
  const SimLoop sl_element,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  );

void
place_Generic_FCC_lattice(
  SimLoop& sl,
  const SimLoop sl_element,
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  int cellsFromXYPlane = 0,
  double a = 1.0*Ao,
  double b = 1.0*Ao,
  double c = 1.0*Ao,
  bool negative = false
  );
}

#endif
