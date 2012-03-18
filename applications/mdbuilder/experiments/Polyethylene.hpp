/*
   Building of polyethylene stuctures

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

#ifndef MDBUILDER_Polyethylene_HPP
#define MDBUILDER_Polyethylene_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_Ethylene(AtomsArray& sl);

void
place_Polyethylene_fold(
  AtomsArray& sl,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  );

void
place_Polyethylene_cell(
  AtomsArray& sl,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  );

void
place_Polyethylene_lattice(
  AtomsArray& sl,
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  int cellsFromXYPlane = 0,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
  );

void
place_Polyethylene_folds(
  AtomsArray& sl,
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  int cellsFromXYPlane = 0,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
  );

void
place_Polyethylene_folded_chains(
  AtomsArray& sl,
  const mdtk::SimLoop& sl_with_chain,
  int a_num = 8,
  int b_num = 12,
  double a = 7.417*Ao,
  double b = 4.945*Ao
  );

SimLoop
build_Polyethylene_lattice_without_folds(
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
  );

SimLoop
build_Polyethylene_lattice_with_folds(
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
  );
}

#endif
