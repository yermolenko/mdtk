/*
   Building of H2 molecule

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

#ifndef MDBUILDER_H2_HPP
#define MDBUILDER_H2_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

inline
void
place_H2_simple(mdtk::SimLoop& sl)
{
  place(H_EL,sl,Vector3D(0,0,0));
  place(H_EL,sl,Vector3D(1.0*mdtk::Ao,0,0));
}

inline
void
place_H2(mdtk::SimLoop& sl)
{
  glPushMatrix();
  place(H_EL,sl);
  glTranslated(1.0*mdtk::Ao,0,0);
  place(H_EL,sl);
  glPopMatrix();
}

}

#endif