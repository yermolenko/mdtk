/*
   Building of graphite stuctures

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

#include "Graphite.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_Graphite_cell(
  AtomsArray& sl,
  double a,
  double b,
  double c,
  double gamma
  )
{
  glPushMatrix();

  {
    glPushMatrix();
    glTranslated(0,0,0);
    place(C_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(a*cos(60.0*Deg),a*sin(60.0*Deg)/3.0,0);
    place(C_EL,sl);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(a*cos(60.0*Deg)*3.0,a*sin(60.0*Deg),c/2.0);
    place(C_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(a*cos(60.0*Deg)*2.0,a*sin(60.0*Deg)/3.0*2.0,c/2.0);
    place(C_EL,sl);
    glPopMatrix();
  }

  glPopMatrix();
}

void
place_Graphite_lattice(
  AtomsArray& sl,
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  double a,
  double b,
  double c,
  double gamma
  )
{
  glPushMatrix();

  Vector3D va = Vector3D(a, 0, 0);
  Vector3D vb = Vector3D(b*cos(gamma), b*sin(gamma), 0);
  Vector3D vc = Vector3D(0, 0, c);

  for(int ia = 0; ia < a_num; ia++)
    for(int ib = 0; ib < b_num; ib++)
      for(int ic = 0; ic < c_num; ic++)
      {
        glPushMatrix();

        glTranslated(
          (va*ia+vb*ib+vc*ic).x,
          (va*ia+vb*ib+vc*ic).y,
          (va*ia+vb*ib+vc*ic).z
          );

        glPushMatrix();
        place_Graphite_cell(sl,a,b,c,gamma);
        glPopMatrix();

        if (fixBottomCellLayer)
        {
          if (ic == c_num-1)
          {
            for(size_t i = 1; i <= 4; i++)
            {
              Atom& a = sl[sl.size()-i];
              a.fix();
            }
          }
        }

        glPopMatrix();
      }

  glPopMatrix();
}

SimLoop
build_Graphite_lattice(
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  double a,
  double b,
  double c,
  double gamma
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  place_Graphite_lattice(sl.atoms,a_num,b_num,c_num,fixBottomCellLayer,a,b,c);

  sl.setPBC(Vector3D(a*a_num, b*b_num*sin(gamma), NO_PBC.z));
  sl.thermalBath.zMin = c*c_num-c*3-c/4.0;
  sl.thermalBath.dBoundary = 3.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,1.0*K);

  sl.atoms.removeMomentum();

  return sl;
}

}
