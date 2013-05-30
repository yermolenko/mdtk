/*
   Building of FCC stuctures

   Copyright (C) 2008, 2011, 2012, 2013 Oleksandr Yermolenko
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

#include "FCC.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_FCC_cell(
  AtomsArray& sl,
  ElementID el,
  double a,
  double b,
  double c
  )
{
  glPushMatrix();

  glPushMatrix();
  glTranslated(0,0,0);
  place(el,sl);
  glPopMatrix();

  glPushMatrix();
  glTranslated(a/2.0,a/2.0,0);
  place(el,sl);
  glPopMatrix();

  glPushMatrix();
  glTranslated(0,a/2.0,a/2.0);
  place(el,sl);
  glPopMatrix();

  glPushMatrix();
  glTranslated(a/2.0,0,a/2.0);
  place(el,sl);
  glPopMatrix();

  glPopMatrix();
}

void
place_FCC_lattice(
  AtomsArray& sl,
  int a_num,
  int b_num,
  int c_num,
  ElementID el,
  bool fixBottomLayer,
  double a,
  double b,
  double c
  )
{
  glPushMatrix();

  Vector3D va = Vector3D(a, 0, 0);
  Vector3D vb = Vector3D(0, b, 0);
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
        place_FCC_cell(sl,el,a,b,c);
        glPopMatrix();

        if (fixBottomLayer)
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
build_FCC_lattice(
  int a_num,
  int b_num,
  int c_num,
  ElementID el,
  bool fixBottomLayer,
  double a,
  double b,
  double c
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  place_FCC_lattice(sl.atoms,a_num,b_num,c_num,el,fixBottomLayer,a,b,c);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_BOX;
  sl.thermalBathGeomBox.zMin = c*c_num-c*3-c/4.0;
  sl.thermalBathGeomBox.dBoundary = 3.0*Ao;
  sl.thermalBathGeomBox.zMinOfFreeZone = -5.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,1.0*K);

  sl.atoms.removeMomentum();

  return sl;
}

void
place_Generic_FCC_cell(
  AtomsArray& sl,
  const AtomsArray sl_element,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  )
{
  glPushMatrix();

  {
    glPushMatrix();
    glTranslated(0,0,0);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(((va+vb)/2.0).x,((va+vb)/2.0).y,0.0);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(0.0,((vb+vc)/2.0).y,((vb+vc)/2.0).z);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(((va+vc)/2.0).x,0.0,((va+vc)/2.0).z);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  glPopMatrix();
}

void
place_Generic_NegFCC_cell(
  AtomsArray& sl,
  const AtomsArray sl_element,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  )
{
  glPushMatrix();

  {
    glPushMatrix();
    glTranslated(((va+vb+vc)/2.0).x,((va+vb+vc)/2.0).y,((va+vb+vc)/2.0).z);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated((va/2.0).x,0.0,0.0);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(0.0,(vb/2.0).y,0.0);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  {
    glPushMatrix();
    glTranslated(0.0,0.0,(vc/2.0).z);
    place_Cluster(sl,sl_element);
    glPopMatrix();
  }

  glPopMatrix();
}

void
place_Generic_FCC_lattice(
  AtomsArray& sl,
  const AtomsArray sl_element,
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  int cellsFromXYPlane,
  double a,
  double b,
  double c,
  bool negative
  )
{
  glPushMatrix();

  Vector3D va = Vector3D(a,0,0);
  Vector3D vb = Vector3D(0,b,0);
  Vector3D vc = Vector3D(0,0,c);

  for(int ia = 0; ia < a_num; ia++)
    for(int ib = 0; ib < b_num; ib++)
      for(int ic = cellsFromXYPlane; ic < c_num; ic++)
      {
        glPushMatrix();

        glTranslated(
          (va*ia+vb*ib+vc*ic).x,
          (va*ia+vb*ib+vc*ic).y,
          (va*ia+vb*ib+vc*ic).z
          );

        glPushMatrix();
        if (negative)
          place_Generic_NegFCC_cell(sl,sl_element,va,vb,vc);
        else
          place_Generic_FCC_cell(sl,sl_element,va,vb,vc);
        glPopMatrix();

        if (fixBottomCellLayer)
        {
          if (ic == c_num-1)
          {
            for(size_t i = 1; i <= 4*sl_element.size(); i++)
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

}
