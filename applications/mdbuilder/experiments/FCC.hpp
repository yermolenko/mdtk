/*
   Building of FCC stuctures

   Copyright (C) 2008, 2011 Oleksandr Yermolenko
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

#ifndef MDBUILDER_Cu_HPP
#define MDBUILDER_Cu_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

inline
void
place_FCC_cell(mdtk::SimLoop& sl, 
               ElementID el = Cu_EL,
               double a = 3.615*Ao,
               double b = 3.615*Ao,
               double c = 3.615*Ao
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

inline
void
place_FCC_lattice(mdtk::SimLoop& sl, 
                  int a_num = 14,
                  int b_num = 14,
                  int c_num = 7,
                  ElementID el = Cu_EL,
                  bool fixBottomLayer = true,
                  double a = 3.615*Ao,
                  double b = 3.615*Ao,
                  double c = 3.615*Ao
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
              Atom& a = *(sl.atoms[sl.atoms.size()-i]);
              a.tag |= ATOMTAG_FIXED;
              if (a.tag & ATOMTAG_FIXED)
              {
                a.M = INFINITE_MASS;a.V=0.0;a.an=0.0;a.an_no_tb=0.0;
              }
            }
          }
        }

        glPopMatrix();
      }

  glPopMatrix();
}

inline
mdtk::SimLoop
build_FCC_lattice(int a_num = 14,
                  int b_num = 14,
                  int c_num = 7,
                  ElementID el = Cu_EL,
                  bool fixBottomLayer = true,
                  double a = 3.615*Ao,
                  double b = 3.615*Ao,
                  double c = 3.615*Ao
                  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  place_FCC_lattice(sl,a_num,b_num,c_num,el,fixBottomLayer,a,b,c);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBath.zMin = c*c_num-c*3-c/4.0;
  sl.thermalBath.dBoundary = 3.0*Ao;

  quench(sl,0.01*ps,0.01*ps,"_tmp-FCC");

  removeMomentum(sl.atoms);

  return sl;
}

}

#endif
