/*
   Building of polyethylene stuctures

   Copyright (C) 2007, 2008, 2009, 2011 Oleksandr Yermolenko
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

inline
void
place_Ethylene(
  mdtk::SimLoop& sl)
{
  Float CC_dist = 1.53*Ao;
  Float CH_dist = 1.07*Ao;
  Float CCC_angle = 112.0*Deg;
  Float HCH_angle = 107.0*Deg;

  glPushMatrix();

  {
    glPushMatrix();

    glTranslated(CC_dist*cos(CCC_angle/2.0)/2.0,0.0,0.0);
    place(C_EL,sl);

    glPushMatrix();
    glTranslated(CH_dist*cos(HCH_angle/2.0),CH_dist*sin(HCH_angle/2.0),0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(CH_dist*cos(HCH_angle/2.0),-CH_dist*sin(HCH_angle/2.0),0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPopMatrix();

    glTranslated(0.0,0.0,CC_dist*sin(CCC_angle/2.0));

    glPushMatrix();

    glTranslated(-CC_dist*cos(CCC_angle/2.0)/2.0,0.0,0.0);
    place(C_EL,sl);

    glPushMatrix();
    glTranslated(-CH_dist*cos(HCH_angle/2.0),CH_dist*sin(HCH_angle/2.0),0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-CH_dist*cos(HCH_angle/2.0),-CH_dist*sin(HCH_angle/2.0),0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPopMatrix();
  }

  glPopMatrix();
}

inline
void
place_Polyethylene_cell(
  mdtk::SimLoop& sl, 
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  )
{
  glPushMatrix();
  
  {
    glPushMatrix();
    glRotated(45.0,0.0,0.0,1.0);
    place_Ethylene(sl);
    glPopMatrix();
  }
    
  {
    glPushMatrix();
    glTranslated(((va+vb)/2.0).x,((va+vb)/2.0).y,0.0);
    glRotated(-90.0,0.0,0.0,1.0);
    glRotated(45.0,0.0,0.0,1.0);
    place_Ethylene(sl);
    glPopMatrix();
  }

  glPopMatrix();
}

inline
void
place_Polyethylene_lattice(
  mdtk::SimLoop& sl, 
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  int cellsFromXYPlane = 0,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
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
        place_Polyethylene_cell(sl,va,vb,vc);
        glPopMatrix();
        
        if (fixBottomCellLayer)
        {     
          if (ic == c_num-1)
          {
            for(size_t i = 1; i <= 12; i++)
            {
              Atom& a = *(sl.atoms[sl.atoms.size()-i]);
              a.fix();
            }
          }
        }

        glPopMatrix();
      }        

  glPopMatrix();
}

inline
mdtk::SimLoop
build_Polyethylene_lattice_without_folds(
  int a_num = 8,
  int b_num = 12,
  int c_num = 10,
  bool fixBottomCellLayer = true,
  double a = 7.417*Ao,
  double b = 4.945*Ao,
  double c = 2.547*Ao
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  place_Polyethylene_lattice(sl,a_num,b_num,c_num,fixBottomCellLayer,0,a,b,c);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBath.zMin = (c_num > 3)?(c*(c_num-3)-0.5*Ao):(0.0);
  sl.thermalBath.dBoundary = 3.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,1.0*K);

  removeMomentum(sl.atoms);

  return sl;
}

}

#endif
