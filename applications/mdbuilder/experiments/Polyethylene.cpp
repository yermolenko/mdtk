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

#include "Polyethylene.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_Ethylene(SimLoop& sl)
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

void
place_Polyethylene_fold(
  SimLoop& sl,
  Vector3D va,
  Vector3D vb,
  Vector3D vc
  )
{
  glTranslated(0.0,0.0,+vc.z*0.3);

  glPushMatrix();
  {
    glTranslated(((va+vb)/2.0).x *(+0.1),((va+vb)/2.0).y *(+0.1),0.0);

    glRotated(45.0,0.0,0.0,1.0);

    place(C_EL,sl);

    glPushMatrix();
    glTranslated(+0.3*Ao,-1.0*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(-0.3*Ao,1.0*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glTranslated(0.0,0.0,-vc.z*0.7);
//    glTranslated(((va+vb)/2.0).x *(+0.3),((va+vb)/2.0).y *(+0.1),0.0);

    place(C_EL,sl);

    glPushMatrix();
    glTranslated(0.0,-1.0*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated( 0.0,1.0*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();
  }
  glPopMatrix();

  glPushMatrix();
  {
    glTranslated(((va+vb)/2.0).x,((va+vb)/2.0).y,0.0);
    glTranslated(((va+vb)/2.0).x *(-0.1),((va+vb)/2.0).y *(-0.1),0.0);

    glRotated(-90.0,0.0,0.0,1.0);
    glRotated( 45.0,0.0,0.0,1.0);

    place(C_EL,sl);

    glPushMatrix();
    glTranslated(-1.0*Ao,+0.3*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated( 1.0*Ao,-0.3*Ao,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glTranslated(0.0,0.0,-vc.z*0.7);
//    glTranslated(((va+vb)/2.0).x *(-0.3),((va+vb)/2.0).y *(-0.1),0.0);

    place(C_EL,sl);

    glPushMatrix();
    glTranslated(-1.0*Ao,0.0,0.0);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated( 1.0*Ao,0.0,0.0);
    place(H_EL,sl);
    glPopMatrix();
  }
  glPopMatrix();

  glTranslated(0.0,0.0,-vc.z*0.9);

  glPushMatrix();
  {
    glTranslated(((va+vb)/2.0).x*(+0.5),((va+vb)/2.0).y*(+0.5),/*+c*0.7*/0.0);

    place(C_EL,sl);

    glPushMatrix();
    glTranslated(-0.7*Ao,0.7*Ao,-0.5*Ao);
    place(H_EL,sl);
    glPopMatrix();

    glPushMatrix();
    glTranslated(0.7*Ao,-0.7*Ao,-0.5*Ao);
    place(H_EL,sl);
    glPopMatrix();
  }
  glPopMatrix();
}

void
place_Polyethylene_cell(
  SimLoop& sl,
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

void
place_Polyethylene_lattice(
  SimLoop& sl,
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  int cellsFromXYPlane,
  double a,
  double b,
  double c
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

void
place_Polyethylene_folds(
  SimLoop& sl,
  int a_num,
  int b_num,
  int c_num,
  int cellsFromXYPlane,
  double a,
  double b,
  double c
  )
{
  glPushMatrix();

  Vector3D va = Vector3D(a,0,0);
  Vector3D vb = Vector3D(0,b,0);
  Vector3D vc = Vector3D(0,0,c);

  for(int ia = 0; ia < a_num; ia++)
    for(int ib = 0; ib < b_num; ib++)
      for(int ic = cellsFromXYPlane-1; ic < cellsFromXYPlane; ic++)
      {
        glPushMatrix();

        glTranslated(
          (va*ia+vb*ib+vc*ic).x,
          (va*ia+vb*ib+vc*ic).y,
          (va*ia+vb*ib+vc*ic).z
          );

        glPushMatrix();
        place_Polyethylene_fold(sl,va,vb,vc);
        glPopMatrix();

        glPopMatrix();
      }

  glPopMatrix();
}

void
place_Polyethylene_folded_chains(
  SimLoop& sl,
  const SimLoop& sl_with_chain,
  int a_num,
  int b_num,
  double a,
  double b
  )
{
  glPushMatrix();

  Vector3D va = Vector3D(a,0,0);
  Vector3D vb = Vector3D(0,b,0);

  for(int ia = 0; ia < a_num; ia++)
    for(int ib = 0; ib < b_num; ib++)
    {
      glPushMatrix();

      glTranslated(
        (va*ia+vb*ib).x,
        (va*ia+vb*ib).y,
        (va*ia+vb*ib).z
        );

      glPushMatrix();

      for(size_t ai = 0; ai < sl_with_chain.atoms.size(); ai++)
      {
        const mdtk::Atom& atom = *sl_with_chain.atoms[ai];
        REQUIRE(atom.ID == C_EL || atom.ID == H_EL);
        place_and_inherit(sl,atom,getPosition()+atom.coords);
      };

      glPopMatrix();

      glPopMatrix();
    }

  glPopMatrix();
}

SimLoop
build_Polyethylene_lattice_without_folds(
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  double a,
  double b,
  double c
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

SimLoop
build_Polyethylene_lattice_with_folds(
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  double a,
  double b,
  double c
  )
{
  SimLoopDump sl_with_chain;
  initialize_simloop(sl_with_chain);

  {
    SimLoopDump sl_rebo;
    initialize_simloop_REBO_only(sl_rebo);

    place_Polyethylene_lattice(sl_rebo,1,1,c_num,fixBottomCellLayer,2,a,b,c);

    std::vector<size_t> fixedAtoms =
      fixNotFixedAtoms(sl_rebo.atoms,0,sl_rebo.atoms.size());
    {
      place_Polyethylene_folds(sl_rebo,1,1,c_num,2,a,b,c);

      sl_rebo.enableDump();

      sl_rebo.dumpConst(0.95);
      relax(sl_rebo,0.05*ps);

      sl_rebo.dumpConst(0.97);
      relax(sl_rebo,0.05*ps);

      sl_rebo.dumpConst(0.99);
      relax(sl_rebo,0.05*ps);

      sl_rebo.disableDump();

      quench(sl_rebo,0.01*K);

      {
        SimLoopDump sl_airebo(sl_rebo);
        initialize_simloop(sl_airebo);

        sl_airebo.enableDump();

        sl_airebo.dumpConst(0.95);
        relax/*_flush*/(sl_airebo,0.05*ps,"_tmp-X-relax_flush-folds-airebo");

        sl_airebo.dumpConst(0.97);
        relax(sl_airebo,0.05*ps);

        sl_airebo.dumpConst(0.99);
        relax(sl_airebo,0.05*ps);

        sl_airebo.disableDump();

        quench(sl_airebo,0.01*K);

        sl_with_chain = sl_airebo;
      }
    }
    unfixAtoms(sl_with_chain.atoms,fixedAtoms);

    {
      yaatk::text_ofstream fomde("_tmp-X-relax_flush-folds-airebo-unfixed.mde");
      sl_with_chain.saveToMDE(fomde);
      fomde.close();
    }

    sl_with_chain.enableDump();

    sl_with_chain.dumpConst(0.95);
    relax/*_flush*/(sl_with_chain,0.05*ps,"_tmp-X-relax_flush-folds-airebo-unfixed");

    sl_with_chain.dumpConst(0.97);
    relax(sl_with_chain,0.05*ps);

    sl_with_chain.dumpConst(0.99);
    relax(sl_with_chain,0.05*ps);

    sl_with_chain.disableDump();

    quench(sl_with_chain,0.01*K);
  }

  removeMomentum(sl_with_chain.atoms);

  SimLoopDump sl;
  initialize_simloop(sl);

  place_Polyethylene_folded_chains(sl,sl_with_chain,a_num,b_num,a,b);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));

  {
    yaatk::text_ofstream fomde("_tmp-X-relax_flush-FINAL.mde");
    sl.saveToMDE(fomde);
    fomde.close();
  }

  sl.enableDump();

  sl.dumpConst(0.95);
  relax/*_flush*/(sl,0.05*ps,"_tmp-X-relax_flush-FINAL");

  sl.dumpConst(0.97);
  relax(sl,0.05*ps);

  sl.dumpConst(0.99);
  relax(sl,0.05*ps);

  sl.disableDump();

  quench(sl,0.01*K);

  sl.thermalBath.zMin = (c_num > 3)?(c*(c_num-3)-0.5*Ao):(0.0);
  sl.thermalBath.dBoundary = 3.0*Ao;
  sl.thermalBath.zMinOfFreeZone = -2.0*Ao;

  removeMomentum(sl.atoms);

  return sl;
}

}
