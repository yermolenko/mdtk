/*
   Building of polyethylene stuctures

   Copyright (C) 2007, 2008, 2009, 2011, 2012, 2013 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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
place_Ethylene(AtomsArray& sl)
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
  AtomsArray& sl,
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
  AtomsArray& sl,
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
  AtomsArray& sl,
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
              Atom& a = sl[sl.size()-i];
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
  AtomsArray& sl,
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
  AtomsArray& sl,
  const AtomsArray& sl_with_chain,
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

      for(size_t ai = 0; ai < sl_with_chain.size(); ai++)
      {
        const mdtk::Atom& atom = sl_with_chain[ai];
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

  place_Polyethylene_lattice(sl.atoms,a_num,b_num,c_num,fixBottomCellLayer,0,a,b,c);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_BOX;
  sl.thermalBathGeomBox.zMinOfFreeZone = -5.0*Ao;
  sl.thermalBathGeomBox.zMin = (c_num > 3)?(c*(c_num-3)-0.5*Ao):(0.0);
  sl.thermalBathGeomBox.dBoundary = 3.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,0.1*K);

  sl.atoms.removeMomentum();

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
  yaatk::ChDir cd("_build_PE");

  const gsl_rng_type* T;
  gsl_rng* r;

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc(T);
  REQUIRE(r != NULL);

  gsl_rng_set(r, 123);

  REQUIRE(gsl_rng_min(r) == 0);
  REQUIRE(gsl_rng_max(r) > 1000);

  AtomsArray bulkAtoms;

  {
    SimLoop sl_bulk_PE;
    initialize_simloop_REBO_only(sl_bulk_PE);

    place_Polyethylene_lattice(sl_bulk_PE.atoms,a_num,b_num,c_num,false,2,a,b,c);
    sl_bulk_PE.heatUpEveryAtom(0.001*eV, r);
    sl_bulk_PE.thermalBathGeomType = SimLoop::TB_GEOM_UNIVERSE;
    sl_bulk_PE.atoms.PBC(Vector3D(a*a_num, b*b_num, c*(c_num-2)));

    relax(sl_bulk_PE,1.0*ps,"000-relax-bulk-PE-REBO");
    quench(sl_bulk_PE,0.01*K, 200*ps, 0.01*ps,"001-cooling-bulk-PE-REBO");

    sl_bulk_PE.atoms.PBC(NO_PBC);
    bulkAtoms = sl_bulk_PE.atoms;
  }

  {
    SimLoop sl_bulk_PE;
    initialize_simloop(sl_bulk_PE);

    sl_bulk_PE.atoms = bulkAtoms;
    sl_bulk_PE.thermalBathGeomType = SimLoop::TB_GEOM_UNIVERSE;
    sl_bulk_PE.atoms.PBC(Vector3D(a*a_num, b*b_num, c*(c_num-2)));

    relax(sl_bulk_PE,1.0*ps,"002-relax-bulk-PE-AIREBO");
    quench(sl_bulk_PE,0.01*K, 200*ps, 0.01*ps,"003-cooling-bulk-PE-AIREBO");

    sl_bulk_PE.atoms.PBC(NO_PBC);
    bulkAtoms = sl_bulk_PE.atoms;
  }

  AtomsArray chainCap;

  {
    SimLoopDump sl_with_chain;
    initialize_simloop(sl_with_chain);

    size_t atomsInChainCap;

    {
      SimLoopDump sl_rebo;
      initialize_simloop_REBO_only(sl_rebo);

      place_Polyethylene_lattice(sl_rebo.atoms,1,1,c_num,fixBottomCellLayer,2,a,b,c);

      std::vector<size_t> fixedAtoms =
        sl_rebo.atoms.fixNotFixedAtoms(0,sl_rebo.atoms.size());
      {
        size_t atomsCount_wo_caps = sl_rebo.atoms.size();
        place_Polyethylene_folds(sl_rebo.atoms,1,1,c_num,2,a,b,c);
        size_t atomsCount_with_caps = sl_rebo.atoms.size();
        atomsInChainCap = atomsCount_with_caps - atomsCount_wo_caps;

        sl_rebo.enableDump();

        sl_rebo.dumpConst(0.95);
        relax(sl_rebo,0.05*ps,"010-folds-relax-REBO");

        sl_rebo.dumpConst(0.97);
        relax(sl_rebo,0.05*ps,"011-folds-relax-REBO");

        sl_rebo.dumpConst(0.99);
        relax(sl_rebo,0.05*ps,"012-folds-relax-REBO");

        sl_rebo.disableDump();

        quench(sl_rebo,0.01*K, 200*ps, 0.01*ps, "013-folds-quench-REBO");

        {
          SimLoopDump sl_airebo(sl_rebo);
          initialize_simloop(sl_airebo);

          sl_airebo.enableDump();

          sl_airebo.dumpConst(0.95);
          relax(sl_airebo,0.05*ps,"020-folds-relax-AIREBO");

          sl_airebo.dumpConst(0.97);
          relax(sl_airebo,0.05*ps,"021-folds-relax-AIREBO");

          sl_airebo.dumpConst(0.99);
          relax(sl_airebo,0.05*ps,"022-folds-relax-AIREBO");

          sl_airebo.disableDump();

          quench(sl_airebo,0.01*K, 200*ps, 0.01*ps,"023-folds-quench-AIREBO");

          sl_with_chain = sl_airebo;
        }
      }
      sl_with_chain.atoms.unfixAtoms(fixedAtoms);

      {
        yaatk::text_ofstream fomdloop("023-folds-quench-AIREBO.mdloop");
        sl_with_chain.saveToStream(fomdloop);
        fomdloop.close();
      }

      sl_with_chain.enableDump();

      sl_with_chain.dumpConst(0.95);
      relax(sl_with_chain,0.05*ps,"030-chain-relax-REBO");

      sl_with_chain.dumpConst(0.97);
      relax(sl_with_chain,0.05*ps,"031-chain-relax-REBO");

      sl_with_chain.dumpConst(0.99);
      relax(sl_with_chain,0.05*ps,"032-chain-relax-REBO");

      sl_with_chain.disableDump();

      quench(sl_with_chain,0.01*K, 200*ps, 0.01*ps,"033-chain-quench-REBO");
    }

    sl_with_chain.atoms.removeMomentum();

    REQUIRE(atomsInChainCap > 0 && atomsInChainCap < 20);

    for(size_t i = sl_with_chain.atoms.size() - atomsInChainCap;
        i < sl_with_chain.atoms.size();
        ++i)
      chainCap.push_back(sl_with_chain.atoms[i]);

    REQUIRE(chainCap.size() == atomsInChainCap);
  }

  SimLoopDump sl;
  initialize_simloop(sl);

  sl.atoms = bulkAtoms;

  place_Polyethylene_folded_chains(sl.atoms,chainCap,a_num,b_num,a,b);
  sl.heatUpEveryAtom(0.001*eV, r);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));

  {
    yaatk::text_ofstream fomdloop("033-crystal-unrelaxed.mdloop");
    sl.saveToStream(fomdloop);
    fomdloop.close();
  }

  sl.enableDump();

  sl.dumpConst(0.95);
  relax(sl,0.05*ps,"040-crystal-relax-AIREBO");

  sl.dumpConst(0.97);
  relax(sl,0.05*ps,"041-crystal-relax-AIREBO");

  sl.dumpConst(0.99);
  relax(sl,0.05*ps,"042-crystal-relax-AIREBO");

  sl.disableDump();

  quench(sl,0.01*K, 200*ps, 0.01*ps,"043-crystal-quench-AIREBO");

  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_BOX;
  sl.thermalBathGeomBox.zMin = (c_num > 3)?(c*(c_num-3)-0.5*Ao):(0.0);
  sl.thermalBathGeomBox.dBoundary = 3.0*Ao;
  sl.thermalBathGeomBox.zMinOfFreeZone = 0.0*Ao;

  sl.atoms.removeMomentum();

  gsl_rng_free(r);

  return sl;
}

}
