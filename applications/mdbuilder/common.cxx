/*
   Common routines for mdbuilder

   Copyright (C) 2010, 2011, 2012, 2013 Oleksandr Yermolenko
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

#include "common.hpp"

namespace mdbuilder
{

void
quench(
  SimLoop& sl,
  Float forTemp,
  Float forTime,
  Float checkTime,
  std::string tmpDir
  )
{
  MDBUILDER_DRY_RUN_HOOK;
  yaatk::ChDir cd(tmpDir);

  bool preventFileOutput_backup = sl.preventFileOutput;
  sl.preventFileOutput = true;

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.0;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  mdtk::SimLoop::TB_GEOM_TYPE tb_type_bak = sl.thermalBathGeomType;
  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;
  while (1)
  {
    sl.simTimeFinal += checkTime;
    sl.writestate();
    {
      sl.preventFileOutput = false;
      sl.writetrajXVA();
      sl.preventFileOutput = true;
    }
    sl.execute();
    Float Temp = sl.energyKin()/(3.0/2.0*kb*sl.atoms.size());
    if (Temp < forTemp) break;
  }
  sl.thermalBathGeomType = tb_type_bak;

  {
    sl.preventFileOutput = false;
    sl.writestate();
    sl.preventFileOutput = true;
  }

  sl.preventFileOutput = preventFileOutput_backup;
}

void
quench(
  AtomsArray& atoms,
  Float forTemp,
  Float forTime,
  Float checkTime,
  std::string tmpDir
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  sl.atoms = atoms;

  quench(sl,forTemp,forTime,checkTime,tmpDir);

  atoms = sl.atoms;
}

void
heatUp(
  SimLoop& sl,
  Float forTemp,
  bool uniformHeatUp,
  Float forTime,
  Float checkTime,
  std::string tmpDir
  )
{
  MDBUILDER_DRY_RUN_HOOK;
  yaatk::ChDir cd(tmpDir);

  bool preventFileOutput_backup = sl.preventFileOutput;
  sl.preventFileOutput = true;

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.0;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  mdtk::SimLoop::TB_GEOM_TYPE tb_type_bak = sl.thermalBathGeomType;
  if (uniformHeatUp)
    sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;
  const int steps = 20;
  sl.thermalBathCommon.To = 0.0;
  while (1)
  {
    if (sl.thermalBathCommon.To < forTemp)
      sl.thermalBathCommon.To += forTemp/steps;

    sl.simTimeFinal += checkTime;
    sl.writestate();
    {
      sl.preventFileOutput = false;
      sl.writetrajXVA();
      sl.preventFileOutput = true;
    }
    sl.execute();
    if (sl.simTimeFinal > checkTime*(steps+3)) break;
  }
  sl.thermalBathGeomType = tb_type_bak;

  {
    sl.preventFileOutput = false;
    sl.writestate();
    sl.preventFileOutput = true;
  }

  sl.preventFileOutput = preventFileOutput_backup;
}

void
heatUp(
  AtomsArray& atoms,
  Float forTemp,
  bool uniformHeatUp,
  Float forTime,
  Float checkTime,
  std::string tmpDir
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  sl.atoms = atoms;

  heatUp(sl, forTemp, uniformHeatUp, forTime, checkTime, tmpDir);

  atoms = sl.atoms;
}

void
relax(
  SimLoop& sl,
  Float forTime,
  std::string tmpDir
  )
{
  MDBUILDER_DRY_RUN_HOOK;
  yaatk::ChDir cd(tmpDir);
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.execute();
  sl.writestate();
}

void
relax(
  AtomsArray& atoms,
  Float forTime,
  std::string tmpDir
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  sl.atoms = atoms;

  relax(sl, forTime, tmpDir);

  atoms = sl.atoms;
}

void
relax_flush(
  SimLoop& sl,
  Float forTime,
  std::string tmpDir
  )
{
  MDBUILDER_DRY_RUN_HOOK;
  yaatk::ChDir cd(tmpDir);

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.001*ps;
  sl.writestate();
  sl.execute();
}

void
relax_flush(
  AtomsArray& atoms,
  Float forTime,
  std::string tmpDir
  )
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  sl.atoms = atoms;

  relax_flush(sl, forTime, tmpDir);

  atoms = sl.atoms;
}

void
randomizeVelocities(SimLoop& sl, Float energyPerAtom)
{
  const gsl_rng_type* T;
  gsl_rng* r;

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc(T);
  REQUIRE(r != NULL);

  gsl_rng_set(r, 123);

  REQUIRE(gsl_rng_min(r) == 0);
  REQUIRE(gsl_rng_max(r) > 1000);

  sl.heatUpEveryAtom(energyPerAtom, r);

  gsl_rng_free(r);
}

void
initialize_simloop(mdtk::SimLoop& sl)
{
  setupPotentials(sl);
}

void
initialize_simloop_REBO_only(SimLoop& sl)
{
  using namespace mdtk;

  REQUIRE(sl.fpot.potentials.size()==0);

  mdtk::FGeneral* pot = NULL;

  pot = new mdtk::REBO();
  sl.fpot.addPotential(pot);
/*
  pot = new mdtk::AIREBO((CREBO*)pot);
  simloop.fpot.addPotential(pot);

  pot = new mdtk::ETors();
  simloop.fpot.addPotential(pot);
*/
/*
  pot = new mdtk::Brenner(Brenner::POTENTIAL2);
  simloop.fpot.addPotential(pot);
*/
}

void
place_Cluster(
  AtomsArray& sl,
  const AtomsArray sl_element
  )
{
  glPushMatrix();

  for(size_t i = 0; i < sl_element.size(); i++)
  {
    const Atom& a = sl_element[i];
    place_and_inherit(sl,a,getPosition()+a.coords);
  }

  glPopMatrix();
}

void
setupSpherical(
  mdtk::SimLoop& sl,
  const mdtk::Float centerHeightAboveSurface
  )
{
//  sl.atoms.unfoldPBC(); // already called in atoms.PBC(NO_PBC)
  Vector3D PBC = sl.atoms.PBC();
  sl.atoms.PBC(NO_PBC); // essential

  mdtk::AtomsArray::Dimensions dim = sl.atoms.dimensions();
  Vector3D sphereCenter = dim.center();
  sphereCenter.z = dim.z_min - centerHeightAboveSurface;

  Float minLateralSize1 = min2(PBC.x,PBC.y);
  Float minLateralSize2 = min2(dim.x_len,dim.y_len);
  Float minLateralSize = min2(minLateralSize1,minLateralSize2);
  Float radius = sqrt(SQR(minLateralSize/2.0) + SQR(centerHeightAboveSurface));

  {
    TRACE(centerHeightAboveSurface/Ao);
    TRACE(dim.z_len/Ao);
    TRACE(radius/Ao);
    REQUIRE(centerHeightAboveSurface + dim.z_len >= radius);
  }

#define CHECK_RADIUS(axis)                                   \
  {                                                          \
    Vector3D edgeCenter = dim.center();                      \
    edgeCenter.axis = dim.axis##_max;                        \
    edgeCenter.z = dim.z_min;                                \
    TRACE("dim max check");                                  \
    PRINT(#axis"\n");                                        \
    TRACE((sphereCenter - edgeCenter).module()/Ao);          \
    TRACE(radius/Ao);                                        \
    REQUIRE((sphereCenter - edgeCenter).module() >= radius); \
  }

  CHECK_RADIUS(x);
  CHECK_RADIUS(y);

#define CHECK_RADIUS_USING_PBC(axis)                         \
  {                                                          \
    Vector3D edgeCenter = PBC/2.0;                           \
    edgeCenter.axis += PBC.axis/2.0;                         \
    edgeCenter.z = dim.z_min;                                \
    TRACE("dim PBC check");                                  \
    PRINT(#axis"\n");                                        \
    TRACE((sphereCenter - edgeCenter).module()/Ao);          \
    TRACE(radius/Ao);                                        \
    REQUIRE((sphereCenter - edgeCenter).module() >= radius); \
  }

  CHECK_RADIUS_USING_PBC(x);
  CHECK_RADIUS_USING_PBC(y);

  Float fixedSphereRadius = radius - 5.5*Ao;
  REQUIRE(fixedSphereRadius > 0);
  AtomsArray allAtoms = sl.atoms;
  sl.atoms.clear();
  for(size_t i = 0; i < allAtoms.size(); ++i)
  {
    if ((allAtoms[i].coords - sphereCenter).module() < radius)
    {
      sl.atoms.push_back(allAtoms[i]);
      if ((sl.atoms.back().coords - sphereCenter).module() >= fixedSphereRadius)
        sl.atoms.back().fix();
    }
  }

  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_SPHERE;
  sl.thermalBathGeomSphere.center = sphereCenter;
  sl.thermalBathGeomSphere.radius = fixedSphereRadius - 3.0*Ao;

  {
    size_t fixedAtomsCount = 0;
    for(size_t i = 0; i < sl.atoms.size(); ++i)
      if (sl.atoms[i].isFixed())
        ++fixedAtomsCount;
    TRACE(fixedAtomsCount);
  }
}

}
