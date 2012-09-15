/*
   Common routines for mdbuilder

   Copyright (C) 2010, 2011, 2012 Oleksandr Yermolenko
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
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
  yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.0;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  Float tb_zMin_bak = sl.thermalBath.zMin;
  sl.thermalBath.zMin = -100000.0*Ao;
  while (1)
  {
    sl.simTimeFinal += checkTime;
    sl.writestate();
    sl.execute();
    Float Temp = sl.energyKin()/(3.0/2.0*kb*sl.atoms.size());
    if (Temp < forTemp) break;
  }
  sl.thermalBath.zMin = tb_zMin_bak;
  yaatk::chdir("..");
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
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());

  yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
  yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.0;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  Float tb_zMin_bak = sl.thermalBath.zMin;
  if (uniformHeatUp)
    sl.thermalBath.zMin = -100000.0*Ao;
  const int steps = 20;
  sl.thermalBath.To = 0.0;
  while (1)
  {
    if (sl.thermalBath.To < forTemp)
      sl.thermalBath.To += forTemp/steps;

    sl.simTimeFinal += checkTime;
    sl.writestate();
    sl.execute();
    if (sl.simTimeFinal > checkTime*(steps+3)) break;
  }
  sl.thermalBath.zMin = tb_zMin_bak;
  yaatk::chdir("..");
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
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.execute();
  yaatk::chdir("..");
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
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
  yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.001*ps;
  sl.writestate();
  sl.execute();
  yaatk::chdir("..");
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
rotate(
  AtomsArray& atoms,
  Vector3D rotVector,
  Float rotAngle,
  bool aroundMassCenter)
{
  if (atoms.size() < 2 && aroundMassCenter)
    return;

  AtomsArray atoms_rotated;

  Vector3D massCenterOrig = atoms.massCenter();
  if (aroundMassCenter)
    atoms.shiftToOrigin();

  glPushMatrix();

  TRACE(rotAngle/M_PI*180.0);
  TRACE(rotVector);

  glLoadIdentity();

  glRotated(rotAngle/M_PI*180.0,rotVector.x,rotVector.y,rotVector.z);

  for(size_t i = 0; i < atoms.size(); i++)
  {
    const Atom& a = atoms[i];
    glPushMatrix();
    glTranslated(a.coords.x,a.coords.y,a.coords.z);
    place_and_inherit(atoms_rotated,a);
    glPopMatrix();
  }

  glPopMatrix();

  if (aroundMassCenter)
    atoms_rotated.shiftToPosition(massCenterOrig);

  if (aroundMassCenter)
    REQUIRE((atoms_rotated.massCenter()-massCenterOrig).module() < 0.000001*Ao);
  if (atoms.size() > 0 && rotAngle != 0.0)
  {
    TRACE(atoms_rotated[0].coords/Ao);
    TRACE(atoms[0].coords/Ao);
    REQUIRE(atoms_rotated[0].coords != atoms[0].coords);
  }

  atoms = atoms_rotated;
}

void
rotate(
  AtomsArray& atoms,
  Vector3D vecBegin,
  Vector3D vecEnd,
  bool aroundMassCenter)
{
  Float rot_angle = acos(scalarmul(vecBegin,vecEnd)/(vecBegin.module()*vecEnd.module()));
  Vector3D rot_vec = vectormul(vecBegin,vecEnd);
  TRACE(rot_angle);
  TRACE(rot_vec);
  TRACE(rot_angle == M_PI);
  TRACE(rot_vec == 0.0);
  rotate(atoms,rot_vec,rot_angle,aroundMassCenter);
}

}
