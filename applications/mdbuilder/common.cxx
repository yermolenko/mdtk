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

  sl.initialize();
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

  sl.initialize();
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
relax(
  SimLoop& sl,
  Float forTime,
  std::string tmpDir
  )
{
  MDBUILDER_DRY_RUN_HOOK;
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.execute();
  yaatk::chdir("..");
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

  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.001*ps;
  sl.writestate();
  sl.execute();
  yaatk::chdir("..");
}

Float
mass(const AtomsContainer& atoms)
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    moleculeMass += atom.M;
  }
  return moleculeMass;
}

Vector3D
velocity(const AtomsContainer& atoms)
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    if (atom.isFixed()) continue;
    sumOfM += atom.M;
    sumOfP += atom.V*atom.M;
  };
  return sumOfP/sumOfM;
}

Vector3D
massCenter(const AtomsContainer& atoms)
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.coords*atom.M;
  };
  return sumOfP/sumOfM;
}

void
removeMomentum(AtomsContainer& atoms)
{
  Vector3D v = velocity(atoms);
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    mdtk::Atom& atom = *atoms[ai];
    if (atom.isFixed()) continue;
    atom.V -= v;
  };
}

void
addTranslationalEnergy(AtomsContainer& atoms, Float energy, Vector3D direction)
{
  direction.normalize();
  Vector3D v = velocity(atoms);
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    mdtk::Atom& a = *atoms[ai];
    if (a.isFixed()) continue;
    a.V += sqrt(2.0*energy/(mass(atoms)))*direction;
  };
}

Vector3D
geomCenter(AtomsContainer& atoms)
{
  Float clusterRadius = 0.0;
  Vector3D clusterCenter(0,0,0);
  
  for(size_t i = 0; i < atoms.size(); i++)
    clusterCenter += atoms[i]->coords;
  clusterCenter /= atoms.size();

  return clusterCenter;
}

Float
maxDistanceFrom(AtomsContainer& atoms, Vector3D point)
{
  Float clusterRadius = 0.0;
  
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Float currentDist = (atoms[i]->coords-point).module();
    clusterRadius = (currentDist>clusterRadius)?currentDist:clusterRadius;
  }

  return clusterRadius;
}

Float
radius(AtomsContainer& atoms)
{
  return maxDistanceFrom(atoms,geomCenter(atoms));
}

void
shiftToOrigin(AtomsContainer& atoms)
{
  if (atoms.size() == 0) return;

  Vector3D clusterCenter = massCenter(atoms);

  for(size_t i = 0; i < atoms.size(); i++)
     atoms[i]->coords -= clusterCenter;
}

void
shiftToPosition(AtomsContainer& atoms, Vector3D v)
{
  shiftToOrigin(atoms);
  for(size_t i = 0; i < atoms.size(); i++)
     atoms[i]->coords += v;
}

Dimensions
dimensions(AtomsContainer& atoms)
{
  Dimensions d;

  const Atom& a = *atoms[0];

  Float x_max = a.coords.x;
  Float x_min = a.coords.x;
  Float y_max = a.coords.y;
  Float y_min = a.coords.y;
  Float z_max = a.coords.z;
  Float z_min = a.coords.z;

  for(size_t i = 0; i < atoms.size(); i++)
  {
    const Atom& a = *atoms[i];

    if (a.coords.x > d.x_max)
      d.x_max = a.coords.x;
    if (a.coords.x < d.x_min)
      d.x_min = a.coords.x;

    if (a.coords.y > d.y_max)
      d.y_max = a.coords.y;
    if (a.coords.y < d.y_min)
      d.y_min = a.coords.y;

    if (a.coords.z > d.z_max)
      d.z_max = a.coords.z;
    if (a.coords.z < d.z_min)
      d.z_min = a.coords.z;
  }

  d.x_len = d.x_max - d.x_min;
  d.y_len = d.y_max - d.y_min;
  d.z_len = d.z_max - d.z_min;

  return d;
}

void
initialize_simloop(mdtk::SimLoop& sl) 
{
  setupPotentials(sl);
  sl.initialize();
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
  sl.initialize();
}

std::vector<size_t>
fixNotFixedAtoms(mdtk::AtomsContainer& atoms,
                    const size_t begin, const size_t end)
{
  std::vector<size_t> fixated;
  for(size_t i = 0; i < end; i++)
    if (!atoms[i]->isFixed())
    {
      atoms[i]->fix();
      fixated.push_back(i);
    }
  return fixated;
}

std::vector<size_t>
unfixFixedAtoms(mdtk::AtomsContainer& atoms,
                const size_t begin, const size_t end)
{
  std::vector<size_t> unfixated;
  for(size_t i = 0; i < end; i++)
    if (atoms[i]->isFixed())
    {
      atoms[i]->unfix();
      unfixated.push_back(i);
    }
  return unfixated;
}

std::vector<size_t>
fixUnfixedCHAtoms(mdtk::AtomsContainer& atoms,
           const size_t begin, const size_t end)
{
  std::vector<size_t> fixated;
  for(size_t i = 0; i < end; i++)
    if (!atoms[i]->isFixed())
      if (atoms[i]->ID == C_EL || atoms[i]->ID == H_EL)
      {
        atoms[i]->fix();
        fixated.push_back(i);
      }
  return fixated;
}

void
unfixAtoms(mdtk::AtomsContainer& atoms,
              const std::vector<size_t> fixedAtoms)
{
  for(size_t i = 0; i < fixedAtoms.size(); i++)
    atoms[fixedAtoms[i]]->unfix();
}

void
fixAtoms(mdtk::AtomsContainer& atoms,
         const std::vector<size_t> atomsToFix)
{
  for(size_t i = 0; i < atomsToFix.size(); i++)
    atoms[atomsToFix[i]]->fix();
}

void
place_Cluster(
  mdtk::SimLoop& sl,
  const mdtk::SimLoop sl_element
  )
{
  glPushMatrix();

  for(size_t i = 0; i < sl_element.atoms.size(); i++)
  {
    Atom& a = *(sl_element.atoms[i]);
    place_and_inherit(sl,a,getPosition()+a.coords);
  }

  glPopMatrix();
}

static int mdbuilder_common_dummy = 1;

}
