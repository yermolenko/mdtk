/*
   Common routines for mdbuilder (header)

   Copyright (C) 2008, 2010, 2011, 2012 Oleksandr Yermolenko
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

#ifndef MDBUILDER_COMMON_HPP
#define MDBUILDER_COMMON_HPP

#include "../mdtrajview/VisBox.hpp"

#include <mdtk/tools.hpp>
#include <mdtk/SimLoop.hpp>

#include "../common.h"

namespace mdbuilder
{

using namespace mdtk;

#define ATOMTAG_FIXED 1<<0
#define ATOMTAG_SUBSTRATE 1<<1
#define ATOMTAG_CLUSTER   1<<2
#define ATOMTAG_NOTAG 0

inline
mdtk::Vector3D getPosition()
{
  GLdouble m[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,m);
  return Vector3D(m[3*4+0],m[3*4+1],m[3*4+2]);
}

inline
Atom*
place(ElementID id, mdtk::SimLoop& sl, Vector3D pos = getPosition())
{
  Atom* a = new Atom;
  a->ID = id;
  a->setAttributesByElementID();
  a->coords = pos;
  sl.atoms.push_back(a);
  return a;
}

inline
Atom*
place_and_inherit(mdtk::SimLoop& sl, const Atom& base,
             Vector3D pos = getPosition())
{
  Atom* a = new Atom(base);
  a->coords = pos;
  sl.atoms.push_back(a);
  return a;
}

inline
void
quench(mdtk::SimLoop& sl, 
       Float forTemp = 1.0*K,
       Float forTime = 200*ps,
       Float checkTime = 0.01*ps,
       std::string tmpDir = "_tmp-X")
{
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

inline
void
heatUp(mdtk::SimLoop& sl,
       Float forTemp = 300.0*K,
       bool uniformHeatUp = true,
       Float forTime = 200*ps,
       Float checkTime = 2.0*ps,
       std::string tmpDir = "_tmp-heatUp")
{
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

inline
void
relax(mdtk::SimLoop& sl, 
      Float forTime = 0.2*ps,
      std::string tmpDir = "_tmp-X")
{
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.execute();
  yaatk::chdir("..");
}

inline
void
relax_flush(mdtk::SimLoop& sl, 
      Float forTime = 0.2*ps,
      std::string tmpDir = "_tmp-X")
{
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

class SimLoopDump : public SimLoop
{
  Float dumpConstant;
  bool dumpEnabled;
public:
  bool isDumpEnabled() {return dumpEnabled;}
  void disableDump() {dumpEnabled = false;}
  void enableDump() {dumpEnabled = true;}
  Float dumpConst() {return dumpConstant;}
  void  dumpConst(Float dc) {dumpConstant = dc;}
  SimLoopDump(Float dc = 0.95)
    :SimLoop(),
     dumpConstant(dc),
     dumpEnabled(true)
    {
    }
  virtual void doBeforeIteration()
    {
      if (dumpEnabled)
        for(size_t i = 0; i < atoms.size(); ++i)
          atoms[i]->V *= dumpConstant;
    };
};

inline
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

inline
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

inline
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

inline
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

inline
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

inline
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

inline
Float
radius(AtomsContainer& atoms)
{
  return maxDistanceFrom(atoms,geomCenter(atoms));
}

inline
void
shiftToOrigin(AtomsContainer& atoms)
{
  if (atoms.size() == 0) return;

  Vector3D clusterCenter = massCenter(atoms);

  for(size_t i = 0; i < atoms.size(); i++)
     atoms[i]->coords -= clusterCenter;
}

inline
void
initialize_simloop(mdtk::SimLoop& sl) 
{
  setupPotentials(sl);
  sl.initialize();
}

inline
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

inline
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

inline
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

inline
void
unfixAtoms(mdtk::AtomsContainer& atoms,
              const std::vector<size_t> fixedAtoms)
{
  for(size_t i = 0; i < fixedAtoms.size(); i++)
    atoms[fixedAtoms[i]]->unfix();
}

inline
void
fixAtoms(mdtk::AtomsContainer& atoms,
         const std::vector<size_t> atomsToFix)
{
  for(size_t i = 0; i < atomsToFix.size(); i++)
    atoms[atomsToFix[i]]->fix();
}

inline
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

}

#endif
