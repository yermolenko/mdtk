/*
   Common routines for mdbuilder (header)

   Copyright (C) 2008, 2010, 2011, 2012, 2013 Oleksandr Yermolenko
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
#include <mdtk/SimLoopSaver.hpp>

#include "../common.h"

//#define MDBUILDER_DRY_RUN_HOOK return;
//#define MDBUILDER_DRY_RUN_HOOK {forTime = 0.05*ps;}
#define MDBUILDER_DRY_RUN_HOOK ;

namespace mdbuilder
{

using namespace mdtk;

inline
mdtk::Vector3D getPosition()
{
  GLdouble m[16];
  glGetDoublev(GL_MODELVIEW_MATRIX,m);
  return Vector3D(m[3*4+0],m[3*4+1],m[3*4+2]);
}

inline
Atom
place(ElementID id, mdtk::AtomsArray& atoms, Vector3D pos = getPosition())
{
  Atom a;
  a.ID = id;
  a.setAttributesByElementID();
  a.coords = pos;
  atoms.push_back(a);
  return a;
}

inline
Atom
place_and_inherit(mdtk::AtomsArray& atoms, const Atom& base, Vector3D pos = getPosition())
{
  Atom a(base);
  a.coords = pos;
  atoms.push_back(a);
  return a;
}

void
quench(
  mdtk::SimLoop& sl,
  Float forTemp = 0.1*K,
  Float forTime = 200*ps,
  Float checkTime = 0.01*ps,
  std::string tmpDir = "_tmp-X"
  );

void
quench(
  mdtk::AtomsArray& atoms,
  Float forTemp = 0.1*K,
  Float forTime = 200*ps,
  Float checkTime = 0.01*ps,
  std::string tmpDir = "_tmp-X"
  );

void
heatUp(
  mdtk::SimLoop& sl,
  Float forTemp = 300.0*K,
  bool uniformHeatUp = true,
  Float forTime = 200*ps,
  Float checkTime = 2.0*ps,
  std::string tmpDir = "_tmp-heatUp"
  );

void
heatUp(
  AtomsArray& atoms,
  Float forTemp = 300.0*K,
  bool uniformHeatUp = true,
  Float forTime = 200*ps,
  Float checkTime = 2.0*ps,
  std::string tmpDir = "_tmp-heatUp"
  );

void
relax(
  mdtk::SimLoop& sl,
  Float forTime = 0.2*ps,
  std::string tmpDir = "_tmp-X"
  );

void
relax(
  AtomsArray& atoms,
  Float forTime = 0.2*ps,
  std::string tmpDir = "_tmp-X"
  );

void
relax_flush(
  mdtk::SimLoop& sl,
  Float forTime = 0.2*ps,
  std::string tmpDir = "_tmp-X"
);

void
relax_flush(
  AtomsArray& atoms,
  Float forTime = 0.2*ps,
  std::string tmpDir = "_tmp-X"
);

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
          atoms[i].V *= dumpConstant;
    };
};

class SimLoopNoMomentums : public SimLoop
{
  void removeMomentums()
    {
      atoms.removeMomentum();
      atoms.removeAngularMomentum();
    }
public:
  SimLoopNoMomentums()
    :SimLoop()
    {
    }
  SimLoopNoMomentums(const SimLoop &c)
  :SimLoop(c)
    {
    }
  virtual void doBeforeIteration()
    {
      removeMomentums();
    }
  virtual void doAfterIteration()
    {
//      removeMomentums();
    }
};

void
initialize_simloop(mdtk::SimLoop& sl);

void
initialize_simloop_REBO_only(SimLoop& sl);

void
place_Cluster(
  mdtk::AtomsArray& sl,
  const mdtk::AtomsArray sl_element
  );

void
setupSpherical(
  mdtk::SimLoop& sl,
  const mdtk::Float centerHeightAboveSurface
  );

}

#endif
