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

//#define MDBUILDER_DRY_RUN_HOOK return;
#define MDBUILDER_DRY_RUN_HOOK ;

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

void
quench(
  mdtk::SimLoop& sl,
  Float forTemp = 1.0*K,
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
relax(
  mdtk::SimLoop& sl,
  Float forTime = 0.2*ps,
  std::string tmpDir = "_tmp-X"
  );

void
relax_flush(
  mdtk::SimLoop& sl,
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
          atoms[i]->V *= dumpConstant;
    };
};

Float
mass(const AtomsContainer& atoms);

Vector3D
velocity(const AtomsContainer& atoms);

Vector3D
massCenter(const AtomsContainer& atoms);

void
removeMomentum(AtomsContainer& atoms);

void
addTranslationalEnergy(
  AtomsContainer& atoms,
  Float energy, Vector3D direction
  );

Vector3D
geomCenter(AtomsContainer& atoms);

Float
maxDistanceFrom(AtomsContainer& atoms, Vector3D point);

Float
radius(AtomsContainer& atoms);

void
shiftToOrigin(AtomsContainer& atoms);

void
shiftToPosition(AtomsContainer& atoms, Vector3D v);

struct Dimensions
{
  Float x_max;
  Float x_min;
  Float x_len;

  Float y_max;
  Float y_min;
  Float y_len;

  Float z_max;
  Float z_min;
  Float z_len;
};

Dimensions
dimensions(AtomsContainer& atoms);

void
initialize_simloop(mdtk::SimLoop& sl);

void
initialize_simloop_REBO_only(SimLoop& sl);

std::vector<size_t>
fixNotFixedAtoms(
  mdtk::AtomsContainer& atoms,
  const size_t begin, const size_t end
  );

std::vector<size_t>
unfixFixedAtoms(
  mdtk::AtomsContainer& atoms,
  const size_t begin, const size_t end
  );

std::vector<size_t>
fixUnfixedCHAtoms(
  mdtk::AtomsContainer& atoms,
  const size_t begin, const size_t end
  );

void
unfixAtoms(
  mdtk::AtomsContainer& atoms,
  const std::vector<size_t> fixedAtoms
  );

void
fixAtoms(
  mdtk::AtomsContainer& atoms,
  const std::vector<size_t> atomsToFix
  );

void
place_Cluster(
  mdtk::SimLoop& sl,
  const mdtk::SimLoop sl_element
  );

}

#endif
