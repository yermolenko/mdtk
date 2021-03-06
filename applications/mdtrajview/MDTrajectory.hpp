/*
   MDTrajectory data structure (header file)

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

#ifndef	mde_CollisionTree_hpp
#define	mde_CollisionTree_hpp

#include "mdtk/Atom.hpp"
#include "mdtk/SnapshotList.hpp"
#include <map>
#include <vector>

namespace xmde
{

using namespace mdtk;

class MDSnapshot
{
public:
  AtomsArray atoms;
  std::vector<bool> upToDate;
  std::vector<bool> accurate;
  Float time;
  std::string name;
  MDSnapshot();
  virtual ~MDSnapshot();
  MDSnapshot(SimLoop ml);
  MDSnapshot(SimLoop ml, const std::string xva);
  MDSnapshot(SimLoop ml, const SnapshotList& snapshots, size_t index);
  void setCustomName(std::string srcName);
};

typedef std::map<Float, MDSnapshot> MDTrajectory;

SimLoop MDTrajectory_read(
  MDTrajectory& mdt,
  const std::string basefile,
  const std::vector<std::string>& xvas
  );

SimLoop MDTrajectory_read_from_SnapshotList(
  MDTrajectory& mdt,
  const std::string basefile
  );

void MDTrajectory_read_from_basefiles(
  MDTrajectory& mdt
  );

void MDTrajectory_read_from_mdloop_states(
  MDTrajectory& mdt,
  const std::vector<std::string>& mdloopStates
  );

SimLoop MDTrajectory_read_ng(
  MDTrajectory& mdt,
  const std::vector<std::string>& states,
  const std::vector<std::string>& xvas,
  bool loadPartialSnapshots = false
  );

void MDTrajectory_add_from_simulation(
  MDTrajectory& mdt,
  SimLoop slInit,
  bool quench
  );

}

#endif
