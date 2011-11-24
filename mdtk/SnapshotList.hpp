/*
   SnapshotList class header file.

   Copyright (C) 2011 Oleksandr Yermolenko
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

#ifndef	mdtk_SnapshotList_hpp
#define	mdtk_SnapshotList_hpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>

#include <mdtk/Vector3D.hpp>
#include <mdtk/SimLoop.hpp>

namespace mdtk
{

struct SnapshotList
{
  struct AtomSnapshot
  {
    IntVector3D PBC_count;
    Vector3D pos;
    Vector3D v;
    AtomSnapshot(const Atom& a):
      PBC_count(a.PBC_count),
      pos(a.coords),
      v(a.V)
      {
      }
    AtomSnapshot():
      PBC_count(),
      pos(),
      v()
      {
      }
  };
  typedef std::vector<AtomSnapshot> SelectedAtomSnapshotList;
  typedef std::pair<Float,SelectedAtomSnapshotList> TimeSnapshot;
  std::vector<size_t> atomsSelectedForSaving;
  std::vector<TimeSnapshot> snapshots;
  void initSelectedAtomList(const SimLoop& sl)
    {
      atomsSelectedForSaving.clear();
      for(size_t ai = 0; ai < sl.atoms.size(); ++ai)
        if (sl.atoms[ai]->ID == Cu_EL ||
            sl.atoms[ai]->ID == Au_EL ||
            sl.atoms[ai]->ID == Ar_EL ||
            sl.atoms[ai]->ID == Xe_EL)
          atomsSelectedForSaving.push_back(ai);
    }
  SnapshotList():
    atomsSelectedForSaving(),
    snapshots()
    {
    }
  void getSnapshot(const SimLoop& sl)
    {
      size_t prevCount = atomsSelectedForSaving.size();
      initSelectedAtomList(sl);
      REQUIRE(prevCount == 0 || prevCount == atomsSelectedForSaving.size());
      SelectedAtomSnapshotList alist;
      for(size_t i = 0; i < atomsSelectedForSaving.size(); ++i)
      {
        REQUIRE(atomsSelectedForSaving[i] < sl.atoms.size());
        alist.push_back(AtomSnapshot(*sl.atoms[atomsSelectedForSaving[i]]));
      }
      bool alreadyAccounted = false;
      for(size_t i = 0; i < snapshots.size(); ++i)
        if (snapshots[i].first == sl.simTime)
          alreadyAccounted = true;
      if (!alreadyAccounted)
        snapshots.push_back(TimeSnapshot(sl.simTime,alist));
    }
  void writestate()
    {
      yaatk::binary_ofstream state("snapshots.conf");
      size_t size;

      size = atomsSelectedForSaving.size();
      YAATK_BIN_WRITE(state,size);
      for(size_t i = 0; i < atomsSelectedForSaving.size(); ++i)
      {
        YAATK_BIN_WRITE(state,atomsSelectedForSaving[i]);
      }

      size = snapshots.size();
      YAATK_BIN_WRITE(state,size);
      for(size_t shotIndex = 0; shotIndex < snapshots.size(); ++shotIndex)
      {
        YAATK_BIN_WRITE(state,snapshots[shotIndex].first);
        size = snapshots[shotIndex].second.size();
        YAATK_BIN_WRITE(state,size);
        for(size_t ai = 0; ai < snapshots[shotIndex].second.size(); ++ai)
        {
          YAATK_BIN_WRITE(state,snapshots[shotIndex].second[ai]);
        }
      }

      state.close();
    }
  void loadstate()
    {
      yaatk::binary_ifstream state("snapshots.conf");
      size_t size;

      YAATK_BIN_READ(state,size);
      atomsSelectedForSaving.resize(size);
      for(size_t i = 0; i < atomsSelectedForSaving.size(); ++i)
      {
        YAATK_BIN_READ(state,atomsSelectedForSaving[i]);
      }

      YAATK_BIN_READ(state,size);
      snapshots.resize(size);
      for(size_t shotIndex = 0; shotIndex < snapshots.size(); ++shotIndex)
      {
        YAATK_BIN_READ(state,snapshots[shotIndex].first);
        YAATK_BIN_READ(state,size);
        snapshots[shotIndex].second.resize(size);
        for(size_t ai = 0; ai < snapshots[shotIndex].second.size(); ++ai)
        {
          YAATK_BIN_READ(state,snapshots[shotIndex].second[ai]);
        }
      }

      state.close();
    }
  void saveInText(std::string fname)
    {
      std::ofstream fo(fname.c_str());

      for(size_t i = 0; i < atomsSelectedForSaving.size(); ++i)
      {
        fo << atomsSelectedForSaving[i] << " ";
      }
      fo << "\n";

      for(size_t shotIndex = 0; shotIndex < snapshots.size(); ++shotIndex)
      {
        fo << snapshots[shotIndex].first << "\n";
        for(size_t ai = 0; ai < snapshots[shotIndex].second.size(); ++ai)
        {
          const AtomSnapshot& as = snapshots[shotIndex].second[ai];
          fo << as.PBC_count << "\n";
          fo << as.pos << "\n";
          fo << as.v << "\n";
        }
      }

      fo.close();
    }
};

}

#endif

