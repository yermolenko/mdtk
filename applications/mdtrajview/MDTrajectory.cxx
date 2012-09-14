/*
   MDTrajectory data structure (implementation)

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

#include "MDTrajectory.hpp"
#include <mdtk/SimLoop.hpp>

#include <FL/Fl.H>
#include "MainWindow.hpp"
#include "applications/common.h"

namespace xmde
{

MDSnapshot::MDSnapshot()
  :atoms(),upToDate(),accurate(),time(),name()
{
}

MDSnapshot::~MDSnapshot()
{
}

void
MDSnapshot::setCustomName(std::string srcName)
{
  std::ostringstream slabel;
  slabel << std::fixed << std::setprecision(5)
         << time/ps << " ps : "
         << srcName;

  name = slabel.str();
}

MDSnapshot::MDSnapshot(SimLoop ml)
  :atoms(),upToDate(),accurate(),time(),name()
{
  upToDate.assign(ml.atoms.size(),true);
  accurate.assign(ml.atoms.size(),true);

  atoms.clear();
  for(size_t i = 0; i < ml.atoms.size(); ++i)
    atoms.push_back(ml.atoms[i]);

  time = ml.simTime;

  std::ostringstream slabel;
  slabel << std::fixed << std::setprecision(5)
         << time/ps << " ps : "
         << "complete snapshot";

  name = slabel.str();
}

MDSnapshot::MDSnapshot(SimLoop ml, const std::string xva)
  :atoms(),upToDate(),accurate(),time(),name()
{
  ml.allowPartialLoading = true; // hack, disables essential checks

  upToDate.assign(ml.atoms.size(),true);
  accurate.assign(ml.atoms.size(),false);

  yaatk::text_ifstream fixva(xva.c_str());
  ml.loadFromStreamXVA(fixva);
  fixva.close();

  atoms.clear();
  for(size_t i = 0; i < ml.atoms.size(); ++i)
    atoms.push_back(ml.atoms[i]);

  time = ml.simTime;

  std::ostringstream slabel;
  slabel << std::fixed << std::setprecision(5)
         << time/ps << " ps : "
         << xva;

  name = slabel.str();
}

MDSnapshot::MDSnapshot(SimLoop ml, const SnapshotList& snapshots, size_t index)
  :atoms(),upToDate(),accurate(),time(),name()
{
  ml.allowPartialLoading = true; // hack, disables essential checks

  upToDate.assign(ml.atoms.size(),false);
  accurate.assign(ml.atoms.size(),false);

  ml.simTime = snapshots.snapshots[index].first;
  for(size_t ai = 0; ai < snapshots.snapshots[index].second.size(); ++ai)
  {
    const SnapshotList::AtomSnapshot& as =
      snapshots.snapshots[index].second[ai];
    size_t atomIndex = snapshots.atomsSelectedForSaving[ai];
    mdtk::Atom& a = ml.atoms[atomIndex];
    as.restoreToAtom(a);
    upToDate[atomIndex] = true;
    accurate[atomIndex] = true;
  }

  atoms.clear();
  for(size_t i = 0; i < ml.atoms.size(); ++i)
    atoms.push_back(ml.atoms[i]);

  time = ml.simTime;

  std::ostringstream slabel;
  slabel << std::fixed << std::setprecision(5)
         << ml.simTime/ps << " ps : "
         << "partial snapshot #" << std::setprecision(20) << index;

  name = slabel.str();
}

/*
void
MDSnapshot::updateFromSimLoop(
    const SimLoop& ml,
    const std::string& name)
{
  REQUIRE(atoms.size()==ml.atoms.size());
  for(size_t i = 0; i < ml.atoms.size(); ++i)
  {
    atoms[i] = ml.atoms[i];
//    atoms[i].container = NULL;
  }
}
*/

SimLoop MDTrajectory_read(
  MDTrajectory& mdt,
  const std::string basefile,
  const std::vector<std::string>& xvas
  )
{
  SimLoop ml;
  if (basefile.find("simloop.conf") != std::string::npos) 
  {
    ml.loadstate();
  }
  else
  {
    yaatk::text_ifstream fi(basefile.c_str()); 

    ml.initNLafterLoading = false;

    if (basefile.find(".xva") == std::string::npos)
      ml.loadFromStream(fi);
    else
    {
      ml.loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
      ml.allowPartialLoading = true; // hack, disables essential checks
      ml.atoms.prepareForSimulatation();
    }
    fi.close(); 
  }

  MDSnapshot s_base(ml);
  mdt[s_base.time] = s_base;

  for(size_t i = 0; i < xvas.size(); ++i)
  {
    TRACE(xvas[i]);

    MDSnapshot s_xva(ml,xvas[i]);
    mdt[s_xva.time] = s_xva;
  }

  return ml;
}

SimLoop MDTrajectory_read_from_SnapshotList(
  MDTrajectory& mdt,
  const std::string basefile
  )
{
  SimLoop ml;
  if (basefile.find("simloop.conf") != std::string::npos) 
  {
    ml.loadstate();
  }
  else
  {
    yaatk::text_ifstream fi(basefile.c_str()); 

    ml.initNLafterLoading = false;

    if (basefile.find(".xva") == std::string::npos)
      ml.loadFromStream(fi);
    else
    {
      ml.loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
      ml.allowPartialLoading = true; // hack, disables essential checks
      ml.atoms.prepareForSimulatation();
    }
    fi.close(); 
  }

  MDSnapshot s_base(ml);
  mdt[s_base.time] = s_base;

  SnapshotList shots;
  shots.loadstate();

  for(size_t i = 0; i < shots.snapshots.size(); ++i)
  {
    TRACE(i);

    MDSnapshot s_shots(ml,shots,i);
    mdt[s_shots.time] = s_shots;
  }

  return ml;
}

void MDTrajectory_read_from_basefiles(
  MDTrajectory& mdt
  )
{
  std::string basename;
  basename = "in.mde";
  if (yaatk::exists(basename.c_str()))
  {
    SimLoop ml;
    yaatk::text_ifstream fi(basename.c_str());
    ml.loadFromMDE(fi);
    ml.atoms.prepareForSimulatation();
    fi.close();
    MDSnapshot s(ml);
    s.setCustomName(basename.c_str());
    mdt[s.time] = s;
  }
  basename = "mde_final";
  if (yaatk::exists(basename.c_str()))
  {
    SimLoop ml;
    yaatk::text_ifstream fi(basename.c_str());
    ml.loadFromStream(fi);
    fi.close();
    MDSnapshot s(ml);
    s.setCustomName(basename.c_str());
    mdt[s.time] = s;
  }
  basename = "mde_init";
  if (yaatk::exists(basename.c_str()))
  {
    SimLoop ml;
    yaatk::text_ifstream fi(basename.c_str());
    ml.loadFromStream(fi);
    fi.close();
    MDSnapshot s(ml);
    s.setCustomName(basename.c_str());
    mdt[s.time] = s;
  }
  basename = "simloop.conf";
  if (yaatk::exists(basename.c_str()))
  {
    SimLoop ml;
    ml.loadstate();
    MDSnapshot s(ml);
    s.setCustomName(basename.c_str());
    mdt[s.time] = s;
  }
}

class InteractiveSimLoop : public SimLoop
{
  bool
  isItTimeToSave(Float interval)
    {
      return (simTime == 0.0 ||
              int(simTime/interval) != int((simTime - dt_prev)/interval));
    }
public:
  InteractiveSimLoop():SimLoop()
    {
//      verboseTrace = false;
      preventFileOutput = true;
    }
  InteractiveSimLoop(const SimLoop &c):SimLoop(c)
    {
//      verboseTrace = false;
      preventFileOutput = true;
    }

  void doBeforeIteration() {}
  void doAfterIteration()
    {
      if (isItTimeToSave(5e-16/**5*/))
      {
        TRACE("***Saving simulated state");
        MDSnapshot s(*this);
        {
          std::ostringstream slabel;
          slabel << std::fixed << std::setprecision(5)
                 << s.time/ps << " ps : "
                 << "simulated";
          s.name = slabel.str();
        }
        MainWindow_GlobalPtr->addMDSnapshot(s);
      }
      Fl::check();
      if (!MainWindow_GlobalPtr->btn_simulate->value())
        breakSimLoop = true;
    }
};

void MDTrajectory_add_from_simulation(
  MDTrajectory& mdt,
  SimLoop slInit,
  bool quench
  )
{
  {
    MDTrajectory::iterator t = mdt.begin();
    while(t != mdt.end())
    {
      if (t->first > slInit.simTime)
        mdt.erase(t++);
      else
        ++t;
    }
  }

  InteractiveSimLoop sl(slInit);
  setupPotentials(sl);
  if (quench)
    sl.thermalBath.zMin = -100000.0*Ao;
  if (sl.simTimeFinal < sl.simTime + 4.0*ps)
    sl.simTimeFinal = sl.simTime + 4.0*ps;
  sl.execute();
  if (sl.simTimeFinal <= sl.simTime)
    MainWindow_GlobalPtr->btn_simulate->value(0);
}

}
