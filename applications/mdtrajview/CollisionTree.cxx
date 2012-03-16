/*
   Collision Tree class (implementation)

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

#include "CollisionTree.hpp"
#include <mdtk/SimLoop.hpp>

namespace xmde
{

void getAtomsFromSimLoop(const SimLoop& ml, std::vector<Atom>& atoms)
{
  atoms.clear();
  for(size_t i = 0; i < ml.atoms.size(); ++i)
  {
    atoms.push_back(*(ml.atoms[i]));
//    atoms.back().container = NULL;
  }
}

void updateAtomsFromSimLoop(const SimLoop& ml, std::vector<Atom>& atoms)
{
  REQUIRE(atoms.size()==ml.atoms.size());
  for(size_t i = 0; i < ml.atoms.size(); ++i)
  {
    atoms[i] = *(ml.atoms[i]);
//    atoms[i].container = NULL;
  }
}

void MDTrajectory_read(MDTrajectory& mdt,
		       const std::string basefile, 
		       const std::vector<std::string>& xvas)
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

    if (basefile.find("mde_init") != std::string::npos)
      ml.loadFromStream(fi);
    else
    {
      ml.loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
      ml.allowPartialLoading = true; // hack, disables essential checks
      ml.updateGlobalIndexes();
    }
    fi.close(); 
  }

  std::vector<Atom> atoms;

  for(size_t i = 0; i < xvas.size(); ++i)
  {
    TRACE(xvas[i]);
    yaatk::text_ifstream fixva(xvas[i].c_str()); 
    ml.loadFromStreamXVA(fixva);
    fixva.close(); 
    /*
      yaatk::binary_ifstream fixva(file.c_str()); 
      ml_->loadFromStreamXVA_bin(fixva);
      fixva.close(); 
    */
    if (i == 0)
      getAtomsFromSimLoop(ml,atoms);
    else
      updateAtomsFromSimLoop(ml,atoms);
    mdt[ml.simTime] = atoms;
    TRACE(ml.simTime);
    TRACE(atoms.size());
  }  
}

void MDTrajectory_read_from_SnapshotList(
  MDTrajectory& mdt,
  const std::string basefile)
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

    if (basefile.find("mde_init") != std::string::npos)
      ml.loadFromStream(fi);
    else
    {
      ml.loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
      ml.allowPartialLoading = true; // hack, disables essential checks
      ml.updateGlobalIndexes();
    }
    fi.close(); 
  }

  std::vector<Atom> atoms;
  SnapshotList shots;
  shots.loadstate();

  for(size_t i = 0; i < shots.snapshots.size(); ++i)
  {
    ml.simTime = shots.snapshots[i].first;
    for(size_t ai = 0; ai < shots.snapshots[i].second.size(); ++ai)
    {
      const SnapshotList::AtomSnapshot& as =
        shots.snapshots[i].second[ai];
      Atom& a = *ml.atoms[shots.atomsSelectedForSaving[ai]];
      as.restoreToAtom(a);
    }

    if (i == 0)
      getAtomsFromSimLoop(ml,atoms);
    else
      updateAtomsFromSimLoop(ml,atoms);
    mdt[ml.simTime] = atoms;
    TRACE(ml.simTime);
    TRACE(atoms.size());
  }  
}

Atom getNearestAtom(const Atom& a,const std::vector<Atom>& atoms)
{
  Atom an;
  Float dmin = 100000.0*Ao;

  int atomFound = 0;
  
  for(size_t i = 0; i < atoms.size(); ++i)
  {
    const Atom& ai = atoms[i];
    if (ai.globalIndex == a.globalIndex) atomFound++;
    if (ai.globalIndex == a.globalIndex) continue;
    Float d = depos(ai,a).module();
    if (d < dmin)
    {
      an = ai;
      dmin = d;
    }
  }

  REQUIRE(atomFound == 1);
  return an;
}

CollisionTree::CollisionTree(const Atom& atom, 
			     MDTrajectory::const_iterator time, 
			     const MDTrajectory& mdt)
  :a(atom),t(time->first),t1(NULL),t2(NULL)
{
  if (a.globalIndex != 10695) return;
  TRACE("*******BEGIN*******");
  TRACE(time->first/ps);
  TRACE(a.globalIndex);
  if (time == mdt.end()) return;
  Float distance = 10000.0*Ao;
  Atom an;
  MDTrajectory::const_iterator t = time;
  while (t != mdt.end())
  {
    a = t->second[a.globalIndex];
    an = getNearestAtom(a,t->second);
    ++t;
//    TRACE(an.globalIndex);
//    TRACE(distance/Ao);
//    TRACE(t->first/ps);
//    TRACE(depos(a,an).module()/Ao);
//    TRACE(distance/Ao);
    
//    TRACE(depos(a,an).module() > distance && distance < 3.0*Ao);
    if (depos(a,an).module() > distance && distance < 3.0*Ao)
      break;
    distance = depos(a,an).module();
  }
  TRACE(t != mdt.end());
  TRACE(t->first/ps);
  TRACE(an.globalIndex);
  TRACE(depos(a,an).module()/Ao);

  if (t != mdt.end())
  {
    t1 = new CollisionTree(a, t,mdt);
    t2 = new CollisionTree(an,t,mdt);
  }
  TRACE("*******END*******");
}

}
