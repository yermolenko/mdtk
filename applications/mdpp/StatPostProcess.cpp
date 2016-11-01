/* 
   Molecular dynamics postprocessor, main classes

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015,
   2016 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include "StatPostProcess.hpp"
#include "NeighbourList.hpp"

#include <algorithm>

//using namespace std;
#include <fstream>
#include <mdtk/tools.hpp>
#include <mdtk/SnapshotList.hpp>


namespace mdepp
{

std::string initialWorkingDirectory = yaatk::getcwd() + DIR_DELIMIT_STR;

bool
StatPostProcess::ProcessAll(const ClassicMolecule&)
{
  return true;
}

bool
StatPostProcess::ProcessProjectile(const ClassicMolecule& mol)
{
  return mol.hasProjectileAtoms();
}

bool
StatPostProcess::ProcessCluster(const ClassicMolecule& mol)
{
  return mol.hasClusterAtoms();
}

bool
StatPostProcess::ProcessFullerene(const ClassicMolecule& mol)
{
  return mol.hasFullereneAtoms();
}

bool
StatPostProcess::ProcessSubstrate(const ClassicMolecule& mol)
{
  return mol.hasSubstrateAtoms();
}

bool
StatPostProcess::ProcessClusterAndSubstrate(const ClassicMolecule& mol)
{
//  return mol.hasSubstrateAtoms() || mol.hasClusterAtoms();
  return mol.hasOnlySubstrateOrClusterAtoms();
}

std::string
StatPostProcess::FProcessClassicMoleculeToString(FProcessClassicMolecule fpm)
{
  if (fpm == ProcessAll) return "All";
  if (fpm == ProcessProjectile) return "Projectile";
  if (fpm == ProcessCluster) return "Cluster";
  if (fpm == ProcessFullerene) return "Fullerene";
  if (fpm == ProcessSubstrate) return "Substrate";
  if (fpm == ProcessClusterAndSubstrate) return "ClusterAndSubstrate";

  throw Exception("Unknown FProcessClassicMolecule filter");
}

StatPostProcess::TrajData::TrajData() :
  trajDir(),
  molecules(),
  trajProjectile(),
  trajCluster(),
  PBC(),
  SPOTTED_DISTANCE(-5.0*mdtk::Ao)
{
}

void
StatPostProcess::TrajData::saveToStream(std::ostream& os) const
{
  os << trajDir << "\n";
  os << molecules.size() << "\n";
  for(size_t i = 0; i < molecules.size(); i++)
    molecules[i].saveToStream(os);

  {
    os << trajProjectile.size() << "\n";
    std::map< Float, AtomGroup >::const_iterator i;
    for( i = trajProjectile.begin(); i != trajProjectile.end() ; ++i )
      os << i->first << "\t " << i->second << "\n";
  }

  {
    os << trajCluster.size() << "\n";
    std::map< Float, AtomGroup >::const_iterator i;
    for( i = trajCluster.begin(); i != trajCluster.end() ; ++i )
      os << i->first << "\t " << i->second << "\n";
  }

  os << PBC << "\n";

  os << SPOTTED_DISTANCE << "\n";
}

void
StatPostProcess::TrajData::loadFromStream(std::istream& is)
{
  is >> trajDir;
  size_t sz, i;
  is >> sz;
  molecules.resize(sz);
  for(i = 0; i < molecules.size(); i++)
    molecules[i].loadFromStream(is);

  is >> sz;
  for(i = 0; i < sz; ++i)
  {
    Float t;
    AtomGroup f;
    is >> t >> f;
    trajProjectile[t] = f;
  }

  is >> sz;
  for(i = 0; i < sz; ++i)
  {
    Float t;
    AtomGroup f;
    is >> t >> f;
    trajCluster[t] = f;
  }

  is >> PBC;

  is >> SPOTTED_DISTANCE;
}

void
StatPostProcess::TrajData::execute(mdtk::SimLoop& state)
{
    yaatk::ChDir cd(initialWorkingDirectory +
                    DIR_DELIMIT_STR +
                    trajDir, false);

    {
      state.initNLafterLoading = false;

      TRACE(yaatk::getcwd());

      mdtk::SimLoopSaver mds(state);
      REQUIRE(mds.loadIterationLatest() & mdtk::SimLoopSaver::LOADED_R);

      state.atoms.prepareForSimulatation();
      setTags(state);
      {
        REQUIRE(mds.listIterations().size() > 0);
        cout << "State " << *mds.listIterations().rbegin() << " loaded." << std::endl;
      }
    }
    {
      SimLoop::Check& check = state.check;
      Float Eo_plus_Eb = check.initialEnergy + check.energyTransferredFromBath;
      Float dE = check.currentEnergy - Eo_plus_Eb;
      Float dE_by_Eo_plus_Eb = dE/Eo_plus_Eb;
      if (fabs(dE_by_Eo_plus_Eb) > 0.01)
      {
        cerr << "Trajectory " << trajDir
             << " has bad energy conservation."
             << " Handle it separately." << endl;
        throw mdepp::BadTrajectoryException();
      }
    }
    {
      if (state.simTimeFinal > state.simTime)
      {
        cerr << "Trajectory " << trajDir
             << " seems to be unfinished."
             << " Handle it separately." << endl;
        throw mdepp::BadTrajectoryException();
      }
    }

    NeighbourList nl(state.atoms);

    Molecule projectile;
    {
      projectile.buildByTag(state,ATOMTAG_PROJECTILE);
      trajProjectile[state.simTime] = projectile;

      projectile.update(state);
      trajProjectile[state.simTime] = projectile;
    }

    Molecule cluster;
    {
      cluster.buildByTag(state,ATOMTAG_CLUSTER);
      trajCluster[state.simTime] = cluster;

      cluster.update(state);
      trajCluster[state.simTime] = cluster;
    }

    PBC = state.atoms.PBC();
    TRACE(PBC/mdtk::Ao);

//    TRACE(getAboveSpottedHeight(*state));

    buildSputteredClassicMolecules(state,STATE_FINAL,nl);
//    buildClusterDynamics(*state,trajIndex,STATE_FINAL,nl);
//    buildProjectileDynamics(*state,trajIndex,STATE_FINAL);

//  if (molecules.size() > 0)

    if (0)
    {
      SnapshotList snapshots;
      snapshots.loadstate();

      {
        if (state.simTimeFinal -
            snapshots.snapshots[snapshots.snapshots.size()-1].first >
            (5e-16*50*4)*2)
        {
          cerr << "Trajectory " << trajDir
               << " seems to have incomplete partial snapshot info."
               << " Handle it separately." << endl;
          throw mdepp::BadTrajectoryException();
        }

        if (snapshots.snapshots[0].first > (5e-16*5*4)*2)
        {
          cerr << "Trajectory " << trajDir
               << " seems to have incomplete partial snapshot info."
               << " Handle it separately." << endl;
          throw mdepp::BadTrajectoryException();
        }
      }

      for(size_t index = 0; index < snapshots.snapshots.size(); ++index)
      {
        state.simTime = snapshots.snapshots[index].first;
        for(size_t ai = 0; ai < snapshots.snapshots[index].second.size(); ++ai)
        {
          const SnapshotList::AtomSnapshot& as =
            snapshots.snapshots[index].second[ai];
          size_t atomIndex = snapshots.atomsSelectedForSaving[ai];
          mdtk::Atom& a = state.atoms[atomIndex];
          as.restoreToAtom(a);
          // upToDate[atomIndex] = true;
          // accurate[atomIndex] = true;
        }
        {
          projectile.update(state);
          trajProjectile[state.simTime] = projectile;
        }
        {
          cluster.update(state);
          trajCluster[state.simTime] = cluster;
        }
      }
    }

#if 0

//    if (0)
    {
      mdtk::SimLoop mde_init;
      mde_init.allowToFreePotentials = true;
//      mde_init.allowToFreeAtoms = true;
      setupPotentials(mde_init);
      std::string mde_init_filename = trajDir+"mde_init";
      cout << "Loading state " << mde_init_filename << std::endl;
      yaatk::text_ifstream fi(mde_init_filename.c_str());
      mde_init.initNLafterLoading = false;
      mde_init.loadFromStream(fi);
      mde_init.atoms.prepareForSimulatation();
      setTags(mde_init);
      NeighbourList nl(mde_init.atoms);
      fi.close();

      cout << "State " << mde_init_filename << " loaded." << std::endl;

      buildSputteredClassicMolecules(mde_init,trajIndex,STATE_INIT,nl);
//      buildClusterDynamics(mde_init,trajIndex,STATE_INIT,nl);
//      buildProjectileDynamics(mde_init,trajIndex,STATE_INIT);

      {
        projectile.update(mde_init);
        trajProjectile[mde_init.simTime] = projectile;
      }

      {
        cluster.update(mde_init);
        trajCluster[mde_init.simTime] = cluster;
      }
    }

//  if (molecules.size() > 0)
    {
      std::vector<std::string> interStates;
      findIntermediateStates(trajDir,interStates);

      mdtk::SimLoop mde_inter;
      mde_inter.allowToFreePotentials = true;
//      mde_inter.allowToFreeAtoms = true;
      setupPotentials(mde_inter);

      std::string trajFinalName = trajDir+"mde_init";
      cout << "Loading state " << trajFinalName << std::endl;
      yaatk::text_ifstream fi(trajFinalName.c_str());
      mde_inter.initNLafterLoading = false;
      mde_inter.loadFromStream(fi);
      fi.close();

      for(int stateIndex = interStates.size()-1; stateIndex >= 0; stateIndex--)
      {
//	if (stateIndex != interStates.size()-1) continue;
	std::string mde_inter_filename = trajDir+interStates[stateIndex];

	{
	  //to remove .GZ simply resize
	  mde_inter_filename.resize(mde_inter_filename.size()-3);
	}

	TRACE(mde_inter_filename);
	yaatk::text_ifstream fi(mde_inter_filename.c_str());
	mde_inter.loadFromStreamXVA(fi);
	mde_inter.atoms.prepareForSimulatation();
	setTags(mde_inter);
	NeighbourList nl(mde_inter.atoms);
	fi.close();

	buildSputteredClassicMolecules(mde_inter,trajIndex,STATE_INTER,nl);
//	buildClusterDynamics(mde_inter,trajIndex,STATE_INTER,nl);
//	buildProjectileDynamics(mde_inter,trajIndex,STATE_INTER);

//        if (bombardingWithFullerene)
        {
          projectile.update(mde_inter);
          trajProjectile[mde_inter.simTime] = projectile;
        }
        {
          cluster.update(mde_inter);
          trajCluster[mde_inter.simTime] = cluster;
        }
      }
    }

#endif
}

void
StatPostProcess::TrajData::buildSputteredClassicMolecules(
  mdtk::SimLoop& state,
  StatPostProcess::TrajData::StateType s,
  NeighbourList& nl)
{
  if (s == STATE_FINAL)
  {
    cout << "Building molecules for state ..." << "\n";
    for(size_t atomIndex = 0; atomIndex < state.atoms.size(); atomIndex++)
    {
      mdtk::Atom &atom = state.atoms[atomIndex];
      if (atom.coords.z < SPOTTED_DISTANCE)
      {
	bool account_atom = true;
	for(size_t mi = 0; mi < molecules.size(); mi++)
	{
	  if (molecules[mi].hasAtom(atom))
	  {
	    account_atom = false;
	    break;
	  }
	}
	if (account_atom)
	{
	  ClassicMolecule molecule;
	  molecule.buildFromAtom(atom,nl,SPOTTED_DISTANCE);
	  if (molecule.atoms.size() > 0 && molecule.getVelocity().z < 0.0)
	  {
	    cout << "Adding molecule." << "\n";
	    molecules.push_back(molecule);
	  }
	}
      }
    }
    std::sort(molecules.begin(),molecules.end());
  }
  else if (s == STATE_INIT)
  {
    mdtk::SimLoop* mde_init = &state;
    for(size_t mi = 0; mi < molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < molecules[mi].atoms.size(); ai++)
      {
        const mdtk::Atom& atom = molecules[mi].atoms[ai];
        const mdtk::Atom& atom_init = mde_init->atoms[atom.globalIndex];
        molecules[mi].atoms_init.push_back(atom_init);
        REQUIRE(molecules[mi].atoms[ai].globalIndex ==
		molecules[mi].atoms_init[ai].globalIndex);
      }
      REQUIRE(molecules[mi].atoms.size() ==
	      molecules[mi].atoms_init.size());
    }
  }
  else if (s == STATE_INTER)
  {
    mdtk::SimLoop* mde_inter = &state;
    for(size_t mi = 0; mi < molecules.size(); mi++)
    {
      bool formedNowOrEarlier = true;
      bool escapedNowOrEarlier = true;
      ClassicMolecule& molecule = molecules[mi];
      for(size_t ai = 0; ai < molecule.atoms.size(); ai++)
      {
	mdtk::Atom& atom_i = molecule.atoms[ai];
	mdtk::Atom& atom_i_inter = mde_inter->atoms[atom_i.globalIndex];
	for(size_t aj = 0; aj < molecule.atoms.size(); aj++)
	  if (ai != aj)
	  {
	    mdtk::Atom& atom_j = molecule.atoms[aj];
	    mdtk::Atom& atom_j_inter = mde_inter->atoms[atom_j.globalIndex];
	    Float distance_ij =
	      depos(atom_i,atom_j).module();
	    Float distance_ij_inter =
	      depos(atom_i_inter,atom_j_inter).module();

	    if (
	      (
		molecule.Rc(atom_i      ,atom_j      ) >= distance_ij &&
		molecule.Rc(atom_i_inter,atom_j_inter) <  distance_ij_inter
		)||
	      (
		molecule.Rc(atom_i      ,atom_j      ) <  distance_ij &&
		molecule.Rc(atom_i_inter,atom_j_inter) >= distance_ij_inter
		)
	      )
	    {
	      formedNowOrEarlier = false;
	      break;
	    }
	  }
      }

//        if (formedNowOrEarlier)
      for(size_t ai = 0; ai < molecule.atoms.size(); ai++)
      {
	const mdtk::Atom& atom_i = molecule.atoms[ai];
	ClassicMolecule molecule_inter;
	molecule_inter.
	  buildFromAtom(mde_inter->atoms[atom_i.globalIndex],
			nl,SPOTTED_DISTANCE);

	if (molecule_inter.atoms.size() == 0)
	{
	  escapedNowOrEarlier = false;
	  break;
	}
	if (!(molecule_inter.getVelocity().z < 0.0))
	{
	  escapedNowOrEarlier = false;
	  break;
	}
	for(size_t ak_inter = 0; ak_inter < molecule_inter.atoms.size(); ak_inter++)
	{
	  const mdtk::Atom& atom_k_inter = molecule_inter.atoms[ak_inter];
	  if (atom_k_inter.coords.z >= SPOTTED_DISTANCE)
	  {
	    escapedNowOrEarlier = false;
	    break;
	  }
	}
      }
      if (formedNowOrEarlier)
	molecule.formationTime  =
	  (mde_inter->simTime>0.01*mdtk::ps)?mde_inter->simTime:0.0*mdtk::ps;
      if (escapedNowOrEarlier)
	molecule.escapeTime =
	  (mde_inter->simTime>0.01*mdtk::ps)?mde_inter->simTime:0.0*mdtk::ps;
    }
  }
}

void
StatPostProcess::TrajData::buildDummyDynamics(
  mdtk::SimLoop& state,
  StatPostProcess::TrajData::StateType s)
{
  if (s == STATE_FINAL)
  {
  }
  else if (s == STATE_INIT)
  {
    mdtk::SimLoop* mde_init = &state;
  }
  else if (s == STATE_INTER)
  {
    mdtk::SimLoop* mde_inter = &state;
  }
}

std::map <Float,AtomGroup>
StatPostProcess::TrajData::getTrajWithPartialSnapshots(
      const std::map<Float,AtomGroup>& trajectorySeed) const
{
  std::map<Float,AtomGroup> trajectory = trajectorySeed;

  {
    yaatk::ChDir cd(initialWorkingDirectory +
                    DIR_DELIMIT_STR +
                    trajDir, false);

    SnapshotList snapshots;
    snapshots.loadstate();

    {
      if (/*state.simTimeFinal*/ 10.0*mdtk::ps -
          snapshots.snapshots[snapshots.snapshots.size()-1].first >
          (5e-16*50*4)*2)
      {
        cerr << "Trajectory " << trajDir
             << " seems to have incomplete partial snapshot info."
             << " Handle it separately." << endl;
        throw mdepp::BadTrajectoryException();
      }

      if (snapshots.snapshots[0].first > (5e-16*5*4)*2)
      {
        cerr << "Trajectory " << trajDir
             << " seems to have incomplete partial snapshot info."
             << " Handle it separately." << endl;
        throw mdepp::BadTrajectoryException();
      }
    }

    for(size_t index = 0; index < snapshots.snapshots.size(); ++index)
    {
      mdtk::Float simTime = snapshots.snapshots[index].first;
      AtomGroup m = trajectorySeed.begin()->second;
      {
        m.update(snapshots.snapshots[index].second,snapshots);
        trajectory[simTime] = m;
      }
    }
  }

  return trajectory;
}

void
StatPostProcess::TrajData::setSpottedDistanceFromInit()
{
  using mdtk::Exception;

  // REQUIRE(trajData.size() > 0);

  SimLoop mdbase;
  mdtk::SimLoopSaver mds(mdbase);

//  setupPotentials(mdloop);
  REQUIRE(mds.load("base") & mdtk::SimLoopSaver::LOADED_R);

  setTags(mdbase);
  mdtk::AtomsArray& mde_init_atoms = mdbase.atoms;
  Float minInitZ = 1000000.0*mdtk::Ao;
  for(size_t ai = 0; ai < mde_init_atoms.size(); ai++)
  {
    const mdtk::Atom& init_atom = mde_init_atoms[ai];
    REQUIRE(!isProjectileAtom(init_atom) || ai == mde_init_atoms.size()-1);
    if (init_atom.coords.z < minInitZ)
      minInitZ = init_atom.coords.z;
  }
  SPOTTED_DISTANCE = minInitZ - 0.05*mdtk::Ao;
  SPOTTED_DISTANCE = -4.55*mdtk::Ao;
  SPOTTED_DISTANCE = -1.0*mdtk::Ao;
  TRACE(SPOTTED_DISTANCE/mdtk::Ao);
}

int
StatPostProcess::TrajData::getAboveSpottedHeight(mdtk::SimLoop& state) const
{
  int spotted;
  {
    spotted = 0;
    for(size_t atomIndex = 0; atomIndex < state.atoms.size(); atomIndex++)
    {
      mdtk::Atom &atom = state.atoms[atomIndex];
      if (atom.coords.z < SPOTTED_DISTANCE)
      {
        spotted++;
      }
    }
  }
  return spotted;
}

bool
StatPostProcess::TrajData::hasIntactClusterSputtering() const
{
  Id id(yaatk::extractItemFromEnd(trajDir,2));

  bool hasClusterAmongSputtered = false;

  for(size_t mi = 0; mi < molecules.size(); mi++)
  {
    if (molecules[mi].atoms.size() < id.clusterSize)
      continue;
    size_t clusterAtomsCount = 0;
    for(size_t ai = 0; ai < molecules[mi].atoms.size(); ai++)
    {
      if (molecules[mi].atoms[ai].ID == id.clusterElement)
        ++clusterAtomsCount;
    }
    REQUIRE(clusterAtomsCount <= id.clusterSize);
    if (clusterAtomsCount == id.clusterSize)
      hasClusterAmongSputtered = true;
  }

  return hasClusterAmongSputtered;
}

std::map <Float,Float>
StatPostProcess::TrajData::plot_Ekin_t(const std::map<Float,AtomGroup>& trajectory)
{
  std::map <Float,Float> plot;

  std::map<Float,AtomGroup>::const_iterator i;
  for(i = trajectory.begin(); i != trajectory.end(); ++i)
    plot[i->first] = i->second.kineticEnergy();

  return plot;
}

int
StatPostProcess::TrajData::yield(FProcessClassicMolecule fpm) const
{
  int sputtered = 0;
  for(size_t mi = 0; mi < molecules.size(); mi++)
  {
    const ClassicMolecule& mol = molecules[mi];
    if (!fpm(mol)) continue;
    sputtered += mol.atoms.size();
  }
  return sputtered;
}

bool
StatPostProcess::TrajFilterProcessAll(const TrajData&)
{
  return true;
}

bool
StatPostProcess::TrajFilterProcessIntactClusterOnly(const TrajData& td)
{
  return td.hasIntactClusterSputtering();
}

StatPostProcess::TrajFilter
StatPostProcess::trajFilter =
  StatPostProcess::TrajFilterProcessAll;

StatPostProcess::Id::Id(std::string s)
  :str(s),
   clusterElement(),
   clusterSize(),
   ionElement(),
   ionEnergy()
{
  {
    std::string clusterElementString = s.substr(0,2);
    bool elementRecognized = false;
    if (clusterElementString == "Cu")
    {
      clusterElement = Cu_EL; elementRecognized = true;
    }
    if (clusterElementString == "Au")
    {
      clusterElement = Au_EL; elementRecognized = true;
    }
    REQUIRE(elementRecognized);
  }

  {
    istringstream is(s.substr(2,3));
    is >> clusterSize;
  }

  {
    size_t istart = s.find("_by_");
    REQUIRE(istart != std::string::npos);
    istart += 4;

    std::string ionElementString = s.substr(istart,2);
    bool elementRecognized = false;
    if (ionElementString == "Ar")
    {
      ionElement = Ar_EL; elementRecognized = true;
    }
    if (ionElementString == "Xe")
    {
      ionElement = Xe_EL; elementRecognized = true;
    }
    REQUIRE(elementRecognized);
  }

  {
    istringstream is(s.substr(s.find("_by_")+7,4));
    is >> ionEnergy;
  }

  yaatk::VerboseOutput suppressVerboseOutput(false);
  TRACE(str);
  TRACE(ElementIDtoString(clusterElement));
  TRACE(clusterSize);
  TRACE(ElementIDtoString(ionElement));
  TRACE(ionEnergy);
  REQUIRE(str.size()>1);
  REQUIRE(*str.begin()=='C' || *str.begin()=='A');
  REQUIRE(*(str.end()-1)=='V');
}

StatPostProcess::Id::Id() :
  str(),
  clusterElement(),
  clusterSize(),
  ionElement(),
  ionEnergy()
{
}

void
StatPostProcess::Id::saveToStream(std::ostream& os) const
{
  os << str << "\n";
  os << clusterElement << "\n";
  os << clusterSize << "\n";
  os << ionElement << "\n";
  os << ionEnergy << "\n";
}

void
StatPostProcess::Id::loadFromStream(std::istream& is)
{
  is >> str;
  int ID;
  is >> ID; clusterElement = ElementID(ID);
  is >> clusterSize;
  is >> ID; ionElement = ElementID(ID);
  is >> ionEnergy;
}

std::string StatPostProcess::cacheDir = initialWorkingDirectory + DIR_DELIMIT_STR + "cache";

std::string
StatPostProcess::getCacheFilename(std::string trajDir)
{
  std::string trajDirID;
  for(size_t itemFromEnd = 0; itemFromEnd < 3; ++itemFromEnd)
  {
    std::string pathItem = yaatk::extractItemFromEnd(trajDir,itemFromEnd);
    REQUIRE(pathItem != ".");
    REQUIRE(pathItem != "..");
    REQUIRE(pathItem != "...");
    REQUIRE(pathItem.find(DIR_DELIMIT_STR) == std::string::npos);
    REQUIRE(pathItem.find("/") == std::string::npos);
    REQUIRE(pathItem.find("\\") == std::string::npos);
    trajDirID = pathItem + "_" + trajDirID;
  }
  return cacheDir + DIR_DELIMIT_STR + trajDirID;
}

std::map<std::string,StatPostProcess::TrajData> StatPostProcess::ramCache;

mdtk::SimLoop StatPostProcess::stateTemplate;

StatPostProcess::StatPostProcess(const std::vector<std::string> trajdirs)
  :testProcessClassicMolecule(&ProcessAll),
   trajData(),
   id(yaatk::extractItemFromEnd(trajdirs[0],1))
{
  yaatk::VerboseOutput suppressVerboseOutput(false);

  instanceCounter++;
  REQUIRE(instanceCounter <= 1);

#ifdef __linux__
  if (memUsage() > ramSize()/2.0)
  {
    FTRACE(memUsage());
    FTRACE(ramSize());
    FPRINT("Clearing RAM cache.\n");
    ramCache.clear();
  }
#else
  if (ramCache.size() > 10000)
    ramCache.clear();
#endif

  stateTemplate.atoms.clear();

  REQUIRE(trajdirs.size() > 0);
  REQUIRE(trajdirs.size() < 5);
  for(size_t tdi = 0; tdi < trajdirs.size(); ++tdi)
  {
    std::string trajsetDir = trajdirs[tdi];

    if (tdi > 0)
    {
      TRACE(trajdirs[0]);
      TRACE(trajsetDir);
      Id id0(yaatk::extractItemFromEnd(trajdirs[0],1));
      Id idi(yaatk::extractItemFromEnd(trajsetDir,1));
      REQUIRE(id0.clusterElement == idi.clusterElement);
      REQUIRE(id0.clusterSize == idi.clusterSize);
      REQUIRE(id0.ionElement == idi.ionElement);
      REQUIRE(id0.ionEnergy == idi.ionEnergy);
    }

    TRACE(yaatk::getcwd());
    TRACE(trajsetDir);

    using mdtk::Exception;

    std::vector<std::string> savedStateNames;
    mdepp::FProcessTrajectory fpt = mdepp::trajProcess_Custom2;
    {
      yaatk::ChDir cd(
        initialWorkingDirectory,
        false);
      mdepp::addTrajDirNames(savedStateNames,trajsetDir.c_str(),fpt);
    }
    // std::sort(savedStateNames.begin(),savedStateNames.end());

    // std::vector<_SavedStateSortStruct> sorted;
    // for(size_t i = 0; i < savedStateNames.size(); i++)
    //   sorted.push_back(savedStateNames[i]);

    // sort(sorted.begin(), sorted.end());

    // for(size_t i = 0; i < sorted.size(); i++)
    for(size_t i = 0; i < savedStateNames.size(); i++)
    {
      std::string trajDir = savedStateNames[i];

      if (ramCache.find(trajDir) != ramCache.end())
      {
        TRACE("**** getting trajData from cache");
        TrajData& td = ramCache.find(trajDir)->second;
        if (trajFilter(td))
          trajData.push_back(&td);
      }
      else
      {
        try
        {
          ramCache[trajDir].trajDir = trajDir;
          TrajData& td = ramCache[trajDir];

          bool alreadyProcessed =
            yaatk::exists(getCacheFilename(trajDir));
          if (!alreadyProcessed)
          {
            TRACE("**** getting trajData from execute()");
            TRACE("**** execute()");

            cout << "execute() started : " << trajDir << std::endl;
            cerr << "execute() started : " << trajDir << std::endl;

            {
              yaatk::ChDir cd(
                initialWorkingDirectory + DIR_DELIMIT_STR +
                trajsetDir +
                DIR_DELIMIT_STR + ".." + DIR_DELIMIT_STR,
                false);
              td.setSpottedDistanceFromInit();
            }

            td.execute(stateTemplate);

            cout << "execute() done : " << trajDir << std::endl;
            cerr << "execute() done : " << trajDir << std::endl;

            TRACE("**** ~execute()");

            {
              TRACE(yaatk::getcwd());
              TRACE(getCacheFilename(trajDir));

              yaatk::text_ofstream fo(getCacheFilename(trajDir));
              td.saveToStream(fo);
              fo.close();
            }
          }
          else
          {
            TRACE("**** getting trajData from disk");
            yaatk::text_ifstream fi(getCacheFilename(trajDir));
            td.loadFromStream(fi);
            fi.close();
          }

          if (trajFilter(td))
            trajData.push_back(&td);
        }
        catch (mdepp::BadTrajectoryException& e)
        {
          EPRINT("Trajectory is bad. Excluding it from postprocess.");
          ramCache.erase(trajDir);
        }
      }
    }
  }
}

// StatPostProcess::StatPostProcess()
//   :testProcessClassicMolecule(&ProcessAll),
//    trajData(),
//    id()
// {
//   instanceCounter++;
//   REQUIRE(instanceCounter <= 1);
// }

StatPostProcess::~StatPostProcess()
{
  REQUIRE(instanceCounter > 0);
  instanceCounter--;
}

size_t
StatPostProcess::instanceCounter = 0;

std::map<ClassicMolecule, size_t>
StatPostProcess::buildMassSpectrum(FProcessClassicMolecule fpm) const
{
  std::map<ClassicMolecule, size_t> massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = *trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      if (!fpm(td.molecules[mi])) continue;
      massSpectrum[td.molecules[mi]]++;
    }
  }

  return massSpectrum;
}

void
StatPostProcess::printClassicMoleculesTotal() const
{
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj]->molecules.size() > 0)
  {
    cout << "ClassicMolecules for trajectory " << traj <<
     " ("  << trajData[traj]->trajDir << ") " << " :\n";
    printClassicMolecules(traj);
    cout << "\n";
  }
}

void
StatPostProcess::printClassicMolecules(size_t trajIndex) const
{
  const TrajData& td = *trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    cout << "ClassicMolecule #" << mi << " : ";
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      cout << "( Z = " << atom.Z/mdtk::e << " , gi = " << atom.globalIndex << " ) ";
    }
    cout << "\n            formation time = " << td.molecules[mi].formationTime/mdtk::ps << " ps"
         << ", escape time = " << td.molecules[mi].escapeTime/mdtk::ps << " ps";
    if (td.molecules[mi].hasSubstrateAtoms() && td.molecules[mi].hasClusterAtoms())
      cout << "\nThis is HETEROGENEOUS molecule!!!";
    cout << "\n";
  }
}

void
StatPostProcess::printFullereneInfo() const
{
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj]->trajProjectile.size() > 0)
  {
    cout << "Fullerene for trajectory " << traj << 
     " ("  << trajData[traj]->trajDir << ") " << " :\n";
    printFullereneInfo(traj);
    cout << "\n";
  }
}

void
StatPostProcess::printFullereneInfo(size_t trajIndex) const
{
  using namespace mdtk;
  const TrajData& td = *trajData[trajIndex];
  std::map< Float, AtomGroup >::const_iterator i;
//  REQUIRE(fabs(td.trajProjectile.begin()->first-0.0*ps)<0.05*ps);
  REQUIRE(fabs(td.trajProjectile.rbegin()->first-6.0*ps)<0.05*ps);
  for( i = td.trajProjectile.begin(); i != td.trajProjectile.end() ; ++i )
  {
    std::cout << "Time : " << i->first/mdtk::ps << "\n";
    AtomGroup f = i->second;
//    std::cout << "\t "; TRACE(getVelocity(f.atoms));
    std::cout << "\t "; TRACE(f.maxMolecule().size());
//    std::cout << "\t "; TRACE(f.minDistanceFromMassCenter()/Ao);
//    std::cout << "\t "; TRACE(f.maxDistanceFromMassCenter()/Ao);
//    std::cout << "\t "; TRACE(f.isUnparted());
//    std::cout << "\t "; TRACE(f.isIntegral());
    std::cout << "\t "; TRACE(f.massCenter().z/Ao);
//    std::cout << "\t "; TRACE(f.isEndoFullerene());
//    std::cout << "\t "; TRACE(f.cluster.maxMolecule().size());
  }
}

bool
isAmongSputtered(const Atom& a, const std::vector<ClassicMolecule>& molecules)
{
  for(size_t mi = 0; mi < molecules.size(); mi++)
  {
    for(size_t ai = 0; ai < molecules[mi].atoms.size(); ai++)
    {
      const mdtk::Atom& atom = molecules[mi].atoms[ai];
      if (atom.globalIndex == a.globalIndex)
        return true;
    }
  }
  return false;
}

bool
isAmongSputtered(const AtomGroup& atoms, const std::vector<ClassicMolecule>& molecules)
{
  for(size_t pi = 0; pi < atoms.atoms.size(); pi++)
  {
    const mdtk::Atom& atom = atoms.atoms[pi];
    if (!isAmongSputtered(atom,molecules))
      return false;
  }
  return true;
}

bool consistsMostlyOfTargetAtoms(const AtomGroup& atoms, double threshold = 0.5)
{
  size_t count = 0;
  for(size_t pi = 0; pi < atoms.atoms.size(); pi++)
  {
    const mdtk::Atom& atom = atoms.atoms[pi];
    if (atom.hasTag(ATOMTAG_TARGET))
      count++;
  }
  return double(count)/atoms.atoms.size() >= threshold;
}

void
StatPostProcess::printCoefficients() const
{
  std::ofstream fo("Coefficients.txt");

  size_t stickedProjectiles = 0;
  size_t backscatteredProjectiles = 0;

  size_t stickedIntegralProjectiles = 0;
  size_t backscatteredIntegralProjectiles = 0;

  size_t stickedProjectileAtoms = 0;
  size_t backscatteredProjectileAtoms = 0;

  size_t sputteredTargetAtoms = 0;
  size_t sputteredTargetMolecules = 0;
  size_t sputteredIntegralTargetMolecules = 0;

  size_t sputteredMolecules = 0;
  size_t sputteredIntegralMolecules = 0;

  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    const TrajData& td = *trajData[traj];
    REQUIRE(td.trajProjectile.size() > 0);
//    REQUIRE(fabs(td.trajProjectile.begin()->first-0.0*ps)<0.05*ps);
    REQUIRE(fabs(td.trajProjectile.rbegin()->first-6.0*ps)<0.05*ps);

    std::map< Float, AtomGroup >::const_iterator i = td.trajProjectile.end();
    i--;
    const AtomGroup& projectile = i->second;
    REQUIRE(fabs(i->first-6.0*ps)<0.05*ps);

    if (projectile.isMolecule())
    {
      if (isAmongSputtered(projectile,td.molecules))
      {
        backscatteredProjectiles++;
        if (projectile.isMonomer())
          backscatteredIntegralProjectiles++;
        if (projectile.isFullerene() && Fullerene(projectile).isIntegral())
          backscatteredIntegralProjectiles++;
      }
      else
      {
        stickedProjectiles++;
        if (projectile.isMonomer())
          stickedIntegralProjectiles++;
        if (projectile.isFullerene() && Fullerene(projectile).isIntegral())
          stickedIntegralProjectiles++;
      }
    }

    for(size_t pi = 0; pi < projectile.size(); pi++)
    {
      const mdtk::Atom& patom = projectile.atoms[pi];
      if (!isAmongSputtered(patom,td.molecules))
        stickedProjectileAtoms++;
    }

    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
      {
        const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
        if (atom.hasTag(ATOMTAG_PROJECTILE))
          backscatteredProjectileAtoms++;
        if (atom.hasTag(ATOMTAG_TARGET))
          sputteredTargetAtoms++;
      }
      const AtomGroup m(td.molecules[mi]);
      if (m.isMonomer() && m.isMetalCluster())
      {
        sputteredIntegralMolecules++;
        sputteredMolecules++;
        if (consistsMostlyOfTargetAtoms(m))
        {
          sputteredIntegralTargetMolecules++;
          sputteredTargetMolecules++;
        }
      }
      if (m.isFullerene())
      {
        sputteredMolecules++;
        if (consistsMostlyOfTargetAtoms(m))
          sputteredTargetMolecules++;
        if (Fullerene(m).isIntegral())
        {
          sputteredIntegralMolecules++;
          if (consistsMostlyOfTargetAtoms(m))
            sputteredIntegralTargetMolecules++;
        }
      }
    }
  }

  Float trajCount = trajData.size();
  fo TRACESS(trajCount);

  fo << "\n";

  fo TRACESS(stickedProjectiles/trajCount);
  fo TRACESS(backscatteredProjectiles/trajCount);

  fo << "\n";

  fo TRACESS(stickedIntegralProjectiles/trajCount);
  fo TRACESS(backscatteredIntegralProjectiles/trajCount);

  fo << "\n";

  fo TRACESS(stickedProjectileAtoms/trajCount);
  fo TRACESS(backscatteredProjectileAtoms/trajCount);

  fo << "\n";

  fo TRACESS(sputteredTargetAtoms/trajCount);
  fo TRACESS(sputteredTargetMolecules/trajCount);
  fo TRACESS(sputteredIntegralTargetMolecules/trajCount);

  fo << "\n";

  fo TRACESS(sputteredMolecules/trajCount);
  fo TRACESS(sputteredIntegralMolecules/trajCount);

  REQUIRE(stickedIntegralProjectiles <= stickedProjectiles);
  REQUIRE(backscatteredIntegralProjectiles <= backscatteredProjectiles);
  REQUIRE(sputteredIntegralMolecules <= sputteredMolecules);
  REQUIRE(sputteredIntegralTargetMolecules <= sputteredTargetMolecules);

  fo.close();
}

int
StatPostProcess::getYieldSum( FProcessClassicMolecule fpm) const
{
  int totalSpotted = 0;
  for(size_t i = 0; i < trajData.size();i++)
  {
    totalSpotted += getYield(i,fpm);
  }
  return totalSpotted;
}

int
StatPostProcess::getYield(size_t trajIndex, FProcessClassicMolecule fpm) const
{
  return trajData[trajIndex]->yield(fpm);
}

Float
StatPostProcess::getYieldAverage( FProcessClassicMolecule fpm) const
{
  return Float(getYieldSum(fpm))/trajData.size();
}

Float
StatPostProcess::getYieldNormalizedByClusterSize(size_t trajIndex, FProcessClassicMolecule fpm) const
{
  return Float(getYield(trajIndex,fpm))/id.clusterSize;
}

Float
StatPostProcess::getYieldAverageProgress( FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"average_Yield.dat");
  {
    std::ofstream fo((std::string(ofilename)+".plt").c_str());
    fo << "reset\n#set yrange [0:*]\nplot '" << ofilename << "' with impulses";
    fo.close();
  }

  std::ofstream fo(ofilename);
  Float avg = 0.0;
  size_t i;
  for(i = 0; i < trajData.size();i++)
  {
    avg += getYield(i,fpm);
    fo << i+1 << " " << avg/(i+1) << "\n";
  }
  fo.close();
  return avg/i;
}

int
StatPostProcess::getYieldFragmentsCount(size_t trajIndex, FProcessClassicMolecule fpm) const
{
  int spotted = 0;
  const TrajData& td = *trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const ClassicMolecule& mol = td.molecules[mi];
    if (!fpm(mol)) continue;
    spotted++;
  }
  return spotted;
}

Float
StatPostProcess::getYieldFragmentsCountNormalizedByClusterSize(size_t trajIndex, FProcessClassicMolecule fpm) const
{
  return Float(getYieldFragmentsCount(trajIndex,fpm))/id.clusterSize;
}

Float
StatPostProcess::getTotalEnergyOfSputtered( FProcessClassicMolecule fpm) const
{
  Float totalSpotted = 0;
  for(size_t i = 0; i < trajData.size();i++)
  {
    totalSpotted += getEnergyOfSputtered(i,fpm);
  }
  return totalSpotted;
}

Float
StatPostProcess::getEnergyOfSputtered(size_t trajIndex, FProcessClassicMolecule fpm) const
{
  Float spotted = 0;
  const TrajData& td = *trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const ClassicMolecule& mol = td.molecules[mi];
    if (!fpm(mol)) continue;
    for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)    {
      const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
      spotted += SQR(atom.V.module())*atom.M/2.0;
    }
  }
  return spotted;
}

Float
StatPostProcess::getAverageEnergyOfSputtered( FProcessClassicMolecule fpm) const
{
  return Float(getTotalEnergyOfSputtered(fpm))/trajData.size();
}

void
depthHist2file(const char* filename, std::vector<Float>& depth, Float scale = 1.0)
{
  using mdtk::Exception;

/*
  const Float minDepth = -10.0;
  const Float maxDepth =  20.0;
  const int n = (maxDepth-minDepth)/(2.547);// for polyethylene // prev was (0.50);
*/

/*
  // Graphite
  const Float minDepth_desired   = -40.0;
  const Float maxDepth_desired   =  40.0;
  const Float matchPoint         =   0.0;
  const Float c = 6.708;
  const Float histStep           =   c/6.0;//2.547; // for polyethylene // prev was (0.50);
*/
/*
  // Cu
  const Float c = 3.61;
  const Float minDepth_desired   = -c/2.0*15.123;
  const Float maxDepth_desired   =  c/2.0*15.123;
  const Float matchPoint         = -c/4.0;
  const Float histStep           =   c/2.0;//2.547; // for polyethylene // prev was (0.50);
*/

  // PE
  const Float c = 2.547;
  const Float minDepth_desired   = -100.0;
  const Float maxDepth_desired   =  100.0;
  const Float matchPoint         =    0.0;
  const Float histStep           =   c/3.0;//2.547; // for polyethylene // prev was (0.50);

  REQUIRE(minDepth_desired < matchPoint);
  REQUIRE(maxDepth_desired > matchPoint);
  const int n_below_matchPoint = int( (matchPoint       - minDepth_desired)/histStep ) +1;
  const int n_above_matchPoint = int( (maxDepth_desired - matchPoint      )/histStep ) +1;

  const int n = n_above_matchPoint + n_below_matchPoint;
  const Float minDepth   = matchPoint - n_below_matchPoint*histStep;
  const Float maxDepth   = matchPoint + n_above_matchPoint*histStep;

  {
    std::ofstream fo((std::string(filename)+".plt").c_str());
    fo << "reset\n#set yrange [*:0]\nset yrange [0:*]\nset xrange [-3:25]\nset format x \"%.1f\"\nset xtics " << c/2.0 << "\nset grid xtics\nplot \'" << filename << "\' with boxes\n";
    fo.close();
  }
  std::ofstream fo(filename);

  fo << "# min depth = " << minDepth << " Ao" << "\n"
     << "# max depth = " << maxDepth << " Ao" << "\n"
     << "# number of bins = " << n << "\n";

  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, minDepth, maxDepth);

  for(size_t i = 0; i < depth.size(); i++)
    gsl_histogram_increment (h, depth[i]);

  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    fo << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i)*scale << "\n";
  }

  gsl_histogram_free (h);

  fo.close();
}

void
saveDepth(std::vector<Float>& depth, const char *filename)
{
  std::ofstream fo(filename);
  fo << depth.size() << "\n";
  for(size_t i = 0; i < depth.size(); i++)
    fo << depth[i] << "\n";
  fo.close();
}

#define MDEPP_BYDEPTH_DIR "_by_depth"

void
StatPostProcess::spottedByDepth() const
{
  yaatk::mkdir(MDEPP_BYDEPTH_DIR);
  yaatk::chdir(MDEPP_BYDEPTH_DIR);

  std::vector<Float> depth;
  std::vector<Float> depth_H;
  std::vector<Float> depth_C;
  for(size_t traj = 0; traj < trajData.size(); traj++)
  {
    if (trajData[traj]->molecules.size() == 0)
      continue;
    const TrajData& td = *trajData[traj];

    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
      {
        const mdtk::Atom& atom_init = td.molecules[mi].atoms_init[ai];
        depth.push_back(atom_init.coords.z/mdtk::Ao);
        if (atom_init.ID == mdtk::H_EL)
          depth_H.push_back(atom_init.coords.z/mdtk::Ao);
        if (atom_init.ID == mdtk::C_EL)
          depth_C.push_back(atom_init.coords.z/mdtk::Ao);
      }
    }
  }
  {
    saveDepth(depth,"spotted_depths_unsorted.dat");
    sort(depth.begin(),depth.end());
    saveDepth(depth,"spotted_depths.dat");
  }
  {
    saveDepth(depth_C,"spotted_depths_C_unsorted.dat");
    sort(depth_C.begin(),depth_C.end());
    saveDepth(depth_C,"spotted_depths_C.dat");
  }
  {
    saveDepth(depth_H,"spotted_depths_H_unsorted.dat");
    sort(depth_H.begin(),depth_H.end());
    saveDepth(depth_H,"spotted_depths_H.dat");
  }
  {
    depthHist2file("spots_by_depth.dat",depth);
    depthHist2file("spots_by_depth_H.dat",depth_H);
    depthHist2file("spots_by_depth_C.dat",depth_C);
  }

  yaatk::chdir("..");
}

double
StatPostProcess::moleculeEnergy(const ClassicMolecule& mol)
{
  return 0.5*SQR(mol.getVelocity().module())*mol.getMass();
}

double
StatPostProcess::moleculeEnergyInEV(const ClassicMolecule& mol)
{
  return 0.5*SQR(mol.getVelocity().module())*mol.getMass()/mdtk::eV;
}

double
StatPostProcess::moleculeMass(const ClassicMolecule& mol)
{
  return mol.getMass();
}

double
StatPostProcess::moleculeMassInAMU(const ClassicMolecule& mol)
{
  return mol.getAMUMass();
}

double
StatPostProcess::moleculeCount(const ClassicMolecule& mol)
{
  return 1.0;
}

double
StatPostProcess::moleculeAtomsCount(const ClassicMolecule& mol)
{
  return mol.atoms.size();
}

double
StatPostProcess::moleculeEnergyByAtom(const ClassicMolecule& mol)
{
  return moleculeEnergy(mol)/moleculeAtomsCount(mol);
}

std::string
StatPostProcess::FMoleculeAttributeToString(FMoleculeAttribute fma)
{
  if (fma == moleculeEnergy) return "moleculeEnergy";
  if (fma == moleculeEnergyInEV) return "moleculeEnergy";
  if (fma == moleculeMass) return "moleculeMass";
  if (fma == moleculeMassInAMU) return "moleculeMassInAMU";
  if (fma == moleculeCount) return "moleculeCount";
  if (fma == moleculeAtomsCount) return "moleculeAtomsCount";
  if (fma == moleculeEnergyByAtom) return "moleculeEnergyByAtom";

  throw Exception("Unknown FMoleculeAttribute");
}

std::map<Float, Float>
StatPostProcess::distByAngle(
  AngleType angleType,
  const int n,
  FMoleculeAttribute fma,
  FProcessClassicMolecule fpm) const
{
  gsl_histogram * h = gsl_histogram_alloc (n);

  const Float minAngle = (angleType == ANGLE_POLAR) ?  0.0 : -180.0;
  const Float maxAngle = (angleType == ANGLE_POLAR) ? 90.0 : +180.0;

  gsl_histogram_set_ranges_uniform (h, minAngle, maxAngle);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = *trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
      Float azimuthal = atan2(v.y,v.x)/mdtk::Deg;
      gsl_histogram_accumulate(
        h,
        (angleType == ANGLE_POLAR) ? polar : azimuthal,
        fma(mol)/Float(trajData.size()));
    }
  }

  std::map<Float, Float> histData;

  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    histData[(lower+upper)/2.0] = gsl_histogram_get(h,i);
  }

  gsl_histogram_free (h);

  return histData;;
}


Float
StatPostProcess::maxMoleculeAttribute(
  FMoleculeAttribute fma, FProcessClassicMolecule moleculeFilter) const
{
  Float maxVal = 0.0;

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const StatPostProcess::TrajData& td = *trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!moleculeFilter(mol)) continue;
      Float val = fma(mol);
      if (val > maxVal)
        maxVal = val;
    }
  }

  return maxVal;
}

Float
StatPostProcess::suggestedBinWidth(
  FMoleculeAttribute fma, FProcessClassicMolecule moleculeFilter,
  ElementID ionElement, size_t clusterSize, ElementID clusterElement)
{
  Float binWidth = 0.0;

  if (fma == moleculeMass)
  {
    ClassicMolecule mol;
    mol.atoms.push_back(Atom(H_EL));
    binWidth = mol.getMass();
    if (*moleculeFilter == StatPostProcess::ProcessCluster ||
        *moleculeFilter == StatPostProcess::ProcessAll)
    {
      mol.atoms.clear();
      mol.atoms.push_back(Atom(clusterElement));
      binWidth = mol.getMass();
    }
    if (*moleculeFilter == StatPostProcess::ProcessProjectile)
    {
      mol.atoms.clear();
      mol.atoms.push_back(Atom(ionElement));
      binWidth = mol.getMass();
    }
    if (*moleculeFilter == StatPostProcess::ProcessSubstrate)
    {
      mol.atoms.clear();
      mol.atoms.push_back(Atom(C_EL));
      binWidth = mol.getMass();
    }
  }

  if (fma == moleculeEnergy)
  {
    binWidth = 10.0*eV;
  }

  REQUIRE(binWidth > 0.0);

  return binWidth;
}

std::map<Float, Float>
StatPostProcess::distBy(
  FMoleculeAttribute histFunc,
  Float binWidth,
  Float histMin,
  Float histMax,
  FMoleculeAttribute fma,
  FProcessClassicMolecule fpm
  ) const
{
  const int n = (histMax - histMin)/binWidth + 1;
  REQUIRE(n > 0);
  REQUIRE(n < 1000000);
  double range[n+1];
  range[0] = 0.0;
  for(size_t binIndex = 0; binIndex < n; ++binIndex)
    range[binIndex + 1] = range[binIndex] + binWidth;

  gsl_histogram* h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges(h, range, n+1);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = *trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      gsl_histogram_accumulate(
        h,
        histFunc(mol),
        fma(mol)/Float(trajData.size()));
    }
  }

  std::map<Float, Float> histData;

  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    histData[(lower+upper)/2.0] = gsl_histogram_get(h,i);
  }

  gsl_histogram_free (h);

  return histData;
}

std::map<Float, Float>
StatPostProcess::divideHistograms(
  std::map<Float, Float>& h1,
  std::map<Float, Float>& h2)
{
  REQUIRE(h1.size() == h2.size());

  std::map<Float, Float> h;
  std::map<Float, Float>::iterator i;
  for(i = h1.begin(); i != h1.end(); ++i)
    h[i->first] = h1[i->first] / h2[i->first];

  return h;
}

void
saveHistogram(gsl_histogram *h, const char *datFileName)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());
  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << "\n";
  }
  foByEscapeHist.close();
  foByEscapeHistPlt << "#reset\nset yrange [0:*]\nplot \'" << byEscapeTimeDatHist << "\' with boxes fs solid 1.0\n";
  foByEscapeHistPlt.close();
}

void
saveHistogram_new(gsl_histogram *h, const char *datFileName)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());
  for(size_t i = 0; i < gsl_histogram_bins(h); i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << "\n";
  }
  foByEscapeHist.close();
  foByEscapeHistPlt << "reset\nset yrange [0:*]\nset xrange [-180:180]\nplot \'" << byEscapeTimeDatHist << "\' with boxes\n";
  foByEscapeHistPlt << "pause -1 \"Press Enter\"\n";
  foByEscapeHistPlt.close();
}

void
saveHistogram_polar(gsl_histogram *h, const char *datFileName/*, bool halfshift = false*/)
{
  std::string byEscapeTimeDatHist(datFileName);
  std::ofstream foByEscapeHist(byEscapeTimeDatHist.c_str());
  std::ofstream foByEscapeHistPlt((byEscapeTimeDatHist+".plt").c_str());

  int n;
  n = gsl_histogram_bins(h);

  for(size_t i = 0; i < gsl_histogram_bins(h)+1; i++)
  {
    double lower, upper;
    double index = i;
    if (i == gsl_histogram_bins(h)) index = 0;
    gsl_histogram_get_range (h, index, &lower, &upper);
    Float ang = (lower+upper)/2.0;//-(halfshift?((360.0/n)/2.0):0.0);
//    if (ang < -180) ang += ;
    foByEscapeHist << ang << " " << gsl_histogram_get(h,index) << "\n";
  }
  foByEscapeHist.close();
  foByEscapeHistPlt << "reset\nset style fill pattern 1\nset polar\nset angles degrees\nset size ratio -1\n\
\nset grid polar\n\nplot \'" << byEscapeTimeDatHist << "\' with filledcurves lw 2 notitle,\\\n \'" << byEscapeTimeDatHist << "\' with impulses lw 2 notitle\n";
//  foByEscapeHistPlt << "pause -1 \"Press Enter\"\n";
  foByEscapeHistPlt.close();
}

#define MDEPP_BYTIME_DIR "_by_time"

void
StatPostProcess::buildByTime(FProcessClassicMolecule fpm) const
{
  yaatk::mkdir(MDEPP_BYTIME_DIR);
  yaatk::chdir(MDEPP_BYTIME_DIR);

  char ofilename[1024];

  sprintf(ofilename,"mass_by_formation_time.dat");

  std::string byFormationTimeDat(ofilename);
  std::ofstream foByFormation(byFormationTimeDat.c_str());
  std::ofstream foByFormationPlt((byFormationTimeDat+".plt").c_str());
  foByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << byFormationTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_formation_time.dat");

  std::string EnergyByFormationTimeDat(ofilename);
  std::ofstream foEnergyByFormation(EnergyByFormationTimeDat.c_str());
  std::ofstream foEnergyByFormationPlt((EnergyByFormationTimeDat+".plt").c_str());
  foEnergyByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByFormationTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_atom_by_formation_time.dat");

  std::string EnergyByAtomByFormationTimeDat(ofilename);
  std::ofstream foEnergyByAtomByFormation(EnergyByAtomByFormationTimeDat.c_str());
  std::ofstream foEnergyByAtomByFormationPlt((EnergyByAtomByFormationTimeDat+".plt").c_str());
  foEnergyByAtomByFormationPlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByAtomByFormationTimeDat << "\' with points\n";

  const Float minTime =    -0.05 /* *mdtk::ps*/;
  const Float maxTime =     6.05 /* *mdtk::ps*/;
  const int n = ceil((maxTime-minTime)/(0.1 /* *mdtk::ps*/));

  gsl_histogram * h_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_by_formation, minTime, maxTime);
  gsl_histogram * h_count_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_count_by_formation, minTime, maxTime);
  gsl_histogram * h_energy_by_formation = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_energy_by_formation, minTime, maxTime);

  sprintf(ofilename,"mass_by_escape_time.dat");

  std::string byEscapeTimeDat(ofilename);
  std::ofstream foByEscape(byEscapeTimeDat.c_str());
  std::ofstream foByEscapePlt((byEscapeTimeDat+".plt").c_str());
  foByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << byEscapeTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_escape_time.dat");

  std::string EnergyByEscapeTimeDat(ofilename);
  std::ofstream foEnergyByEscape(EnergyByEscapeTimeDat.c_str());
  std::ofstream foEnergyByEscapePlt((EnergyByEscapeTimeDat+".plt").c_str());
  foEnergyByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByEscapeTimeDat << "\' with points\n";

  sprintf(ofilename,"energy_by_atom_by_escape_time.dat");

  std::string EnergyByAtomByEscapeTimeDat(ofilename);
  std::ofstream foEnergyByAtomByEscape(EnergyByAtomByEscapeTimeDat.c_str());
  std::ofstream foEnergyByAtomByEscapePlt((EnergyByAtomByEscapeTimeDat+".plt").c_str());
  foEnergyByAtomByEscapePlt << "reset\n#set yrange [0:*]\nplot \'" << EnergyByAtomByEscapeTimeDat << "\' with points\n";


  gsl_histogram * h_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_by_escape, minTime, maxTime);
  gsl_histogram * h_count_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_count_by_escape, minTime, maxTime);
  gsl_histogram * h_energy_by_escape = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h_energy_by_escape, minTime, maxTime);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = *trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& m = td.molecules[mi];
//      if (m.hasProjectileAtoms()) continue;
      if (!fpm(m)) continue;

      Float atomsCount = m.atoms.size();
      Float moleculeMass = m.getAMUMass();
      Float formationTime = m.formationTime /mdtk::ps;
      Float escapeTime = m.escapeTime /mdtk::ps;
      mdtk::Vector3D velocity = m.getVelocity();
      Float kineticEnergy = 0.5*moleculeMass*mdtk::amu*SQR(velocity.module())/mdtk::eV;
      foByFormation << formationTime << " " << moleculeMass << "\n";
      foByEscape    << escapeTime    << " " << moleculeMass << "\n";
      foEnergyByFormation << formationTime << " " << kineticEnergy << "\n";
      foEnergyByEscape    << escapeTime    << " " << kineticEnergy << "\n";
      foEnergyByAtomByFormation << formationTime << " " << kineticEnergy/atomsCount << "\n";
      foEnergyByAtomByEscape    << escapeTime    << " " << kineticEnergy/atomsCount << "\n";
      gsl_histogram_accumulate(h_by_formation, formationTime, moleculeMass);
      gsl_histogram_accumulate(h_by_escape,       escapeTime, moleculeMass);
      gsl_histogram_accumulate(h_count_by_formation, formationTime, atomsCount);
      gsl_histogram_accumulate(h_count_by_escape,       escapeTime, atomsCount);
      gsl_histogram_accumulate(h_energy_by_formation, formationTime, kineticEnergy);
      gsl_histogram_accumulate(h_energy_by_escape,       escapeTime, kineticEnergy);
    }
  }

  foByEscape.close();
  foByEscapePlt.close();
  foEnergyByEscape.close();
  foEnergyByEscapePlt.close();
  foEnergyByAtomByEscape.close();
  foEnergyByAtomByEscapePlt.close();

  foByFormation.close();
  foByFormationPlt.close();
  foEnergyByFormation.close();
  foEnergyByFormationPlt.close();
  foEnergyByAtomByFormation.close();
  foEnergyByAtomByFormationPlt.close();

  sprintf(ofilename,"mass_by_escape_time_HIST.dat");
  saveHistogram(h_by_escape,ofilename);
  gsl_histogram_free(h_by_escape);

  sprintf(ofilename,"mass_by_formation_time_HIST.dat");
  saveHistogram(h_by_formation,ofilename);
  gsl_histogram_free(h_by_formation);

  sprintf(ofilename,"atomscount_by_escape_time_HIST.dat");
  saveHistogram(h_count_by_escape,ofilename);
  gsl_histogram_free(h_count_by_escape);

  sprintf(ofilename,"energy_by_escape_time_HIST.dat");
  saveHistogram(h_energy_by_escape,ofilename);
  gsl_histogram_free(h_energy_by_escape);

  sprintf(ofilename,"atomscount_by_formation_time_HIST.dat");
  saveHistogram(h_count_by_formation,ofilename);
  gsl_histogram_free(h_count_by_formation);

  sprintf(ofilename,"energy_by_formation_time_HIST.dat");
  saveHistogram(h_energy_by_formation,ofilename);
  gsl_histogram_free(h_energy_by_formation);

  yaatk::chdir("..");
}

}
