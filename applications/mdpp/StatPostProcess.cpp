/* 
   Molecular dynamics postprocessor, main classes

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014
   Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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


namespace mdepp
{

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

  TRACE(str);
  TRACE(ElementIDtoString(clusterElement));
  TRACE(clusterSize);
  TRACE(ElementIDtoString(ionElement));
  TRACE(ionEnergy);
  REQUIRE(str.size()>1);
  REQUIRE(*str.begin()=='C' || *str.begin()=='A');
  REQUIRE(*(str.end()-1)=='V');
}

void
StatPostProcess::buildSputteredClassicMolecules(mdtk::SimLoop& state,size_t trajIndex,
  StatPostProcess::StateType s, NeighbourList& nl)
{
  TrajData& td = trajData[trajIndex];
  if (s == STATE_FINAL)
  {
    cout << "Building molecules for state ..." << std::endl;
    for(size_t atomIndex = 0; atomIndex < state.atoms.size(); atomIndex++)
    {
      mdtk::Atom &atom = state.atoms[atomIndex];
      if (atom.coords.z < SPOTTED_DISTANCE)
      {
	bool account_atom = true;
	for(size_t mi = 0; mi < td.molecules.size(); mi++)
	{
	  if (td.molecules[mi].hasAtom(atom))
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
	    cout << "Adding molecule." << std::endl;
	    td.molecules.push_back(molecule);
	  }
	}
      }
    }
    std::sort(td.molecules.begin(),td.molecules.end());
  }
  else if (s == STATE_INIT)
  {
    mdtk::SimLoop* mde_init = &state;
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      for(size_t ai = 0; ai < td.molecules[mi].atoms.size(); ai++)
      {
        const mdtk::Atom& atom = td.molecules[mi].atoms[ai];
        const mdtk::Atom& atom_init = mde_init->atoms[atom.globalIndex];
        td.molecules[mi].atoms_init.push_back(atom_init);
        REQUIRE(td.molecules[mi].atoms[ai].globalIndex ==
		td.molecules[mi].atoms_init[ai].globalIndex);
      }
      REQUIRE(td.molecules[mi].atoms.size() ==
	      td.molecules[mi].atoms_init.size());
    }
  }
  else if (s == STATE_INTER)
  {
    mdtk::SimLoop* mde_inter = &state;
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      bool formedNowOrEarlier = true;
      bool escapedNowOrEarlier = true;
      ClassicMolecule& molecule = td.molecules[mi];
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
StatPostProcess::buildDummyDynamics(mdtk::SimLoop& state,size_t trajIndex,
  StatPostProcess::StateType s)
{
  TrajData& td = trajData[trajIndex];
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

void
StatPostProcess::removeBadTrajectories()
{
  std::vector<size_t> badTrajIndices;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    TrajData& td = trajData[trajIndex];

    yaatk::ChDir cd(td.trajDir, false);

    mdtk::SimLoop state;
    state.allowToFreePotentials = true;
    setupPotentials(state);

    state.initNLafterLoading = false;

    mdtk::SimLoopSaver mds(state);
    REQUIRE(mds.loadIterationLatest() & mdtk::SimLoopSaver::LOADED_R);

    state.atoms.prepareForSimulatation();

    {
      SimLoop::Check& check = state.check;
      Float Eo_plus_Eb = check.initialEnergy + check.energyTransferredFromBath;
      Float dE = check.currentEnergy - Eo_plus_Eb;
      Float dE_by_Eo_plus_Eb = dE/Eo_plus_Eb;
      if (fabs(dE_by_Eo_plus_Eb) > 0.01)
      {
        cerr << "Trajectory " << trajData[trajIndex].trajDir << " has bad energy conservation. Handle it separately." << endl;
        badTrajIndices.push_back(trajIndex);
      }
      else
      {
        if (state.simTimeFinal > state.simTime)
        {
          cerr << "Trajectory " << trajData[trajIndex].trajDir << " seems to be unfinished. Handle it separately." << endl;
          badTrajIndices.push_back(trajIndex);
        }
      }
    }
  }

  size_t oldSize = trajData.size();

  for(int i = badTrajIndices.size()-1; i >= 0; i--)
  {
    cerr << "Reoving trajectory from postprocess: "
         << (trajData.begin()+badTrajIndices[i])->trajDir << endl;
    trajData.erase(trajData.begin()+badTrajIndices[i]);
  }

  cerr << "Removed " << oldSize - trajData.size() << " trajectories." << endl;
}

void
StatPostProcess::setSpottedDistanceFromInit()
{
  using mdtk::Exception;

  REQUIRE(trajData.size() > 0);

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

void
StatPostProcess::execute()
{
  using mdtk::Exception;

  removeBadTrajectories();

  if (0)
  {
    REQUIRE(trajData.size() >= 1);
    trajData.resize(1);
  }

  cout << "PostProcess::execute() started." << std::endl;
  cerr << "PostProcess::execute() started." << std::endl;

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    TrajData& td = trajData[trajIndex];

    yaatk::ChDir cd(td.trajDir, false);

    mdtk::SimLoop state;
    state.allowToFreePotentials = true;
    setupPotentials(state);

    {
      state.initNLafterLoading = false;

      mdtk::SimLoopSaver mds(state);
      REQUIRE(mds.loadIterationLatest() & mdtk::SimLoopSaver::LOADED_R);

      state.atoms.prepareForSimulatation();
      setTags(state);
      {
        REQUIRE(mds.listIterations().size() > 0);
        cout << "State " << *mds.listIterations().rbegin() << " loaded." << std::endl;
      }
    }
    NeighbourList nl(state.atoms);

    Molecule projectile;
    {
      projectile.buildByTag(state,ATOMTAG_PROJECTILE);
      td.trajProjectile[state.simTime] = projectile;

      projectile.update(state);
      td.trajProjectile[state.simTime] = projectile;
    }

    Molecule cluster;
    {
      cluster.buildByTag(state,ATOMTAG_CLUSTER);
      td.trajCluster[state.simTime] = cluster;

      cluster.update(state);
      td.trajCluster[state.simTime] = cluster;
    }

    td.PBC = state.atoms.PBC();
    TRACE(td.PBC/mdtk::Ao);

//    TRACE(getAboveSpottedHeight(*state));

    buildSputteredClassicMolecules(state,trajIndex,STATE_FINAL,nl);
//    buildClusterDynamics(*state,trajIndex,STATE_FINAL,nl);
//    buildProjectileDynamics(*state,trajIndex,STATE_FINAL);

//  if (td.molecules.size() > 0)

#if 0

//    if (0)
    {
      mdtk::SimLoop mde_init;
      mde_init.allowToFreePotentials = true;
//      mde_init.allowToFreeAtoms = true;
      setupPotentials(mde_init);
      std::string mde_init_filename = td.trajDir+"mde_init";
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
        td.trajProjectile[mde_init.simTime] = projectile;
      }

      {
        cluster.update(mde_init);
        td.trajCluster[mde_init.simTime] = cluster;
      }
    }

//  if (td.molecules.size() > 0)
    {
      std::vector<std::string> interStates;
      findIntermediateStates(td.trajDir,interStates);

      mdtk::SimLoop mde_inter;
      mde_inter.allowToFreePotentials = true;
//      mde_inter.allowToFreeAtoms = true;
      setupPotentials(mde_inter);

      std::string trajFinalName = td.trajDir+"mde_init";
      cout << "Loading state " << trajFinalName << std::endl;
      yaatk::text_ifstream fi(trajFinalName.c_str());
      mde_inter.initNLafterLoading = false;
      mde_inter.loadFromStream(fi);
      fi.close();

      for(int stateIndex = interStates.size()-1; stateIndex >= 0; stateIndex--)
      {
//	if (stateIndex != interStates.size()-1) continue;
	std::string mde_inter_filename = td.trajDir+interStates[stateIndex];

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
          td.trajProjectile[mde_inter.simTime] = projectile;
        }
        {
          cluster.update(mde_inter);
          td.trajCluster[mde_inter.simTime] = cluster;
        }
      }
    }

#endif

    cout << "Building molecules for state done." << std::endl;
  }

  cout << "PostProcess::execute() done." << std::endl;
  cerr << "PostProcess::execute() done." << std::endl;
}

void
StatPostProcess::addHalo(const StatPostProcess& pp)
{
  REQUIRE(id.str == pp.id.str);
  REQUIRE(SPOTTED_DISTANCE == pp.SPOTTED_DISTANCE);
  for(size_t i = 0; i < pp.trajData.size(); i++)
    trajData.push_back(pp.trajData[i]);
}

void
StatPostProcess::buildMassSpectrum(FProcessClassicMolecule fpm) const
{
  std::ofstream fo("massSpectrum.txt");

  std::map<ClassicMolecule, size_t> massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      if (!fpm(td.molecules[mi])) continue;
      massSpectrum[td.molecules[mi]]++;
    }
  }
/*
  size_t specieYield = 0;
  size_t specieIndex = 0;
*/
  std::map<ClassicMolecule, size_t>::iterator i = massSpectrum.begin();
  while (i != massSpectrum.end())
  {
    fo << i->first.formula()
       << " (" << i->first.getAMUMass() << " amu) : "
       << i->second << "/" << trajData.size()
       << " = " << double(i->second)/trajData.size() << "\n";
    i++;
  }

  fo.close();
}

void
StatPostProcess::printClassicMoleculesTotal() const
{
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].molecules.size() > 0)
  {
    cout << "ClassicMolecules for trajectory " << traj <<
     " ("  << trajData[traj].trajDir << ") " << " :\n";
    printClassicMolecules(traj);
    cout << std::endl;
  }
}

void
StatPostProcess::printClassicMolecules(size_t trajIndex) const
{
  const TrajData& td = trajData[trajIndex];
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
    cout << std::endl;
  }
}

void
StatPostProcess::printFullereneInfo() const
{
  for(size_t traj = 0; traj < trajData.size(); traj++)
  if (trajData[traj].trajProjectile.size() > 0)
  {
    cout << "Fullerene for trajectory " << traj << 
     " ("  << trajData[traj].trajDir << ") " << " :\n";
    printFullereneInfo(traj);
    cout << std::endl;
  }
}

void
StatPostProcess::printFullereneInfo(size_t trajIndex) const
{
  using namespace mdtk;
  const TrajData& td = trajData[trajIndex];
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
    const TrajData& td = trajData[traj];
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
StatPostProcess::getAboveSpottedHeight(mdtk::SimLoop& state) const
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
  int spotted = 0;
  const TrajData& td = trajData[trajIndex];
  for(size_t mi = 0; mi < td.molecules.size(); mi++)
  {
    const ClassicMolecule& mol = td.molecules[mi];
    if (!fpm(mol)) continue;
    spotted += mol.atoms.size();
  }
  return spotted;
}

Float
StatPostProcess::getAverageYield( FProcessClassicMolecule fpm) const
{
  return Float(getYieldSum(fpm))/trajData.size();
}

Float
StatPostProcess::getAverageYieldProgress( FProcessClassicMolecule fpm) const
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
  const TrajData& td = trajData[trajIndex];
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

  fo << "# min depth = " << minDepth << " Ao" << std::endl
     << "# max depth = " << maxDepth << " Ao" << std::endl
     << "# number of bins = " << n << std::endl;

  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, minDepth, maxDepth);

  for(size_t i = 0; i < depth.size(); i++)
    gsl_histogram_increment (h, depth[i]);

  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    fo << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i)*scale << std::endl;
  }

  gsl_histogram_free (h);

  fo.close();
}

void
saveDepth(std::vector<Float>& depth, const char *filename)
{
  std::ofstream fo(filename);
  fo << depth.size() << std::endl;
  for(size_t i = 0; i < depth.size(); i++)
    fo << depth[i] << std::endl;
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
    if (trajData[traj].molecules.size() == 0)
      continue;
    const TrajData& td = trajData[traj];

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
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << std::endl;
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
    foByEscapeHist << (lower+upper)/2.0 << " " << gsl_histogram_get(h,i) << std::endl;
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
    foByEscapeHist << ang << " " << gsl_histogram_get(h,index) << std::endl;
  }
  foByEscapeHist.close();
  foByEscapeHistPlt << "reset\nset style fill pattern 1\nset polar\nset angles degrees\nset size ratio -1\n\
\nset grid polar\n\nplot \'" << byEscapeTimeDatHist << "\' with filledcurves lw 2 notitle,\\\n \'" << byEscapeTimeDatHist << "\' with impulses lw 2 notitle\n";
//  foByEscapeHistPlt << "pause -1 \"Press Enter\"\n";
  foByEscapeHistPlt.close();
}

std::string
StatPostProcess::buildAtomByEnergy(const Float energyStep, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atom_by_energy_%.2f.dat", energyStep);

  const Float minEnergy =     0.0;
  const Float maxEnergy =  1000.0;
  const int n = (maxEnergy-minEnergy)/(energyStep);

  gsl_histogram * h = gsl_histogram_alloc (n);
  gsl_histogram_set_ranges_uniform (h, minEnergy, maxEnergy);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
      gsl_histogram_accumulate (h, energy/mol.atoms.size(), mol.atoms.size());
    }
  }

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

void
StatPostProcess::histEnergyByPolar(gsl_histogram* h, bool byAtom, FProcessClassicMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
      if (!byAtom)
        gsl_histogram_accumulate(h, polar, energy);
      else
        gsl_histogram_accumulate(h, polar, energy/mol.atoms.size());
    }
  }
}

void
StatPostProcess::histAtomsCountByPolar(gsl_histogram* h, FProcessClassicMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
/*
  Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
*/
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
      gsl_histogram_accumulate(h, polar, mol.atoms.size()/Float(trajData.size()));
    }
  }
}

void
StatPostProcess::histEnergyByPolarByAtomsInRange(gsl_histogram* h, FProcessClassicMolecule fpm) const
{
  const Float minPolar =   0.0;
  const Float maxPolar =  90.0;

  int n = gsl_histogram_bins(h);

  gsl_histogram * hEnergy = gsl_histogram_alloc (n);
  gsl_histogram * hCount  = gsl_histogram_alloc (n);

  gsl_histogram_set_ranges_uniform (hEnergy, minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (hCount, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(sqrt(SQR(v.x)+SQR(v.y)),-v.z)/mdtk::Deg;
      gsl_histogram_accumulate (hEnergy, polar, energy /* /mol.atoms.size()*/);
      gsl_histogram_accumulate (hCount, polar, mol.atoms.size());
    }
  }

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);
  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    if (gsl_histogram_get(hCount,i) > 0.0)
    gsl_histogram_accumulate(h, (lower+upper)/2.0,
                             gsl_histogram_get(hEnergy,i)/gsl_histogram_get(hCount,i));
  }

  gsl_histogram_free (hEnergy);
  gsl_histogram_free (hCount);
}

std::string
StatPostProcess::buildEnergyByPolar(const int n, bool byAtom, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_polar_%s%05d.dat", byAtom?"by_atom_":"",n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByPolar(h, byAtom, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

std::string
StatPostProcess::buildAtomsCountByPolar(const int n, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atomsCount_by_polar_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histAtomsCountByPolar(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

std::string
StatPostProcess::buildEnergyByPolarByAtomsInRange(const int n, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_polar_by_atomsInRange_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByPolarByAtomsInRange(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

void
StatPostProcess::histEnergyByAzimuth(gsl_histogram *h, bool byAtom, FProcessClassicMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(v.y,v.x)/mdtk::Deg;
      if (!byAtom)
        gsl_histogram_accumulate(h, polar, energy);
      else
        gsl_histogram_accumulate(h, polar, energy/mol.atoms.size());
    }
  }
}

void
StatPostProcess::histAtomsCountByAzimuth(gsl_histogram *h,  FProcessClassicMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
/*
  Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
*/
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(v.y,v.x)/mdtk::Deg;
      gsl_histogram_accumulate(h, polar, mol.atoms.size() / Float(trajData.size()));
    }
  }
}

void
StatPostProcess::histEnergyByAzimuthByAtomsInRange(gsl_histogram *h, FProcessClassicMolecule fpm) const
{
  const Float minPolar = -180.0;
  const Float maxPolar = +180.0;

  int n = gsl_histogram_bins(h);

  gsl_histogram * hEnergy = gsl_histogram_alloc (n);
  gsl_histogram * hCount  = gsl_histogram_alloc (n);

  gsl_histogram_set_ranges_uniform (hEnergy, minPolar, maxPolar);
  gsl_histogram_set_ranges_uniform (hCount, minPolar, maxPolar);

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
    {
      const ClassicMolecule& mol = td.molecules[mi];
      if (!fpm(mol)) continue;
      Float energy = 0.5*SQR(mol.getVelocity().module())*mol.getAMUMass()*mdtk::amu/mdtk::eV;
      mdtk::Vector3D v = mol.getVelocity();
      Float polar = atan2(v.y,v.x)/mdtk::Deg;
      gsl_histogram_accumulate (hEnergy, polar, energy /* /mol.atoms.size()*/);
      gsl_histogram_accumulate (hCount, polar, mol.atoms.size());
    }
  }

  gsl_histogram_set_ranges_uniform (h, minPolar, maxPolar);
  for(int i = 0; i < n; i++)
  {
    double lower, upper;
    gsl_histogram_get_range (h, i, &lower, &upper);
    if (gsl_histogram_get(hCount,i) > 0.0)
    gsl_histogram_accumulate(h, (lower+upper)/2.0,
                             gsl_histogram_get(hEnergy,i)/gsl_histogram_get(hCount,i));
  }

  gsl_histogram_free (hEnergy);
  gsl_histogram_free (hCount);

}

std::string
StatPostProcess::buildEnergyByAzimuth(const int n, bool byAtom, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_azimuth_%s%05d.dat",
    byAtom?"by_atom_":"",n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByAzimuth(h, byAtom, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

std::string
StatPostProcess::buildAtomsCountByAzimuth(const int n, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"atomsCount_by_azimuth_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histAtomsCountByAzimuth(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

std::string
StatPostProcess::buildEnergyByAzimuthByAtomsInRange(const int n, FProcessClassicMolecule fpm) const
{
  char ofilename[1024];
  sprintf(ofilename,"energy_by_azimuth_by_atomsInRange_%05d.dat",
    n);

  gsl_histogram * h = gsl_histogram_alloc (n);

  histEnergyByAzimuthByAtomsInRange(h, fpm);

  saveHistogram(h, ofilename);
  gsl_histogram_free (h);

  return std::string(ofilename);
}

void
appendFileToStream(std::ofstream& fo, std::string& filename)
{
  char tempStr[10000];
  std::ifstream fi(filename.c_str());
  while(fi.getline(tempStr, 10000-1, '\n')) {fo << tempStr << "\n";}
  fi.close();
}

void
StatPostProcess::buildAngular(FProcessClassicMolecule fpm) const
{
  int n_to_leave_pol = 3;
  int n_to_leave_az  = 12;
  std::string subdir = "_angular";

  yaatk::mkdir(subdir.c_str());
  yaatk::chdir(subdir.c_str());

  {
    std::string datFileName;
    char pltFileName[1000];
    sprintf(pltFileName,"!!!AtomByEnergy.plt");
    std::ofstream fo(pltFileName);
    fo << "\
reset\n\
set xrange [0:10]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
//    for(Float energyStep = 0.05; energyStep <= 10.0; energyStep += 0.05)
    for(Float energyStep = 0.1; energyStep <= 10.0; energyStep += 0.1)
    {
      datFileName = "-";
      fo << "\
reset\n\
set xrange [0:10]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildAtomByEnergy(energyStep,fpm);
      fo << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo, datFileName);
      std::remove(datFileName.c_str());
      std::remove((datFileName+".plt").c_str());
      fo << "e\n";
      fo << "pause -1 \"Press Return\" \n";
    }
    fo.close();
  }

  {
    std::string datFileName;
    char pltFileName[1000];
    sprintf(pltFileName,"!!!EnergyByPolar.plt");
    std::ofstream fo(pltFileName);
    sprintf(pltFileName,"!!!AtomsCountByPolar.plt");
    std::ofstream fo_count(pltFileName);
    sprintf(pltFileName,"!!!EnergyByPolarByAtom.plt");
    std::ofstream fo_by_atom(pltFileName);
    sprintf(pltFileName,"!!!EnergyByPolarByAtomsInRange.plt");
    std::ofstream fo_by_range(pltFileName);
    fo << "\
reset\n\
set xrange [0:90]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_count << "\
reset\n\
set xrange [0:90]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_by_atom << "\
reset\n\
set xrange [0:90]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_by_range << "\
reset\n\
set xrange [0:90]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
//    for(int n = 1; n <= 180; n++)
    for(int n = 1; n <= 90; n++)
    {
      datFileName = "-";
      fo << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByPolar(n,false,fpm);
      fo << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo, datFileName);
      if (n!=n_to_leave_pol)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo << "e\n";
      fo << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_count << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildAtomsCountByPolar(n,fpm);
      fo_count << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_count, datFileName);
      if (n!=n_to_leave_pol)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_count << "e\n";
      fo_count << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_by_atom << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByPolar(n,true,fpm);
      fo_by_atom << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_by_atom, datFileName);
      if (n!=n_to_leave_pol)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_by_atom << "e\n";
      fo_by_atom << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_by_range << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByPolarByAtomsInRange(n,fpm);
      fo_by_range << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_by_range, datFileName);
      if (n!=n_to_leave_pol)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_by_range << "e\n";
      fo_by_range << "pause -1 \"Press Return\" \n";
    }
    fo.close();
    fo_count.close();
    fo_by_atom.close();
    fo_by_range.close();
  }

  {
    std::string datFileName;
    char pltFileName[1000];
    sprintf(pltFileName,"!!!EnergyByAzimuth.plt");
    std::ofstream fo(pltFileName);
    sprintf(pltFileName,"!!!AtomsCountByAzimuth.plt");
    std::ofstream fo_count(pltFileName);
    sprintf(pltFileName,"!!!EnergyByAzimuthByAtom.plt");
    std::ofstream fo_by_atom(pltFileName);
    sprintf(pltFileName,"!!!EnergyByAzimuthByAtomsInRange.plt");
    std::ofstream fo_by_range(pltFileName);
    fo << "\
reset\n\
set xrange [-180:180]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_count << "\
reset\n\
set xrange [-180:180]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_by_atom << "\
reset\n\
set xrange [-180:180]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
    fo_by_range << "\
reset\n\
set xrange [-180:180]\n\
set yrange [0:*]\n\
set format y \"%10g\"\n\
";
//    for(int n = 1; n <= 180; n++)
    for(int n = 1; n <= 360; n++)
//    for(int n = 1; n <= 360*2; n++)
    {
      datFileName = "-";
      fo << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByAzimuth(n,false,fpm);
      fo << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo, datFileName);
      if (n!=n_to_leave_az)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo << "e\n";
      fo << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_count << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildAtomsCountByAzimuth(n,fpm);
      fo_count << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_count, datFileName);
      if (n!=n_to_leave_az)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_count << "e\n";
      fo_count << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_by_atom << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByAzimuth(n,true,fpm);
      fo_by_atom << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_by_atom, datFileName);
      if (n!=n_to_leave_az)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_by_atom << "e\n";
      fo_by_atom << "pause -1 \"Press Return\" \n";

      datFileName = "-";
      fo_by_range << "plot \'" << datFileName << "\' with boxes fs solid 1.0";
      datFileName =
        buildEnergyByAzimuthByAtomsInRange(n,fpm);
      fo_by_range << " title \"" << datFileName << "\"\n";
      appendFileToStream(fo_by_range, datFileName);
      if (n!=n_to_leave_az)
      {
        std::remove(datFileName.c_str());
        std::remove((datFileName+".plt").c_str());
      }
      fo_by_range << "e\n";
      fo_by_range << "pause -1 \"Press Return\" \n";
    }

    fo.close();
    fo_count.close();
    fo_by_atom.close();
    fo_by_range.close();
  }

  yaatk::chdir("..");
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
    const TrajData& td = trajData[trajIndex];
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
