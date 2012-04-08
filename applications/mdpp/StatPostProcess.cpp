/* 
   Molecular dynamics postprocessor, main classes

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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
   projectile(),
   target(),
   projectileEnergy()
{
  REQUIRE(s.substr(0,4) == "bomb");
  {
    size_t istart = 5;
    size_t iend = s.find("_",istart);
    REQUIRE(iend != std::string::npos);

    target = s.substr(istart,iend-istart);
    bool targetRecognized = false;
    if (target == "Fullerite")
      targetRecognized = true;
    if (target == "Cu")
      targetRecognized = true;
    REQUIRE(targetRecognized);
  }

  {
    size_t istart = s.find("_with_");
    REQUIRE(istart != std::string::npos);
    istart += 6;
    size_t iend = s.find("_",istart);
    REQUIRE(iend != std::string::npos);

    projectile = s.substr(istart,iend-istart);
    bool projectileRecognized = false;
    if (projectile == "C60")
      projectileRecognized = true;
    if (projectile == "Cu")
      projectileRecognized = true;
    REQUIRE(projectileRecognized);
  }

  {
    istringstream is(s.substr(s.size()-6,4));
    is >> projectileEnergy;
  }

  TRACE(str);
  TRACE(projectile);
  TRACE(target);
  TRACE(projectileEnergy);
  REQUIRE(s.size()>1);
  REQUIRE(*(s.end()-1)=='V');
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

    mdtk::SimLoop* state = new mdtk::SimLoop();
    state->allowToFreePotentials = true;
    setupPotentials(*state);

    std::string trajFinalName = trajData[trajIndex].trajDir+"mde_final";
    cout << "Loading state " << trajFinalName << std::endl;
    yaatk::text_ifstream fi(trajFinalName.c_str());
    state->initNLafterLoading = false;
    state->loadFromStream(fi);
    state->atoms.prepareForSimulatation();
    fi.close();

    {
      SimLoop::Check& check = state->check;
      Float Eo_plus_Eb = check.initialEnergy + check.energyTransferredFromBath;
      Float dE = check.currentEnergy - Eo_plus_Eb;
      Float dE_by_Eo_plus_Eb = dE/Eo_plus_Eb;
      if (fabs(dE_by_Eo_plus_Eb) > 0.01)
      {
        cerr << "Trajectory " << trajData[trajIndex].trajDir << " has bad energy conservation. Handle it separately." << endl;
        badTrajIndices.push_back(trajIndex);
      }
    }

    delete state;
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

  std::string mde_init_filename = trajData[0].trajDir+"mde_init";

  mdtk::SimLoop* mde_init = new mdtk::SimLoop();
  mde_init->allowToFreePotentials = true;
//  mde_init->allowToFreeAtoms = true;
  setupPotentials(*mde_init);
  yaatk::text_ifstream fi(mde_init_filename.c_str());
  mde_init->initNLafterLoading = false;
  mde_init->loadFromStream(fi/*,false*/);
  setTags(mde_init);
  mde_init->initNLafterLoading = true;
  fi.close();
  mdtk::AtomsArray& mde_init_atoms = mde_init->atoms;
  Float minInitZ = 1000000.0*mdtk::Ao;
  for(size_t ai = 0; ai < mde_init_atoms.size(); ai++)
  {
    const mdtk::Atom& init_atom = mde_init_atoms[ai];
    if (init_atom.coords.z < minInitZ && !isProjectileAtom(init_atom))
      minInitZ = init_atom.coords.z;
  }
  SPOTTED_DISTANCE = minInitZ - 0.05*mdtk::Ao;
  SPOTTED_DISTANCE = -4.55*mdtk::Ao;
  SPOTTED_DISTANCE = -1.0*mdtk::Ao;
  TRACE(SPOTTED_DISTANCE/mdtk::Ao);
  delete mde_init;
}

void
StatPostProcess::execute()
{
  using mdtk::Exception;

  removeBadTrajectories();

  cout << "PostProcess::execute() started." << std::endl;
  cerr << "PostProcess::execute() started." << std::endl;

  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    TrajData& td = trajData[trajIndex];

    mdtk::SimLoop* state = new mdtk::SimLoop();
    state->allowToFreePotentials = true;
    setupPotentials(*state);

    {
      std::string trajFinalName = trajData[trajIndex].trajDir+"mde_final";
      cout << "Loading state " << trajFinalName << std::endl;
      {
        yaatk::text_ifstream fi(trajFinalName.c_str());
        state->initNLafterLoading = false;
        state->loadFromStream(fi);
        fi.close();
      }
      state->atoms.prepareForSimulatation();
      setTags(state);
      cout << "State " << trajFinalName << " loaded." << std::endl;
    }
    NeighbourList nl(state->atoms);

    Molecule projectile;
    {
      projectile.buildByTag(*state,ATOMTAG_PROJECTILE);
      td.trajProjectile[state->simTime] = projectile;

      projectile.update(*state);
      td.trajProjectile[state->simTime] = projectile;
    }

    td.PBC = state->atoms.PBC();
    TRACE(td.PBC/mdtk::Ao);

//    TRACE(getAboveSpottedHeight(*state));

    buildSputteredClassicMolecules(*state,trajIndex,STATE_FINAL,nl);
//    buildClusterDynamics(*state,trajIndex,STATE_FINAL,nl);
//    buildProjectileDynamics(*state,trajIndex,STATE_FINAL);

//  if (td.molecules.size() > 0)

#if 0

//    if (0)
    {
      mdtk::SimLoop* mde_init = new mdtk::SimLoop();
      mde_init->allowToFreePotentials = true;
//      mde_init->allowToFreeAtoms = true;
      setupPotentials(*mde_init);
      std::string mde_init_filename = td.trajDir+"mde_init";
      cout << "Loading state " << mde_init_filename << std::endl;
      yaatk::text_ifstream fi(mde_init_filename.c_str());
      mde_init->initNLafterLoading = false;
      mde_init->loadFromStream(fi);
      mde_init->atoms.prepareForSimulatation();
      setTags(mde_init);
      NeighbourList nl(mde_init->atoms);
      fi.close();

      cout << "State " << mde_init_filename << " loaded." << std::endl;

      buildSputteredClassicMolecules(*mde_init,trajIndex,STATE_INIT,nl);
//      buildClusterDynamics(*mde_init,trajIndex,STATE_INIT,nl);
//      buildProjectileDynamics(*mde_init,trajIndex,STATE_INIT);

      {
        projectile.update(*mde_init);
        td.trajProjectile[mde_init->simTime] = projectile;
      }

      delete mde_init;
    }

//  if (td.molecules.size() > 0)
    {
      std::vector<std::string> interStates;
      findIntermediateStates(td.trajDir,interStates);

      mdtk::SimLoop* mde_inter = new mdtk::SimLoop();
      mde_inter->allowToFreePotentials = true;
//      mde_inter->allowToFreeAtoms = true;
      setupPotentials(*mde_inter);

      std::string trajFinalName = td.trajDir+"mde_init";
      cout << "Loading state " << trajFinalName << std::endl;
      yaatk::text_ifstream fi(trajFinalName.c_str());
      mde_inter->initNLafterLoading = false;
      mde_inter->loadFromStream(fi);
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
	mde_inter->loadFromStreamXVA(fi);
	mde_inter->atoms.prepareForSimulatation();
	setTags(mde_inter);
	NeighbourList nl(mde_inter->atoms);
	fi.close();

	buildSputteredClassicMolecules(*mde_inter,trajIndex,STATE_INTER,nl);
//	buildClusterDynamics(*mde_inter,trajIndex,STATE_INTER,nl);
//	buildProjectileDynamics(*mde_inter,trajIndex,STATE_INTER);

//        if (bombardingWithFullerene)
        {
          projectile.update(*mde_inter);
          td.trajProjectile[mde_inter->simTime] = projectile;
        }
      }
      delete mde_inter;
    }

#endif

    cout << "Building molecules for state done." << std::endl;

    delete state;
  }

  cout << "PostProcess::execute() done." << std::endl;
  cerr << "PostProcess::execute() done." << std::endl;
}

void
StatPostProcess::buildMassSpectrum() const
{
  std::ofstream fo("massSpectrum.txt");

  std::map<ClassicMolecule, size_t> massSpectrum;
  for(size_t trajIndex = 0; trajIndex < trajData.size(); trajIndex++)
  {
    const TrajData& td = trajData[trajIndex];
    for(size_t mi = 0; mi < td.molecules.size(); mi++)
      massSpectrum[td.molecules[mi]]++;
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

}
