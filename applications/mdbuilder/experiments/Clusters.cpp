/*
   Building of various clusters

   Copyright (C) 2007, 2008, 2010, 2011, 2012, 2013, 2015 Oleksandr
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

#include "Clusters.hpp"
#include "FCC.hpp"
#include "Graphite.hpp"
#include "Polyethylene.hpp"
#include "Fullerite.hpp"
#include "Fulleride.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_C60(AtomsArray& sl)
{
  yaatk::text_ifstream fcoords("C60.coords");
  if (!fcoords.isOpened())
    cerr << "Can't find file with C60 configuration." << endl;
  REQUIRE(fcoords.isOpened());
  size_t atomsCount;
  fcoords >> atomsCount;

  for(size_t i = 0; i < atomsCount; i++)
  {
    int index, tmp1, tmp2;
    Float x,y,z;
    fcoords >> index >> tmp1 >> tmp2 >> x >> y >> z;
    place(C_EL,sl,Vector3D(x*Ao,y*Ao,z*Ao));
  }

  fcoords.close();
}

AtomsArray
C60()
{
  AtomsArray atoms;

  glLoadIdentity();
  place_C60(atoms);

  quench(atoms,0.1*K);

  atoms.shiftToOrigin();

  return atoms;
}

void
buildFromAtom(const mdtk::Atom& a, const AtomsArray& atoms, std::set<size_t>& molecule)
{
  if (molecule.find(a.globalIndex) != molecule.end())
    return;

  molecule.insert(a.globalIndex);

  for(size_t i = 0; i < atoms.size(); i++)
  {
    const mdtk::Atom& nb_a = atoms[i];
    if (depos(a,nb_a).module() <= 7.0*Ao)
      buildFromAtom(nb_a,atoms,molecule);
  }
}

bool areUnparted(const AtomsArray& atoms)
{
  std::set<size_t> molecule;

  buildFromAtom(atoms[0],atoms,molecule);

  TRACE(atoms.size());
  TRACE(molecule.size());
  TRACE(molecule.size() == atoms.size());

  return molecule.size() != atoms.size();
}

struct OptiSnapshot
{
  AtomsArray atoms;
  Float T;
  OptiSnapshot(const AtomsArray &atomsArray = AtomsArray(),
               const Float temperature = 0.0*K)
    :atoms(atomsArray), T(temperature)
    {
    }
};

Float
optimize_single(SimLoop& simloop, gsl_rng* rng)
{
  yaatk::VerboseOutput vo(false);

  SimLoopNoMomentums mdloop(simloop);
  initialize_simloop(mdloop);

  const Float ENERGYINF = 1e6*eV;

  Float minPotEnergy = ENERGYINF;

  mdloop.simTime = 0.0*ps;
  mdloop.simTimeFinal = 0.0*ps;
  mdloop.iteration = 0;
  mdloop.dt = 1e-20;

  mdloop.simTimeSaveTrajInterval = 1000.0*ps;
  mdloop.iterationFlushStateInterval = 1000000;

  mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;

  bool stopHeating = false;
  std::vector<OptiSnapshot> snapshots;

  snapshots.push_back(OptiSnapshot(mdloop.atoms,0));

  {
    yaatk::ChDir cd("heating");

    mdloop.preventFileOutput = true;

    while (!stopHeating && mdloop.thermalBathCommon.To <= 10000.0*K)
    {
      Float T = mdloop.temperature();
      cerr << "To( " << mdloop.simTime/ps << " ps ) = "
           << mdloop.thermalBathCommon.To << " K" << endl;
      cerr << "T ( " << mdloop.simTime/ps << " ps ) = "
           << T << " K" << endl;

      Float dTo = 0.5*K;
      mdloop.thermalBathCommon.To = (dTo)/(1.0*ps)*mdloop.simTime;

      Float ToSnapshotInterval = 10.0*K;
      if (int(mdloop.thermalBathCommon.To/ToSnapshotInterval) !=
          int((mdloop.thermalBathCommon.To - dTo)/ToSnapshotInterval))
      {
        REQUIRE(mdloop.atoms.size() > 0);
        ElementID id = mdloop.atoms[0].ID;
        {
          snapshots.push_back(OptiSnapshot(mdloop.atoms,T));
          cerr << "Snapshot saved." << std::endl;
        }
      }

      mdloop.simTimeFinal += 1.0*ps;

      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();
      if (mdloop.execute())
        return ENERGYINF;
      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();

      if (areUnparted(mdloop.atoms))
        stopHeating = true;
    }

    mdloop.preventFileOutput = false;
    mdloop.writestate();
  }

  VETRACE(mdloop.thermalBathCommon.To);
  VETRACE(snapshots.size());
  REQUIRE(snapshots.size() > 0);

  size_t minPotEnergySnapshotIndex = 0;

  for(size_t i = 0; i < snapshots.size(); ++i)
  {
    std::ostringstream dirname;
    dirname << "snapshot";
    PRINT2STREAM_FW(dirname, i, '0', 10);
    dirname << "-";
    PRINT2STREAM_FWP(dirname, snapshots[i].T, '0', 8, 4);
    dirname << "K";
    yaatk::ChDir cd(dirname.str());

    mdloop.atoms = snapshots[i].atoms;
    mdloop.thermalBathCommon.To = snapshots[i].T;

    mdloop.simTime = 0.0*ps;
    mdloop.simTimeFinal = 0.0*ps;
    mdloop.iteration = 0;
    mdloop.dt = 1e-20;

    mdloop.simTimeSaveTrajInterval = 1000.0*ps;
    mdloop.iterationFlushStateInterval = 1000000;

    mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;

    mdloop.writetrajXVA();
    mdloop.preventFileOutput = true;

    bool stopCooling = false;

    while (!stopCooling && mdloop.thermalBathCommon.To > 2.0*K)
    {
      Float T = mdloop.energyKin()/(3.0/2.0*kb*mdloop.atoms.size());
      cerr << "To( " << mdloop.simTime/ps << " ps ) = "
           << mdloop.thermalBathCommon.To << " K" << endl;
      cerr << "T ( " << mdloop.simTime/ps << " ps ) = "
           << T << " K" << endl;

      mdloop.thermalBathCommon.To -= 10.0*K;
      if (mdloop.thermalBathCommon.To < 0)
        mdloop.thermalBathCommon.To = 0;

      mdloop.simTimeFinal += 1.0*ps;

      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();
      if (mdloop.execute())
        return ENERGYINF;
      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();

      if (areUnparted(mdloop.atoms))
        stopCooling = true;
    }

    mdloop.thermalBathCommon.To = 0.0;
    mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;
    mdloop.dt = 1e-20;

    while (!stopCooling && mdloop.temperature() > 0.001*K)
    {
      if (areUnparted(mdloop.atoms))
        stopCooling = true;

      mdloop.simTimeFinal += 0.05*ps;

      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();
      if (mdloop.execute())
        return ENERGYINF;
      mdloop.atoms.removeMomentum();
      mdloop.atoms.removeAngularMomentum();
    }

    mdloop.preventFileOutput = false;
    mdloop.writestate();

    std::ofstream fo("mdloop.opti");
    mdloop.saveToStream(fo);
    fo.close();

    SimLoopSaver mds(mdloop);
    mds.write("optimal");
    mds.removeAttributesButPosVel("optimal");

    Float energyValue = mdloop.energyPot();
    if (energyValue < minPotEnergy)
    {
      simloop.atoms = mdloop.atoms;
      minPotEnergy = energyValue;
      minPotEnergySnapshotIndex = i;
    }

    std::ofstream fo_energyValue("mdloop.opti.energy");
    fo_energyValue << energyValue/eV << " "
                   << energyValue/eV/mdloop.atoms.size() << std::endl;
    fo_energyValue.close();

    {
      std::ostringstream fname;
      fname << "../mdloop.opti-";
      PRINT2STREAM_FW(fname, i, '0', 10);
      fname << "-";
      PRINT2STREAM_P(fname, energyValue/eV/mdloop.atoms.size(), 10);
      fname << "_eV_per_atom";

      std::ofstream fo(fname.str().c_str());
      mdloop.saveToStream(fo);
      fo.close();
    }
  }

  FETRACE(minPotEnergySnapshotIndex);

  {
    std::ofstream fo("mdloop.opti.last");
    mdloop.saveToStream(fo);
    fo.close();
  }
  {
    std::ofstream fo("mdloop.opti.best");
    simloop.saveToStream(fo);
    fo.close();
  }

  return minPotEnergy;
}

void
add1atomInit(AtomsArray& atoms, ElementID id)
{
  Atom newAtom;
  newAtom.ID = id;
  newAtom.setAttributesByElementID();
  atoms.push_back(newAtom);
}

void
findMinMax(AtomsArray& atoms, Float xm[3][2], size_t xmi[3][2])
{
  size_t xi, xs;

  for(xi = 0; xi < 3; xi++)
  {
    for(xs = 0; xs < 2; xs++)
    {
      Float& xM = xm[xi][xs];
      size_t& xMIndex = xmi[xi][xs];
      xM = atoms[0].coords.X(xi);
      xMIndex = 0;
      for(size_t i = 0; i < atoms.size(); i++)
      {
        if (xs == 0)
          if (atoms[i].coords.X(xi) < xM)
          {
            xM = atoms[i].coords.X(xi);xMIndex = i;
          }
        if (xs == 1)
          if (atoms[i].coords.X(xi) > xM)
          {
            xM = atoms[i].coords.X(xi);xMIndex = i;
          }
      }
    }
  }
}

void
add1atom(AtomsArray& atoms, ElementID id, gsl_rng* r)
{
  Float  xm[3][2];
  size_t xmi[3][2];

  findMinMax(atoms, xm, xmi);

  size_t xi, xs;

  for(size_t i = 0; i < atoms.size(); i++)
  {
    atoms[i].coords.x -= (xm[0][0]+xm[0][1])/2.0;
    atoms[i].coords.y -= (xm[1][0]+xm[1][1])/2.0;
    atoms[i].coords.z -= (xm[2][0]+xm[2][1])/2.0;
  }

  findMinMax(atoms, xm, xmi);

  for(xi = 0; xi < 3; xi++)
  {
    for(xs = 0; xs < 2; xs++)
    {
      fprintf(stderr,"xm[%lu][%lu]=%f\n",xi,xs,xm[xi][xs]/Ao);
    }
  }

  xi = gsl_rng_get(r)%3;
  xs = gsl_rng_get(r)%2;

  Float sign = (xs==0)?-1:+1;

  Float& xM = xm[xi][xs];
  size_t& xMIndex = xmi[xi][xs];

  if (mdtk::verboseTrace)
  {
    ETRACE(xi);
    ETRACE(xs);
    ETRACE(sign);
    ETRACE(xM/Ao);
    ETRACE(xMIndex);
  }

  Atom newAtom;
  newAtom.ID  = id;
  newAtom.setAttributesByElementID();
  Vector3D vn(gsl_rng_uniform(r),
              gsl_rng_uniform(r),
              gsl_rng_uniform(r));
  vn.normalize();

  Vector3D dv(0,0,0);
  dv.X(xi) = 2.29*Ao*sign;
  if (mdtk::verboseTrace)
    ETRACE(dv);
  dv += 0.03*Ao*vn;
  if (mdtk::verboseTrace)
    ETRACE(dv);

  newAtom.coords = atoms[xMIndex].coords + dv;

  if (mdtk::verboseTrace)
    ETRACE(newAtom.coords/Ao);

  atoms.push_back(newAtom);
}

AtomsArray
cluster(ElementID id, int clusterSize)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc (T);
  REQUIRE(r != NULL);

  gsl_rng_set(r, 123);

  REQUIRE(gsl_rng_min(r) == 0);
  REQUIRE(gsl_rng_max(r) > 1000);

  SimLoop sl;
  initialize_simloop(sl);

  if (clusterSize == 0) return sl.atoms;

  yaatk::ChDir cd("_tmp-Cu-clusters");

  std::ofstream foGlobal("energy.min.all",std::ios::app);

  if (sl.atoms.size() > 0)
    add1atom(sl.atoms,id,r);
  else
    add1atomInit(sl.atoms, id);

  for(int atomsCount = sl.atoms.size();
      atomsCount <= clusterSize;
      atomsCount++)
  {
    if (mdtk::verboseTrace)
      ETRACE(sl.atoms.size());

    sl.executeDryRun();

    char dirname[100];
    sprintf(dirname,"%03lu",sl.atoms.size());
    yaatk::ChDir cd(dirname);

    Float minPotEnergyOf = optimize_single(sl, r);

    foGlobal << std::setw (10) << sl.atoms.size() << " "
             << std::setw (20) << minPotEnergyOf/eV << " "
             << std::setw (20) << minPotEnergyOf/eV/sl.atoms.size() << std::endl;

    if (atomsCount < clusterSize)
      add1atom(sl.atoms, id, r);
  }

  foGlobal.close();

  gsl_rng_free (r);

  sl.atoms.shiftToOrigin();

  return sl.atoms;
}

AtomsArray
clusterFromCrystal(const AtomsArray& atoms, int clusterSize, Vector3D c)
{
  REQUIRE(atoms.size() >= clusterSize);

  const gsl_rng_type * T;
  gsl_rng * r;

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc (T);
  REQUIRE(r != NULL);

  gsl_rng_set(r, 123);

  REQUIRE(gsl_rng_min(r) == 0);
  REQUIRE(gsl_rng_max(r) > 1000);

  SimLoop sl;
  REQUIRE(!sl.atoms.PBCEnabled());
  initialize_simloop(sl);

  if (clusterSize == 0) return sl.atoms;

  {
    if (c == NO_PBC)
      c = atoms.geomCenter();

    Float cutRadius = 0.01*Ao;
    std::vector<bool> alreadyAdded(atoms.size());
    while (1)
    {
      for(size_t i = 0; i < atoms.size() && sl.atoms.size() < clusterSize; ++i)
        if ((atoms[i].coords - c).module() <= cutRadius && !alreadyAdded[i])
        {
          sl.atoms.push_back(atoms[i]);
          alreadyAdded[i] = true;
        }
      if (sl.atoms.size() == clusterSize)
        break;
      cutRadius += 0.1*Ao;
      REQUIRE(cutRadius <= 1000.0*Ao);
    }
  }

  yaatk::ChDir cd("_tmp-clusters");

  std::ofstream foGlobal("energy.min.all",std::ios::app);

  {
    if (mdtk::verboseTrace)
      ETRACE(sl.atoms.size());

    sl.executeDryRun();

    char dirname[100];
    sprintf(dirname,"%03lu",sl.atoms.size());
    yaatk::ChDir cd(dirname);

    Float minPotEnergyOf = optimize_single(sl, r);

    foGlobal << std::setw (10) << sl.atoms.size() << " "
             << std::setw (20) << minPotEnergyOf/eV << " "
             << std::setw (20) << minPotEnergyOf/eV/sl.atoms.size() << std::endl;
  }

  foGlobal.close();

  gsl_rng_free (r);

  sl.atoms.shiftToOrigin();

  return sl.atoms;
}

AtomsArray
clusterFromFCCCrystal(ElementID id, int clusterSize)
{
  int num = ceil(pow(clusterSize/4.0,1.0/3.0));
  REQUIRE(num >= 1 && num <= 100);
  num++;
  if (num%2 != 0)
    num++;
  double a;
  switch (id)
  {
  case Cu_EL : a = 3.615*Ao; break;
  case Au_EL : a = 4.0781*Ao; break;
  case Ag_EL : a = 4.086*Ao; break;
  default: throw Exception("Unknown FCC element.");
  }

  AtomsArray atoms;
  place_FCC_lattice(atoms,num,num,num,id,false,a,a,a);

  return clusterFromCrystal(atoms,clusterSize,Vector3D(num/2*a,num/2*a,num/2*a));
}

AtomsArray
embed(AtomsArray cluster, AtomsArray shell)
{
  cluster.shiftToOrigin();
  shell.shiftToOrigin();

  AtomsArray atoms = cluster;

//shrinking
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& a = atoms[i];
    a.coords *= 0.8;
  }
//~shrinking

  atoms.addAtoms(shell);

  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& a = atoms[i];
    a.V = 0.0;
  }

  quench(atoms,0.1*K);

  atoms.removeMomentum();

  return atoms;
}

void
add_rotational_motion(
  AtomsArray& atoms,
  Float totalRotEnergy,
  Vector3D rotAxis
  )
{
  Vector3D clusterCenterOfM = atoms.massCenter();

  Float sumMVSQR = 0.0;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& a = atoms[i];
    Vector3D r = a.coords - clusterCenterOfM;
    sumMVSQR += a.M*SQR(vectormul(r,rotAxis).module());
  }

  Float maxVrot = sqrt(2.0*totalRotEnergy/sumMVSQR);
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& a = atoms[i];
    Vector3D r = a.coords - clusterCenterOfM;
    a.V = vectormul(r,rotAxis)*maxVrot;
  }

  SimLoop sl;
  initialize_simloop(sl);

  sl.atoms = atoms;

  TRACE(sl.energyKin()/eV);
  sl.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_NONE;
  relax/*_flush*/(sl,0.2*ps);
  TRACE(sl.energyKin()/eV);
  sl.atoms.removeMomentum();
  TRACE(sl.energyKin()/eV);

  atoms = sl.atoms;
}

SimLoop
build_target_by_cluster_bombardment(
  const SimLoop& sl_target,
  AtomsArray cluster,
  Float clusterEnergy,
  Float interactionDistance
  )
{
  cluster.shiftToOrigin();

  Float clusterRadius
    = cluster.maxDistanceFrom(cluster.geomCenter());
  TRACE(clusterRadius/Ao);

  Vector3D dCluster;
  if (sl_target.atoms.PBCEnabled())
  {
    dCluster = sl_target.atoms.PBC()/2.0;
  }
  else
  {
    if (sl_target.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_SPHERE)
      dCluster = sl_target.thermalBathGeomSphere.center;
    else
      dCluster = sl_target.atoms.dimensions().center();
  }
  dCluster.z = -(10.0*Ao+clusterRadius);
  TRACE(dCluster/Ao);

  for(size_t i = 0; i < cluster.size(); i++)
  {
    Atom& a = cluster[i];
    a.coords += dCluster;
  }

  for(int it = 0;;it = 1)
  {
    TRACE("Landing fullerene ...");
    bool areInteracting = false;
    for(size_t cli = 0; !areInteracting && cli < cluster.size(); cli++)
      for(size_t surfi = 0; !areInteracting && surfi < sl_target.atoms.size(); surfi++)
      {
        Atom& clusterAtom = cluster[cli];
        const Atom& surfaceAtom = sl_target.atoms[surfi];
        if (depos(clusterAtom,surfaceAtom).module() < (/*(clusterEnergy==0.0)?3.0*Ao:*/ /*6.0*Ao*/interactionDistance))
        {
          areInteracting = true;
        }
      }
    if (areInteracting)
    {
      if (it == 0)
      {
        TRACE("*****ERROR:Fullerene is already interacting");throw;
      }
      break;
    }
    for(size_t i = 0; i < cluster.size(); i++)
    {
      Atom& a = cluster[i];
      a.coords.z += 0.05*Ao;
    }
    TRACE(cluster[0].coords.z/Ao);
  };
  TRACE("Fullerene landed ok.");

  for(size_t i = 0; i < cluster.size(); i++)
  {
    Atom& a = cluster[i];
    a.V += Vector3D(0,0,sqrt(2.0*clusterEnergy/(cluster.mass())));
  }

  SimLoop sl;
  initialize_simloop(sl);

  sl = sl_target;
  sl.atoms.addAtoms(cluster);

//  removeMomentum(sl.atoms);

  return sl;
}

void
prepare_Cu_by_Cu_at_C60_bobardment()
{
  std::vector<int> cluster_sizes;
  cluster_sizes.push_back(0);
  cluster_sizes.push_back(1);
  cluster_sizes.push_back(6);
  cluster_sizes.push_back(13);

  std::vector<Float> trans_energies;
  for(Float e = 0; e <= 400; e += 50)
    trans_energies.push_back(e*eV);
  trans_energies.push_back(10.0*eV);

  std::vector<Float> rot_energies;
  for(Float e = 0; e <= 100; e += 10)
    rot_energies.push_back(e*eV);

  std::vector<Vector3D> rot_axes;
  rot_axes.push_back(Vector3D(0,0,1));
  rot_axes.push_back(Vector3D(0,1,0));
  rot_axes.push_back(Vector3D(0,1,1));

  SimLoop sl_Cu = build_FCC_lattice(20,20,10,Cu_EL);
  AtomsArray C60 = mdbuilder::C60();

  for(size_t i_cluster_size = 0; i_cluster_size < cluster_sizes.size(); i_cluster_size++)
  {
    int& cluster_size = cluster_sizes[i_cluster_size];
    AtomsArray cluster = mdbuilder::cluster(Cu_EL,cluster_size);

    for(size_t i_rot_axis = 0; i_rot_axis < rot_axes.size(); i_rot_axis++)
    {
      for(size_t i_rot_energy = 0; i_rot_energy < rot_energies.size(); i_rot_energy++)
      {
        AtomsArray endo = embed(cluster,C60);
        Vector3D& rot_axis = rot_axes[i_rot_axis];
        Float& rot_energy = rot_energies[i_rot_energy];
        add_rotational_motion(endo,rot_energy,rot_axis);

        for(size_t i_trans_energy = 0; i_trans_energy < trans_energies.size(); i_trans_energy++)
        {
          Float& trans_energy = trans_energies[i_trans_energy];
          SimLoop sl = build_target_by_cluster_bombardment(sl_Cu,endo,trans_energy);

          TRACE(sl.energyKin()/eV);

          sl.iteration = 0;
          sl.simTime = 0.0*ps;
          sl.simTimeFinal = 10.0*ps;
          sl.simTimeSaveTrajInterval = 0.05*ps;

          char id_string[1000];
          sprintf(id_string,
                  "%s_by_Cu%02d@%s_%s%03deV_(%01d%01d%01d)%03deV",
                  "Cu",
                  cluster_size,
                  "C60",
                  "n",
                  int(trans_energy/eV),
                  int(rot_axis.x),int(rot_axis.y),int(rot_axis.z),
                  int(rot_energy/eV));

          yaatk::ChDir cd(id_string);
          yaatk::text_ofstream fomde("in.mde");
          sl.saveToMDE(fomde);
          fomde.close();
        }
      }
    }
  }
}

void
prepare_Graphite_by_Cu_at_C60_bombardment()
{
  std::vector<int> cluster_sizes;
  cluster_sizes.push_back(0);
//  cluster_sizes.push_back(1);
//  cluster_sizes.push_back(6);
//  cluster_sizes.push_back(13);

  std::vector<Float> trans_energies;
  for(Float e = 0; e <= 400; e += 50)
    trans_energies.push_back(e*eV);
  trans_energies.push_back(10.0*eV);

  std::vector<Float> rot_energies;
  for(Float e = 0; e <= 100; e += 10)
    rot_energies.push_back(e*eV);

  std::vector<Vector3D> rot_axes;
  rot_axes.push_back(Vector3D(0,0,1));
//  rot_axes.push_back(Vector3D(0,1,0));
//  rot_axes.push_back(Vector3D(0,1,1));

  SimLoop sl_Graphite = build_Graphite_lattice(12,14,3);
  AtomsArray C60 = mdbuilder::C60();

  for(size_t i_cluster_size = 0; i_cluster_size < cluster_sizes.size(); i_cluster_size++)
  {
    int& cluster_size = cluster_sizes[i_cluster_size];
    AtomsArray cluster = mdbuilder::cluster(Cu_EL,cluster_size);

    for(size_t i_rot_axis = 0; i_rot_axis < rot_axes.size(); i_rot_axis++)
    {
      for(size_t i_rot_energy = 0; i_rot_energy < rot_energies.size(); i_rot_energy++)
      {
        AtomsArray endo = embed(cluster,C60);
        Vector3D& rot_axis = rot_axes[i_rot_axis];
        Float& rot_energy = rot_energies[i_rot_energy];
        add_rotational_motion(endo,rot_energy,rot_axis);

        for(size_t i_trans_energy = 0; i_trans_energy < trans_energies.size(); i_trans_energy++)
        {
          Float& trans_energy = trans_energies[i_trans_energy];
          SimLoop sl = build_target_by_cluster_bombardment(sl_Graphite,endo,trans_energy);

          TRACE(sl.energyKin()/eV);

          sl.iteration = 0;
          sl.simTime = 0.0*ps;
          sl.simTimeFinal = 10.0*ps;
          sl.simTimeSaveTrajInterval = 0.05*ps;

          char id_string[1000];
          sprintf(id_string,
                  "%s_by_Cu%02d@%s_%s%03deV_(%01d%01d%01d)%03deV",
                  "Cu",
                  cluster_size,
                  "C60",
                  "n",
                  int(trans_energy/eV),
                  int(rot_axis.x),int(rot_axis.y),int(rot_axis.z),
                  int(rot_energy/eV));

          yaatk::ChDir cd(id_string);
          yaatk::text_ofstream fomde("in.mde");
          sl.saveToMDE(fomde);
          fomde.close();
        }
      }
    }
  }
}

SimLoop
build_Cluster_Landed_on_Substrate(
  const SimLoop sl_Substrate,
  AtomsArray cluster,
  bool applyPBCtoCluster
  )
{
  std::ostringstream sbuildSubdir;
  sbuildSubdir << "_build_" << ElementIDtoString(cluster[0].ID) << cluster.size() << "_on_PE";
  yaatk::ChDir cd(sbuildSubdir.str());

  cluster.tag(ATOMTAG_CLUSTER);
  cluster.removeMomentum();

  SimLoop sl;
  initialize_simloop(sl);
  sl = build_target_by_cluster_bombardment(sl_Substrate,cluster,0.0*eV,3.3*Ao);

  TRACE(sl.energyKin()/eV);

  {
    yaatk::text_ofstream fomdloop("000.mdloop");
    sl.saveToStream(fomdloop);
    fomdloop.close();
  }


  // std::vector<size_t> fixedAtoms = sl.atoms.fixUnfixedCHAtoms(0,sl.atoms.size());
  // if (0)
  // {
  //   Float tb_zMin_bak = sl.thermalBath.zMin;
  //   sl.thermalBath.zMin = -1.0*Ao;

  //   relax(sl,5.0*ps,"001-relax");

  //   sl.thermalBath.zMin = tb_zMin_bak;
  // }
  // sl.atoms.unfixAtoms(fixedAtoms);

  relax(sl,40.0*ps,"010-relax");
  quench(sl,0.01*K,200*ps,0.01*ps,"011-quench");

  relax(sl,10.0*ps,"020-relax");
  quench(sl,0.01*K,200*ps,0.01*ps,"021-quench");

  if (!applyPBCtoCluster) // TODO: should be forced if PBC are disabled
  {
    for(size_t i = 0; i < sl.atoms.size(); ++i)
    {
      Atom& a = sl.atoms[i];
      if (a.ID == cluster[0].ID)
      {
        REQUIRE(a.PBC_count.x == 0);
        REQUIRE(a.PBC_count.y == 0);
        REQUIRE(a.PBC_count.z == 0);
        a.PBC = NO_PBC;
      }
    }
  }

  sl.atoms.removeMomentum();

  return sl;
}

void
bomb_Cluster_with_Ions(
  std::string dirname,
  const SimLoop& target,
  std::vector<size_t> clusterAtomIndices,
  ElementID ionElement,
  Float ionEnergy,
  std::set<Float> halos,
  size_t numberOfImpacts
  )
{
  REQUIRE(halos.size() > 0);

  yaatk::ChDir cd(dirname);

  SimLoop sl(target);
  sl.atoms.tag(ATOMTAG_TARGET);

  TRACE(clusterAtomIndices.size());
  REQUIRE(clusterAtomIndices.size() > 0);

  for(size_t i = 0; i < clusterAtomIndices.size(); i++)
  {
    TRACE(i);
    Atom& clusterAtom = sl.atoms[clusterAtomIndices[i]];
  }

  const Atom& clusterAtom = target.atoms[clusterAtomIndices[0]];

  Float clusterXMax = clusterAtom.coords.x;
  Float clusterXMin = clusterAtom.coords.x;
  Float clusterYMax = clusterAtom.coords.y;
  Float clusterYMin = clusterAtom.coords.y;
  Float clusterZMax = clusterAtom.coords.z;
  Float clusterZMin = clusterAtom.coords.z;

  for(size_t i = 0; i < clusterAtomIndices.size(); i++)
  {
    TRACE(i);
    const Atom& clusterAtom = target.atoms[clusterAtomIndices[i]];

    if (clusterAtom.coords.x > clusterXMax)
      clusterXMax = clusterAtom.coords.x;
    if (clusterAtom.coords.x < clusterXMin)
      clusterXMin = clusterAtom.coords.x;
    if (clusterAtom.coords.y > clusterYMax)
      clusterYMax = clusterAtom.coords.y;
    if (clusterAtom.coords.y < clusterYMin)
      clusterYMin = clusterAtom.coords.y;
    if (clusterAtom.coords.z > clusterZMax)
      clusterZMax = clusterAtom.coords.z;
    if (clusterAtom.coords.z < clusterZMin)
      clusterZMin = clusterAtom.coords.z;
  }

  {
    Float halo = *halos.rbegin();
    clusterXMin -= halo;
    clusterXMax += halo;
    clusterYMin -= halo;
    clusterYMax += halo;
  }

  Float a = clusterXMax - clusterXMin;
  Float b = clusterYMax - clusterYMin;

  if (sl.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_BOX)
  {
    REQUIRE(sl.atoms.PBCEnabled());
    REQUIRE(clusterXMin > 0.0 + sl.thermalBathGeomBox.dBoundary && clusterXMax < sl.atoms.PBC().x - sl.thermalBathGeomBox.dBoundary);
    REQUIRE(clusterYMin > 0.0 + sl.thermalBathGeomBox.dBoundary && clusterYMax < sl.atoms.PBC().y - sl.thermalBathGeomBox.dBoundary);
  }
  else
  {
    mdtk::AtomsArray::Dimensions dim = sl.atoms.dimensions();
    REQUIRE(clusterXMin > dim.x_min && clusterXMax < dim.x_max);
    REQUIRE(clusterYMin > dim.y_min && clusterYMax < dim.y_max);
  }

//  yaatk::ChDir cd_dataset("dataset");

  {
    Atom projectile(ionElement,Vector3D(0,0,clusterZMin-5.5*Ao));
    projectile.V = Vector3D(0,0,sqrt(2.0*ionEnergy/(projectile.M)));
    projectile.tag(ATOMTAG_PROJECTILE);
    sl.atoms.push_back(projectile);

    sl.iteration = 0;
    sl.simTime = 0.0*ps;
    sl.simTimeFinal = 10.0*ps;
    sl.simTimeSaveTrajInterval = 0.1*ps;

    SimLoopSaver mds(sl);
    mds.write("base");
  }

  for(std::set<Float>::iterator hi = halos.begin(); hi != halos.end(); ++hi)
  {
    std::ostringstream ossHaloId;
    ossHaloId << "halo-" << *hi/mdtk::Ao << "Ao";
    yaatk::ChDir cd(ossHaloId.str());

    gsl_qrng * coord2d_qrng = gsl_qrng_alloc (/*gsl_qrng_sobol*/ gsl_qrng_niederreiter_2, 2);
    REQUIRE(coord2d_qrng != NULL);

    std::ofstream rngUsedTrace("rng.used.trace");
    REQUIRE(rngUsedTrace != NULL);
    std::ofstream rngSkippedTrace("rng.skipped.trace");
    REQUIRE(rngSkippedTrace != NULL);

    yaatk::binary_ofstream ionpos_bin("ionpos.bin");
    yaatk::text_ofstream ionpos_txt("ionpos.txt");

    for(int trajIndex = 0; trajIndex < numberOfImpacts; trajIndex++)
    {
      Float bombX = 0.0;
      Float bombY = 0.0;
      bool allowBomb = false;
      Float cell_part_x;
      Float cell_part_y;
      do
      {
        double v[2];
        gsl_qrng_get (coord2d_qrng, v);
        cell_part_x = v[0];
        cell_part_y = v[1];

        allowBomb = false;
        Float dMin = 1e6*Ao;
        for(size_t i = 0; i < clusterAtomIndices.size(); i++)
        {
          const Atom& clusterAtom = target.atoms[clusterAtomIndices[i]];

          bombX = clusterXMin + cell_part_x*(a+b);
          bombY = clusterYMin + cell_part_y*(a+b);

          Float d = distance(clusterAtom.coords,
                             Vector3D(bombX,bombY,clusterAtom.coords.z));

          if (d < dMin)
            dMin = d;
        }
        REQUIRE(dMin < 1e5*Ao);

        if (hi == halos.begin())
        {
          if (dMin < *hi)
            allowBomb = true;
        }
        else
        {
          std::set<Float>::iterator hprev = hi;
          --hprev;
          if (dMin < *hi && dMin >= *hprev)
            allowBomb = true;
        }

        if ((cell_part_x >= a/(a+b) || cell_part_y >= b/(a+b)) || !allowBomb)
          rngSkippedTrace << cell_part_x << " " << cell_part_y << "\n";

      }while ( (cell_part_x >= a/(a+b) || cell_part_y >= b/(a+b)) || !allowBomb);
      rngUsedTrace << cell_part_x << " " << cell_part_y << "\n";

      if (sl.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_BOX)
      {
        REQUIRE(bombX > 0.0 + sl.thermalBathGeomBox.dBoundary && bombX < sl.atoms.PBC().x - sl.thermalBathGeomBox.dBoundary);
        REQUIRE(bombY > 0.0 + sl.thermalBathGeomBox.dBoundary && bombY < sl.atoms.PBC().y - sl.thermalBathGeomBox.dBoundary);
      }
      else
      {
        mdtk::AtomsArray::Dimensions dim = sl.atoms.dimensions();
        REQUIRE(bombX > dim.x_min && bombX < dim.x_max);
        REQUIRE(bombY > dim.y_min && bombY < dim.y_max);
      }

      YAATK_BIN_WRITE(ionpos_bin,bombX);
      YAATK_BIN_WRITE(ionpos_bin,bombY);
      ionpos_txt << bombX << " "
                 << bombY << "\n";
    }

    ionpos_txt.close();
    ionpos_bin.close();

    gsl_qrng_free (coord2d_qrng);
    rngUsedTrace.close();
    rngSkippedTrace.close();
  }
}

void
bomb_landedCluster_with_Ions(
  const SimLoop& target,
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  size_t numberOfImpacts
  )
{
  yaatk::VerboseOutput vo(true);

  std::vector<size_t> clusterAtomIndices;
  for(size_t ai = 0; ai < target.atoms.size(); ++ai)
    if (target.atoms[ai].hasTag(ATOMTAG_CLUSTER))
      clusterAtomIndices.push_back(ai);
  REQUIRE(clusterAtomIndices.size() > 0);
  ElementID clusterElement = target.atoms[clusterAtomIndices[0]].ID;


  for(size_t ionEnergyIndex = 0;
      ionEnergyIndex < ionEnergies.size();
      ++ionEnergyIndex)
  {
    Float ionEnergy = ionEnergies[ionEnergyIndex];
    for(size_t ionElementIndex = 0;
        ionElementIndex < ionElements.size();
        ++ionElementIndex)
    {
      ElementID ionElement = ionElements[ionElementIndex];

      char id_string[1000];
      sprintf(id_string,
              "%s%03d_on_PE_by_%s_%04deV",
              ElementIDtoString(clusterElement).c_str(),
              int(clusterAtomIndices.size()),
              ElementIDtoString(ionElement).c_str(),
              int(ionEnergy/eV));
      std::string dirname(id_string);

      std::set<Float> halos;
      Float maxHalo = 6.0*Ao;
      {
        AtomsArray::Dimensions dim = target.atoms.dimensions();
        Float minLateralSize = min2(dim.x_len,dim.y_len);
        maxHalo = minLateralSize/4.0;
      }
      REQUIRE(maxHalo > 5.5*Ao);
      halos.insert(5.5*Ao);
      halos.insert(1.3*Ao);
      TRACE(clusterElement);
      TRACE(ionElement);
      if (clusterElement == Cu_EL && ionElement == Ar_EL)
      {
        halos.clear();
        halos.insert(5.5*Ao);
        halos.insert(1.243*Ao); // ~40 eV
        // halos.insert(1.326*Ao); // ~30 eV
      }
      if (clusterElement == Au_EL && ionElement == Ar_EL)
      {
        halos.clear();
        halos.insert(5.5*Ao);
        halos.insert(1.397*Ao); // ~40 eV
        // halos.insert(1.481*Ao); // ~30 eV
      }
      if (clusterElement == Cu_EL && ionElement == Xe_EL)
      {
        halos.clear();
        halos.insert(5.5*Ao);
        halos.insert(1.425*Ao); // ~40 eV
        // halos.insert(1.511*Ao); // ~30 eV
      }
      if (clusterElement == Au_EL && ionElement == Xe_EL)
      {
        halos.clear();
        halos.insert(5.5*Ao);
        halos.insert(1.582*Ao); // ~40 eV
        // halos.insert(1.668*Ao); // ~30 eV
      }
      for(Float halo = 10.0*Ao; halo <= maxHalo; halo += 10.0*Ao)
        halos.insert(halo);
      for(std::set<Float>::iterator halo = halos.begin(); halo != halos.end(); ++halo)
        TRACE(*halo/Ao);

      bomb_Cluster_with_Ions(dirname,
                             target,
                             clusterAtomIndices,
                             ionElement, ionEnergy,
                             halos,
                             numberOfImpacts);
    }
  }
}

void
bomb_orthorhombic_with_clusters(
  std::string dirname,
  SimLoop cluster,
  const SimLoop target,
  int a_num,
  int b_num,
  double a,
  double b,
  size_t numberOfImpacts
  )
{
  yaatk::ChDir cd(dirname);

  std::ofstream rngSelected("rng.selected.log");
  REQUIRE(rngSelected != NULL);
  std::ofstream rngExcluded("rng.excluded.log");
  REQUIRE(rngExcluded != NULL);
  std::ofstream bombXY("bombXY.log");
  REQUIRE(bombXY != NULL);

  {
    bombXY << 0 << " " << 0 << "\n"
           << 0 << " " << b*b_num/Ao << "\n"
           << a*a_num/Ao << " " << b*b_num/Ao << "\n"
           << a*a_num/Ao << " " << 0 << "\n";
  }

  Float bombX0 = a*(a_num/2.0)-a/2.0;
  Float bombY0 = b*(b_num/2.0)-b/2.0;

  gsl_qrng* qrng_2d_pos = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
  REQUIRE(qrng_2d_pos != NULL);

  yaatk::ChDir cd_dataset("dataset");

  Float allowToBomb = false;

  for(int trajIndex = 0; trajIndex < numberOfImpacts; trajIndex++)
  {
    Float bombX = 0.0;
    Float bombY = 0.0;
    Float d = (a>b)?a:b;
    Float cell_part_x;
    Float cell_part_y;
    bool positionNotFound = true;

    while (positionNotFound)
    {
      double v[2];
      gsl_qrng_get(qrng_2d_pos, v);
      cell_part_x = v[0];
      cell_part_y = v[1];

      bombX = bombX0 + cell_part_x*d;
      bombY = bombY0 + cell_part_y*d;

      allowToBomb = true;

      positionNotFound =
        (cell_part_x >= a/d || cell_part_y >= b/d) || !allowToBomb;

      if (positionNotFound)
        rngExcluded << cell_part_x << " " << cell_part_y << "\n";
    }

    rngSelected << cell_part_x << " " << cell_part_y << "\n";

    Vector3D initialClusterPosition
      (bombX,
       bombY,
       target.atoms.dimensions().z_min - (5.5*Ao + cluster.atoms.radius()));

    cluster.atoms.shiftToPosition(initialClusterPosition);

    bombXY << cluster.atoms.massCenter().x/Ao << " "
           << cluster.atoms.massCenter().y/Ao << "\n";

    SimLoop sl(target);
    sl.add_simloop(cluster);

    for(size_t i = 0; i < sl.atoms.size(); i++)
    {
//      sl.atoms[i]->apply_PBC=true;
      sl.atoms[i].apply_ThermalBath=true;
    }

    char trajDirName[1024];
    sprintf(trajDirName,"%08d",trajIndex);
    yaatk::ChDir cd(trajDirName);
    {
      sl.forgetHistory();
      sl.atoms.prepareForSimulatation();
      sl.simTime = 0.0*ps;
      sl.simTimeFinal = 6.0*ps;
      sl.simTimeSaveTrajInterval = 0.1*ps;

      yaatk::text_ofstream fomde("mde_init");
      sl.saveToStream(fomde);
      fomde.close();
    }
  }

  rngSelected.close();
  rngExcluded.close();
  bombXY.close();
  gsl_qrng_free(qrng_2d_pos);
}

void
build_FCC_metal_bombardment_with_ions(
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  size_t numberOfImpacts,
  int a_num,
  int b_num,
  int c_num,
  double a,
  double b,
  double c,
  ElementID metalElement
  )
{
  SimLoop sl_target =
    build_FCC_lattice(a_num,b_num,c_num,metalElement,true,a,b,c);
  {
    for(size_t ionEnergyIndex = 0;
        ionEnergyIndex < ionEnergies.size();
        ++ionEnergyIndex)
    {
      Float ionEnergy = ionEnergies[ionEnergyIndex];
      for(size_t ionElementIndex = 0;
          ionElementIndex < ionElements.size();
          ++ionElementIndex)
      {
        ElementID ionElement = ionElements[ionElementIndex];

        char id_string[1000];
        sprintf(id_string,
                "bomb_%s_with_%s_%04deV",
                ElementIDtoString(metalElement).c_str(),
                ElementIDtoString(ionElement).c_str(),
                int(ionEnergy/eV));
        std::string dirname(id_string);

        SimLoop sl_ion;
        Atom atom(ionElement);
        atom.V = Vector3D(0,0,sqrt(2.0*ionEnergy/(atom.M)));
        sl_ion.atoms.push_back(atom);

        bomb_orthorhombic_with_clusters(dirname,
                                        sl_ion,
                                        sl_target,
                                        a_num,b_num,
                                        a,b,
                                        numberOfImpacts);
      }
    }
  }
}

void
build_FCC_metal_bombardment_with_C60(
  std::vector<Float> fullereneEnergies,
  size_t numberOfImpacts,
  int a_num,
  int b_num,
  int c_num,
  double a,
  double b,
  double c,
  ElementID metalElement
  )
{
  SimLoop sl_target =
    build_FCC_lattice(a_num,b_num,c_num,metalElement,true,a,b,c);

  AtomsArray fullerene = mdbuilder::C60();
  fullerene.removeMomentum();
  fullerene.shiftToOrigin();

  for(size_t energyIndex = 0;
      energyIndex < fullereneEnergies.size();
      ++energyIndex)
  {
    Float fullereneEnergy = fullereneEnergies[energyIndex];

    char id_string[1000];
    sprintf(id_string,
            "bomb_%s_with_C60_%04deV",
            ElementIDtoString(metalElement).c_str(),
            int(fullereneEnergy/eV));
    std::string dirname(id_string);

    SimLoop sl_energeticFullerene;
    sl_energeticFullerene.atoms = fullerene;
    sl_energeticFullerene.atoms.setTranslationalEnergy(
      fullereneEnergy,
      Vector3D(0,0,1));

    bomb_orthorhombic_with_clusters(dirname,
                                    sl_energeticFullerene,
                                    sl_target,
                                    a_num,b_num,
                                    a,b,
                                    numberOfImpacts);
  }
}

void
build_fullerite_bombardment_with_ions(
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  size_t numberOfImpacts,
  int a_num,
  int b_num,
  int c_num,
  double a
  )
{
  Float b = a;
  SimLoop sl_target =
    build_Fullerite_C60(a_num,b_num,c_num,true,a);
  {
    for(size_t ionEnergyIndex = 0;
        ionEnergyIndex < ionEnergies.size();
        ++ionEnergyIndex)
    {
      Float ionEnergy = ionEnergies[ionEnergyIndex];
      for(size_t ionElementIndex = 0;
          ionElementIndex < ionElements.size();
          ++ionElementIndex)
      {
        ElementID ionElement = ionElements[ionElementIndex];

        char id_string[1000];
        sprintf(id_string,
                "bomb_%s_with_%s_%04deV",
                "Fullerite",
                ElementIDtoString(ionElement).c_str(),
                int(ionEnergy/eV));
        std::string dirname(id_string);

        SimLoop sl_ion;
        Atom atom(ionElement);
        atom.V = Vector3D(0,0,sqrt(2.0*ionEnergy/(atom.M)));
        sl_ion.atoms.push_back(atom);

        bomb_orthorhombic_with_clusters(dirname,
                                        sl_ion,
                                        sl_target,
                                        a_num,b_num,
                                        a,b,
                                        numberOfImpacts);
      }
    }
  }
}

std::string
metal_C60_mixing_id(std::string what,
                    std::string withWhat,
                    Float energy)
{
  char id_string[1000];
  sprintf(id_string,
          "bomb_%s_with_%s_%04deV",
          what.c_str(),
          withWhat.c_str(),
          int(energy/eV));
  return id_string;
}

void
build_metal_C60_mixing(
  std::vector<Float> impactEnergies,
  ElementID metalElement
  )
{
  size_t numberOfImpacts = 128;

  AtomsArray fullerene = mdbuilder::C60();
  fullerene.removeMomentum();
  fullerene.shiftToOrigin();
  fullerene.tag(ATOMTAG_PROJECTILE | ATOMTAG_FULLERENE | ATOMTAG_CLUSTER);

  SimLoop sl_metalAtom;
  sl_metalAtom.atoms.push_back(Atom(metalElement));
  sl_metalAtom.atoms.removeMomentum();
  sl_metalAtom.atoms.shiftToOrigin();
  sl_metalAtom.atoms.tag(ATOMTAG_PROJECTILE | ATOMTAG_MONOMER);

  int a_num_metal;
  int b_num_metal;
  int c_num_metal;
  Float a_metal;
  Float b_metal;
  Float c_metal;
  SimLoop sl_metalCrystal;
  switch (metalElement)
  {
  case Cu_EL :
    a_num_metal = 12;
    b_num_metal = 12;
    c_num_metal = 12;
    a_metal = 3.615*Ao;
    b_metal = 3.615*Ao;
    c_metal = 3.615*Ao;
    sl_metalCrystal = build_FCC_lattice(
      a_num_metal,b_num_metal,c_num_metal,
      metalElement,
      true,
      a_metal,
      b_metal,
      c_metal);
    break;
  case Ag_EL :
    a_num_metal = 12;
    b_num_metal = 12;
    c_num_metal = 12;
    a_metal = 4.086*Ao;
    b_metal = 4.086*Ao;
    c_metal = 4.086*Ao;
    sl_metalCrystal = build_FCC_lattice(
      a_num_metal,b_num_metal,c_num_metal,
      metalElement,
      true,
      a_metal,
      b_metal,
      c_metal);
    break;
  case Au_EL :
    a_num_metal = 12;
    b_num_metal = 12;
    c_num_metal = 12;
    a_metal = 4.078*Ao;
    b_metal = 4.078*Ao;
    c_metal = 4.078*Ao;
    sl_metalCrystal = build_FCC_lattice(
      a_num_metal,b_num_metal,c_num_metal,
      metalElement,
      true,
      a_metal,
      b_metal,
      c_metal);
    break;
  default:
    throw Exception("Unknow metal!");
  }
  sl_metalCrystal.atoms.tag(ATOMTAG_TARGET | ATOMTAG_SUBSTRATE);

  int a_num_fullerite = 3;
  int b_num_fullerite = 3;
  int c_num_fullerite = 3;
  Float a_fullerite = 14.17*Ao;
  Float b_fullerite = 14.17*Ao;
  Float c_fullerite = 14.17*Ao;
  SimLoop sl_fulleriteCrystal =
    build_Fullerite_C60(
      a_num_fullerite,
      b_num_fullerite,
      c_num_fullerite,
      true,
      a_fullerite);
  sl_fulleriteCrystal.atoms.tag(ATOMTAG_TARGET | ATOMTAG_SUBSTRATE);

  for(size_t energyIndex = 0;
      energyIndex < impactEnergies.size();
      ++energyIndex)
  {
    Float impactEnergy = impactEnergies[energyIndex];

    SimLoop sl_energeticFullerene;
    sl_energeticFullerene.atoms = fullerene;
    sl_energeticFullerene.atoms.setTranslationalEnergy(
      impactEnergy,Vector3D(0,0,1));

    SimLoop sl_energeticMetalAtom(sl_metalAtom);
    sl_energeticMetalAtom.atoms.setTranslationalEnergy(
      impactEnergy,
      Vector3D(0,0,1));

    {
      std::string id = metal_C60_mixing_id(
        ElementIDtoString(metalElement),
        ElementIDtoString(metalElement),
        impactEnergy);
      bomb_orthorhombic_with_clusters(
        id,
        sl_energeticMetalAtom,
        sl_metalCrystal,
        a_num_metal,b_num_metal,
        a_metal,b_metal,
        numberOfImpacts);
    }
    {
      std::string id = metal_C60_mixing_id(
        ElementIDtoString(metalElement),
        "C60",
        impactEnergy);
      bomb_orthorhombic_with_clusters(
        id,
        sl_energeticFullerene,
        sl_metalCrystal,
        a_num_metal,b_num_metal,
        a_metal,b_metal,
        numberOfImpacts);
    }
    {
      std::string id = metal_C60_mixing_id(
        "Fullerite",
        ElementIDtoString(metalElement),
        impactEnergy);
      bomb_orthorhombic_with_clusters(
        id,
        sl_energeticMetalAtom,
        sl_fulleriteCrystal,
        a_num_fullerite,b_num_fullerite,
        a_fullerite,b_fullerite,
        numberOfImpacts);
    }
    {
      std::string id = metal_C60_mixing_id(
        "Fullerite",
        "C60",
        impactEnergy);
      bomb_orthorhombic_with_clusters(
        id,
        sl_energeticFullerene,
        sl_fulleriteCrystal,
        a_num_fullerite,b_num_fullerite,
        a_fullerite,b_fullerite,
        numberOfImpacts);
    }
  }
}

}
