/*
   Building of various clusters

   Copyright (C) 2007, 2008, 2010, 2011, 2012 Oleksandr Yermolenko
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

  quench(atoms,1.0*K);

  atoms.shiftToOrigin();

  return atoms;
}

void
checkOnEnergyPotMin(SimLoop& simloop, Float& minPotEnergy)
{
  Float currentPotEnergy = simloop.energyPot();
  TRACE(currentPotEnergy/eV);
  if (currentPotEnergy < minPotEnergy)
  {
    minPotEnergy = currentPotEnergy;
    TRACE(minPotEnergy/eV);
    TRACE(minPotEnergy/eV/simloop.atoms.size());
    {
      std::ofstream fo("in.mde.min");
      simloop.saveToMDE(fo);
      fo.close();
    }
    {
      std::ofstream fo("energy.min",std::ios::app);
      fo << minPotEnergy/eV/simloop.atoms.size() << " "
         << minPotEnergy/eV/simloop.atoms.size() << std::endl;
      fo.close();
    }
    ERRTRACE(minPotEnergy/eV/simloop.atoms.size());
  }
}

Float
optimize_single(SimLoop *modloop, gsl_rng* rng)
{
  const Float ENERGYINF = 1e6*eV;

  Float minPotEnergy = ENERGYINF;
  TRACE(minPotEnergy/eV);

  Float freezingEnergy = 0.0001*2.5*eV*modloop->atoms.size();
  TRACE(freezingEnergy/eV);
  ERRTRACE(freezingEnergy/eV);

  for(Float maxHeatUpEnergy = 0.0001*eV;
      maxHeatUpEnergy <= 0.75*eV;
      maxHeatUpEnergy += 0.05*eV)
  {
    modloop->simTimeSaveTrajInterval = 1000.0*ps;
    modloop->iterationFlushStateInterval = 1000000;

    modloop->thermalBath.zMin = +100000.0*Ao;
    modloop->simTime = 0.0*ps;
    modloop->simTimeFinal = 1.0*ps;
    modloop->dt = 1e-20;
    modloop->iteration = 0;
    cerr << "Heating up every atom to " << maxHeatUpEnergy/eV << " eV" << std::endl;
    modloop->heatUpEveryAtom(maxHeatUpEnergy, rng);
    TRACE(modloop->atoms.PBC());
    int retval;
    cerr << "Releasing..." << std::endl;
    retval = modloop->execute();
    if (retval) return ENERGYINF;
    ERRTRACE(modloop->simTime/ps);
    checkOnEnergyPotMin(*modloop,minPotEnergy);

    cerr << "Cooling..." << std::endl;

    modloop->thermalBath.zMin = -100000.0*Ao;
    modloop->simTime = 0.0*ps;
    modloop->simTimeFinal = 0.0*ps;
    modloop->dt = 1e-20;
    modloop->iteration = 0;
    TRACE(modloop->atoms.PBC());
    while (modloop->simTimeFinal < 4.0*ps)
    {
      modloop->simTimeFinal += 0.05*ps;
      retval = modloop->execute();
      ERRTRACE(modloop->simTime/ps);
      if (modloop->energyKin() < freezingEnergy) break;
    }
    if (retval) return ENERGYINF;
    checkOnEnergyPotMin(*modloop,minPotEnergy);
  }
  TRACE(minPotEnergy/eV);
  TRACE(minPotEnergy/eV/modloop->atoms.size());

  {
    yaatk::text_ifstream fi("in.mde.min");
    modloop->loadFromMDE(fi);
    fi.close();
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

  ERRTRACE(xi);
  ERRTRACE(xs);
  ERRTRACE(sign);
  ERRTRACE(xM/Ao);
  ERRTRACE(xMIndex);

  Atom newAtom;
  newAtom.ID  = id;
  newAtom.setAttributesByElementID();
  Vector3D vn(gsl_rng_uniform(r),
              gsl_rng_uniform(r),
              gsl_rng_uniform(r));
  vn.normalize();

  Vector3D dv(0,0,0);
  dv.X(xi) = 2.29*Ao*sign;
  ERRTRACE(dv);
  dv += 0.03*Ao*vn;
  ERRTRACE(dv);

  newAtom.coords = atoms[xMIndex].coords + dv;

  ERRTRACE(newAtom.coords/Ao);

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

  yaatk::mkdir("_tmp-Cu-clusters");
  yaatk::chdir("_tmp-Cu-clusters");

  std::ofstream foGlobal("energy.min.all",std::ios::app);

  if (sl.atoms.size() > 0)
    add1atom(sl.atoms,id,r);
  else
    add1atomInit(sl.atoms, id);

  for(int atomsCount = sl.atoms.size();
      atomsCount <= clusterSize;
      atomsCount++)
  {
    ERRTRACE(sl.atoms.size());

    sl.initialize();

    char dirname[100];
    sprintf(dirname,"%03lu",sl.atoms.size());
    yaatk::mkdir(dirname);
    yaatk::chdir(dirname);

    Float minPotEnergyOf = optimize_single(&sl, r);

    foGlobal << std::setw (10) << sl.atoms.size() << " "
             << std::setw (20) << minPotEnergyOf/eV << " "
             << std::setw (20) << minPotEnergyOf/eV/sl.atoms.size() << std::endl;

    yaatk::chdir("..");

    if (atomsCount < clusterSize)
      add1atom(sl.atoms, id, r);
  }

  foGlobal.close();

  yaatk::chdir("..");

  gsl_rng_free (r);

  sl.atoms.shiftToOrigin();

  return sl.atoms;
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

  quench(atoms,1.0*K);

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
  sl.thermalBath.zMin = +100000.0*Ao;
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
  Float clusterEnergy
  )
{
  cluster.shiftToOrigin();

  Float clusterRadius
    = cluster.maxDistanceFrom(cluster.geomCenter());
  TRACE(clusterRadius/Ao);

  Vector3D dCluster = sl_target.atoms.PBC()/2.0;
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
        if (depos(clusterAtom,surfaceAtom).module() < (/*(clusterEnergy==0.0)?3.0*Ao:*/6.0*Ao))
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

          yaatk::mkdir(id_string);
          yaatk::chdir(id_string);
          yaatk::text_ofstream fomde("in.mde");
          sl.saveToMDE(fomde);
          fomde.close();
          yaatk::chdir("..");
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

          yaatk::mkdir(id_string);
          yaatk::chdir(id_string);
          yaatk::text_ofstream fomde("in.mde");
          sl.saveToMDE(fomde);
          fomde.close();
          yaatk::chdir("..");
        }
      }
    }
  }
}

SimLoop
build_Cluster_Landed_on_Substrate(
  const SimLoop sl_Substrate,
  ElementID id,
  int clusterSize
  )
{
  AtomsArray Cluster = cluster(id,clusterSize);
  Cluster.removeMomentum();

  SimLoop sl;
  initialize_simloop(sl);
  sl = build_target_by_cluster_bombardment(sl_Substrate,Cluster,0.2*eV*clusterSize);

  TRACE(sl.energyKin()/eV);

  {
    yaatk::text_ofstream fomde("_tmp-X-relax_flush-landing.mde");
    sl.saveToMDE(fomde);
    fomde.close();
  }


  std::vector<size_t> fixedAtoms = sl.atoms.fixUnfixedCHAtoms(0,sl.atoms.size());
  {
    Float tb_zMin_bak = sl.thermalBath.zMin;
    sl.thermalBath.zMin = -1.0*Ao;

    relax/*_flush*/(sl,5.0*ps,"_tmp-X-landing-fixed-CH-relax_flush");

    sl.thermalBath.zMin = tb_zMin_bak;
  }
  sl.atoms.unfixAtoms(fixedAtoms);

  relax/*_flush*/(sl,5.0*ps,"_tmp-X-landing-unfixed-CH-relax_flush");

  quench(sl,0.01*K);

  return sl;
}

void
bomb_Cluster_with_Ions(
  std::string dirname,
  const SimLoop& target,
  std::vector<size_t> clusterAtomIndices,
  ElementID ionElement,
  Float ionEnergy,
  Float halo,
  size_t numberOfImpacts
  )
{
  yaatk::mkdir(dirname.c_str());
  yaatk::chdir(dirname.c_str());

  SimLoop sl(target);

  std::ofstream rngout("rng.out");
  REQUIRE(rngout != NULL);
  std::ofstream rngNOTout("rng.NOT.out");
  REQUIRE(rngNOTout != NULL);

  TRACE(clusterAtomIndices.size());
  REQUIRE(clusterAtomIndices.size() > 0);

  for(size_t i = 0; i < clusterAtomIndices.size(); i++)
  {
    TRACE(i);
    Atom& clusterAtom = sl.atoms[clusterAtomIndices[i]];

//    clusterAtom.apply_PBC = false;
    clusterAtom.apply_ThermalBath = false;
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

  clusterXMin -= halo;
  clusterXMax += halo;
  clusterYMin -= halo;
  clusterYMax += halo;
  Float a = clusterXMax - clusterXMin;
  Float b = clusterYMax - clusterYMin;

  gsl_qrng * coord2d_qrng = gsl_qrng_alloc (/*gsl_qrng_sobol*/ gsl_qrng_niederreiter_2, 2);
  REQUIRE(coord2d_qrng != NULL);

  yaatk::mkdir("dataset");
  yaatk::chdir("dataset");

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
      for(size_t i = 0; i < clusterAtomIndices.size(); i++)
      {
        const Atom& clusterAtom = target.atoms[clusterAtomIndices[i]];

        bombX = clusterXMin + cell_part_x*(a+b);
        bombY = clusterYMin + cell_part_y*(a+b);

        if (distance(clusterAtom.coords,
                     Vector3D(bombX,bombY,clusterAtom.coords.z))
            < halo)
        {
          allowBomb = true; break;
        }
      }

      if ((cell_part_x >= a/(a+b) || cell_part_y >= b/(a+b)) || !allowBomb)
        rngNOTout << cell_part_x << " " << cell_part_y << "\n";

    }while ( (cell_part_x >= a/(a+b) || cell_part_y >= b/(a+b)) || !allowBomb);
    rngout << cell_part_x << " " << cell_part_y << "\n";

    Atom projectile(ionElement,Vector3D(bombX,bombY,clusterZMin-5.5*Ao));

    projectile.V = Vector3D(0,0,sqrt(2.0*ionEnergy/(projectile.M)));
//    projectile->apply_PBC=false;
    projectile.apply_ThermalBath=false;
    sl.atoms.push_back(projectile);

    char trajDirName[1024];
    sprintf(trajDirName,"%08d",trajIndex);
    yaatk::mkdir(trajDirName);
    yaatk::chdir(trajDirName);

    sl.simTime = 0.0*ps;
    sl.simTimeFinal = 10.0*ps;
    sl.simTimeSaveTrajInterval = 0.1*ps;

    yaatk::text_ofstream fomde("in.mde");
    sl.saveToMDE(fomde);
    fomde.close();
    yaatk::chdir("..");

    sl.atoms.resize(sl.atoms.size()-1);
  }

  yaatk::chdir("..");

  rngout.close();
  rngNOTout.close();
  gsl_qrng_free (coord2d_qrng);

  yaatk::chdir("..");
}

void
bomb_MetalCluster_on_Polyethylene_with_Ions(
  int a_num,
  int b_num,
  int c_num,
  std::vector<int> clusterSizes,
  std::vector<ElementID> clusterElements,
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies
  )
{
  SimLoop sl_Polyethylene =
    build_Polyethylene_lattice_with_folds(a_num,b_num,c_num);
  for(size_t clusterElementIndex = 0;
      clusterElementIndex < clusterElements.size();
      ++clusterElementIndex)
  {
    ElementID clusterElement = clusterElements[clusterElementIndex];
    for(size_t sizeIndex = 0;
        sizeIndex < clusterSizes.size();
        ++sizeIndex)
    {
      int clusterSize = clusterSizes[sizeIndex];
      SimLoop sl_Landed =
        build_Cluster_Landed_on_Substrate(sl_Polyethylene,
                                          clusterElement,
                                          clusterSize);
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
                  clusterSize,
                  ElementIDtoString(ionElement).c_str(),
                  int(ionEnergy/eV));
          std::string dirname(id_string);

          Float halo = 5.5*Ao;

          std::vector<size_t> clusterAtomIndices;
          for(size_t ai = 0; ai < sl_Landed.atoms.size(); ++ai)
            if (sl_Landed.atoms[ai].ID == clusterElement)
              clusterAtomIndices.push_back(ai);

          bomb_Cluster_with_Ions(dirname,
                                 sl_Landed,
                                 clusterAtomIndices,
                                 ionElement, ionEnergy,
                                 halo);
        }
      }
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
  yaatk::mkdir(dirname.c_str());
  yaatk::chdir(dirname.c_str());

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

  yaatk::mkdir("dataset");
  yaatk::chdir("dataset");

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
    yaatk::mkdir(trajDirName);
    yaatk::chdir(trajDirName);
    {
      sl.simTime = 0.0*ps;
      sl.simTimeFinal = 6.0*ps;
      sl.simTimeSaveTrajInterval = 0.1*ps;

      yaatk::text_ofstream fomde("mde_init");
      sl.saveToStream(fomde);
      fomde.close();
    }
    yaatk::chdir("..");
  }

  yaatk::chdir("..");

  rngSelected.close();
  rngExcluded.close();
  bombXY.close();
  gsl_qrng_free(qrng_2d_pos);

  yaatk::chdir("..");
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
    sl_energeticFullerene.atoms.addTranslationalEnergy(
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
    sl_energeticFullerene.atoms.addTranslationalEnergy(
      impactEnergy,Vector3D(0,0,1));

    SimLoop sl_energeticMetalAtom(sl_metalAtom);
    sl_energeticMetalAtom.atoms.addTranslationalEnergy(
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
