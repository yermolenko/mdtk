/*
   Building of various clusters

   Copyright (C) 2007, 2008, 2010, 2011 Oleksandr Yermolenko
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

#ifndef MDBUILDER_Clusters_HPP
#define MDBUILDER_Clusters_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

inline
void
place_C60(mdtk::SimLoop& sl)
{
  yaatk::text_ifstream fcoords("C60.coords");
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

inline
mdtk::SimLoop
build_C60_optimized()
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  glLoadIdentity();
  mdbuilder::place_C60(sl);
  
  quench(sl,1.0*K);

  return sl;
}

inline
void
checkOnEnergyPotMin(SimLoop& simloop,Float& minPotEnergy)
{
  Float currentPotEnergy = simloop.energyPot();
  TRACE(currentPotEnergy/mdtk::eV);
  if (currentPotEnergy < minPotEnergy)
  {
    minPotEnergy = currentPotEnergy;
    TRACE(minPotEnergy/mdtk::eV);
    TRACE(minPotEnergy/mdtk::eV/simloop.atoms_.size());
    { 
      std::ofstream fo("in.mde.min"); 
      simloop.saveToMDE(fo); 
      fo.close(); 
    }
    { 
      std::ofstream fo("energy.min",std::ios::app); 
      fo << minPotEnergy/mdtk::eV/simloop.atoms_.size() << " "
         << minPotEnergy/mdtk::eV/simloop.atoms_.size() << std::endl;
      fo.close(); 
    }
    ERRTRACE(minPotEnergy/mdtk::eV/simloop.atoms_.size());
  }
}

inline
Float
optimize_single(SimLoop *modloop, gsl_rng* rng)
{
  const Float ENERGYINF = 1e6*eV;

  Float minPotEnergy = ENERGYINF;
  TRACE(minPotEnergy/mdtk::eV);
  
  Float freezingEnergy = 0.0001*2.5*mdtk::eV*modloop->atoms_.size();
  TRACE(freezingEnergy/mdtk::eV);
  ERRTRACE(freezingEnergy/mdtk::eV);
  
  for(Float maxHeatUpEnergy = 0.0001*eV;
      maxHeatUpEnergy <= 0.75*eV;
      maxHeatUpEnergy += 0.05*eV)
  {
    modloop->simTimeSaveTrajInterval = 1000.0*ps;
    modloop->iterationFlushStateInterval = 1000000;
    
    modloop->thermalBath.zMin = +100000.0*Ao;
    modloop->simTime = 0.0*ps;
    modloop->simTimeFinal = 1.0*ps;
    modloop->dt_ = 1e-20;
    modloop->iteration = 0;
    cerr << "Heating up every atom to " << maxHeatUpEnergy/mdtk::eV << " eV" << std::endl;
    modloop->heatUpEveryAtom(maxHeatUpEnergy, rng);
    TRACE(modloop->getPBC());
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
    modloop->dt_ = 1e-20;
    modloop->iteration = 0;
    TRACE(modloop->getPBC());
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
  TRACE(minPotEnergy/mdtk::eV);
  TRACE(minPotEnergy/mdtk::eV/modloop->atoms_.size());

  { 
    yaatk::text_ifstream fi("in.mde.min");
    modloop->loadFromMDE(fi); 
    fi.close(); 
  }
  
  return minPotEnergy;
}

inline
void
add1atomInit(mdtk::AtomsContainer& atoms, ElementID id)
{
  mdtk::Atom* newAtom = new Atom;
  newAtom->ID = id;
  newAtom->setAttributesByElementID();
  atoms.push_back(newAtom);
}

inline
void
findMinMax(mdtk::AtomsContainer& atoms, Float xm[3][2], size_t xmi[3][2])
{
  size_t xi, xs;
  
  for(xi = 0; xi < 3; xi++)
  {
    for(xs = 0; xs < 2; xs++)
    {
      Float& xM = xm[xi][xs];
      size_t& xMIndex = xmi[xi][xs];
      xM = atoms[0]->coords.X(xi);
      xMIndex = 0;
      for(size_t i = 0; i < atoms.size(); i++)
      {
        if (xs == 0)
          if (atoms[i]->coords.X(xi) < xM)
          {
            xM = atoms[i]->coords.X(xi);xMIndex = i;
          }
        if (xs == 1)
          if (atoms[i]->coords.X(xi) > xM)
          {
            xM = atoms[i]->coords.X(xi);xMIndex = i;
          }
      }
    }
  }
}

inline
void 
add1atom(mdtk::AtomsContainer& atoms, ElementID id, gsl_rng* r)
{
  Float  xm[3][2];
  size_t xmi[3][2];

  findMinMax(atoms, xm, xmi);  

  size_t xi, xs;
  
  for(size_t i = 0; i < atoms.size(); i++)
  {
    atoms[i]->coords.x -= (xm[0][0]+xm[0][1])/2.0;
    atoms[i]->coords.y -= (xm[1][0]+xm[1][1])/2.0;
    atoms[i]->coords.z -= (xm[2][0]+xm[2][1])/2.0;
  }

  findMinMax(atoms, xm, xmi);  

  for(xi = 0; xi < 3; xi++)
  {
    for(xs = 0; xs < 2; xs++)
    {
      fprintf(stderr,"xm[%lu][%lu]=%f\n",xi,xs,xm[xi][xs]/mdtk::Ao);
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
  ERRTRACE(xM/mdtk::Ao);
  ERRTRACE(xMIndex);

  mdtk::Atom* newAtom = new Atom;
  newAtom->ID  = id;
  newAtom->setAttributesByElementID();
  Vector3D vn(gsl_rng_uniform(r),
              gsl_rng_uniform(r),
              gsl_rng_uniform(r));
  vn.normalize();
  
  Vector3D dv(0,0,0);
  dv.X(xi) = 2.29*Ao*sign;
  ERRTRACE(dv);
  dv += 0.03*Ao*vn;
  ERRTRACE(dv);
  
  newAtom->coords = atoms[xMIndex]->coords + dv;
  
  ERRTRACE(newAtom->coords/mdtk::Ao);
  
  atoms.push_back(newAtom);
}

inline
SimLoop
build_cluster(ElementID id, int clusterSize)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  T = gsl_rng_ranlxd2;
  r = gsl_rng_alloc (T);
  REQUIRE(r != NULL);

  gsl_rng_set(r, 123);

  REQUIRE(gsl_rng_min(r) == 0);
  REQUIRE(gsl_rng_max(r) > 1000);

  mdtk::SimLoop sl;
  initialize_simloop(sl);

  if (clusterSize == 0) return sl;

  yaatk::mkdir("_tmp-Cu-clusters");
  yaatk::chdir("_tmp-Cu-clusters");

  std::ofstream foGlobal("energy.min.all",std::ios::app); 

  if (sl.atoms.size() > 0)
    add1atom(sl.atoms_,id,r);
  else
    add1atomInit(sl.atoms_, id);

  for(int atomsCount = sl.atoms_.size(); 
      atomsCount <= clusterSize; 
      atomsCount++)
  {    
    ERRTRACE(sl.atoms_.size());

    sl.initialize();

    char dirname[100];
    sprintf(dirname,"%03lu",sl.atoms_.size());
    yaatk::mkdir(dirname);
    yaatk::chdir(dirname);

    Float minPotEnergyOf = optimize_single(&sl, r);

    foGlobal << std::setw (10) << sl.atoms_.size() << " "
             << std::setw (20) << minPotEnergyOf/mdtk::eV << " "
             << std::setw (20) << minPotEnergyOf/mdtk::eV/sl.atoms_.size() << std::endl;

    yaatk::chdir("..");

    if (atomsCount < clusterSize)
      add1atom(sl.atoms_, id, r);
  }

  foGlobal.close(); 

  yaatk::chdir("..");

  gsl_rng_free (r);

  return sl;
}

inline
SimLoop
build_embed(const mdtk::SimLoop& sl_cluster, 
            const mdtk::SimLoop& sl_shell)
{
  mdtk::SimLoop sl;
  initialize_simloop(sl);

  shiftToOrigin(sl_cluster.atoms);
  shiftToOrigin(sl_shell.atoms);
  
  sl = sl_cluster;

//shrinking
  for(size_t i = 0; i < sl.atoms.size(); i++)
  {
    Atom& a = *(sl.atoms[i]);
    a.coords *= 0.8;
  }
//~shrinking

  sl.add_simloop(sl_shell);

  for(size_t i = 0; i < sl.atoms.size(); i++)
  {
    Atom& a = *(sl.atoms[i]);
    a.V = 0.0;
  }

  quench(sl,1.0*K);

  removeMomentum(sl.atoms);

  return sl;
}

inline
void
add_rotational_motion(mdtk::SimLoop& sl,
                      Float totalRotEnergy = 100.0*eV,
                      Vector3D rotAxis = Vector3D(0.0,1.0,0.0))
{
  Vector3D clusterCenterOfM = massCenter(sl.atoms);

  Float sumMVSQR = 0.0;
  for(size_t i = 0; i < sl.atoms.size(); i++)
  {
    Atom& a = *(sl.atoms[i]);
    Vector3D r = a.coords - clusterCenterOfM;
    sumMVSQR += a.M*SQR(vectormul(r,rotAxis).module());
  }

  Float maxVrot = sqrt(2.0*totalRotEnergy/sumMVSQR);
  for(size_t i = 0; i < sl.atoms.size(); i++)
  {
    Atom& a = *(sl.atoms[i]);
    Vector3D r = a.coords - clusterCenterOfM;
    a.V = vectormul(r,rotAxis)*maxVrot;
  }

  TRACE(sl.energyKin()/eV);
  sl.thermalBath.zMin = +100000.0*Ao;
  relax/*_flush*/(sl,0.2*ps);
  TRACE(sl.energyKin()/eV);
  removeMomentum(sl.atoms);
  TRACE(sl.energyKin()/eV);
}

inline
SimLoop
build_target_by_cluster_bombardment(
  const mdtk::SimLoop& sl_target, 
  const mdtk::SimLoop sl_cluster,
  Float clusterEnergy = 100*eV)
{
  shiftToOrigin(sl_cluster.atoms);

  Float clusterRadius 
    = maxDistanceFrom(sl_cluster.atoms,geomCenter(sl_cluster.atoms));
  TRACE(clusterRadius/mdtk::Ao);

  Vector3D dCluster = sl_target.getPBC()/2.0;
  dCluster.z = -(10.0*Ao+clusterRadius);
  TRACE(dCluster/mdtk::Ao);

  for(size_t i = 0; i < sl_cluster.atoms.size(); i++)
  {
    Atom& a = *(sl_cluster.atoms_[i]);
    a.coords += dCluster;
  }

  for(int it = 0;;it = 1)
  {
    TRACE("Landing fullerene ...");
    bool areInteracting = false;
    for(size_t cli = 0; !areInteracting && cli < sl_cluster.atoms.size(); cli++)
      for(size_t surfi = 0; !areInteracting && surfi < sl_target.atoms.size(); surfi++)
      {
        Atom& clusterAtom = *(sl_cluster.atoms[cli]);
        Atom& surfaceAtom = *(sl_target.atoms[surfi]);
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
    for(size_t i = 0; i < sl_cluster.atoms.size(); i++)
    {
      Atom& a = *(sl_cluster.atoms[i]);
      a.coords.z += 0.05*Ao;
    }
    TRACE(sl_cluster.atoms[0]->coords.z/mdtk::Ao);
  };
  TRACE("Fullerene landed ok.");

  for(size_t i = 0; i < sl_cluster.atoms.size(); i++)
  {
    Atom& a = *(sl_cluster.atoms[i]);
    a.V += Vector3D(0,0,sqrt(2.0*clusterEnergy/(mass(sl_cluster.atoms))));
  }

  mdtk::SimLoop sl;
  initialize_simloop(sl);

  sl = sl_target;
  sl.add_simloop(sl_cluster);

//  removeMomentum(sl.atoms);

  return sl;
}

inline
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

  mdtk::SimLoop sl_Cu = mdbuilder::build_FCC_lattice(20,20,10,Cu_EL);
  mdtk::SimLoop sl_C60 = mdbuilder::build_C60_optimized();

  for(size_t i_cluster_size = 0; i_cluster_size < cluster_sizes.size(); i_cluster_size++)
  {
    int& cluster_size = cluster_sizes[i_cluster_size];
    mdtk::SimLoop sl_cluster = mdbuilder::build_cluster(Cu_EL,cluster_size);

    for(size_t i_rot_axis = 0; i_rot_axis < rot_axes.size(); i_rot_axis++)
    {
      for(size_t i_rot_energy = 0; i_rot_energy < rot_energies.size(); i_rot_energy++)
      {
        mdtk::SimLoop sl_endo = mdbuilder::build_embed(sl_cluster,sl_C60);
        Vector3D& rot_axis = rot_axes[i_rot_axis];
        Float& rot_energy = rot_energies[i_rot_energy];
        mdbuilder::add_rotational_motion(sl_endo,rot_energy,rot_axis);

        for(size_t i_trans_energy = 0; i_trans_energy < trans_energies.size(); i_trans_energy++)
        {
          Float& trans_energy = trans_energies[i_trans_energy];
          mdtk::SimLoop sl = mdbuilder::build_target_by_cluster_bombardment(sl_Cu,sl_endo,trans_energy);
          
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

inline
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

  mdtk::SimLoop sl_Graphite = mdbuilder::build_Graphite_lattice(12,14,3);
  mdtk::SimLoop sl_C60 = mdbuilder::build_C60_optimized();

  for(size_t i_cluster_size = 0; i_cluster_size < cluster_sizes.size(); i_cluster_size++)
  {
    int& cluster_size = cluster_sizes[i_cluster_size];
    mdtk::SimLoop sl_cluster = mdbuilder::build_cluster(Cu_EL,cluster_size);

    for(size_t i_rot_axis = 0; i_rot_axis < rot_axes.size(); i_rot_axis++)
    {
      for(size_t i_rot_energy = 0; i_rot_energy < rot_energies.size(); i_rot_energy++)
      {
        mdtk::SimLoop sl_endo = mdbuilder::build_embed(sl_cluster,sl_C60);
        Vector3D& rot_axis = rot_axes[i_rot_axis];
        Float& rot_energy = rot_energies[i_rot_energy];
        mdbuilder::add_rotational_motion(sl_endo,rot_energy,rot_axis);

        for(size_t i_trans_energy = 0; i_trans_energy < trans_energies.size(); i_trans_energy++)
        {
          Float& trans_energy = trans_energies[i_trans_energy];
          mdtk::SimLoop sl = mdbuilder::build_target_by_cluster_bombardment(sl_Graphite,sl_endo,trans_energy);
          
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

inline
mdtk::SimLoop
build_Cluster_Landed_on_Substrate(
  const mdtk::SimLoop sl_Substrate,
  ElementID id,
  int clusterSize
  )
{
  mdtk::SimLoop sl_Cluster = mdbuilder::build_cluster(id,clusterSize);
  removeMomentum(sl_Cluster.atoms);

  mdtk::SimLoop sl = mdbuilder::build_target_by_cluster_bombardment(sl_Substrate,sl_Cluster,0.2*eV*clusterSize);

  TRACE(sl.energyKin()/eV);

  {
    yaatk::text_ofstream fomde("_tmp-X-relax_flush-landing.mde");
    sl.saveToMDE(fomde);
    fomde.close();
  }


  std::vector<size_t> fixedAtoms = fixUnfixedCHAtoms(sl.atoms,0,sl.atoms.size());
  {
    Float tb_zMin_bak = sl.thermalBath.zMin;
    sl.thermalBath.zMin = -1.0*Ao;

    relax/*_flush*/(sl,5.0*ps,"_tmp-X-landing-fixed-CH-relax_flush");

    sl.thermalBath.zMin = tb_zMin_bak;
  }
  unfixAtoms(sl.atoms,fixedAtoms);

  relax/*_flush*/(sl,5.0*ps,"_tmp-X-landing-unfixed-CH-relax_flush");

  quench(sl,0.01*K);

  return sl;
}

}

#endif
