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
void
quench(mdtk::SimLoop& sl, 
       Float forTime = 0.2*ps,
       std::string tmpDir = "_tmp-X")
{
  Float freezingEnergy = 0.0001*2.5*mdtk::eV*sl.atoms_.size();
  TRACE(freezingEnergy/mdtk::eV);
  ERRTRACE(freezingEnergy/mdtk::eV);

  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  setupPotentials(sl);
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = 0.0;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.thermalBath.zMin = -100000.0*Ao;
  while (sl.simTimeFinal < forTime)
  {
    sl.simTimeFinal += 0.05*ps;
    sl.execute();
    if (sl.energyKin() < freezingEnergy) break;
  }
  yaatk::chdir("..");
}

inline
void
relax(mdtk::SimLoop& sl, 
      Float forTime = 0.2*ps,
      std::string tmpDir = "_tmp-X")
{
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  setupPotentials(sl);
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.05*ps;
  sl.execute();
  yaatk::chdir("..");
}

inline
void
relax_flush(mdtk::SimLoop& sl, 
      Float forTime = 0.2*ps,
      std::string tmpDir = "_tmp-X")
{
  yaatk::mkdir(tmpDir.c_str());
  yaatk::chdir(tmpDir.c_str());
  setupPotentials(sl);
  sl.initialize();
  sl.simTime = 0.0*ps;
  sl.simTimeFinal = forTime;
  sl.simTimeSaveTrajInterval = 0.001*ps;
  sl.execute();
  yaatk::chdir("..");
}

inline
void
build_C60_optimized(mdtk::SimLoop& sl)
{
  glLoadIdentity();
  mdbuilder::place_C60(sl);
  
  quench(sl,0.2*ps,"_tmp-C60-optimized");
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
optimize_single(SimLoop *modloop)
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
    modloop->heatUpEveryAtom(maxHeatUpEnergy);
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
add1atom(mdtk::AtomsContainer& atoms, ElementID id)
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

  xi = rand()%3;
  xs = rand()%2;

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
  Vector3D vn(rand()/Float(RAND_MAX),
              rand()/Float(RAND_MAX),
              rand()/Float(RAND_MAX));    
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
void
build_cluster(mdtk::SimLoop& sl, ElementID id, int clusterSize)
{
  if (clusterSize == 0) return;

  yaatk::mkdir("_tmp-Cu-clusters");
  yaatk::chdir("_tmp-Cu-clusters");

  std::ofstream foGlobal("energy.min.all",std::ios::app); 

  setupPotentials(sl);

  if (sl.atoms.size() > 0)
    add1atom(sl.atoms_,id);
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

    Float minPotEnergyOf = optimize_single(&sl);

    foGlobal << std::setw (10) << sl.atoms_.size() << " "
             << std::setw (20) << minPotEnergyOf/mdtk::eV << " "
             << std::setw (20) << minPotEnergyOf/mdtk::eV/sl.atoms_.size() << std::endl;

    yaatk::chdir("..");

    if (atomsCount < clusterSize)
      add1atom(sl.atoms_, id);
  }

  foGlobal.close(); 

  yaatk::chdir("..");
}

Float
mass(const AtomsContainer& atoms)
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    moleculeMass += atom.M;
  }
  return moleculeMass;
}

Vector3D
velocity(const AtomsContainer& atoms)
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.V*atom.M;
  };
  return sumOfP/sumOfM;
}

Vector3D
massCenter(const AtomsContainer& atoms)
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = *atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.coords*atom.M;
  };
  return sumOfP/sumOfM;
}

void
removeMomentum(AtomsContainer& atoms)
{
  Vector3D v = velocity(atoms);
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    mdtk::Atom& atom = *atoms[ai];
    atom.V -= v;
  };
}

inline
Vector3D
geomCenter(AtomsContainer& atoms)
{
  Float clusterRadius = 0.0;
  Vector3D clusterCenter(0,0,0);
  
  for(size_t i = 0; i < atoms.size(); i++)
    clusterCenter += atoms[i]->coords;
  clusterCenter /= atoms.size();

  return clusterCenter;
}

inline
Float
maxDistanceFrom(AtomsContainer& atoms, Vector3D point)
{
  Float clusterRadius = 0.0;
  
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Float currentDist = (atoms[i]->coords-point).module();
    clusterRadius = (currentDist>clusterRadius)?currentDist:clusterRadius;
  }

  return clusterRadius;
}

inline
Float
shiftToOrigin(AtomsContainer& atoms)
{
  Float clusterRadius = 0.0;
  Vector3D clusterCenter(0,0,0);
  
  for(size_t i = 0; i < atoms.size(); i++)
    clusterCenter += atoms[i]->coords;
  clusterCenter /= atoms.size();
  for(size_t i = 0; i < atoms.size(); i++)
    atoms[i]->coords -= clusterCenter;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Float currentDist = atoms[i]->coords.module();
    clusterRadius = (currentDist>clusterRadius)?currentDist:clusterRadius;
  }

  return clusterRadius;
}

inline
void
build_embed(mdtk::SimLoop& sl_cluster, 
            mdtk::SimLoop& sl_shell,
            mdtk::SimLoop& sl)
{
  shiftToOrigin(sl_cluster.atoms);
  shiftToOrigin(sl_shell.atoms);
  
//shrinking
  for(size_t i = 0; i < sl_cluster.atoms.size(); i++)
  {
    Atom& a = *(sl_cluster.atoms[i]);
    a.coords *= 0.8;
  }
//~shrinking

  setupPotentials(sl);
  sl.thermalBath.zMin = -100000.0*Ao;

  for(size_t i = 0; i < sl_shell.atoms.size(); i++)
  {
    Atom& a = *(sl_shell.atoms[i]);
    sl.atoms.push_back(&a);
  }

  for(size_t i = 0; i < sl_cluster.atoms.size(); i++)
  {
    Atom& a = *(sl_cluster.atoms[i]);
    sl.atoms.push_back(&a);
  }

  for(size_t i = 0; i < sl.atoms.size(); i++)
  {
    Atom& a = *(sl.atoms[i]);
    a.V = 0.0;
  }

  quench(sl,0.5*ps,"_tmp-embed");

  removeMomentum(sl.atoms);
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
  relax_flush(sl,0.2*ps);
  TRACE(sl.energyKin()/eV);
  removeMomentum(sl.atoms);
  TRACE(sl.energyKin()/eV);
}

inline
void
copy_simloop(
  mdtk::SimLoop& sl_dest, 
  mdtk::SimLoop& sl_src)
{
  setupPotentials(sl_dest);
  sl_dest.setPBC(sl_src.getPBC());
  sl_dest.thermalBath = sl_src.thermalBath;
  for(size_t i = 0; i < sl_src.atoms_.size(); i++)
  {
    Atom& a = *(sl_src.atoms_[i]);
    sl_dest.atoms_.push_back(&a);
  }
}

inline
void
add_simloop(
  mdtk::SimLoop& sl, 
  mdtk::SimLoop& sl_addon)
{
  setupPotentials(sl);
  for(size_t i = 0; i < sl_addon.atoms_.size(); i++)
  {
    Atom& a = *(sl_addon.atoms_[i]);
    sl.atoms_.push_back(&a);
  }
}

inline
void
build_target_by_cluster_bombardment(
  mdtk::SimLoop& sl_target, 
  mdtk::SimLoop& sl_cluster,
  mdtk::SimLoop& sl,
  Float clusterEnergy = 100*eV)
{
  copy_simloop(sl,sl_target);
  add_simloop(sl,sl_cluster);

  shiftToOrigin(sl_cluster.atoms);

  Float clusterRadius 
    = maxDistanceFrom(sl_cluster.atoms,geomCenter(sl_cluster.atoms));
  TRACE(clusterRadius/mdtk::Ao);

  Vector3D dCluster = sl.getPBC()/2.0;
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
        if (depos(clusterAtom,surfaceAtom).module() < ((clusterEnergy==0.0)?3.0*Ao:6.0*Ao))
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
    a.V += Vector3D(0,0,sqrt(2.0*clusterEnergy/(a.M)/sl_cluster.atoms.size()));
  }

  removeMomentum(sl.atoms);
}

}

#endif
