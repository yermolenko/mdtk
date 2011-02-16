/*
   The molecular dynamics simulation loop class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011
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

#include "Exception.hpp"
#include "SimLoop.hpp"
#include <fstream>
#include "release_info.hpp"

#include <cmath>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>

namespace mdtk
{

using namespace std;

SimLoop::SimLoop()
  : allowToFreePotentials(true),
    allowToFreeAtoms(true),
    atoms_(),
    atoms(atoms_),
    check(),
    simTime(0.0), 
    breakSimLoop(false),
    iteration(0) ,
    barrier(100,0.0),
    thermalBath(),
    initNLafterLoading(true),
    allowPartialLoading(false),
    fpot()
    ,CPUTimeUsed_prev(0)
    ,CPUTimeUsed_total(0)
{
  verboseTrace = true;
  
  timeaccel_ = 1.0*Ao;

  dt_ = 1e-20; // initial dt_, it changes adaptive during simulation
//  dt_ = 1e-17;
  
  iterationFlushStateInterval = 10;
  simTimeFinal = 4.0*ps;

  simTimeSaveTrajInterval = 0.1*ps;

  check.checkEnergy = true;
  check.checkForce = true;
  check.checkEnergyAfter = 1; // dummy, will be removed soon
}

SimLoop::SimLoop(const SimLoop &c)
  : allowToFreePotentials(true),
    allowToFreeAtoms(true),
    atoms_(),
    atoms(atoms_),
    check(),
    simTime(0.0), 
    breakSimLoop(false),
    iteration(0) ,
    barrier(100,0.0),
    thermalBath(),
    initNLafterLoading(true),
    allowPartialLoading(false),
    fpot()
    ,CPUTimeUsed_prev(0)
    ,CPUTimeUsed_total(0)
{
  setPBC(c.getPBC());
  thermalBath = c.thermalBath;
  for(size_t i = 0; i < c.atoms_.size(); i++)
  {
    Atom& a = *(c.atoms_[i]);
    atoms_.push_back(&(*(new Atom()) = a));
  }
}

SimLoop&
SimLoop::operator =(const SimLoop &c) 
{
  if (this == &c) return *this;

  atoms_.clear();

  setPBC(c.getPBC());
  thermalBath = c.thermalBath;
  for(size_t i = 0; i < c.atoms_.size(); i++)
  {
    Atom& a = *(c.atoms_[i]);
    atoms_.push_back(&(*(new Atom()) = a));
  }

  return *this;
}

void
SimLoop::add_simloop(const SimLoop &sl_addon)
{
  for(size_t i = 0; i < sl_addon.atoms_.size(); i++)
  {
    Atom& a = *(sl_addon.atoms_[i]);
    atoms_.push_back(&(*(new Atom()) = a));
  }
}

SimLoop::~SimLoop()
{
  freePotentials();
  freeAtoms();
}

void
SimLoop::checkOnSpot(Atom& atom) // obsolete
{
  if (!atom.apply_barrier) return;
  if (
      ! (atom.ejected) &&
      atom.coords.z > barrier.z &&
      atom.V.z > 0.0)
  {
    Float Vn;
    Float En;
    Vn = atom.V.z;
    En = atom.M*SQR(Vn)/2.0;
    En += barrier.dE;
    cout << "***Attempt to eject ";
    if (En>0)
    {
      atom.V.z = sqrt(2.0*En/atom.M);
      cout << "is successfull ***" << endl;
      atom.ejected = true;
    }
    else
    {
      atom.V.z = -(atom.V.z);
      atom.coords.z = barrier.z;
      cout << endl;
    }
  }
}  

int
SimLoop::execute()
{
  try
  {
    fpot.diagnose();
    PTRACE(usePBC());
    PTRACE(getPBC()/Ao);
    PTRACE(thermalBath.zMin/Ao);
    PTRACE(thermalBath.dBoundary/Ao);
    PTRACE(thermalBath.zMinOfFreeZone/Ao);
    PTRACE(thermalBath.To/K);
    return execute_wo_checks();
  }
  catch (Exception& e)
  {
    std::cerr << "Caught mdtk Exception: " << e.what() << std::endl;
    std::cerr << "Flushing state.....";
    writestate();
    { 
      std::ofstream fo("in.mde.after_crash"); 
      saveToMDE(fo); 
      fo.close(); 
    }
    std::cerr << "done." << std::endl;
    return -1;
  }
  catch (MPI_Exception& e)
  {
    std::cerr << "Caught MPI Exception: " << e.what() << std::endl;
    std::cerr << "Flushing state.....";
    writestate();
    std::cerr << "done." << std::endl;
    return -1;
  }
}  

bool
SimLoop::checkMIC()
{
  Vector3D PBC;
  PBC = getPBC();
  if (PBC.x <= fpot.getRcutoff()*2.0 || PBC.y <= fpot.getRcutoff()*2.0)
  {
    return false;
  }
  else
  {
    return true;
  }  
}  

bool
SimLoop::fitInCell()
{
  Vector3D PBC;
  PBC = getPBC();
  size_t i;
  for(i = 0; i < atoms_.size(); i++)
  {
    if (!(atoms_[i]->apply_PBC)) continue;
    Vector3D aci = atoms_[i]->coords;
    if (PBC.x < MDTK_MAX_PBC) if (aci.x < 0 || aci.x>=PBC.x) {TRACE(i);TRACE(aci.x);TRACE(PBC.x);return false;}
    if (PBC.y < MDTK_MAX_PBC) if (aci.y < 0 || aci.y>=PBC.y) {TRACE(i);TRACE(aci.y);TRACE(PBC.y);return false;}
    if (PBC.z < MDTK_MAX_PBC) if (aci.z < 0 || aci.z>=PBC.z) {TRACE(i);TRACE(aci.z);TRACE(PBC.z);return false;}
  }
  return true;
}  

void
SimLoop::applyPBC(Atom& a)
{
  Vector3D& ac = a.coords;
  if (!usePBC()) return;
  if (!(a.apply_PBC)) return;
  Vector3D PBC(getPBC());

  if (PBC.x < MDTK_MAX_PBC)
  {
    while(ac.x < 0)      {ac.x += PBC.x;--a.PBC_count.x;}
    while(ac.x >= PBC.x) {ac.x -= PBC.x;++a.PBC_count.x;}
  }  
  if (PBC.y < MDTK_MAX_PBC)
  {
    while(ac.y < 0)      {ac.y += PBC.y;--a.PBC_count.y;}
    while(ac.y >= PBC.y) {ac.y -= PBC.y;++a.PBC_count.y;}
  }
  if (PBC.z < MDTK_MAX_PBC)
  {
    while(ac.z < 0)      {ac.z += PBC.z;--a.PBC_count.z;}
    while(ac.z >= PBC.z) {ac.z -= PBC.z;++a.PBC_count.z;}
  }  
}  

void
SimLoop::initPBC()
{
  size_t j, atoms_count = atoms_.size();
  for(j = 0; j < atoms_count; j++)
  {
     applyPBC(*(atoms_[j]));
  }  
}  

int
SimLoop::execute_wo_checks()
{
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_ranlxd2);
  gsl_rng_set(r, 697860L);

  breakSimLoop = false;
  time_t startWallTime = time(NULL);
  time_t curWallTime = startWallTime;

  procmon::ProcmonTimer pmtTotal;

  if (usePBC() && !checkMIC()) 
  {
    cerr << "Rcutoff is too large for given PBC !" << endl << flush;
    throw Exception("Rcutoff is too large for given PBC !");
  }  

  if (usePBC()) initPBC();

  if (usePBC() && !fitInCell()) 
  {
    cerr << "Atoms do not fit given PBC cell !" << endl << flush;
    throw Exception("Atoms do not fit given PBC cell !");
  }  

  int j,atoms_count;
  Float new_v_max,tmp_v;
  
  atoms_count = atoms_.size();
  Vector3D *new_coords = new Vector3D[atoms_count]; REQUIRE(new_coords != 0);

  while (simTime < simTimeFinal && !breakSimLoop)
  {
     doBeforeIteration();
      if (verboseTrace)
      {
        cout << "---------" << endl;
        cout << "t : " << simTime << endl;
        cout << "dt : " << dt_ << endl;
      }

    if (iteration%iterationFlushStateInterval == 0/* && iteration != 0*/)
    {
      cout << "Writing state ... " ;
      writestate();
      cout << "done. " << endl;
    };

    if (simTime == 0.0 || int(simTime/simTimeSaveTrajInterval) != int((simTime - dt_)/simTimeSaveTrajInterval))
    {
      cout << "Writing trajectory ... " ;
      writetrajXVA();
      cout << "done. " << endl;
    };

    fpot.NL_UpdateIfNeeded(atoms_);

      if (iteration == 0)
        init_check_energy();
      else
        do_check_energy();  

    Float actualThermalBathTemp = actualTemperatureOfThermalBath();
    {
      Float Tb = actualThermalBathTemp;
      PTRACE(Tb);
    }

    if (check.checkForce)
    {
      check.fullForce = 0;
    }    

    new_v_max = 0.0;

    for(j = 0; j < atoms_count; j++)
    {
      Atom& atom = *(atoms_[j]); 
      Vector3D force;

      force = (-fpot.grad(atom,this->atoms_));

      if (check.checkForce)
      {
        check.fullForce += force;
      }  

      if (isWithinThermalBath(atom.coords) && atom.apply_ThermalBath)
      {
//        Float T = check.temperatureCur;
        Float T = actualThermalBathTemp;
        Float To_by_T = (fabs(T)<1e-5)?0:(thermalBath.To/T);
        Float max_To_by_T = 5.0;
        if (To_by_T < -max_To_by_T) To_by_T = -max_To_by_T;
        if (To_by_T > +max_To_by_T) To_by_T = +max_To_by_T;

        Float gamma = 1.0e13;
//        Float gamma = 1.0/(1000.0*dt_);

        Vector3D dforce = -atom.V*atom.M*gamma*(1.0-sqrt(To_by_T));

        {
          Vector3D dv_no_tb,dv;

          Vector3D an_new_no_tb(force         /atom.M);
          if (iteration == 0) atom.an_no_tb = an_new_no_tb;
          dv_no_tb = (an_new_no_tb+atom.an_no_tb)*dt_/2.0;
          atom.an_no_tb = an_new_no_tb;

          Vector3D an_new(      (force+dforce)/atom.M);
          dv       = (an_new      +      atom.an)*dt_/2.0;
          
          Float dEkin = atom.M*SQR((atom.V+dv)      .module())/2.0                      - atom.M*SQR((atom.V+dv_no_tb).module())/2.0;
          check.energyStart += dEkin;        }          force += dforce;
      }  

      {
        Vector3D an_new(force/atom.M);
        if (iteration == 0) atom.an = an_new;
        
        atom.V += (an_new+atom.an)*dt_/2.0;
        new_coords[j] = atom.coords + atom.V*dt_+an_new*dt_*dt_/2.0;
        atom.an = an_new;
      }
        
      tmp_v = atom.V.module();
      if (tmp_v > new_v_max  && fpot.hasNB(atom) && atom.coords.module() < 500.0*Ao)  new_v_max = tmp_v;
    };

    for(j = 0; j < atoms_count; j++)
    {
      fpot.incDisplacement(*(atoms_[j]),new_coords[j]-atoms_[j]->coords);
    }  

    for(j = 0; j < atoms_count; j++)
    {
      atoms_[j]->coords = new_coords[j];
 
      applyPBC(*(atoms_[j]));
      checkOnSpot(*(atoms_[j])); // obsolete
    }

    doAfterIteration();
    
    fpot.NL_checkRequestUpdate(atoms_);

    simTime += dt_;

    if (check.checkForce)
    {
      if (verboseTrace)
      {
        Float ffmod = check.fullForce.module();
        PTRACE(ffmod);
      }

      if (check.fullForce.module() > 1e-8) 
      {
        cerr << "FullForce != 0" << endl << flush;
        cout << "FullForce != 0" << endl << flush;
      }  
    }    

    const Float dt_max = 5e-16;
    const Float dt_min = 1e-20;
    if (new_v_max != 0.0)
      dt_ = 0.05*timeaccel_/new_v_max;
    else
      dt_ = dt_max;    
    if (dt_ > dt_max) dt_ = dt_max;
    if (dt_ < dt_min) dt_ = dt_min;


    if (verboseTrace)
    {
      cout << "it : " << iteration << endl;
    }
    curWallTime = time(NULL);
    if (verboseTrace) cout << "tw : " << (curWallTime-startWallTime) << endl;

    CPUTimeUsed_total = CPUTimeUsed_prev + pmtTotal.getTimeInSeconds();

    if (verboseTrace) 
      cout << "tc : " << CPUTimeUsed_total << endl;

    iteration = (iteration+1)%2000000000L;
    
  };
  cout << "Final Modeling time = " << simTime << endl;
  cout << "Final Time step = " << dt_ << endl;
  cout << "Modeling cycle complete " << endl;

  cout << "--------------------------------------------------------- " << endl;
    {
      cout << "Writing trajectory ... " ;
      writetrajXVA();
      cout << "done. " << endl;
    };
  curWallTime = time(NULL);
  cout << "Wall TIME used = " << (curWallTime-startWallTime) << endl;
  
  delete [] new_coords;

  gsl_rng_free (r);

  return 0;
}

Float
SimLoop::energy()
{
  Float energyPotCur = energyPot();
  Float energyKinCur = energyKin();

{
Float Ep = energyPotCur/mdtk::eV;
Float Ek = energyKinCur/mdtk::eV;
Float Et = (energyPotCur+energyKinCur)/mdtk::eV;

PTRACE(Ep);
PTRACE(Ek);
PTRACE(Et);
}

  return energyPotCur+energyKinCur;
}

Float
SimLoop::energyPot()
{
  return fpot(atoms_);
}

Float
SimLoop::energyKin()
{
  Float energyKinCur = 0.0;
  int j,atoms_count;
  atoms_count = atoms_.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = *atoms_[j];
    energyKinCur += atom.M*SQR(atom.V.module())/2.0;
  };

  return energyKinCur;
}

Float
SimLoop::temperature()
{
//  Float T = energyKin()/(3.0/2.0*kb*atoms_.size());
  Float T = temperatureWithoutFixed();
  PTRACE(T);
  return T;
}

Float
SimLoop::actualTemperatureOfThermalBath()
{
  Float energyKinCur = 0.0;
  int j,atoms_count,atoms_accounted = 0;
  atoms_count = atoms_.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = *atoms_[j];
    if (isWithinThermalBath(atom.coords) && atom.apply_ThermalBath)
    {
      energyKinCur += atom.M*SQR(atom.V.module())/2.0;
      atoms_accounted++;
    }
  };

  return energyKinCur/(3.0/2.0*kb*atoms_accounted);
}

Float
SimLoop::temperatureWithoutFixed()
{
  Float energyKinCur = 0.0;
  int j,atoms_count,atoms_accounted = 0;
  atoms_count = atoms_.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = *atoms_[j];
    if (!atom.isFixed())
    {
      energyKinCur += atom.M*SQR(atom.V.module())/2.0;
      atoms_accounted++;
    }
  };

  return energyKinCur/(3.0/2.0*kb*atoms_accounted);
}

void
SimLoop::loadFromStream(istream& is, YAATK_FSTREAM_MODE smode)
{
  freeAtoms();

  if (smode == YAATK_FSTREAM_TEXT) 
  {
    char s[1024];
    is >> s;
    if (s != id) 
    {
      throw Exception("State file is not compatible");
    }
  }  
  int i,atoms_count;
  YAATK_FSTREAM_READ(is,atoms_count,smode);
  cout << "Reading info about " << atoms_count << " atoms..." << endl;
  atoms_.resize(atoms_count);
  for(i = 0; i < atoms_count; i++)
  {
    atoms_[i] = new Atom; REQUIRE(atoms_[i] != 0);
    YAATK_FSTREAM_READ(is,*(atoms_[i]),smode);
  }
  cout << endl;

  cout << "Reading timing info ... " << endl;

  check.LoadFromStream(is,smode);

  YAATK_FSTREAM_READ(is,simTime,smode);
  YAATK_FSTREAM_READ(is,simTimeSaveTrajInterval,smode);
  YAATK_FSTREAM_READ(is,simTimeFinal,smode);
  YAATK_FSTREAM_READ(is,timeaccel_,smode);
  YAATK_FSTREAM_READ(is,dt_,smode);
  YAATK_FSTREAM_READ(is,iteration,smode); //iteration++;
  YAATK_FSTREAM_READ(is,iterationFlushStateInterval,smode);

  barrier.LoadFromStream(is,smode);
  
  thermalBath.LoadFromStream(is,smode);

  fpot.LoadFromStream(is,smode);

  YAATK_FSTREAM_READ(is,CPUTimeUsed_prev,smode);
  CPUTimeUsed_total = CPUTimeUsed_prev;

  atoms_.LoadFromStream(is,smode);

//TRACE(atoms_.getPBC());

  cout << "Parsing of state file done. " << endl;

  initialize();
  
  allowPartialLoading = true;
}

void
SimLoop::saveToStream(ostream& os, YAATK_FSTREAM_MODE smode)
{
  std::streamsize prevStreamSize = os.precision();
  os.precision(FLOAT_PRECISION);

  if (smode == YAATK_FSTREAM_TEXT) os << id << "\n"; 
  int i,atoms_count = atoms_.size();
  YAATK_FSTREAM_WRITE(os,atoms_count,smode);
  for(i = 0; i < atoms_count; i++)
  {
    YAATK_FSTREAM_WRITE(os,*(atoms_[i]),smode);
  }

  check.SaveToStream(os,smode);

  YAATK_FSTREAM_WRITE(os,simTime,smode);
  YAATK_FSTREAM_WRITE(os,simTimeSaveTrajInterval,smode);
  YAATK_FSTREAM_WRITE(os,simTimeFinal,smode);
  YAATK_FSTREAM_WRITE(os,timeaccel_,smode);
  YAATK_FSTREAM_WRITE(os,dt_,smode);
  YAATK_FSTREAM_WRITE(os,iteration,smode);
  YAATK_FSTREAM_WRITE(os,iterationFlushStateInterval,smode);

  barrier.SaveToStream(os,smode);
  
  thermalBath.SaveToStream(os,smode);

  fpot.SaveToStream(os,smode);
  
  YAATK_FSTREAM_WRITE(os,CPUTimeUsed_total,smode);

  os.precision(prevStreamSize);

  atoms_.SaveToStream(os,smode);
}

//#define DONT_USE_XVASCALE

void
SimLoop::loadFromStreamXVA(istream& is)
{
  REQUIRE(allowPartialLoading == true);


Float XVA_VELOCITY_SCALE = 1.0;
Float XVA_DISTANCE_SCALE = 1.0;

#ifndef DONT_USE_XVASCALE

is >> XVA_VELOCITY_SCALE;
is >> XVA_DISTANCE_SCALE;

#endif

TRACE(XVA_VELOCITY_SCALE);
TRACE(XVA_DISTANCE_SCALE);


  size_t i,atoms_count;
  is >> atoms_count;
  REQUIRE(atoms_count == atoms_.size());


  cout << "Reading XVA info about " << atoms_count << " atoms..." << endl;

  for(i = 0; i < atoms_count; i++)
  {
    is >> atoms_[i]->V;
    is >> atoms_[i]->coords;
    is >> atoms_[i]->PBC_count;

    atoms_[i]->V *= XVA_VELOCITY_SCALE;
    atoms_[i]->coords *= XVA_DISTANCE_SCALE;
  }
  cout << endl;

  check.LoadFromStream(is,YAATK_FSTREAM_TEXT);

  is >> simTime;
  is >> dt_;
  is >> iteration; //iteration++;

  is >> CPUTimeUsed_prev;
  CPUTimeUsed_total = CPUTimeUsed_prev;

  cout << "Parsing of state file done. " << endl;

  initialize();
}

void
SimLoop::saveToStreamXVA(ostream& os)
{

Float XVA_VELOCITY_SCALE = 0.0;//1.0e3;
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    XVA_VELOCITY_SCALE += fabs(atoms_[i]->V.x);
    XVA_VELOCITY_SCALE += fabs(atoms_[i]->V.y);
    XVA_VELOCITY_SCALE += fabs(atoms_[i]->V.z);
  }

XVA_VELOCITY_SCALE /= (3.0*atoms_.size());
if (fabs(XVA_VELOCITY_SCALE) < 1e-30) XVA_VELOCITY_SCALE = 1e-30;
XVA_VELOCITY_SCALE = pow(10.0,ceil(log10(fabs(XVA_VELOCITY_SCALE)))-2);


Float XVA_DISTANCE_SCALE = Ao;

#ifdef DONT_USE_XVASCALE
XVA_VELOCITY_SCALE = 1.0;
XVA_DISTANCE_SCALE = 1.0;
#endif

#ifndef DONT_USE_XVASCALE

os << XVA_VELOCITY_SCALE << "\n";
os << XVA_DISTANCE_SCALE << "\n";

#endif

  std::streamsize prevStreamSize = os.precision();

  int i,atoms_count = atoms_.size();
  os << atoms_count << "\n";
  for(i = 0; i < atoms_count; i++)
  {
#ifndef DONT_USE_XVASCALE
    os << fixed;
#endif
    os.precision(2);
    os << atoms_[i]->V.x/XVA_VELOCITY_SCALE  << " ";
    os << atoms_[i]->V.y/XVA_VELOCITY_SCALE  << " ";
    os << atoms_[i]->V.z/XVA_VELOCITY_SCALE  << "\n";

#ifndef DONT_USE_XVASCALE
    os << fixed;
#endif
    os.precision(2);
    os << atoms_[i]->coords.x/XVA_DISTANCE_SCALE  << " ";
    os << atoms_[i]->coords.y/XVA_DISTANCE_SCALE  << " ";
    os << atoms_[i]->coords.z/XVA_DISTANCE_SCALE  << "\n";

    os << atoms_[i]->PBC_count.x << " ";
    os << atoms_[i]->PBC_count.y << " ";
    os << atoms_[i]->PBC_count.z << "\n";
  }
  os << scientific;

  os.precision(FLOAT_PRECISION);

  check.SaveToStream(os,YAATK_FSTREAM_TEXT);

  os << simTime << "\n";
  os << dt_ << "\n";
  os << iteration << "\n";
  
  os << CPUTimeUsed_total << "\n";

  os.precision(prevStreamSize);
}

void
SimLoop::loadFromStreamXVA_bin(istream& is)
{
  REQUIRE(allowPartialLoading == true);

  size_t i,atoms_count;
  YAATK_BIN_READ(is,atoms_count);
  REQUIRE(atoms_count == atoms_.size());


  cout << "Reading XVA info about " << atoms_count << " atoms..." << endl;

  for(i = 0; i < atoms_count; i++)
  {
    YAATK_BIN_READ(is,atoms_[i]->V);
    YAATK_BIN_READ(is,atoms_[i]->coords);
  }
  cout << endl;

  YAATK_BIN_READ(is,check);

  YAATK_BIN_READ(is,simTime);
  YAATK_BIN_READ(is,dt_);
  YAATK_BIN_READ(is,iteration); //iteration++;

  TRACE(iteration);

  YAATK_BIN_READ(is,CPUTimeUsed_prev);
  CPUTimeUsed_total = CPUTimeUsed_prev;

  cout << "Parsing of state file done. " << endl;

  initialize();
}

void
SimLoop::saveToStreamXVA_bin(ostream& os)
{
  int i,atoms_count = atoms_.size();
  YAATK_BIN_WRITE(os,atoms_count);
  for(i = 0; i < atoms_count; i++)
  {
    YAATK_BIN_WRITE(os,atoms_[i]->V);
    YAATK_BIN_WRITE(os,atoms_[i]->coords);
  }

  YAATK_BIN_WRITE(os,check);

  YAATK_BIN_WRITE(os,simTime);
  YAATK_BIN_WRITE(os,dt_);
  YAATK_BIN_WRITE(os,iteration);
  
  YAATK_BIN_WRITE(os,CPUTimeUsed_total);
}



void
SimLoop::writetraj()
{
  static char s[1024];
  sprintf(s,"mde""%010ld",iteration);
  yaatk::text_ofstream fo1(s);
  saveToStream(fo1);
  fo1.close();
}

void
SimLoop::writetrajXVA()
{
  static char s[1024];
  sprintf(s,"mde""%010ld.xva",iteration);
  yaatk::text_ofstream fo1(s);
  saveToStreamXVA(fo1);
  fo1.close();
}


void
SimLoop::writetrajXVA_bin()
{
  static char s[1024];
  sprintf(s,"mde""%010ld.xva.bin",iteration); //#
  yaatk::binary_ofstream fo1(s);
  saveToStreamXVA_bin(fo1);
  fo1.close();
}

void
SimLoop::writetrajXYZ()
{
  static char s[1024];
  sprintf(s,"mde""%010ld.xyz",iteration);
  yaatk::text_ofstream os(s);

  os << atoms_.size() << "\n";
  os << "Sample\n";
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    os << setw(10) << left << ElementString(*(atoms_[i])) << " ";
    os << fixed << right;
    os.precision(3);
    os << setw(10) << atoms_[i]->coords.x/Ao  << " ";
    os << setw(10) << atoms_[i]->coords.y/Ao  << " ";
    os << setw(10) << atoms_[i]->coords.z/Ao  << "\n";
  }

  os.close();
}

void
SimLoop::writestate()
{
  {
  yaatk::binary_ofstream fo1("simloop.conf.bak");
  saveToStream(fo1,YAATK_FSTREAM_BIN);
  fo1.close();
  }

  {
  yaatk::binary_ofstream fo1("simloop.conf");
  saveToStream(fo1,YAATK_FSTREAM_BIN);
  fo1.close();
  }
}

void
SimLoop::loadstate()
{
try
{
  std::cout << "Reading simloop.conf." << std::endl;
  yaatk::binary_ifstream fo1("simloop.conf");
  loadFromStream(fo1,YAATK_FSTREAM_BIN);
  
  REQUIRE(fo1!=0);
  fo1.close();
}  
catch(mdtk::Exception& e)
{ 
    std::cerr << "simloop.conf is corrupted. Trying simloop.conf.bak" << std::endl;
    yaatk::binary_ifstream fo1bak("simloop.conf.bak");
    loadFromStream(fo1bak,YAATK_FSTREAM_BIN);
    if (fo1bak == 0)
    {
      std::cerr << "simloop.conf.bak is ALSO corrupted. Exiting..." << std::endl;
      exit(1);
    }  
    fo1bak.close();
}
}


const std::string SimLoop::id = 
                       std::string("GENERATED_BY_")
                     + "MDE-B"
                     + "-"
                     + "1.0";

/*
void
SimLoop::fixAtom(Atom& atom)
{
  check.energy0 -= atom.M*SQR(atom.V.module())/2.0;
  atom.M = INFINITE_MASS;
  atom.V = Vector3D(0,0,0);
  atom.an = Vector3D(0,0,0);
  atom.an_no_tb = Vector3D(0,0,0);
}
*/


void SimLoop::init_check_energy()
{
  using std::fabs;

  check.temperatureCur = temperature();
  Float energyPotStart = energyPot();
  Float energyKinStart = energyKin();
  check.energyCur = check.energyStart = energyPotStart+energyKinStart;

      check.energy0 = min3(
                            energyPotStart,energyKinStart,
                           (energyPotStart+energyKinStart)
                          )
                      -
                      max3(
                            fabs(energyPotStart),fabs(energyKinStart),
                            fabs(energyPotStart+energyKinStart)
                          );

  cout << "Eo : " << showpos << check.energy0/mdtk::eV << endl;
  cout << "E0 : " << showpos << (check.energyStart)/eV << endl;
  cout << noshowpos;
}  

void SimLoop::do_check_energy()
{
  Float ediff;
    bool reallyPrintCheck = check.checkEnergy && iteration%check.checkEnergyAfter == 0
        && iteration != 0;

  check.temperatureCur = temperature();

  check.energyCur = energy();
  ediff = (check.energyCur-check.energyStart)
          /
          (/*check.energyStart-*/check.energy0);

      if (reallyPrintCheck)
      {
        cout << "E0+Eb : " << showpos << (check.energyStart)/eV << endl;
//        cout << "Ec : " << showpos << (check.energyCur)/eV << endl;
        cout << "E-(E0+Eb) : " << showpos << (check.energyCur-check.energyStart)/eV << endl;
        cout << "dE/Eo : " << showpos << ediff << endl;
        cout << noshowpos;
      }
}


void SimLoop::updateGlobalIndexes()
{
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    atoms_[i]->globalIndex = i;
    atoms_[i]->container = &atoms_;
  }
}  

void SimLoop::initialize()
{
  updateGlobalIndexes();
  TRACE(initNLafterLoading);
  if (initNLafterLoading)
  {
    fpot.NL_init(atoms_);
    fpot.NL_UpdateIfNeeded(atoms_);
  }
}  


void
SimLoop::saveToMDE(std::ostream& fo)
{
  fo << atoms_.size() << std::endl;
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    Atom& Ro_i = *(atoms_[i]);
    YAATK_FSTREAM_WRITE(fo,Ro_i,YAATK_FSTREAM_TEXT);
/*
    fo << Ro_i.Z/mdtk::e << " " << Ro_i.M/mdtk::amu << " " << Ro_i.ID << " " << Ro_i.tag << " " << Ro_i.fixed << " " << Ro_i.thermostat << std::endl;
    fo << Ro_i.V.x << " " << Ro_i.V.y << " " << Ro_i.V.z <<  std::endl;
    fo << Ro_i.coords.x << " " << Ro_i.coords.y << " " << Ro_i.coords.z <<  std::endl;
*/
  }  
  YAATK_FSTREAM_WRITE(fo,simTime,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_WRITE(fo,simTimeSaveTrajInterval,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_WRITE(fo,simTimeFinal,YAATK_FSTREAM_TEXT);

  barrier.SaveToStream(fo,YAATK_FSTREAM_TEXT);

  fo << getPBC() << std::endl;
  thermalBath.SaveToStream(fo,YAATK_FSTREAM_TEXT);
}


void
SimLoop::saveToNanoHive(std::ostream& fo)
{
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    Atom& Ro_i = *(atoms_[i]);
    char et = 'X';
    switch (Ro_i.ID)
    {
      case C_EL: et = 'C'; break;
      case H_EL: et = 'H'; break;
      case Cu_EL: et = 'X'; break;
      case Ar_EL: et = 'X'; break;
      case Ag_EL: et = 'X'; break;
      case Au_EL: et = 'X'; break;
      case DUMMY_EL: et = 'X'; break;
    };
    fo << "<atom id=\"" << i+1 << "\" elementType=\"" << et
       << "\" x3=\"" << Ro_i.coords.x/Ao << "\" y3=\"" << Ro_i.coords.y/Ao << "\" z3=\"" << Ro_i.coords.z/Ao << "\" />\n";
  }  
}

void
SimLoop::loadFromMDE(std::istream& fi)
{
  freeAtoms();
  AtomsContainer& Ro = atoms_;

  Ro.clear();

  int atoms_count;
  fi >> atoms_count;

  cout << "Reading " << atoms_count << " atoms...\n";

  for(int i = 0; i < atoms_count; i++)
  {
   Atom *new_atom;
   new_atom = new Atom();
   YAATK_FSTREAM_READ(fi,*new_atom,YAATK_FSTREAM_TEXT);
   Ro.push_back(new_atom);
  }  

  YAATK_FSTREAM_READ(fi,simTime,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_READ(fi,simTimeSaveTrajInterval,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_READ(fi,simTimeFinal,YAATK_FSTREAM_TEXT);

  barrier.LoadFromStream(fi,YAATK_FSTREAM_TEXT);

  Vector3D PBC;
  fi >> PBC;
  setPBC(PBC);
  thermalBath.LoadFromStream(fi,YAATK_FSTREAM_TEXT);

  initialize();
}


void
SimLoop::loadFromMDE_OLD(std::istream& fi)
{
  freeAtoms();
  AtomsContainer& Ro = atoms_;

  Ro.clear();

  int atoms_count;

  Float x,y,z;
  Float vx,vy,vz;
  Float Z,M;
  int id;
  int tag;
  bool fixed;
  int  thermostat;

  fi >> atoms_count;
  cout << "Reading " << atoms_count << " atoms...\n";
  for(int i = 0; i < atoms_count; i++)
  {

    fi >> Z >> M >> id >> tag >> fixed >> thermostat;
    fi >> vx >> vy >> vz;
    fi >> x >> y >> z;

   Atom *new_atom;
   new_atom = new Atom(e*Z,amu*M,Vector3D(vx,vy,vz),Vector3D(x,y,z)); REQUIRE(new_atom != 0);
   new_atom->ID  = ElementID(id); //!!!
   new_atom->tag  = tag;
   new_atom->fixed  = fixed;
   Ro.push_back(new_atom);
    
  }  

  Vector3D PBC;
  fi >> PBC;
  setPBC(PBC);
  fi >> thermalBath.zMin >> thermalBath.dBoundary;

  initialize();
}

void
SimLoop::heatUpEveryAtom(Float upEnergy)
{
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    Vector3D vn(rand()/Float(RAND_MAX)-1.0,rand()/Float(RAND_MAX)-1.0,rand()/Float(RAND_MAX)-1.0);    
    vn.normalize();
  
    atoms_[i]->V = vn*sqrt(2.0*upEnergy/atoms_[i]->M);
  } 	
}  

void
SimLoop::displaceEveryAtom(Float dist)
{
  for(size_t i = 0; i < atoms_.size(); i++)
  {
    Vector3D vn(rand()/Float(RAND_MAX)-1.0,rand()/Float(RAND_MAX)-1.0,rand()/Float(RAND_MAX)-1.0);    
    vn.normalize();
  
    TRACE(vn*dist/Ao);
    TRACE((vn*dist).module()/Ao);

    atoms_[i]->coords += vn*dist;
  } 	
}  


}// namespace mdtk



