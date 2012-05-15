/*
   The molecular dynamics simulation loop class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012
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
    atoms(),
    check(),
    simTime(0.0),
    simTimeSaveTrajInterval(0.1*ps),
    simTimeFinal(4.0*ps),
    breakSimLoop(false),
    timeaccel(1.0*Ao),
    dt(1e-20),
    dt_prev(1e-20),
    iteration(0),
    iterationFlushStateInterval(1000),
    thermalBath(),
    initNLafterLoading(true),
    allowPartialLoading(false),
    fpot(),
    CPUTimeUsed_prev(0),
    CPUTimeUsed_total(0)
{
  check.checkEnergy = true;
  check.checkForce = true;
}

SimLoop::SimLoop(const SimLoop &c)
  : allowToFreePotentials(true),
    atoms(c.atoms),
    check(c.check),
    simTime(c.simTime),
    simTimeSaveTrajInterval(c.simTimeSaveTrajInterval),
    simTimeFinal(c.simTimeFinal),
    breakSimLoop(false),
    timeaccel(c.timeaccel),
    dt(1e-20), // check this!
    dt_prev(1e-20), // check this!
    iteration(c.iteration),
    iterationFlushStateInterval(c.iterationFlushStateInterval),
    thermalBath(c.thermalBath),
    initNLafterLoading(true),
    allowPartialLoading(false),
    fpot(),
    CPUTimeUsed_prev(0),
    CPUTimeUsed_total(0)
{
  check.checkEnergy = true;
  check.checkForce = true;
}

SimLoop&
SimLoop::operator =(const SimLoop &c)
{
  if (this == &c) return *this;

  atoms.clear();

  allowToFreePotentials = true;
  atoms = c.atoms;
  check = c.check;
  simTime = c.simTime;
  simTimeSaveTrajInterval = c.simTimeSaveTrajInterval;
  simTimeFinal = c.simTimeFinal;
  breakSimLoop = false;
  timeaccel = c.timeaccel;
  dt = 1e-20; // check this!
  dt_prev = 1e-20; // check this!
  iteration = c.iteration;
  iterationFlushStateInterval = c.iterationFlushStateInterval;
  thermalBath = c.thermalBath;
  initNLafterLoading = true;
  allowPartialLoading = false;
//  fpot();
  CPUTimeUsed_prev = 0;
  CPUTimeUsed_total = 0;

  return *this;
}

void
SimLoop::add_simloop(const SimLoop &sl_addon)
{
  atoms.addAtoms(sl_addon.atoms);
}

SimLoop::~SimLoop()
{
  freePotentials();
}

int
SimLoop::execute()
{
  try
  {
    fpot.diagnose();
    PTRACE(atoms.front().PBCEnabled());
    PTRACE(atoms.front().lateralPBCEnabled());
    PTRACE(atoms.back().PBCEnabled());
    PTRACE(atoms.back().lateralPBCEnabled());
    PTRACE(atoms.front().PBC/Ao);
    PTRACE(atoms.back().PBC/Ao);
    PTRACE(thermalBath.zMin/Ao);
    PTRACE(thermalBath.dBoundary/Ao);
    PTRACE(thermalBath.zMinOfFreeZone/Ao);
    PTRACE(thermalBath.To/K);
    return executeMain();
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
    {
      std::ofstream fo("mde_state.after_crash");
      saveToStream(fo);
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

int SimLoop::executeDryRun()
{
// update global indexes, initialize neighbor lists etc
// without actually running simulation
  atoms.prepareForSimulatation();
//  REQUIRE(atoms.checkMIC(fpot.getRcutoff()*2.0));

  PTRACE(initNLafterLoading);
  if (initNLafterLoading)
  {
    fpot.NL_init(atoms);
    fpot.NL_UpdateIfNeeded(atoms);
  }

  return 0;
}

int
SimLoop::executeMain()
{
  gsl_rng * r;
  r = gsl_rng_alloc (gsl_rng_ranlxd2);
  gsl_rng_set(r, 697860L);

  breakSimLoop = false;
  time_t startWallTime = time(NULL);
  time_t curWallTime = startWallTime;

  procmon::ProcmonTimer pmtTotal;

  atoms.prepareForSimulatation();
  if (!atoms.checkMIC(fpot.getRcutoff()*2.0))
  {
    TRACE(atoms.PBC()/Ao);
    TRACE(fpot.getRcutoff()*2.0);
    cerr << "Rcutoff is too large for given PBC !" << endl << flush;
    throw Exception("Rcutoff is too large for given PBC !");
  }

  dt_prev = dt;

  fpot.NL_init(atoms);
  fpot.NL_UpdateIfNeeded(atoms);

  while (simTime < simTimeFinal && !breakSimLoop)
  {
    doBeforeIteration();

    if (verboseTrace)
    {
      cout << "---------" << endl;
      cout << "t : " << simTime << endl;
      cout << "dt : " << dt << endl;
    }

    if (iteration%iterationFlushStateInterval == 0/* && iteration != 0*/)
    {
      if (verboseTrace) cout << "Writing state ... " ;
      writestate();
      if (verboseTrace) cout << "done. " << endl;
    };


    if (simTime == 0.0 || int(simTime/simTimeSaveTrajInterval) != int((simTime - dt_prev)/simTimeSaveTrajInterval))
    {
      if (verboseTrace) cout << "Writing trajectory ... " ;
      writetrajXVA();
      if (verboseTrace) cout << "done. " << endl;
    };

    Float actualThermalBathTemp = actualTemperatureOfThermalBath();
    {
      Float Tb = actualThermalBathTemp;
      PTRACE(Tb);
    }

    if (check.checkForce)
      check.netForce = 0;

    Float v_max = 0.0;

    for(size_t j = 0; j < atoms.size(); j++)
    {
      Atom& atom = atoms[j];

      atom.grad = 0;

      if (atom.isFixed())
      {
        REQUIRE(atom.an == Vector3D(0.0,0.0,0.0));
        REQUIRE(atom.an_no_tb == Vector3D(0.0,0.0,0.0));
        REQUIRE(atom.V == Vector3D(0.0,0.0,0.0));
      }
    }

    if (iteration == 0)
    {
      initEnergyConservationCheck();

      for(size_t j = 0; j < atoms.size(); j++)
      {
        Atom& atom = atoms[j];

        if (atom.isFixed()) continue;

        Vector3D  force = -atom.grad;
        atom.an = force/atom.M;
        atom.an_no_tb = force/atom.M;

        atom.grad = 0;
      }
    }

    for(size_t j = 0; j < atoms.size(); j++)
    {
      Atom& atom = atoms[j];

      if (atom.isFixed()) continue;

      Vector3D dr = atom.V*dt + atom.an*dt*dt/2.0; // eq 1
      atom.coords += dr;
      fpot.incDisplacement(atoms[j],dr);
    }

    fpot.NL_checkRequestUpdate(atoms);
    fpot.NL_UpdateIfNeeded(atoms);

    doEnergyConservationCheck();

    for(size_t j = 0; j < atoms.size(); j++)
    {
      Atom& atom = atoms[j];

      Vector3D  force = -atom.grad;

      if (check.checkForce)
        check.netForce += force;

      if (atom.isFixed()) continue;

      Vector3D  vdt2 = atom.V + atom.an*dt/2.0; // eq 2

      if (isWithinThermalBath(atom) && atom.apply_ThermalBath)
      {
//        Float T = check.temperatureCur;
        Float T = actualThermalBathTemp;
        Float To_by_T = (fabs(T)<1e-5)?0:(thermalBath.To/T);
        Float max_To_by_T = 5.0;
        if (To_by_T < -max_To_by_T) To_by_T = -max_To_by_T;
        if (To_by_T > +max_To_by_T) To_by_T = +max_To_by_T;

        Float gamma = 1.0e13;
//        Float gamma = 1.0/(1000.0*dt);

        Vector3D dforce = -atom.V*atom.M*gamma*(1.0-sqrt(To_by_T));

        // try to account energy transfered to thermalbath
        // only required to perform energy conservation check
        {
          Vector3D dv_no_tb,dv;

          {
            Vector3D an_no_tb_new = force/atom.M;
            // from eq 2 and eq 4
            dv_no_tb = atom.an_no_tb*dt/2.0 + an_no_tb_new*dt/2.0;
            atom.an_no_tb = an_no_tb_new;
          }

          {
            Vector3D an_new = (force + dforce)/atom.M;
            // from eq 2 and eq 4
            dv = atom.an*dt/2.0 + an_new*dt/2.0;
            atom.an = an_new;
          }

          Float dEkin = atom.M*SQR((atom.V+dv)      .module())/2.0
                      - atom.M*SQR((atom.V+dv_no_tb).module())/2.0;
          check.energyTransferredFromBath += dEkin;
        }

        force += dforce;
      }

      atom.an = force/atom.M; //eq 3

      atom.V  = vdt2 + atom.an*dt/2.0; //eq 4

      Float v = atom.V.module();
      if (v > v_max &&
          fpot.hasNB(atom) &&
          atom.coords.module() < 500.0*Ao)
        v_max = v;
    }

    atoms.applyPBC();

    doAfterIteration();

    simTime += dt;

    dt_prev = dt;

    const Float dt_max = 5e-16;
    const Float dt_min = 1e-20;
    if (v_max != 0.0)
      dt = 0.05*timeaccel/v_max;
    else
      dt = dt_max;
    if (dt > dt_max) dt = dt_max;
    if (dt < dt_min) dt = dt_min;

    if (check.checkForce)
    {
      if (check.netForce.module() > 1e-8)
      {
        if (verboseTrace)
        {
          Float ffmod = check.netForce.module();
          PTRACE(ffmod);
        }
        cerr << "FullForce != 0" << endl << flush;
        cout << "FullForce != 0" << endl << flush;
      }
    }

    if (verboseTrace)
      cout << "it : " << iteration << endl;

    curWallTime = time(NULL);
    if (verboseTrace)
      cout << "tw : " << (curWallTime-startWallTime) << endl;

    CPUTimeUsed_total = CPUTimeUsed_prev + pmtTotal.getTimeInSeconds();

    if (verboseTrace)
      cout << "tc : " << CPUTimeUsed_total << endl;

    iteration = (iteration+1)%2000000000L;
  };

  if (verboseTrace)
  {
    cout << "Final Modeling time = " << simTime << endl;
    cout << "Final Time step = " << dt << endl;
    cout << "Modeling cycle complete " << endl;

    cout << "--------------------------------------------------------- " << endl;
  }
  {
    if (verboseTrace) cout << "Writing trajectory ... ";
    writetrajXVA();
    if (verboseTrace) cout << "done. " << endl;
  };

  curWallTime = time(NULL);
  if (verboseTrace)
    cout << "Wall TIME used = " << (curWallTime-startWallTime) << endl;

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
  return fpot(atoms);
}

Float
SimLoop::energyKin()
{
  Float energyKinCur = 0.0;
  int j,atoms_count;
  atoms_count = atoms.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = atoms[j];
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
  atoms_count = atoms.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = atoms[j];
    if (atom.isFixed()) continue;
    if (isWithinThermalBath(atom) && atom.apply_ThermalBath)
    {
      energyKinCur += atom.M*SQR(atom.V.module())/2.0;
      atoms_accounted++;
    }
  };

  if (atoms_accounted == 0)
    return 0.0;

  return energyKinCur/(3.0/2.0*kb*atoms_accounted);
}

Float
SimLoop::temperatureWithoutFixed()
{
  Float energyKinCur = 0.0;
  int j,atoms_count,atoms_accounted = 0;
  atoms_count = atoms.size();
  for(j = 0; j < atoms_count; j++)
  {
    Atom& atom = atoms[j];
    if (!atom.isFixed())
    {
      energyKinCur += atom.M*SQR(atom.V.module())/2.0;
      atoms_accounted++;
    }
  };

  if (atoms_accounted == 0)
    return 0.0;

  return energyKinCur/(3.0/2.0*kb*atoms_accounted);
}

void
SimLoop::loadFromStream(istream& is, YAATK_FSTREAM_MODE smode)
{
  atoms.loadFromStream(is,smode);

  if (verboseTrace) cout << "Reading timing info ... " << endl;

  check.loadFromStream(is,smode);

  YAATK_FSTREAM_READ(is,simTime,smode);
  YAATK_FSTREAM_READ(is,simTimeSaveTrajInterval,smode);
  YAATK_FSTREAM_READ(is,simTimeFinal,smode);
  YAATK_FSTREAM_READ(is,timeaccel,smode);
  YAATK_FSTREAM_READ(is,dt,smode);
  YAATK_FSTREAM_READ(is,iteration,smode);
  YAATK_FSTREAM_READ(is,iterationFlushStateInterval,smode);

  thermalBath.loadFromStream(is,smode);

//  fpot.LoadFromStream(is,smode);

  YAATK_FSTREAM_READ(is,CPUTimeUsed_prev,smode);
  CPUTimeUsed_total = CPUTimeUsed_prev;

  if (verboseTrace) cout << "Parsing of state file done. " << endl;

  executeDryRun();

  allowPartialLoading = true;
}

void
SimLoop::saveToStream(ostream& os, YAATK_FSTREAM_MODE smode)
{
  std::streamsize prevStreamSize = os.precision();
  os.precision(FLOAT_PRECISION);

  atoms.saveToStream(os,smode);

  check.saveToStream(os,smode);

  YAATK_FSTREAM_WRITE(os,simTime,smode);
  YAATK_FSTREAM_WRITE(os,simTimeSaveTrajInterval,smode);
  YAATK_FSTREAM_WRITE(os,simTimeFinal,smode);
  YAATK_FSTREAM_WRITE(os,timeaccel,smode);
  YAATK_FSTREAM_WRITE(os,dt,smode);
  YAATK_FSTREAM_WRITE(os,iteration,smode);
  YAATK_FSTREAM_WRITE(os,iterationFlushStateInterval,smode);

  thermalBath.saveToStream(os,smode);

//  fpot.SaveToStream(os,smode);

  YAATK_FSTREAM_WRITE(os,CPUTimeUsed_total,smode);

  os.precision(prevStreamSize);
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
  REQUIRE(atoms_count == atoms.size());

  if (verboseTrace) cout << "Reading XVA info about " << atoms_count << " atoms..." << endl;

  for(i = 0; i < atoms_count; i++)
  {
    is >> atoms[i].V;
    is >> atoms[i].coords;
    is >> atoms[i].PBC_count;

    atoms[i].V *= XVA_VELOCITY_SCALE;
    atoms[i].coords *= XVA_DISTANCE_SCALE;
  }
  if (verboseTrace) cout << endl;

  check.loadFromStream(is,YAATK_FSTREAM_TEXT);

  is >> simTime;
  is >> dt;
  is >> iteration; //iteration++;

  is >> CPUTimeUsed_prev;
  CPUTimeUsed_total = CPUTimeUsed_prev;

  if (verboseTrace) cout << "Parsing of state file done. " << endl;

  executeDryRun();
}

void
SimLoop::saveToStreamXVA(ostream& os)
{
  Float XVA_VELOCITY_SCALE = 0.0;//1.0e3;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    XVA_VELOCITY_SCALE += fabs(atoms[i].V.x);
    XVA_VELOCITY_SCALE += fabs(atoms[i].V.y);
    XVA_VELOCITY_SCALE += fabs(atoms[i].V.z);
  }

  XVA_VELOCITY_SCALE /= (3.0*atoms.size());
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

  int i,atoms_count = atoms.size();
  os << atoms_count << "\n";
  for(i = 0; i < atoms_count; i++)
  {
#ifndef DONT_USE_XVASCALE
    os << fixed;
#endif
    os.precision(2);
    os << atoms[i].V.x/XVA_VELOCITY_SCALE  << " ";
    os << atoms[i].V.y/XVA_VELOCITY_SCALE  << " ";
    os << atoms[i].V.z/XVA_VELOCITY_SCALE  << "\n";

#ifndef DONT_USE_XVASCALE
    os << fixed;
#endif
    os.precision(2);
    os << atoms[i].coords.x/XVA_DISTANCE_SCALE  << " ";
    os << atoms[i].coords.y/XVA_DISTANCE_SCALE  << " ";
    os << atoms[i].coords.z/XVA_DISTANCE_SCALE  << "\n";

    os << atoms[i].PBC_count.x << " ";
    os << atoms[i].PBC_count.y << " ";
    os << atoms[i].PBC_count.z << "\n";
  }
  os << scientific;

  os.precision(FLOAT_PRECISION);

  check.saveToStream(os,YAATK_FSTREAM_TEXT);

  os << simTime << "\n";
  os << dt << "\n";
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
  REQUIRE(atoms_count == atoms.size());

  if (verboseTrace) cout << "Reading XVA info about " << atoms_count << " atoms..." << endl;

  for(i = 0; i < atoms_count; i++)
  {
    YAATK_BIN_READ(is,atoms[i].V);
    YAATK_BIN_READ(is,atoms[i].coords);
  }
  if (verboseTrace) cout << endl;

  YAATK_BIN_READ(is,check);

  YAATK_BIN_READ(is,simTime);
  YAATK_BIN_READ(is,dt);
  YAATK_BIN_READ(is,iteration); //iteration++;

  TRACE(iteration);

  YAATK_BIN_READ(is,CPUTimeUsed_prev);
  CPUTimeUsed_total = CPUTimeUsed_prev;

  if (verboseTrace) cout << "Parsing of state file done. " << endl;

  executeDryRun();
}

void
SimLoop::saveToStreamXVA_bin(ostream& os)
{
  int i,atoms_count = atoms.size();
  YAATK_BIN_WRITE(os,atoms_count);
  for(i = 0; i < atoms_count; i++)
  {
    YAATK_BIN_WRITE(os,atoms[i].V);
    YAATK_BIN_WRITE(os,atoms[i].coords);
  }

  YAATK_BIN_WRITE(os,check);

  YAATK_BIN_WRITE(os,simTime);
  YAATK_BIN_WRITE(os,dt);
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
  sprintf(s,"mde""%010ld.xva.bin",iteration);
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

  os << atoms.size() << "\n";
  os << "Sample\n";
  for(size_t i = 0; i < atoms.size(); i++)
  {
    os << setw(10) << left << ElementString(atoms[i]) << " ";
    os << fixed << right;
    os.precision(3);
    os << setw(10) << atoms[i].coords.x/Ao  << " ";
    os << setw(10) << atoms[i].coords.y/Ao  << " ";
    os << setw(10) << atoms[i].coords.z/Ao  << "\n";
  }

  os.close();
}

void
SimLoop::writetrajAccumulated(const std::vector<size_t>& atomIndices)
{
  Float XVA_VELOCITY_SCALE = 0.0;//1.0e3;
  Float XVA_DISTANCE_SCALE = Ao;
  {
    for(size_t i = 0; i < atoms.size(); i++)
    {
      XVA_VELOCITY_SCALE += fabs(atoms[i].V.x);
      XVA_VELOCITY_SCALE += fabs(atoms[i].V.y);
      XVA_VELOCITY_SCALE += fabs(atoms[i].V.z);
    }

    XVA_VELOCITY_SCALE /= (3.0*atoms.size());
    if (fabs(XVA_VELOCITY_SCALE) < 1e-30) XVA_VELOCITY_SCALE = 1e-30;
    XVA_VELOCITY_SCALE = pow(10.0,ceil(log10(fabs(XVA_VELOCITY_SCALE)))-2);
  }

  yaatk::text_ifstream accPrev("acc");
  yaatk::text_ofstream acc("acc_next");

  size_t stateCount = 0;
  if (accPrev.isOpened())
  {
    size_t atomsCount;
    accPrev >> stateCount;
    accPrev >> atomsCount;
    REQUIRE(atomsCount == atomIndices.size());
  }
  acc << stateCount+1 << "\n";
  acc << atomIndices.size() << "\n";

  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    size_t atomIndex = 0;
    if (stateCount > 0)
      accPrev >> atomIndex;
    else
      atomIndex = atomIndices[ai];

    REQUIRE(atomIndex == atomIndices[ai]);

    acc << atomIndex << " ";
  }
  acc << "\n";

  Float tempFloat;
  int   tempInt;

  {
    for(size_t si = 0; si < stateCount; si++)
    {
      accPrev >> tempFloat;
      acc << tempFloat << " ";
    }
    acc << XVA_VELOCITY_SCALE << "\n";
  }

  {
    for(size_t si = 0; si < stateCount; si++)
    {
      accPrev >> tempFloat;
      acc << tempFloat << " ";
    }
    acc << XVA_DISTANCE_SCALE << "\n";
  }

  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        accPrev >> tempInt;
        acc << tempInt << " ";
      }
      acc << atoms[atomIndices[ai]].PBC_count.X(ci) << "\n";
    }
  }
  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        accPrev >> tempFloat;
        acc << fixed << setprecision(2) << tempFloat << " ";
      }
      acc << fixed << setprecision(2) << atoms[atomIndices[ai]].coords.X(ci)/XVA_DISTANCE_SCALE << "\n";
    }
  }
  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        accPrev >> tempFloat;
        acc << tempFloat << " ";
      }
      acc << fixed << setprecision(2) << atoms[atomIndices[ai]].V.X(ci)/XVA_VELOCITY_SCALE << "\n";
    }
  }
  if (verboseTrace) cout << endl;

  accPrev.close();
  acc.close();

  yaatk::remove("acc.xz");
  yaatk::rename("acc_next.xz","acc.xz");
}

void
SimLoop::writetrajAccumulated_bin(const std::vector<size_t>& atomIndices)
{
  yaatk::binary_ifstream accPrev("acc.bin");
  yaatk::binary_ofstream acc("acc_next.bin");

  size_t stateCount = 0;
  if (accPrev.isOpened())
  {
    size_t atomsCount;
    YAATK_BIN_READ(accPrev,stateCount);
    YAATK_BIN_READ(accPrev,atomsCount);
    REQUIRE(atomsCount == atomIndices.size());
  }
  size_t stateCount_newval = stateCount+1;
  YAATK_BIN_WRITE(acc,stateCount_newval);
  size_t atomsCount_val = atomIndices.size();
  YAATK_BIN_WRITE(acc,atomsCount_val);

  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    size_t atomIndex = 0;
    if (stateCount > 0)
      YAATK_BIN_READ(accPrev,atomIndex);
    else
      atomIndex = atomIndices[ai];

    REQUIRE(atomIndex == atomIndices[ai]);

    YAATK_BIN_WRITE(acc,atomIndex);
  }

  Float tempFloat;
  int   tempInt;

  {
    for(size_t si = 0; si < stateCount; si++)
    {
      YAATK_BIN_READ(accPrev,tempFloat);
      YAATK_BIN_WRITE(acc,tempFloat);
    }
    YAATK_BIN_WRITE(acc,simTime);
  }

  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        YAATK_BIN_READ(accPrev,tempInt);
        YAATK_BIN_WRITE(acc,tempInt);
      }
      YAATK_BIN_WRITE(acc,atoms[atomIndices[ai]].PBC_count.X(ci));
    }
  }
  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        YAATK_BIN_READ(accPrev,tempFloat);
        YAATK_BIN_WRITE(acc,tempFloat);
      }
      YAATK_BIN_WRITE(acc,atoms[atomIndices[ai]].coords.X(ci));
    }
  }
  for(size_t ai = 0; ai < atomIndices.size(); ai++)
  {
    for(size_t ci = 0; ci < 3; ci++)
    {
      for(size_t si = 0; si < stateCount; si++)
      {
        YAATK_BIN_READ(accPrev,tempFloat);
        YAATK_BIN_WRITE(acc,tempFloat);
      }
      YAATK_BIN_WRITE(acc,atoms[atomIndices[ai]].V.X(ci));
    }
  }
  if (verboseTrace) cout << endl;

  accPrev.close();
  acc.close();

  yaatk::remove("acc.bin.xz");
  yaatk::rename("acc_next.bin.xz","acc.bin.xz");
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

void
SimLoop::initEnergyConservationCheck()
{
  check.currentTemperature = temperature();
  check.currentEnergy = check.initialEnergy = energy();
  check.energyTransferredFromBath = 0;

  if (verboseTrace)
  {
    cout << "Eo : " << showpos << check.initialEnergy/eV << endl;
    cout << noshowpos;
  }

  if(std::fabs(check.initialEnergy) < 0.001*eV)
    if (verboseTrace)
      cerr << "Total energy is less than 0.001*eV." << endl;
}

void SimLoop::doEnergyConservationCheck()
{
  check.currentTemperature = temperature();
  check.currentEnergy = energy();

  Float Eo_plus_Eb = check.initialEnergy + check.energyTransferredFromBath;
  Float dE = check.currentEnergy - Eo_plus_Eb;
  Float dE_by_Eo_plus_Eb = dE/Eo_plus_Eb;

  bool reallyPrintCheck = check.checkEnergy && iteration != 0;
  if (reallyPrintCheck && verboseTrace)
  {
    cout << "Eb : " << showpos <<
      check.energyTransferredFromBath/eV << endl;
    cout << "dE : " << showpos << dE/eV << endl;
    cout << "dE/(Eo+Eb) : " << showpos << dE_by_Eo_plus_Eb << endl;
    cout << noshowpos;
  }

  if(std::fabs(check.currentEnergy) < 0.001*eV)
    if (verboseTrace)
      cerr << "Total energy is less than 0.001*eV." << endl;
}

void
SimLoop::saveToMDE(std::ostream& fo)
{
  fo << atoms.size() << std::endl;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& Ro_i = atoms[i];
    YAATK_FSTREAM_WRITE(fo,Ro_i,YAATK_FSTREAM_TEXT);
  }
  YAATK_FSTREAM_WRITE(fo,simTime,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_WRITE(fo,simTimeSaveTrajInterval,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_WRITE(fo,simTimeFinal,YAATK_FSTREAM_TEXT);

  REQUIRE(atoms.size() > 0);
  fo << atoms.PBC() << std::endl;
  thermalBath.saveToStream(fo,YAATK_FSTREAM_TEXT);
}

void
SimLoop::saveToNanoHive(std::ostream& fo)
{
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Atom& Ro_i = atoms[i];
    char et = 'X';
    switch (Ro_i.ID)
    {
      case C_EL: et = 'C'; break;
      case H_EL: et = 'H'; break;
      case Cu_EL: et = 'X'; break;
      case Ar_EL: et = 'X'; break;
      case Xe_EL: et = 'X'; break;
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
  atoms.clear();

  int atoms_count;
  fi >> atoms_count;

  if (verboseTrace) cout << "Reading " << atoms_count << " atoms...\n";
  atoms.resize(atoms_count);

  for(int i = 0; i < atoms_count; i++)
  {
    atoms[i] = Atom();
    YAATK_FSTREAM_READ(fi,atoms[i],YAATK_FSTREAM_TEXT);
  }

  YAATK_FSTREAM_READ(fi,simTime,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_READ(fi,simTimeSaveTrajInterval,YAATK_FSTREAM_TEXT);
  YAATK_FSTREAM_READ(fi,simTimeFinal,YAATK_FSTREAM_TEXT);

  Vector3D PBC;
  fi >> PBC;
  setPBC(PBC);
  thermalBath.loadFromStream(fi,YAATK_FSTREAM_TEXT);

  executeDryRun();
}

void
SimLoop::heatUpEveryAtom(Float upEnergy, gsl_rng* rng)
{
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Vector3D vn(gsl_rng_uniform(rng)-1.0,gsl_rng_uniform(rng)-1.0,gsl_rng_uniform(rng)-1.0);
    vn.normalize();

    atoms[i].V = vn*sqrt(2.0*upEnergy/atoms[i].M);
  }
}

void
SimLoop::displaceEveryAtom(Float dist, gsl_rng* rng)
{
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Vector3D vn(gsl_rng_uniform(rng)-1.0,gsl_rng_uniform(rng)-1.0,gsl_rng_uniform(rng)-1.0);
    vn.normalize();

    TRACE(vn*dist/Ao);
    TRACE((vn*dist).module()/Ao);

    atoms[i].coords += vn*dist;
  }
}

SimLoop::Check::Check(bool ce)
  :checkForce(true),netForce(0.0,0.0,0.0),
   checkEnergy(ce),initialEnergy(1.0),currentEnergy(1.0),
   energyTransferredFromBath(0.0),
   currentTemperature(0.0)
{
}

void
SimLoop::Check::saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
{
  YAATK_FSTREAM_WRITE(os,checkForce,smode);
  YAATK_FSTREAM_WRITE(os,netForce,smode);

  YAATK_FSTREAM_WRITE(os,checkEnergy,smode);
  YAATK_FSTREAM_WRITE(os,initialEnergy,smode);
  YAATK_FSTREAM_WRITE(os,currentEnergy,smode);
  YAATK_FSTREAM_WRITE(os,energyTransferredFromBath,smode);

  YAATK_FSTREAM_WRITE(os,currentTemperature,smode);
}

void
SimLoop::Check::loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
{
  YAATK_FSTREAM_READ(is,checkForce,smode);
  YAATK_FSTREAM_READ(is,netForce,smode);

  YAATK_FSTREAM_READ(is,checkEnergy,smode);
  YAATK_FSTREAM_READ(is,initialEnergy,smode);
  YAATK_FSTREAM_READ(is,currentEnergy,smode);
  YAATK_FSTREAM_READ(is,energyTransferredFromBath,smode);

  YAATK_FSTREAM_READ(is,currentTemperature,smode);
}

/*
void
ThermalBath::disableGlobally()
{
  zMin = 1000000.0*Ao;dBoundary = 0.0;zMinOfFreeZone = 0.0;
}
*/

SimLoop::ThermalBath::ThermalBath(Float zMin_, Float dBoundary_, Float zMinOfFreeZone_)
  :zMin(zMin_), dBoundary(dBoundary_), zMinOfFreeZone(zMinOfFreeZone_), To(0.0)
{
}

void
SimLoop::ThermalBath::saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
{
  YAATK_FSTREAM_WRITE(os,zMin,smode);
  YAATK_FSTREAM_WRITE(os,dBoundary,smode);
  YAATK_FSTREAM_WRITE(os,zMinOfFreeZone,smode);
  YAATK_FSTREAM_WRITE(os,To,smode);
}

void
SimLoop::ThermalBath::loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
{
  YAATK_FSTREAM_READ(is,zMin,smode);
  YAATK_FSTREAM_READ(is,dBoundary,smode);
  YAATK_FSTREAM_READ(is,zMinOfFreeZone,smode);
  YAATK_FSTREAM_READ(is,To,smode);
}

bool
SimLoop::isWithinThermalBath(const Atom& a)
{
  return (a.coords.z > thermalBath.zMin) ||
    (a.lateralPBCEnabled() &&
     (
       (a.coords.x < 0.0 + thermalBath.dBoundary) ||
       (a.coords.x > a.PBC.x - thermalBath.dBoundary) ||
       (a.coords.y < 0.0 + thermalBath.dBoundary) ||
       (a.coords.y > a.PBC.y - thermalBath.dBoundary)
       ) && a.coords.z > thermalBath.zMinOfFreeZone
      );
}

}// namespace mdtk
