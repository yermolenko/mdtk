/*
   The molecular dynamics simulation loop class (header file).

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

#ifndef mdtk_SimLoop_hpp
#define mdtk_SimLoop_hpp

#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>

#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>
#include <mdtk/Vector3D.hpp>
#include <mdtk/Atom.hpp>

#include <mdtk/potentials/FProxy.hpp>

#include <string>

#include <mdtk/procmon.hpp>

namespace mdtk
{

class SimLoop
{
public:
  virtual void doBeforeIteration() {};
  virtual void doAfterIteration() {};
  bool allowToFreePotentials;
  bool preventFileOutput;
  void freePotentials()
  {
    if (!allowToFreePotentials) return;
    fpot.freePotentials();
  }
  AtomsArray atoms;
public:
  struct Check
  {
    bool checkForce;
    Vector3D netForce;

    bool checkEnergy;
    Float initialEnergy;
    Float currentEnergy;
    Float energyTransferredFromBath;

    Float currentTemperature;

    Check(bool ce = true);
    void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode);
    void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode);
  }check;
  void forgetHistory() { iteration = 0; }
public:
  Float simTime;
  Float simTimeSaveTrajInterval;
  Float simTimeFinal;
  bool  breakSimLoop;
protected:
  Float timeaccel;
public:
  Float dt;
  Float dt_prev;
  unsigned long iteration;
  unsigned long iterationFlushStateInterval;
public:
  Float energy();
  Float energyPot();
  Float energyKin();
  Float temperature();
  Float temperatureWithoutFixed();
protected:
  void doEnergyConservationCheck();
  void initEnergyConservationCheck();
public:
  struct ThermalBath
  {
    Float zMin;
    Float dBoundary;
    Float zMinOfFreeZone;
    Float To;
    Float gamma;
    ThermalBath(Float zMin_ = 1000000.0*Ao, Float dBoundary_ = 0.0, Float zMinOfFreeZone_ = -5.0*mdtk::Ao);
    void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode);
    void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode);
  }thermalBath;
  Float actualTemperatureOfThermalBath();
  bool isWithinThermalBath(const Atom& a);
public:
  bool initNLafterLoading;

  void writetraj();
  void writestate();
  void loadstate();

  SimLoop();
  virtual ~SimLoop();

  SimLoop(const SimLoop &c);
  SimLoop& operator=(const SimLoop &c);
  void add_simloop(const SimLoop &sl_addon);

  int executeDryRun();
  int executeMain();
  int execute();

  void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE = YAATK_FSTREAM_TEXT);
  void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE = YAATK_FSTREAM_TEXT);

  bool allowPartialLoading;
  void saveToStreamXVA(std::ostream& os);
  void loadFromStreamXVA(std::istream& is);
  void saveToStreamXVA_bin(std::ostream& os);
  void loadFromStreamXVA_bin(std::istream& is);
  void writetrajXVA();
  void writetrajXVA_bin();
  void writetrajXYZ();
  void writetrajAccumulated(const std::vector<size_t>& atomIndices);
  void writetrajAccumulated_bin(const std::vector<size_t>& atomIndices);
  void saveToMDE(std::ostream& fo);
  void loadFromMDE(std::istream& fi);
  void loadFromMDE_OLD(std::istream& fi);
  void saveToNanoHive(std::ostream& fo);
public:
  FProxy fpot;
public:
  double getRcutoff() {return fpot.getRcutoff();}
  void heatUpEveryAtom(Float upEnergy, gsl_rng* rng);
  void displaceEveryAtom(Float dist, gsl_rng* rng);
private:
  double CPUTimeUsed_prev;
  double CPUTimeUsed_total;
public:
  // these functions are obsolete
  void setPBC(Vector3D PBC_){atoms.PBC(PBC_);}
};

} // namespace mdtk

#endif

