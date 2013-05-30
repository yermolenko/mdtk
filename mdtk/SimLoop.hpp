/*
   The molecular dynamics simulation loop class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013
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
  struct LegacyThermalBathStruct
  {
    Float zMin;
    Float dBoundary;
    Float zMinOfFreeZone;
    Float To;
    Float gamma;
    LegacyThermalBathStruct(Float zMin_ = 1000000.0*Ao, Float dBoundary_ = 0.0, Float zMinOfFreeZone_ = -5.0*mdtk::Ao);
    void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode);
    void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode);
    bool isInsideThermalBath(const Atom& a);
  }legacyThermalBathStruct;
  void updateThermalBathFromLegacyStruct()
    {
      thermalBathCommon.To = legacyThermalBathStruct.To;
      thermalBathCommon.gamma = legacyThermalBathStruct.gamma;
      thermalBathGeomBox.zMin = legacyThermalBathStruct.zMin;
      thermalBathGeomBox.dBoundary = legacyThermalBathStruct.dBoundary;
      thermalBathGeomBox.zMinOfFreeZone = legacyThermalBathStruct.zMinOfFreeZone;
      thermalBathGeomType = TB_GEOM_BOX;
    }

  struct ThermalBathCommon
  {
    Float To;
    Float gamma;
    ThermalBathCommon()
      : To(0.0),
        gamma(1.0e13)
      {
      }
  }thermalBathCommon;
  struct ThermalBathGeomBox
  {
    Float zMin;
    Float dBoundary;
    Float zMinOfFreeZone;
    ThermalBathGeomBox()
      : zMin(1000000.0*Ao),
        dBoundary(0.0),
        zMinOfFreeZone(-5.0*Ao)
      {
      }
    bool isInsideThermalBath(const Atom& a)
      {
        REQUIRE(a.lateralPBCEnabled());
        return (a.coords.z > zMin) ||
          (a.lateralPBCEnabled() &&
           (
             (a.coords.x < 0.0 + dBoundary) ||
             (a.coords.x > a.PBC.x - dBoundary) ||
             (a.coords.y < 0.0 + dBoundary) ||
             (a.coords.y > a.PBC.y - dBoundary)
             ) && a.coords.z > zMinOfFreeZone
            );
      }
  }thermalBathGeomBox;
  struct ThermalBathGeomSphere
  {
    Vector3D center;
    Float radius;
    Float zMinOfFreeZone;
    ThermalBathGeomSphere()
      : center(Vector3D(0,0,0)),
        radius(1000000.0*Ao),
        zMinOfFreeZone(-5.0*Ao)
      {
//        REQUIRE(radius - dBoundary > 0);
      }
    bool isInsideThermalBath(const Atom& a)
      {
        REQUIRE(a.PBC == NO_PBC);
        return a.coords.z > zMinOfFreeZone &&
          (center - a.coords).module() > radius;
      }
  }thermalBathGeomSphere;
  enum TB_GEOM_TYPE{TB_GEOM_NONE, TB_GEOM_UNIVERSE, TB_GEOM_BOX, TB_GEOM_SPHERE};
  TB_GEOM_TYPE thermalBathGeomType;
  bool isInsideThermalBath(const Atom& a)
    {
      if (thermalBathGeomType == TB_GEOM_UNIVERSE)
        return true;
      if (thermalBathGeomType == TB_GEOM_BOX &&
          thermalBathGeomBox.isInsideThermalBath(a))
        return true;
      if (thermalBathGeomType == TB_GEOM_SPHERE &&
          thermalBathGeomSphere.isInsideThermalBath(a))
        return true;
      return false;
    }
  bool thermalBathShouldBeApplied(const Atom& atom)
    {
      return atom.apply_ThermalBath &&
        isInsideThermalBath(atom);
    }
  Float actualTemperatureOfThermalBath();
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

