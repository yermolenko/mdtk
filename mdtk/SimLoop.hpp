/*
   The molecular dynamics simulation loop class (header file).

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

#ifndef mdtk_SimLoop_hpp
#define mdtk_SimLoop_hpp

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
  static const std::string id;
public:
  virtual void doBeforeIteration() {};
  virtual void doAfterIteration() {};
  bool allowToFreePotentials;
  bool allowToFreeAtoms;
  void freePotentials()
  {
    if (!allowToFreePotentials) return;
    fpot.freePotentials();
  }
  void freeAtoms()
  {
    if (!allowToFreeAtoms) return;
    for ( size_t it = 0; it < atoms_.size(); it++ )
      delete (atoms_[it]);
  }
  AtomsContainer atoms_;
  AtomsContainer& atoms; // added for compatibility
  void updateGlobalIndexes();
  void initialize();
protected:
  struct Check
  {
    Vector3D fullForce;
    Float energyStart;
    Float energyCur;

    Float energy0;
    
    Float temperatureCur;

    unsigned long checkEnergyAfter;
    bool checkEnergy;
    bool checkForce;
    
    Check(bool ce = true)
     :fullForce(0.0,0.0,0.0),energyStart(1.0),energyCur(1.0),
      energy0(0.0), temperatureCur(0.0),
      checkEnergyAfter(1),checkEnergy(ce),checkForce(true)
    {
    }   
    void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_WRITE(os,fullForce,smode);
      YAATK_FSTREAM_WRITE(os,energyStart,smode);
      YAATK_FSTREAM_WRITE(os,energyCur,smode);
      YAATK_FSTREAM_WRITE(os,energy0,smode);
      YAATK_FSTREAM_WRITE(os,temperatureCur,smode);
      
      YAATK_FSTREAM_WRITE(os,checkEnergyAfter,smode);
      YAATK_FSTREAM_WRITE(os,checkEnergy,smode);
      YAATK_FSTREAM_WRITE(os,checkForce,smode);
    }  
    void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_READ(is,fullForce,smode);
      YAATK_FSTREAM_READ(is,energyStart,smode);
      YAATK_FSTREAM_READ(is,energyCur,smode);
      YAATK_FSTREAM_READ(is,energy0,smode);
      YAATK_FSTREAM_READ(is,temperatureCur,smode);
      
      YAATK_FSTREAM_READ(is,checkEnergyAfter,smode);
      YAATK_FSTREAM_READ(is,checkEnergy,smode);
      YAATK_FSTREAM_READ(is,checkForce,smode);
    }  
  }check;                    
public:

  Float simTime;
  Float simTimeSaveTrajInterval;
  Float simTimeFinal;
  bool  breakSimLoop;
protected:    
  Float timeaccel_;
public:
  Float dt_;
  unsigned long iteration;
  unsigned long iterationFlushStateInterval;
public:
  Float getEnergyCur(){return check.energyCur;}
public: 
  Float energy();
  Float energyPot();
  Float energyKin();
  Float temperature();
  Float temperatureWithoutFixed();
protected:  
  void do_check_energy();
  void init_check_energy();

protected:
  bool checkMIC(); // check Minimum Image Criteria  
  bool fitInCell();          
  void applyPBC(Atom&);
//  void fixAtom(Atom&);
//  bool isFixed(Atom&);
public:
  void initPBC();
public:
  struct Bar // obsolete, will be removed soon
  {
    Float z;
    Float dE;
    Bar(Float z_, Float dE_): z(z_), dE(dE_)
    {
    }
    void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_WRITE(os,z,smode);
      YAATK_FSTREAM_WRITE(os,dE,smode);
    }  
    void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_READ(is,z,smode);
      YAATK_FSTREAM_READ(is,dE,smode);
    }  
  }barrier;                
  void checkOnSpot(Atom&);
  struct ThermalBath
  {
    Float zMin;
    Float dBoundary;
    Float zMinOfFreeZone;
    Float To;
//    void disableGlobally() {zMin = 1000000.0*Ao;dBoundary = 0.0;zMinOfFreeZone = 0.0;}
    ThermalBath(Float zMin_ = 1000000.0*Ao, Float dBoundary_ = 0.0, Float zMinOfFreeZone_ = -5.0*mdtk::Ao)
    : zMin(zMin_), dBoundary(dBoundary_), zMinOfFreeZone(zMinOfFreeZone_), To(0.0)
    {
    }
    void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_WRITE(os,zMin,smode);
      YAATK_FSTREAM_WRITE(os,dBoundary,smode);
      YAATK_FSTREAM_WRITE(os,zMinOfFreeZone,smode);
      YAATK_FSTREAM_WRITE(os,To,smode);
    }  
    void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_READ(is,zMin,smode);
      YAATK_FSTREAM_READ(is,dBoundary,smode);
      YAATK_FSTREAM_READ(is,zMinOfFreeZone,smode);
      YAATK_FSTREAM_READ(is,To,smode);
    }  
  }thermalBath;                
  Float actualTemperatureOfThermalBath();
  bool isWithinThermalBath(Vector3D& c)
  {
    return (c.z > thermalBath.zMin) ||
           ( usePBC() &&
             (
               (c.x < 0.0 + thermalBath.dBoundary) ||
               (c.x > getPBC().x - thermalBath.dBoundary) ||
               (c.y < 0.0 + thermalBath.dBoundary) ||
               (c.y > getPBC().y - thermalBath.dBoundary)
             ) && c.z > thermalBath.zMinOfFreeZone
           )
           ;
  }

public:
  bool initNLafterLoading;


  void writetraj();
  void writestate();
  void loadstate();
  
  SimLoop();
  virtual void init(Float /*px*/,Float /*py*/) {};//= 0;
  virtual ~SimLoop();

  SimLoop(const SimLoop &c);
  SimLoop& operator=(const SimLoop &c);
  void add_simloop(const SimLoop &sl_addon);

  int execute_wo_checks();
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
  void writetrajAccumulated();
  void saveToMDE(std::ostream& fo);
  void loadFromMDE(std::istream& fi);
  void loadFromMDE_OLD(std::istream& fi);
  void saveToNanoHive(std::ostream& fo);
public:
  FProxy fpot;
public:
  double getRcutoff() {return fpot.getRcutoff();}
  void heatUpEveryAtom(Float upEnergy);
  void displaceEveryAtom(Float dist);
private:
private:
  double CPUTimeUsed_prev;
  double CPUTimeUsed_total;
public:
  void setPBC(Vector3D PBC_){atoms_.setPBC(PBC_);}
  Vector3D getPBC() const {return atoms_.getPBC();}
  bool usePBC()const{return atoms_.usePBC();};
};

#define MDTK_MAX_PBC 1e-5

} // namespace mdtk

#endif

