/*
   SimLoopSaver class file.

   Copyright (C) 2013 Oleksandr Yermolenko
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

#include <stdint.h>
#include <locale>
#include "SimLoopSaver.hpp"

namespace mdtk
{

SimLoopSaver::SimLoopSaver(SimLoop& mdloopInstance)
  :mdloop(mdloopInstance),
   idPrefix("md")
{
}

struct uint8_t_saver
{
  uint8_t val;
  size_t binarySize() { return sizeof(val); }
  uint8_t_saver(ElementID elid = H_EL):val(uint8_t(elid)) {}
  operator ElementID() { return ElementID(val); }
  uint8_t_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,val); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,val); }
};

struct uint32_t_saver
{
  uint32_t val;
  size_t binarySize() { return sizeof(val); }
  uint32_t_saver(size_t x = 0):val(uint32_t(x)) {}
  operator size_t() { return size_t(val);}
  uint32_t_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,val); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,val); }
};

struct double_saver
{
  double val;
  size_t binarySize() { return sizeof(val); }
  double_saver(double x = 0.0):val(x) {}
  operator double() { return double(val); }
  double_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,val); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,val); }
};

struct bool_saver
{
  uint8_t val;
  size_t binarySize() { return sizeof(val); }
  bool_saver(bool x = false):val(x) {}
  operator bool() { return bool(val); }
  bool_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,val); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,val); }
};

struct Vector3D_double_saver
{
  double x, y, z;
  size_t binarySize() { return 3*sizeof(x); }
  Vector3D_double_saver(Vector3D v = Vector3D(0.0,0.0,0.0)): x(v.x),y(v.y),z(v.z) {}
  operator Vector3D() { return Vector3D(x,y,z); }
  Vector3D_double_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,x); YAATK_BIN_WRITE(os,y); YAATK_BIN_WRITE(os,z); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,x); YAATK_BIN_READ(is,y); YAATK_BIN_READ(is,z); }
};

struct Vector3D_int32_t_saver
{
  int32_t x, y, z;
  size_t binarySize() { return 3*sizeof(x); }
  Vector3D_int32_t_saver(IntVector3D v = IntVector3D(0,0,0)): x(v.x),y(v.y),z(v.z) {}
  operator IntVector3D() { return IntVector3D(int(x),int(y),int(z)); }
  Vector3D_int32_t_saver(yaatk::binary_ifstream& is) { load(is); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,x); YAATK_BIN_WRITE(os,y); YAATK_BIN_WRITE(os,z); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,x); YAATK_BIN_READ(is,y); YAATK_BIN_READ(is,z); }
};

struct Vector3D_bool_saver
{
  uint8_t x_, y_, z_;
  size_t binarySize() { return 3*sizeof(x_); }
  Vector3D_bool_saver(bool x_val = false, bool y_val = false, bool z_val = false): x_(x_val),y_(y_val),z_(z_val) {}
  Vector3D_bool_saver(yaatk::binary_ifstream& is) { load(is); }
  bool x() { return (x_ == 0)?false:true; }
  bool y() { return (y_ == 0)?false:true; }
  bool z() { return (z_ == 0)?false:true; }
  void x(bool v) { (v == false)?(x_ = 0):(x_ = 1); }
  void y(bool v) { (v == false)?(y_ = 0):(y_ = 1); }
  void z(bool v) { (v == false)?(z_ = 0):(z_ = 1); }
  void write(yaatk::binary_ofstream& os) { YAATK_BIN_WRITE(os,x_); YAATK_BIN_WRITE(os,y_); YAATK_BIN_WRITE(os,z_); }
  void load(yaatk::binary_ifstream& is) { YAATK_BIN_READ(is,x_); YAATK_BIN_READ(is,y_); YAATK_BIN_READ(is,z_); }
};

int
SimLoopSaver::write(std::string id)
{
  int retval = 0;

  try
  {
    yaatk::binary_ofstream z_(id + ".z");
    yaatk::binary_ofstream r_(id + ".r");
    yaatk::binary_ofstream v_(id + ".v");
    yaatk::binary_ofstream pbc_rect_(id + ".pbc_rect");
    {
      Vector3D_double_saver(mdloop.atoms.PBC()).write(pbc_rect_);
    }
    yaatk::binary_ofstream pbc_rect_count_(id + ".pbc_rect.count");
    yaatk::binary_ofstream pbc_rect_enabled_(id + ".pbc_rect.enabled");
    yaatk::binary_ofstream a_(id + ".a");

    yaatk::binary_ofstream indexes_(id + ".indices");

    yaatk::binary_ofstream tag_thermal_bath_applicable_(id + ".tag.thermal_bath_applicable");
    yaatk::binary_ofstream tag_fixed_(id + ".tag.fixed");
    yaatk::binary_ofstream tag_target_(id + ".tag.target");
    yaatk::binary_ofstream tag_projectile_(id + ".tag.projectile");
    yaatk::binary_ofstream tag_substrate_(id + ".tag.substrate");
    yaatk::binary_ofstream tag_monomer_(id + ".tag.monomer");
    yaatk::binary_ofstream tag_cluster_(id + ".tag.cluster");
    yaatk::binary_ofstream tag_fullerene_(id + ".tag.fullerene");

    yaatk::binary_ofstream thermal_bath_(id + ".thermal_bath");
    {
      double_saver(mdloop.thermalBath.To).write(thermal_bath_);
      double_saver(mdloop.thermalBath.zMin).write(thermal_bath_);
      double_saver(mdloop.thermalBath.dBoundary).write(thermal_bath_);
      double_saver(mdloop.thermalBath.zMinOfFreeZone).write(thermal_bath_);
      double_saver(mdloop.thermalBath.gamma).write(thermal_bath_);
    }
    yaatk::binary_ofstream time_(id + ".time");
    {
      double_saver(mdloop.simTime).write(time_);
      double_saver(mdloop.simTimeFinal).write(time_);
      double_saver(mdloop.dt).write(time_);
      double_saver(mdloop.dt_prev).write(time_);
    }
    yaatk::binary_ofstream check_(id + ".check");
    {
      bool_saver(mdloop.check.checkForce).write(check_);
      Vector3D_double_saver(mdloop.check.netForce).write(check_);

      bool_saver(mdloop.check.checkEnergy).write(check_);
      double_saver(mdloop.check.initialEnergy).write(check_);
      double_saver(mdloop.check.currentEnergy).write(check_);
      double_saver(mdloop.check.energyTransferredFromBath).write(check_);

      double_saver(mdloop.check.currentTemperature).write(check_);
    }

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      const Atom& atom = mdloop.atoms[i];

      uint8_t_saver(atom.ID).write(z_);
      Vector3D_double_saver(atom.coords).write(r_);
      Vector3D_double_saver(atom.V).write(v_);
      Vector3D_int32_t_saver(atom.PBC_count).write(pbc_rect_count_);
      Vector3D_bool_saver(atom.PBC.x != NO_PBC.x,
                          atom.PBC.y != NO_PBC.y,
                          atom.PBC.z != NO_PBC.z).write(pbc_rect_enabled_);
      Vector3D_double_saver(atom.an).write(a_);
      uint32_t_saver(atom.globalIndex).write(indexes_);

      bool_saver(atom.apply_ThermalBath).write(tag_thermal_bath_applicable_);
      bool_saver(atom.fixed).write(tag_fixed_);
      bool_saver(atom.hasTag(ATOMTAG_TARGET)).write(tag_target_);
      bool_saver(atom.hasTag(ATOMTAG_PROJECTILE)).write(tag_projectile_);
      bool_saver(atom.hasTag(ATOMTAG_SUBSTRATE)).write(tag_substrate_);
      bool_saver(atom.hasTag(ATOMTAG_MONOMER)).write(tag_monomer_);
      bool_saver(atom.hasTag(ATOMTAG_CLUSTER)).write(tag_cluster_);
      bool_saver(atom.hasTag(ATOMTAG_FULLERENE)).write(tag_fullerene_);
    }

    retval = 0;
  }
  catch (...)
  {
    retval = -1;
  }

  return retval;
}

void
SimLoopSaver::prepareForAttributeReading(yaatk::binary_ifstream& is, size_t attributeSize)
{
  int dataLength = is.getDataLength();

  REQUIRE(dataLength % attributeSize == 0);
  size_t atomsCount = dataLength/attributeSize;

  if (mdloop.atoms.size() != 0)
  {
    REQUIRE(atomsCount == mdloop.atoms.size());
  }
  else
  {
    mdloop.atoms.resize(atomsCount);
  }
}

#define MDTK_LOAD_ATOM_ATTRIBUTE(stream,attribute,type)                 \
  {                                                                     \
    prepareForAttributeReading(stream,type##_saver().binarySize());     \
    for(size_t i = 0; i < mdloop.atoms.size(); ++i)                     \
      mdloop.atoms[i].##attribute = type##_saver(stream);               \
  }

#define MDTK_LOAD_ATOM_TAG(stream,tagmask)                              \
  {                                                                     \
    prepareForAttributeReading(stream,bool_saver().binarySize());       \
    for(size_t i = 0; i < mdloop.atoms.size(); ++i)                     \
    {                                                                   \
      bool_saver t(stream);                                             \
      if (bool(t))                                                      \
      {                                                                 \
        mdloop.atoms[i].tag(tagmask);                                   \
      }                                                                 \
      else                                                              \
      {                                                                 \
        mdloop.atoms[i].untag(tagmask);                                 \
      }                                                                 \
    }                                                                   \
  }

int
SimLoopSaver::load(std::string id)
{
  int retval = 0;

  try
  {
    yaatk::binary_ifstream stream(id + ".z");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,ID,uint8_t);
    retval |= LOADED_Z;
  }
  catch (...)
  {
    VEPRINT("Error reading Z numbers.\n");
  }

  mdloop.atoms.setAttributesByElementID();

  try
  {
    yaatk::binary_ifstream stream(id + ".r");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,coords,Vector3D_double);
    retval |= LOADED_R;
  }
  catch (...)
  {
    VEPRINT("Error reading positions.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".v");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,V,Vector3D_double);
    retval |= LOADED_V;
  }
  catch (...)
  {
    VEPRINT("Error reading velocities.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect");
    REQUIRE(stream.getDataLength() == Vector3D_double_saver().binarySize());
    mdloop.atoms.PBC(Vector3D_double_saver(stream));
    retval |= LOADED_PBC_RECT;
  }
  catch (...)
  {
    VEPRINT("Error reading PBC.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect.count");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,PBC_count,Vector3D_int32_t);
    retval |= LOADED_PBC_COUNT;
  }
  catch (...)
  {
    VEPRINT("Error reading PBC counts.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect.enabled");
    {
      prepareForAttributeReading(stream,Vector3D_bool_saver().binarySize());
      for(size_t i = 0; i < mdloop.atoms.size(); ++i)
      {
        Vector3D_bool_saver saver(stream);
        if (saver.x() == false)
          mdloop.atoms[i].PBC.x = NO_PBC.x;
        if (saver.y() == false)
          mdloop.atoms[i].PBC.y = NO_PBC.y;
        if (saver.z() == false)
          mdloop.atoms[i].PBC.z = NO_PBC.z;
      }
    }
    retval |= LOADED_PBC_ENABLED;
  }
  catch (...)
  {
    VEPRINT("Error reading 'PBC-enabled' properties.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".a");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,an,Vector3D_double);
    retval |= LOADED_A;
  }
  catch (...)
  {
    VEPRINT("Error reading accelerations.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".indices");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,globalIndex,uint32_t);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading atoms indices.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.thermal_bath_applicable");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,apply_ThermalBath,bool);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'thermal-bath-applicable' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.fixed");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,fixed,bool);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'fixed' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.target");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_TARGET);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'target' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.projectile");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_PROJECTILE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'projectile' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.substrate");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_SUBSTRATE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'substrate' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.monomer");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_MONOMER);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'monomer' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.cluster");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_CLUSTER);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'cluster' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.fullerene");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_FULLERENE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading 'fullerene' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".thermal_bath");
    REQUIRE(stream.getDataLength() == double_saver().binarySize()*5);
    mdloop.thermalBath.To = double_saver(stream);
    mdloop.thermalBath.zMin = double_saver(stream);
    mdloop.thermalBath.dBoundary = double_saver(stream);
    mdloop.thermalBath.zMinOfFreeZone = double_saver(stream);
    mdloop.thermalBath.gamma = double_saver(stream);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading thermal bath parameters.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".time");
    REQUIRE(stream.getDataLength() == double_saver().binarySize()*4);
    mdloop.simTime = double_saver(stream);
    mdloop.simTimeFinal = double_saver(stream);
    mdloop.dt = double_saver(stream);
    mdloop.dt_prev = double_saver(stream);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading time parameters.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".check");
    REQUIRE(stream.getDataLength() == bool_saver().binarySize()*2 + Vector3D_double_saver().binarySize() + double_saver().binarySize()*4);
    mdloop.check.checkForce = bool_saver(stream);
    mdloop.check.netForce = Vector3D_double_saver(stream);

    mdloop.check.checkEnergy = bool_saver(stream);
    mdloop.check.initialEnergy = double_saver(stream);
    mdloop.check.currentEnergy = double_saver(stream);
    mdloop.check.energyTransferredFromBath = double_saver(stream);

    mdloop.check.currentTemperature = double_saver(stream);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VEPRINT("Error reading check values.\n");
  }

/*
  try
  {
    yaatk::binary_ifstream stream(id + ".tag.");
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_);
//    retval |= LOADED_;
  }
  catch (...) {}

*/

/*
  try
  {
    yaatk::binary_ifstream stream(id + ".");
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,,);
//    retval |= LOADED_;
  }
  catch (...) {}

*/

/*
  try
  {
    yaatk::binary_ifstream stream(id + ".");
    REQUIRE(stream.getDataLength() == _saver().binarySize()*);
    mdloop. = _saver(stream);
//    retval |= LOADED_;
  }
  catch (...) {}

*/

  return retval;
}

std::string
SimLoopSaver::generateId(unsigned long iteration)
{
  std::ostringstream oss;
  oss << idPrefix;
  PRINT2STREAM_FW(oss,iteration,'0',10);
  return oss.str();
}

bool
SimLoopSaver::mayContainData(std::string filename)
{
  return
    filename.find(".z") != std::string::npos ||
    filename.find(".r") != std::string::npos ||
    filename.find(".v") != std::string::npos;
}

std::vector<std::string>
SimLoopSaver::listIds()
{
  std::vector<std::string> ids;

  std::set<std::string> idsSet;

  std::vector<std::string> filenames = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < filenames.size(); ++i)
  {
    std::string filename = filenames[i];

    if (mayContainData(filename))
    {
      size_t idEnd = filename.find(".");

      if (idEnd != std::string::npos)
      {
        std::string id = filename.substr(0,idEnd);
//        TRACE(id);
        idsSet.insert(id);
      }
    }
  }

  {
    std::set<std::string>::iterator it;

    for (it=idsSet.begin(); it!=idsSet.end(); ++it)
    {
      ids.push_back(*it);
//      TRACE(*it);
    }
  }

  return ids;
}

std::vector<unsigned long>
SimLoopSaver::listIterations()
{
  std::vector<unsigned long> iterations;

  std::set<unsigned long> iterationsSet;

  std::vector<std::string> filenames = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < filenames.size(); ++i)
  {
    std::string filename = filenames[i];

    std::locale loc;

    if (filename.find(idPrefix) == 0 &&
        filename.size() > idPrefix.size() &&
        std::isdigit(filename[idPrefix.size()],loc))
    {
      std::istringstream iss(filename.substr(idPrefix.size()));

//      TRACE(filename.substr(idPrefix.size()));

      unsigned long iter;
      iss >> iter;

      iterationsSet.insert(iter);
    }
  }

  {
    std::set<unsigned long>::iterator it;

    for (it=iterationsSet.begin(); it!=iterationsSet.end(); ++it)
    {
      iterations.push_back(*it);
//      TRACE(*it);
    }
  }

  return iterations;
}

int
SimLoopSaver::loadIteration(unsigned long iteration)
{
  mdloop.iteration = iteration;
  return load(generateId(iteration));
}

int
SimLoopSaver::loadIterationLatest()
{
  std::vector<unsigned long> iterations = listIterations();
  REQUIRE(iterations.size() > 0);
  return loadIteration(iterations[iterations.size()-1]);
}

int
SimLoopSaver::write()
{
  return write(generateId(mdloop.iteration));
}

}

