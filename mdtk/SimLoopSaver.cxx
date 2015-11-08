/*
   SimLoopSaver class file.

   Copyright (C) 2013, 2015 Oleksandr Yermolenko
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

  yaatk::DataState ds;

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

    if (mdloop.thermalBathGeomType != mdtk::SimLoop::TB_GEOM_NONE)
    {
      yaatk::binary_ofstream thermal_bath_common_(id + ".thermal_bath.common");
      {
        double_saver(mdloop.thermalBathCommon.To).write(thermal_bath_common_);
        double_saver(mdloop.thermalBathCommon.gamma).write(thermal_bath_common_);
      }
      if (mdloop.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_UNIVERSE)
      {
        yaatk::binary_ofstream thermal_bath_box_(id + ".thermal_bath.univserse");
      }
      if (mdloop.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_BOX)
      {
        yaatk::binary_ofstream thermal_bath_box_(id + ".thermal_bath.box");
        {
          double_saver(mdloop.thermalBathGeomBox.zMin).write(thermal_bath_box_);
          double_saver(mdloop.thermalBathGeomBox.dBoundary).write(thermal_bath_box_);
          double_saver(mdloop.thermalBathGeomBox.zMinOfFreeZone).write(thermal_bath_box_);
        }
      }
      if (mdloop.thermalBathGeomType == mdtk::SimLoop::TB_GEOM_SPHERE)
      {
        yaatk::binary_ofstream thermal_bath_sphere_(id + ".thermal_bath.sphere");
        {
          Vector3D_double_saver(mdloop.thermalBathGeomSphere.center).write(thermal_bath_sphere_);
          double_saver(mdloop.thermalBathGeomSphere.radius).write(thermal_bath_sphere_);
          double_saver(mdloop.thermalBathGeomSphere.zMinOfFreeZone).write(thermal_bath_sphere_);
        }
      }
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
      mdloop.atoms[i].attribute = type##_saver(stream);                 \
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
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,ID,uint8_t);
    retval |= LOADED_Z;
  }
  catch (...)
  {
    VEPRINT("Cannot read Z numbers.\n");
  }

  mdloop.atoms.setAttributesByElementID();

  try
  {
    yaatk::binary_ifstream stream(id + ".r");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,coords,Vector3D_double);
    retval |= LOADED_R;
  }
  catch (...)
  {
    VEPRINT("Cannot read positions.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".v");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,V,Vector3D_double);
    retval |= LOADED_V;
  }
  catch (...)
  {
    VPRINT("Cannot read velocities.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == Vector3D_double_saver().binarySize());
    mdloop.atoms.arrayPBC = Vector3D_double_saver(stream);
    retval |= LOADED_PBC_RECT;
  }
  catch (...)
  {
    VPRINT("Cannot read PBC.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect.count");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,PBC_count,Vector3D_int32_t);
    retval |= LOADED_PBC_COUNT;
  }
  catch (...)
  {
    VPRINT("Cannot read PBC counts.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".pbc_rect.enabled");
    REQUIRE_SILENT(stream.isOpened());
    {
      prepareForAttributeReading(stream,Vector3D_bool_saver().binarySize());
      for(size_t i = 0; i < mdloop.atoms.size(); ++i)
      {
        Vector3D_bool_saver saver(stream);
        if (saver.x() == false)
          mdloop.atoms[i].PBC.x = NO_PBC.x;
        else
          mdloop.atoms[i].PBC.x = mdloop.atoms.arrayPBC.x;
        if (saver.y() == false)
          mdloop.atoms[i].PBC.y = NO_PBC.y;
        else
          mdloop.atoms[i].PBC.y = mdloop.atoms.arrayPBC.y;
        if (saver.z() == false)
          mdloop.atoms[i].PBC.z = NO_PBC.z;
        else
          mdloop.atoms[i].PBC.z = mdloop.atoms.arrayPBC.z;
      }
    }
    retval |= LOADED_PBC_ENABLED;
  }
  catch (...)
  {
    VPRINT("Cannot read 'PBC-enabled' properties.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".a");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,an,Vector3D_double);
    retval |= LOADED_A;
  }
  catch (...)
  {
    VPRINT("Cannot read accelerations.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".indices");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,globalIndex,uint32_t);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read atoms indices.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.thermal_bath_applicable");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,apply_ThermalBath,bool);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'thermal-bath-applicable' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.fixed");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_ATTRIBUTE(stream,fixed,bool);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'fixed' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.target");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_TARGET);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'target' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.projectile");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_PROJECTILE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'projectile' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.substrate");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_SUBSTRATE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'substrate' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.monomer");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_MONOMER);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'monomer' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.cluster");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_CLUSTER);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'cluster' tags.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".tag.fullerene");
    REQUIRE_SILENT(stream.isOpened());
    MDTK_LOAD_ATOM_TAG(stream,ATOMTAG_FULLERENE);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read 'fullerene' tags.\n");
  }

  mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_NONE;

  try
  {
    yaatk::binary_ifstream stream(id + ".thermal_bath.common");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == double_saver().binarySize()*2);
    mdloop.thermalBathCommon.To = double_saver(stream);
    mdloop.thermalBathCommon.gamma = double_saver(stream);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read common thermal bath parameters.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".thermal_bath.univserse");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == 0);
    mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_UNIVERSE;
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read thermal bath parameters for 'univserse' geometry.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".thermal_bath.box");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == double_saver().binarySize()*3);
    mdloop.thermalBathGeomBox.zMin = double_saver(stream);
    mdloop.thermalBathGeomBox.dBoundary = double_saver(stream);
    mdloop.thermalBathGeomBox.zMinOfFreeZone = double_saver(stream);
    mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_BOX;
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read thermal bath parameters for 'box' geometry.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".thermal_bath.sphere");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == Vector3D_double_saver().binarySize() + double_saver().binarySize()*2);
    mdloop.thermalBathGeomSphere.center = Vector3D_double_saver(stream);
    mdloop.thermalBathGeomSphere.radius = double_saver(stream);
    mdloop.thermalBathGeomSphere.zMinOfFreeZone = double_saver(stream);
    mdloop.thermalBathGeomType = mdtk::SimLoop::TB_GEOM_SPHERE;
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read thermal bath parameters for 'sphere' geometry.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".time");
    REQUIRE_SILENT(stream.isOpened());
    REQUIRE(stream.getDataLength() == double_saver().binarySize()*4);
    mdloop.simTime = double_saver(stream);
    mdloop.simTimeFinal = double_saver(stream);
    mdloop.dt = double_saver(stream);
    mdloop.dt_prev = double_saver(stream);
//    retval |= LOADED_;
  }
  catch (...)
  {
    VPRINT("Cannot read time parameters.\n");
  }

  try
  {
    yaatk::binary_ifstream stream(id + ".check");
    REQUIRE_SILENT(stream.isOpened());
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
    VPRINT("Cannot read check values.\n");
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

std::string
SimLoopSaver::extractAttributeName(const std::string filename)
{
  size_t idEnd = filename.find(".");
  REQUIRE(idEnd != std::string::npos);

  size_t attrEnd = filename.rfind(".");
  REQUIRE(attrEnd != std::string::npos);

  REQUIRE(attrEnd - idEnd > 0 && attrEnd - idEnd < 50);

  std::string attr = filename.substr(idEnd + 1, attrEnd - idEnd - 1);

  return attr;
}

std::string
SimLoopSaver::extractId(const std::string filename)
{
  size_t dirEnd = filename.rfind(DIR_DELIMIT_STR);
  if(dirEnd == std::string::npos)
    dirEnd = 0;
  else
    dirEnd++;

  size_t idEnd = filename.find(".",dirEnd);
  REQUIRE(idEnd != 0);

  std::string id = filename.substr(0, idEnd);

  return id;
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

void
SimLoopSaver::removeIterations(const std::vector<unsigned long>& its)
{
  std::vector<std::string> filenames = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < filenames.size(); ++i)
  {
    std::string filename = filenames[i];

    if (filename.find(idPrefix) == 0)
    {
      for(size_t idi = 0; idi < its.size(); ++idi)
      {
        std::string id = generateId(its[idi]);

        if (filename.find(id) == 0 && filename.find(id) != std::string::npos)
        {
          yaatk::remove(filename);
        }
      }
    }
  }
}

void
SimLoopSaver::removeIterations(bool keepFirst, bool keepLast)
{
  std::vector<unsigned long> its = listIterations();

  if (its.size() >= 1 && keepLast)
    its.erase(its.end()-1);
  if (its.size() >= 1 && keepFirst)
    its.erase(its.begin());

  removeIterations(its);
}

void
SimLoopSaver::removeAttributes(const std::string id, const std::set<std::string>& protectedAttributes)
{
  std::vector<std::string> filenames = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < filenames.size(); ++i)
  {
    std::string filename = filenames[i];

    if (filename.find(id) == 0 && filename.find(id) != std::string::npos)
    {
      std::string attribute = extractAttributeName(filename);
      if (protectedAttributes.find(attribute) == protectedAttributes.end())
      {
        yaatk::remove(filename);
      }
    }
  }
}

void
SimLoopSaver::removeAttributesButPos(const std::string id)
{
  std::set<std::string> pa;
  pa.insert("z");
  pa.insert("r");
  removeAttributes(id,pa);
}

void
SimLoopSaver::removeAttributesButPosVel(const std::string id)
{
  std::set<std::string> pa;
  pa.insert("z");
  pa.insert("r");
  pa.insert("v");
  removeAttributes(id,pa);
}

void
SimLoopSaver::removeAttributesButPosVelPBC(const std::string id)
{
  std::set<std::string> pa;
  pa.insert("z");
  pa.insert("r");
  pa.insert("v");
  pa.insert("pbc_rect");
  pa.insert("pbc_rect.count");
  pa.insert("pbc_rect.enabled");
  removeAttributes(id,pa);
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

