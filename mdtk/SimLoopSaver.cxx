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
   filenamePrefix("md")
{
}

#define MDTK_WRITE_ATOM_ATTRIBUTE_SCALAR(stream,attribute,type)       \
  {                                                                   \
    type x = atom.##attribute;                                        \
    YAATK_BIN_WRITE(stream,x);                                        \
  }

#define MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(stream,attribute,type)       \
  {                                                                   \
    type x = atom.##attribute##.x;                                    \
    YAATK_BIN_WRITE(stream,x);                                        \
                                                                      \
    type y = atom.##attribute##.y;                                    \
    YAATK_BIN_WRITE(stream,y);                                        \
                                                                      \
    type z = atom.##attribute##.z;                                    \
    YAATK_BIN_WRITE(stream,z);                                        \
  }

int
SimLoopSaver::write(std::string filenameBase)
{
  int retval = 0;

  try
  {
    yaatk::binary_ofstream z_stream(filenameBase + ".z");
    yaatk::binary_ofstream r_stream(filenameBase + ".r");
    yaatk::binary_ofstream v_stream(filenameBase + ".v");
    yaatk::binary_ofstream pbc_count_stream(filenameBase + ".pbc_count");
    yaatk::binary_ofstream pbc_rect_stream(filenameBase + ".pbc_rect");
    yaatk::binary_ofstream a_stream(filenameBase + ".a");

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      const Atom& atom = mdloop.atoms[i];

      MDTK_WRITE_ATOM_ATTRIBUTE_SCALAR(z_stream,ID,uint8_t);
      MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(r_stream,coords,double);
      MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(v_stream,V,double);
      MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(pbc_count_stream,PBC_count,int32_t);
      MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(pbc_rect_stream,PBC,double);
      MDTK_WRITE_ATOM_ATTRIBUTE_VECTOR(a_stream,an,double);
    }

    retval = 0;
  }
  catch (...)
  {
    retval = -1;
  }

  return retval;
}

#define MDTK_LOAD_ATOM_ATTRIBUTE_SCALAR_Z(stream,attribute,type)      \
  {                                                                   \
    int dataLength = stream.getDataLength();                          \
                                                                      \
    type z;                                                           \
    REQUIRE(dataLength % sizeof(z) == 0);                             \
    size_t atomsCount = dataLength/sizeof(z);                         \
                                                                      \
    if (mdloop.atoms.size() != 0)                                     \
      REQUIRE(atomsCount == mdloop.atoms.size())                      \
    else                                                              \
      mdloop.atoms.resize(atomsCount);                                \
                                                                      \
    for(size_t i = 0; i < mdloop.atoms.size(); ++i)                   \
    {                                                                 \
      Atom& atom = mdloop.atoms[i];                                   \
                                                                      \
      YAATK_BIN_READ(stream,z);                                       \
      atom.##attribute = ElementID(z);                                \
    }                                                                 \
  }

#define MDTK_LOAD_ATOM_ATTRIBUTE_SCALAR(stream,attribute,type)        \
  {                                                                   \
    int dataLength = stream.getDataLength();                          \
                                                                      \
    type z;                                                           \
    REQUIRE(dataLength % sizeof(z) == 0);                             \
    size_t atomsCount = dataLength/sizeof(z);                         \
                                                                      \
    if (mdloop.atoms.size() != 0)                                     \
      REQUIRE(atomsCount == mdloop.atoms.size())                      \
    else                                                              \
      mdloop.atoms.resize(atomsCount);                                \
                                                                      \
    for(size_t i = 0; i < mdloop.atoms.size(); ++i)                   \
    {                                                                 \
      Atom& atom = mdloop.atoms[i];                                   \
                                                                      \
      YAATK_BIN_READ(stream,z);                                       \
      atom.##attribute = z;                                           \
    }                                                                 \
  }

#define MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(stream,attribute,type)        \
  {                                                                   \
    int dataLength = stream.getDataLength();                          \
                                                                      \
    type x;                                                           \
    REQUIRE(dataLength % (sizeof(x)*3) == 0);                         \
    size_t atomsCount = dataLength/(sizeof(x)*3);                     \
                                                                      \
    if (mdloop.atoms.size() != 0)                                     \
      REQUIRE(atomsCount == mdloop.atoms.size())                      \
    else                                                              \
      mdloop.atoms.resize(atomsCount);                                \
                                                                      \
    for(size_t i = 0; i < mdloop.atoms.size(); ++i)                   \
    {                                                                 \
      Atom& atom = mdloop.atoms[i];                                   \
                                                                      \
      YAATK_BIN_READ(stream,x);                                       \
      atom.##attribute##.x = x;                                       \
                                                                      \
      YAATK_BIN_READ(stream,x);                                     \
      atom.##attribute##.y = x;                                       \
                                                                      \
      YAATK_BIN_READ(stream,x);                                     \
      atom.##attribute##.z = x;                                       \
    }                                                                 \
  }

int
SimLoopSaver::load(std::string filenameBase)
{
  int retval = 0;

  try
  {
    yaatk::binary_ifstream z_stream(filenameBase + ".z");
    MDTK_LOAD_ATOM_ATTRIBUTE_SCALAR_Z(z_stream,ID,uint8_t);
    retval |= LOADED_Z;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream r_stream(filenameBase + ".r");
    MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(r_stream,coords,double);
    retval |= LOADED_R;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream v_stream(filenameBase + ".v");
    MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(v_stream,V,double);
    retval |= LOADED_V;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream pbc_count_stream(filenameBase + ".pbc_count");
    MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(pbc_count_stream,PBC_count,int32_t);
    retval |= LOADED_PBC_COUNT;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream pbc_rect_stream(filenameBase + ".pbc_rect");
    MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(pbc_rect_stream,PBC,double);
    retval |= LOADED_PBC_RECT;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream a_stream(filenameBase + ".a");
    MDTK_LOAD_ATOM_ATTRIBUTE_VECTOR(a_stream,an,double);
    retval |= LOADED_A;
  }
  catch (...)
  {
  }

  return retval;
}

std::string
SimLoopSaver::generateFilenameBase(unsigned long iteration)
{
  std::ostringstream oss;
  oss << filenamePrefix;
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
SimLoopSaver::listFilenameBases()
{
  std::vector<std::string> filenameBases;

  std::set<std::string> filenameBasesSet;

  std::vector<std::string> filenames = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < filenames.size(); ++i)
  {
    std::string filename = filenames[i];

    if (mayContainData(filename))
    {
      size_t filenameBaseEnd = filename.find(".");

      if (filenameBaseEnd != std::string::npos)
      {
        std::string filenameBase = filename.substr(0,filenameBaseEnd);
        TRACE(filenameBase);
        filenameBasesSet.insert(filenameBase);
      }
    }
  }

  {
    std::set<std::string>::iterator it;

    for (it=filenameBasesSet.begin(); it!=filenameBasesSet.end(); ++it)
    {
      filenameBases.push_back(*it);
      TRACE(*it);
    }
  }

  return filenameBases;
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

    if (filename.find(filenamePrefix) == 0 &&
        filename.size() > filenamePrefix.size() &&
        std::isdigit(filename[filenamePrefix.size()],loc))
    {
      std::istringstream iss(filename.substr(filenamePrefix.size()));

      TRACE(filename.substr(filenamePrefix.size()));

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
      TRACE(*it);
    }
  }

  return iterations;
}

int
SimLoopSaver::loadIteration(unsigned long iteration)
{
  return load(generateFilenameBase(iteration));
}

int
SimLoopSaver::loadIterationLatest()
{
  std::vector<unsigned long> iterations = listIterations();
  REQUIRE(iterations.size() > 0);
  return loadIteration(iterations.size()-1);
}

int
SimLoopSaver::write()
{
  return write(generateFilenameBase(mdloop.iteration));
}

}

