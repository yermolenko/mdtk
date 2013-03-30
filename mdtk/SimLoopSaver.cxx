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

int
SimLoopSaver::write(std::string filenameBase)
{
  int retval = 0;

  try
  {
    yaatk::binary_ofstream z_stream(filenameBase + ".z");
    yaatk::binary_ofstream r_stream(filenameBase + ".r");
    yaatk::binary_ofstream v_stream(filenameBase + ".v");
    yaatk::binary_ofstream pbc_stream(filenameBase + ".pbc");

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      const Atom& atom = mdloop.atoms[i];

      {
        uint8_t z = atom.ID;
        YAATK_BIN_WRITE(z_stream,z);
      }

      {
        double x = atom.coords.x;
        YAATK_BIN_WRITE(r_stream,x);

        double y = atom.coords.y;
        YAATK_BIN_WRITE(r_stream,y);

        double z = atom.coords.z;
        YAATK_BIN_WRITE(r_stream,z);
      }

      {
        double x = atom.V.x;
        YAATK_BIN_WRITE(v_stream,x);

        double y = atom.V.y;
        YAATK_BIN_WRITE(v_stream,y);

        double z = atom.V.z;
        YAATK_BIN_WRITE(v_stream,z);
      }

      {
        int32_t x = atom.PBC_count.x;
        YAATK_BIN_WRITE(pbc_stream,x);

        int32_t y = atom.PBC_count.y;
        YAATK_BIN_WRITE(pbc_stream,y);

        int32_t z = atom.PBC_count.z;
        YAATK_BIN_WRITE(pbc_stream,z);
      }
    }

    retval = 0;
  }
  catch (...)
  {
    retval = -1;
  }

  return retval;
}

int
SimLoopSaver::load(std::string filenameBase)
{
  int retval = 0;

  size_t atomsCount = 0;

  try
  {
    yaatk::binary_ifstream z_stream(filenameBase + ".z");

    int dataLength = z_stream.getDataLength();

    uint8_t z;
    REQUIRE(dataLength % sizeof(z) == 0);
    atomsCount = dataLength/sizeof(z);

    mdloop.atoms.resize(atomsCount);

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      Atom& atom = mdloop.atoms[i];

      YAATK_BIN_READ(z_stream,z);
      atom.ID = ElementID(z);
    }

    retval |= LOADED_Z;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream x_stream(filenameBase + ".r");

    int dataLength = x_stream.getDataLength();

    double x;
    REQUIRE(dataLength % (sizeof(x)*3) == 0);
    size_t atomsCount_recalc = dataLength/(sizeof(x)*3);

    REQUIRE(atomsCount_recalc == mdloop.atoms.size());

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      Atom& atom = mdloop.atoms[i];

      YAATK_BIN_READ(x_stream,x);
      atom.coords.x = x;

      YAATK_BIN_READ(x_stream,x);
      atom.coords.y = x;

      YAATK_BIN_READ(x_stream,x);
      atom.coords.z = x;
    }

    retval |= LOADED_R;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream v_stream(filenameBase + ".v");

    int dataLength = v_stream.getDataLength();

    double x;
    REQUIRE(dataLength % (sizeof(x)*3) == 0);
    size_t atomsCount_recalc = dataLength/(sizeof(x)*3);

    REQUIRE(atomsCount_recalc == mdloop.atoms.size());

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      Atom& atom = mdloop.atoms[i];

      YAATK_BIN_READ(v_stream,x);
      atom.V.x = x;

      YAATK_BIN_READ(v_stream,x);
      atom.V.y = x;

      YAATK_BIN_READ(v_stream,x);
      atom.V.z = x;
    }

    retval |= LOADED_V;
  }
  catch (...)
  {
  }

  try
  {
    yaatk::binary_ifstream pbc_stream(filenameBase + ".pbc");

    int dataLength = pbc_stream.getDataLength();

    int32_t x;
    REQUIRE(dataLength % (sizeof(x)*3) == 0);
    size_t atomsCount_recalc = dataLength/(sizeof(x)*3);

    REQUIRE(atomsCount_recalc == mdloop.atoms.size());

    for(size_t i = 0; i < mdloop.atoms.size(); ++i)
    {
      Atom& atom = mdloop.atoms[i];

      YAATK_BIN_READ(pbc_stream,x);
      atom.PBC_count.x = x;

      YAATK_BIN_READ(pbc_stream,x);
      atom.PBC_count.y = x;

      YAATK_BIN_READ(pbc_stream,x);
      atom.PBC_count.z = x;
    }

    retval |= LOADED_PBC;
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
    filename.find(".v") != std::string::npos ||
    filename.find(".pbc") != std::string::npos;
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

