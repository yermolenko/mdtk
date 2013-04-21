/*
   SimLoopSaver class header file.

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

#ifndef	mdtk_SimLoopSaver_hpp
#define	mdtk_SimLoopSaver_hpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>

#include <mdtk/Vector3D.hpp>
#include <mdtk/SimLoop.hpp>

namespace mdtk
{

class SimLoopSaver
{
  SimLoop& mdloop;
  const std::string idPrefix;
  std::string generateId(unsigned long iteration);
  std::string extractAttributeName(const std::string filename);
  void prepareForAttributeReading(yaatk::binary_ifstream& is, size_t attributeSize);
public:
  SimLoopSaver(SimLoop& mdloopInstance);
  virtual ~SimLoopSaver() {}

  int write(std::string id);
  int write();

  enum {LOADED_Z = (1<<0)};
  enum {LOADED_R = (1<<1)};
  enum {LOADED_V = (1<<2)};
  enum {LOADED_ZRV = LOADED_Z | LOADED_R | LOADED_V};

  enum {LOADED_PBC_RECT = (1<<3)};
  enum {LOADED_PBC_COUNT = (1<<4)};
  enum {LOADED_PBC_ENABLED = (1<<5)};
  enum {LOADED_PBC = LOADED_PBC_RECT | LOADED_PBC_COUNT | LOADED_PBC_ENABLED};

  enum {LOADED_ZRV_PBC = LOADED_ZRV | LOADED_PBC};

  enum {LOADED_A = (1<<6)};

  int load(std::string id);
  int loadIteration(unsigned long iteration);
  int loadIterationLatest();

  static bool mayContainData(std::string filename);
  static std::vector<std::string> listIds();
  std::vector<unsigned long> listIterations();

  void removeIterations(const std::vector<unsigned long>&);
  void removeIterations(bool keepFirst = true, bool keepLast = true);

  void removeAttributes(const std::string id, const std::set<std::string>& protectedAttributes);
  void removeAttributesButPos(const std::string id);
  void removeAttributesButPos() {removeAttributesButPos(generateId(mdloop.iteration));}
  void removeAttributesButPosVel(const std::string id);
  void removeAttributesButPosVel() {removeAttributesButPosVel(generateId(mdloop.iteration));}
  void removeAttributesButPosVelPBC(const std::string id);
  void removeAttributesButPosVelPBC() {removeAttributesButPosVelPBC(generateId(mdloop.iteration));}
};

}

#endif

