/*
   The Cluster class (header file).

   Copyright (C) 2010, 2012 Oleksandr Yermolenko
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

#ifndef mdpp_Cluster_hpp
#define mdpp_Cluster_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>
#include <vector>
#include "AtomGroup.hpp"

namespace mdepp
{
  using namespace mdtk;

class Cluster : public AtomGroup
{
public:
  Cluster();
  ~Cluster();
  Cluster(const Cluster &c);

  Cluster& operator =(const Cluster &c);

  virtual void get(std::istream& is);
  virtual void put(std::ostream& os) const;

  void build(const mdtk::SimLoop& ml);
};

}

#endif
