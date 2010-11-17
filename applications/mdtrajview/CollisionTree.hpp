/*
   Collision Tree class (header file)

   Copyright (C) 2010 Oleksandr Yermolenko
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

#ifndef	mde_CollisionTree_hpp
#define	mde_CollisionTree_hpp

#include "mdtk/Atom.hpp"
#include <map>
#include <vector>

namespace xmde
{

using namespace mdtk;

typedef std::map<Float,std::vector<Atom> > MDTrajectory;

void MDTrajectory_read(MDTrajectory& mdt,
		       const std::string basefile, 
		       const std::vector<std::string>& xvas);

class CollisionTree
{
public:
  Atom a;
  Float t;
  CollisionTree *t1,*t2;

  CollisionTree(const Atom& atom, 
		MDTrajectory::const_iterator time, 
		const MDTrajectory& mdt);
};

}

#endif
