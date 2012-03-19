/*
   The Atom class.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2012 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include <mdtk/Atom.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>

namespace mdtk
{
using namespace std;

Atom::Atom(ElementID id, Vector3D Cx, Vector3D Vx)
  :ID(id),
   Z(1.0*e),
   M(1.0*amu),
   coords(Cx),
   PBC_count(),
   V(Vx),
   an(0.0,0.0,0.0),
   an_no_tb(0.0,0.0,0.0),
   apply_ThermalBath(true),
   globalIndex(0),
   fixed(false),
   PBC(NO_PBC),
   tagbits(0)
{
  setAttributesByElementID();
}

Atom::Atom(const Atom &C)
{
  ID = C.ID;
  Z = C.Z;
  M = C.M;
  coords = C.coords;
  PBC_count = C.PBC_count;
  V = C.V;
  an = C.an;
  an_no_tb = C.an_no_tb;
  apply_ThermalBath = C.apply_ThermalBath;
  globalIndex = C.globalIndex;
  fixed = C.fixed;
  PBC = C.PBC;
  tagbits = C.tagbits;
}

Atom&
Atom::operator=(const Atom &C)
{
  if (this == &C) return *this;

  ID = C.ID;
  Z = C.Z;
  M = C.M;
  coords = C.coords;
  PBC_count = C.PBC_count;
  V = C.V;
  an = C.an;
  an_no_tb = C.an_no_tb;
  apply_ThermalBath = C.apply_ThermalBath;
  globalIndex = C.globalIndex;
  fixed = C.fixed;
  PBC = C.PBC;
  tagbits = C.tagbits;

  return *this;
}

void
Atom::applyPBC()
{
  Vector3D& c = coords;

  if (PBC.x != NO_PBC.x)
  {
    while(c.x < 0)      {c.x += PBC.x;--PBC_count.x;}
    while(c.x >= PBC.x) {c.x -= PBC.x;++PBC_count.x;}
  }
  if (PBC.y != NO_PBC.y)
  {
    while(c.y < 0)      {c.y += PBC.y;--PBC_count.y;}
    while(c.y >= PBC.y) {c.y -= PBC.y;++PBC_count.y;}
  }
  if (PBC.z != NO_PBC.z)
  {
    while(c.z < 0)      {c.z += PBC.z;--PBC_count.z;}
    while(c.z >= PBC.z) {c.z -= PBC.z;++PBC_count.z;}
  }
}

void
Atom::unfoldPBC()
{
  Vector3D& c = coords;

  if (PBC.x != NO_PBC.x)
  {
    c.x += PBC.x*PBC_count.x;
    PBC_count.x = 0;
  }
  if (PBC.y != NO_PBC.y)
  {
    coords.y += PBC.y*PBC_count.y;
    PBC_count.y = 0;
  }
  if (PBC.z != NO_PBC.z)
  {
    coords.z += PBC.z*PBC_count.z;
    PBC_count.z = 0;
  }
}

istream&
operator>>(istream& is, Atom& a)
{
  /*ElementID*/ int ID;
  Float Z;
  Float M;

  Vector3D coords;
  IntVector3D PBC_count;

  Vector3D V;

  Vector3D an;
  Vector3D an_no_tb;

  bool apply_ThermalBath;

  size_t globalIndex;
  bool fixed;

  Vector3D PBC;

  unsigned int tagbits;

  is >> ID
     >> Z
     >> M
     >> coords
     >> PBC_count
     >> V
     >> an
     >> an_no_tb
     >> apply_ThermalBath
     >> globalIndex
     >> fixed
     >> PBC
     >> tagbits;

  if (is)
  {
    a.ID = ElementID(ID);
    a.Z = Z;
    a.M = M;
    a.coords = coords;
    a.PBC_count = PBC_count;
    a.V = V;
    a.an = an;
    a.an_no_tb = an_no_tb;
    a.apply_ThermalBath = apply_ThermalBath;
    a.globalIndex = globalIndex;
    a.fixed = fixed;
    a.PBC = PBC;
    a.tagbits = tagbits;
  }
  else
  {
    cerr << " Error in reading Atom " << endl << flush;
  }
  return is;
}

ostream&
operator<<(ostream& os, const Atom& a)
{
  os << a.ID << "\n"
     << a.Z << "\n"
     << a.M << "\n"
     << a.coords << "\n"
     << a.PBC_count << "\n"
     << a.V << "\n"
     << a.an << "\n"
     << a.an_no_tb << "\n"
     << a.apply_ThermalBath << "\n"
     << a.globalIndex << "\n"
     << a.fixed << "\n"
     << a.PBC << "\n"
     << a.tagbits << "\n";

  return os;
}

}


