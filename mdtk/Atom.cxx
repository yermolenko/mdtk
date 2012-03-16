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
   apply_PBC(true),
   apply_ThermalBath(true),
   globalIndex(0),
   fixed(false),
   PBC(NO_PBC),
   tag(0)
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
  apply_PBC = C.apply_PBC;
  apply_ThermalBath = C.apply_ThermalBath;
  globalIndex = C.globalIndex;
  fixed = C.fixed;
  PBC = C.PBC;
  tag = C.tag;
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
  apply_PBC = C.apply_PBC;
  apply_ThermalBath = C.apply_ThermalBath;
  globalIndex = C.globalIndex;
  fixed = C.fixed;
  PBC = C.PBC;
  tag = C.tag;

  return *this;
}

Vector3D
depos(const Atom &a1, const Atom &a2)
{
  Vector3D r(a1.coords - a2.coords);
//if (a1.globalIndex == 1806) {TRACE("before getPBC");TRACE(a1.container);}
  if (a1.PBC == NULL) return r;
  Vector3D PBC = a1.PBC;
//if (a1.globalIndex == 1806) {TRACE("after getPBC");TRACE(a1.container);}
  if (PBC == NO_PBC) return r;

  if (!(a1.apply_PBC && a2.apply_PBC)) return r;

  if(fabs(r.x) > PBC.x*0.5) {r.x += (r.x > 0)?(-PBC.x):(PBC.x);}
  if(fabs(r.y) > PBC.y*0.5) {r.y += (r.y > 0)?(-PBC.y):(PBC.y);}
  if(fabs(r.z) > PBC.z*0.5) {r.z += (r.z > 0)?(-PBC.z):(PBC.z);}

  return r;
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

  bool apply_PBC;
  bool apply_ThermalBath;

  size_t globalIndex;
  bool fixed;

  Vector3D PBC;

  unsigned int tag;

  is >> ID
     >> Z
     >> M
     >> coords
     >> PBC_count
     >> V
     >> an
     >> an_no_tb
     >> apply_PBC
     >> apply_ThermalBath
     >> globalIndex
     >> fixed
     >> PBC
     >> tag;

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
    a.apply_PBC = apply_PBC;
    a.apply_ThermalBath = apply_ThermalBath;
    a.globalIndex = globalIndex;
    a.fixed = fixed;
    a.PBC = PBC;
    a.tag = tag;
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
     << a.apply_PBC << "\n"
     << a.apply_ThermalBath << "\n"
     << a.globalIndex << "\n"
     << a.fixed << "\n"
     << a.PBC << "\n"
     << a.tag << "\n";

  return os;
}

}


