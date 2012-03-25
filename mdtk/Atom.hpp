/*
   The Atom class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012
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

#ifndef mdtk_Atom_hpp
#define mdtk_Atom_hpp

#include <mdtk/Vector3D.hpp>
#include <mdtk/consts.hpp>

namespace mdtk
{

class Atom
{
public:
  ElementID ID;

  Float Z;
  Float M;
  void setAttributesByElementID();

  Vector3D coords;
  IntVector3D PBC_count;
  void applyPBC();
  void unfoldPBC();

  Vector3D V;

  Vector3D an;
  Vector3D an_no_tb;

  bool apply_ThermalBath;

  size_t globalIndex;

  bool fixed;
  bool isFixed() const { return fixed; }
  void fix() { fixed = true; V=0.0; an=0.0; an_no_tb=0.0; }
  void unfix() { fixed = false; V=0.0; an=0.0; an_no_tb=0.0; }

  Vector3D PBC;
  bool PBCEnabled()const{return PBC != NO_PBC;};
  bool lateralPBCEnabled()const{return PBC.x != NO_PBC && PBC.y != NO_PBC;};

  Atom(ElementID id=H_EL, Vector3D Cx=Vector3D(0,0,0), Vector3D Vx=Vector3D(0,0,0));

  Atom(const Atom &C);
  Atom& operator=(const Atom &C);

  friend Vector3D depos(const Atom &a1, const Atom &a2);

  friend std::istream& operator>>(std::istream& is, Atom& vec);
  friend std::ostream& operator<<(std::ostream& os, const Atom& vec);

  friend bool operator==(const Atom& v1, const Atom& v2);
  friend bool operator!=(const Atom& v1, const Atom& v2);
  friend bool operator>(const Atom& v1, const Atom& v2);
  friend bool operator<(const Atom& v1, const Atom& v2);
  friend bool operator>=(const Atom& v1, const Atom& v2);
  friend bool operator<=(const Atom& v1, const Atom& v2);

  unsigned int tagbits;

#define ATOMTAG_EXAMPLE (1<<0)
#define ATOMTAG_TARGET (1<<1)
#define ATOMTAG_PROJECTILE (1<<2)
#define ATOMTAG_SUBSTRATE (1<<3)
#define ATOMTAG_MONOMER (1<<4)
#define ATOMTAG_CLUSTER (1<<5)
#define ATOMTAG_FULLERENE (1<<6)

  bool hasTag(unsigned int tagMask) const { return (tagbits & tagMask); }
  void tag(unsigned int tagMask) { tagbits |= tagMask; }
  void untag(unsigned int tagMask) { tagbits &= ~(tagMask); }
  void clearTags() {tagbits = 0;};
};

inline
void
Atom::setAttributesByElementID()
{
  switch (ID)
  {
    case H_EL  : Z =   1.0*e; M =   1.0*amu; break;
    case C_EL  : Z =   6.0*e; M =  12.0*amu; break;
    case Cu_EL : Z =  29.0*e; M =  63.546*amu; break;
    case Ag_EL : Z =  47.0*e; M =  107.868*amu; break;
    case Au_EL : Z =  79.0*e; M =  196.967*amu; break;
    case Ar_EL : Z =  18.0*e; M =  39.948*amu; break;
    case Xe_EL : Z =  54.0*e; M =  131.293*amu; break;
    case DUMMY_EL : break;
  }
}

inline
Vector3D
depos(const Atom &a1, const Atom &a2)
{
  Vector3D r(a1.coords - a2.coords);

//  if (!(a1.PBCEnabled() && a2.PBCEnabled())) return r;
  REQUIRE(a1.PBC.x == a2.PBC.x && a1.PBC.y == a2.PBC.y && a1.PBC.z == a2.PBC.z);
  REQUIRE(a1.PBC != NO_PBC);
  if (a1.PBC.x == a2.PBC.x && a1.PBC.y == a2.PBC.y && a1.PBC.z == a2.PBC.z)
  {
    if(fabs(r.x) > a1.PBC.x*0.5) {r.x += (r.x > 0)?(-a1.PBC.x):(a1.PBC.x);}
    if(fabs(r.y) > a1.PBC.y*0.5) {r.y += (r.y > 0)?(-a1.PBC.y):(a1.PBC.y);}
    if(fabs(r.z) > a1.PBC.z*0.5) {r.z += (r.z > 0)?(-a1.PBC.z):(a1.PBC.z);}
  }

  return r;
}

inline
bool
operator==(const Atom& a1, const Atom& a2)
{
  return a1.globalIndex == a2.globalIndex;
}

inline
bool
operator!=(const Atom& a1, const Atom& a2)
{
  return !(operator==(a1,a2));
}

inline
bool
operator<(const Atom& a1, const Atom& a2)
{
  return a1.globalIndex < a2.globalIndex;
}

inline
bool
operator>(const Atom& a1, const Atom& a2)
{
  return a1.globalIndex > a2.globalIndex;
}

inline
bool
operator<=(const Atom& a1, const Atom& a2)
{
  return a1.globalIndex <= a2.globalIndex;
}

inline
bool
operator>=(const Atom& a1, const Atom& a2)
{
  return a1.globalIndex >= a2.globalIndex;
}

inline
Float
wDeb(Atom &/*atom*/)
{
  return 5.0e13;
//  return 1.0/(0.1*ps);
//  return 1.0/(0.4*ps);
}

inline
std::string
ElementString(Atom &atom)
{
  return ElementIDtoString(atom.ID);
}

typedef std::pair<int,int> AtomPair;

}  // namespace mdtk

#endif


