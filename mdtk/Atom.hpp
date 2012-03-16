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

  Vector3D V;

  Vector3D an;
  Vector3D an_no_tb;

  bool apply_PBC;
  bool apply_ThermalBath;

  size_t globalIndex;

  bool fixed;
  bool isFixed() const { return fixed; }
  void fix() { fixed = true; V=0.0; an=0.0; an_no_tb=0.0; }
  void unfix() { fixed = false; V=0.0; an=0.0; an_no_tb=0.0; }

  Vector3D PBC;
  Vector3D getPBC() const {return PBC;}
  bool usePBC()const{return PBC != NO_PBC;};

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

  unsigned int tag;

#define ATOMTAG_EXAMPLE 1<<0
  bool isTaggedExample() const { return (tag & ATOMTAG_EXAMPLE); }
  void tagExample() { tag |= ATOMTAG_EXAMPLE; }
  void untagExample() { tag &= ~((unsigned int)ATOMTAG_EXAMPLE); }
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


