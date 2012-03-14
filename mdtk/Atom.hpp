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

typedef unsigned int Tag;

class Atom;

typedef std::pair<int,int> AtomPair;

enum ElementID 
{H_EL = 1, C_EL = 12, Cu_EL = 64, Ag_EL = 108, Au_EL = 197, Ar_EL = 40, Xe_EL = 131, DUMMY_EL = -1};

inline
std::string
ElementIDtoString(ElementID id)
{
  std::string str;
  switch (id)
  {
  case H_EL  : str = "H"; break;
  case C_EL  : str = "C"; break;
  case Cu_EL : str = "Cu"; break;
  case Ag_EL : str = "Ag"; break;
  case Au_EL : str = "Au"; break;
  case Ar_EL : str = "Ar"; break;
  case Xe_EL : str = "Xe"; break;
  case DUMMY_EL :
  default : Exception("Unknown element");
  }
  return str;
}

#define INFINITE_MASS 1.0e100*mdtk::amu

const int EL_ID_size = 200;

class FGeneral;
class FManybody;
class AtomsContainer;

class Atom
{
public:
  ElementID ID;

  Float Z;
  Float M;
  Vector3D V;
  Vector3D coords;
  IntVector3D PBC_count;

  Vector3D an;
  Vector3D an_no_tb;

  bool apply_barrier;
  bool apply_PBC;
  bool apply_ThermalBath;
  bool ejected;
  
  Tag tag;

  size_t globalIndex;
  
  bool fixed; // dummy, should be removed!!!

  Atom(Float Zx=1,Float Mx=1,
    Vector3D Vx=Vector3D(0,0,0),
    Vector3D Cx=Vector3D(0,0,0));

  Atom(ElementID id,
    Vector3D Cx=Vector3D(0,0,0), Vector3D Vx=Vector3D(0,0,0));

  Atom( const Atom &C );
  Atom& operator =(const Atom &C);

  friend Vector3D    depos(const Atom &a1, const Atom &a2);

  friend std::istream&  operator>> (std::istream& is, Atom& vec);
  friend std::ostream&  operator<< (std::ostream& os, const Atom& vec);

  friend bool     operator==(const Atom& v1, const Atom& v2);
  friend bool     operator!=(const Atom& v1, const Atom& v2);
  friend bool     operator>(const Atom& v1, const Atom& v2);
  friend bool     operator<(const Atom& v1, const Atom& v2);
  friend bool     operator>=(const Atom& v1, const Atom& v2);
  friend bool     operator<=(const Atom& v1, const Atom& v2);


  AtomsContainer* container;
  void setAttributesByElementID();

#define ATOMTAG_FIXED 1<<0
  bool isFixed() const
    {
      return (tag & ATOMTAG_FIXED);
    }
  void fix()
    {
      tag |= ATOMTAG_FIXED;
      M = INFINITE_MASS;V=0.0;an=0.0;an_no_tb=0.0;
    }
  void unfix()
    {
      tag &= ~((unsigned int)ATOMTAG_FIXED);
      setAttributesByElementID();
      V=0.0;an=0.0;an_no_tb=0.0;
    }
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
{
  return a1.globalIndex >= a2.globalIndex;
}  
  
}  

const Vector3D NO_PBC = Vector3D(10000.0*Ao, 10000.0*Ao, 10000.0*Ao);

class AtomsContainer:public std::vector<Atom*>
{
  Vector3D PBC;
public:
  void setPBC(Vector3D newPBC)
    {
      for (size_t i = 0; i < size(); i++)
      {
        Atom& a = *(at(i));

        if (newPBC.x == NO_PBC.x)
        {
          a.coords.x += PBC.x*a.PBC_count.x;
          a.PBC_count.x = 0;
        }
        if (newPBC.y == NO_PBC.y)
        {
          a.coords.y += PBC.y*a.PBC_count.y;
          a.PBC_count.y = 0;
        }
        if (newPBC.z == NO_PBC.z)
        {
          a.coords.z += PBC.z*a.PBC_count.z;
          a.PBC_count.z = 0;
        }
      }
      PBC = newPBC;
    }
  Vector3D getPBC()const {return PBC;}
  bool usePBC()const{return PBC != NO_PBC;};
  AtomsContainer():std::vector<Atom*>(),PBC(NO_PBC){}
    void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_WRITE(os,PBC,smode);
    }  
    void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
    {
      YAATK_FSTREAM_READ(is,PBC,smode);
    }  
  void setAttributesByElementID()
  {
    size_t i;
    for(i = 0; i < size(); i++)
      operator[](i)->setAttributesByElementID();
  }
void 
normalize()
{
  size_t i;

  Float msum = 0.0;
  Vector3D mvsum = 0.0;
  for (i = 0; i < size(); i++)
  {
    Atom& a = *(at(i));
    if (a.isFixed()) {a.V=0.0;a.an_no_tb=0.0;a.an=0.0;continue;}
    mvsum += a.V*a.M;
    msum += a.M;
  }
  for (i = 0; i < size(); i++)
  {
    Atom& a = *(at(i));
    if (a.isFixed()) {continue;}
    a.V -= mvsum/msum;
  }
}
    void unfoldPBC()
    {
      for (size_t i = 0; i < size(); i++)
      {
        Atom& a = *(at(i));
        if (PBC.x != NO_PBC.x)
        {
          a.coords.x += PBC.x*a.PBC_count.x;
          a.PBC_count.x = 0;
        }
        if (PBC.y != NO_PBC.y)
        {
          a.coords.y += PBC.y*a.PBC_count.y;
          a.PBC_count.y = 0;
        }
        if (PBC.z != NO_PBC.z)
        {
          a.coords.z += PBC.z*a.PBC_count.z;
          a.PBC_count.z = 0;
        }
      }
    }
};


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
  switch (atom.ID)
  {
    case C_EL:  return "C";
    case H_EL:  return "H";
    case Cu_EL:  return "Cu";    case Ar_EL:  return "Ar";    default:    return "Unhandled";
  }  
}  

void
saveToMDE(std::ofstream& fo, AtomsContainer& Ro, Vector3D PBC, Float thermalBath_z);

void
saveToBRENNERMD(std::ofstream& fo, AtomsContainer& Ro, Vector3D PBC);

void
loadFromMDE(std::ifstream& fi, AtomsContainer& Ro, Vector3D& PBC, Float &thermalBath_z);



using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::sqrt;
using std::acos;
using std::FILE;


inline
Atom::Atom(Float Zx,Float Mx,
  Vector3D Vx,Vector3D Cx)
  :
   ID(H_EL),
   Z(Zx),M(Mx),V(Vx),coords(Cx),
   PBC_count(),
   an(0.0,0.0,0.0),
   an_no_tb(0.0,0.0,0.0),
   apply_barrier(false),apply_PBC(true),apply_ThermalBath(true),ejected(false),tag(0),
   globalIndex(0)
   ,fixed(false)
   ,container(NULL)
{
}

inline
Atom::Atom(ElementID id,
  Vector3D Cx,Vector3D Vx)
  :
   ID(id),
   Z(1.0*e),M(1.0*amu),V(Vx),coords(Cx),
   PBC_count(),
   an(0.0,0.0,0.0),
   an_no_tb(0.0,0.0,0.0),
   apply_barrier(false),apply_PBC(true),apply_ThermalBath(true),ejected(false),tag(0),
   globalIndex(0)
   ,fixed(false)
   ,container(NULL)
{
  setAttributesByElementID();
}

inline
Atom::Atom( const Atom &C )
{
  ID = C.ID;
  tag = 0;
  Z = C.Z;
  M = C.M;
  V = C.V;
  coords = C.coords;
  PBC_count = C.PBC_count;
  an = C.an;
  an_no_tb = C.an_no_tb;
  apply_barrier = C.apply_barrier;
  apply_PBC = C.apply_PBC;
  apply_ThermalBath = C.apply_ThermalBath;
  ejected = C.ejected;
  
  tag = C.tag;
  globalIndex = C.globalIndex;
  fixed = C.fixed;

  container = NULL;//C.container;
//  container = C.container;
}

inline
Atom&
Atom::operator =(const Atom &C) 
{
  if (this == &C) return *this;
    ID = C.ID;
    Z = C.Z;
    M = C.M;
    V = C.V;
    coords = C.coords;
    PBC_count = C.PBC_count;
    an = C.an;
    an_no_tb = C.an_no_tb;
    apply_barrier = C.apply_barrier;
    apply_PBC = C.apply_PBC;
    apply_ThermalBath = C.apply_ThermalBath;
    ejected = C.ejected;

    tag = C.tag;
    globalIndex = C.globalIndex;

    container = NULL;//C.container;
//    container = C.container;

  return *this;
}

inline
Vector3D
depos(const Atom &a1, const Atom &a2)
{
  Vector3D r(a1.coords - a2.coords);
//if (a1.globalIndex == 1806) {TRACE("before getPBC");TRACE(a1.container);}
  if (a1.container == NULL) return r;
  Vector3D PBC = a1.container->getPBC();
//if (a1.globalIndex == 1806) {TRACE("after getPBC");TRACE(a1.container);}
  if (PBC == NO_PBC) return r;

  if (!(a1.apply_PBC && a2.apply_PBC)) return r;

  if(fabs(r.x) > PBC.x*0.5) {r.x += (r.x > 0)?(-PBC.x):(PBC.x);}
  if(fabs(r.y) > PBC.y*0.5) {r.y += (r.y > 0)?(-PBC.y):(PBC.y);}
  if(fabs(r.z) > PBC.z*0.5) {r.z += (r.z > 0)?(-PBC.z):(PBC.z);}

  return r;
}


inline
istream&
operator>>(istream& is,  Atom& vec)
{
  /*ElementID*/ int ID;
  Float Z = 0, M = 0;
  int t;
  Vector3D V,C,an, an_no_tb;
  IntVector3D PBC_count;
  bool apply_barrier;
  bool apply_PBC;
  bool apply_ThermalBath;
  bool ejected;
  size_t gi;
  bool fixed;
  is >> ID >> Z >> M >> V >> C >> PBC_count >> an >> an_no_tb >> apply_barrier >> apply_PBC >> apply_ThermalBath >> ejected >> t >> gi >> fixed;
  if (is)
  {
    vec = Atom(Z,M,V,C);
    vec.PBC_count = PBC_count;
    vec.ID = ElementID(ID);
    vec.an = an;
    vec.an_no_tb = an_no_tb;
    vec.apply_barrier = apply_barrier;
    vec.apply_PBC = apply_PBC;
    vec.apply_ThermalBath = apply_ThermalBath;
    vec.ejected = ejected;
    vec.tag = t;
    vec.globalIndex = gi;
  }
  else
  {
    cerr << " Error in reading Atom " << endl << flush;
  }
  return is;
}

inline
ostream&
operator<<(ostream& os, const Atom& vec)
{
  os << vec.ID << "\n"
     << vec.Z << "\n"
     << vec.M << "\n"
     << vec.V << "\n"
     << vec.coords << "\n"
     << vec.PBC_count << "\n"
     << vec.an << "\n"
     << vec.an_no_tb << "\n"
     << vec.apply_barrier << "\n"
     << vec.apply_PBC << "\n"
     << vec.apply_ThermalBath << "\n"
     << vec.ejected << "\n"
     << vec.tag << "\n"
     << vec.globalIndex << "\n"
     << vec.fixed << "\n"
;
  return os;
}

}  // namespace mdtk


#endif


