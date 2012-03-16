/*
   The AtomsContainer class (header file).

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

#ifndef mdtk_AtomContainer_hpp
#define mdtk_AtomContainer_hpp

#include <mdtk/Atom.hpp>

namespace mdtk
{

class AtomsContainer:public std::vector<Atom*>
{
//  Vector3D PBC;
public:
  void setPBC(Vector3D newPBC)
    {
      for (size_t i = 0; i < size(); i++)
      {
        Atom& a = *(at(i));

        if (newPBC.x == NO_PBC.x)
        {
          a.coords.x += a.PBC.x*a.PBC_count.x;
          a.PBC_count.x = 0;
        }
        if (newPBC.y == NO_PBC.y)
        {
          a.coords.y += a.PBC.y*a.PBC_count.y;
          a.PBC_count.y = 0;
        }
        if (newPBC.z == NO_PBC.z)
        {
          a.coords.z += a.PBC.z*a.PBC_count.z;
          a.PBC_count.z = 0;
        }
        a.PBC = newPBC;
      }
    }
//  Vector3D getPBC()const {return PBC;}
//  bool usePBC()const{return PBC != NO_PBC;};
  AtomsContainer():std::vector<Atom*>(){}
    void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
    {
//      YAATK_FSTREAM_WRITE(os,PBC,smode);
    }  
    void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
    {
//      YAATK_FSTREAM_READ(is,PBC,smode);
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
        if (a.PBC.x != NO_PBC.x)
        {
          a.coords.x += a.PBC.x*a.PBC_count.x;
          a.PBC_count.x = 0;
        }
        if (a.PBC.y != NO_PBC.y)
        {
          a.coords.y += a.PBC.y*a.PBC_count.y;
          a.PBC_count.y = 0;
        }
        if (a.PBC.z != NO_PBC.z)
        {
          a.coords.z += a.PBC.z*a.PBC_count.z;
          a.PBC_count.z = 0;
        }
      }
    }
};

}  // namespace mdtk


#endif


