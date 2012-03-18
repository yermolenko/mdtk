/*
   The AtomsContainer class.

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

#include <mdtk/AtomsContainer.hpp>

namespace mdtk
{

void
AtomsContainer::applyPBC()
{
  for(size_t i = 0; i < size(); i++)
    at(i)->applyPBC();
}

void
AtomsContainer::unfoldPBC()
{
  for (size_t i = 0; i < size(); i++)
    at(i)->unfoldPBC();
}

void
AtomsContainer::PBC(Vector3D newPBC)
{
  for (size_t i = 0; i < size(); i++)
  {
    Atom& a = *(at(i));

    a.unfoldPBC();
    a.PBC = newPBC;
    a.applyPBC();
  }
}

Vector3D
AtomsContainer::PBC() const
{
  REQUIRE(size() > 0);
  return front()->PBC;
}

bool
AtomsContainer::PBCEnabled() const
{
  return PBC() != NO_PBC;
}

bool
AtomsContainer::checkMIC(Float Rc) const
{
  if (PBC().x <= Rc)
    return false;
  if (PBC().y <= Rc)
    return false;
  if (PBC().z <= Rc)
    return false;
  return true;
}

bool
AtomsContainer::fitInPBC() const
{
  for(size_t i = 0; i < size(); i++)
  {
    Vector3D PBC(at(i)->PBC);
    Vector3D aci = at(i)->coords;
    if (PBC.x != NO_PBC.x)
      if (aci.x < 0 || aci.x>=PBC.x)
      {
        TRACE(i);TRACE(aci.x);TRACE(PBC.x);return false;
      }
    if (PBC.y != NO_PBC.y)
      if (aci.y < 0 || aci.y>=PBC.y)
      {
        TRACE(i);TRACE(aci.y);TRACE(PBC.y);return false;
      }
    if (PBC.z != NO_PBC.z)
      if (aci.z < 0 || aci.z>=PBC.z)
      {
        TRACE(i);TRACE(aci.z);TRACE(PBC.z);return false;
      }
  }
  return true;
}

void
AtomsContainer::prepareForSimulatation()
{
  REQUIRE(size() > 0);
  for(size_t i = 0; i < size(); i++)
  {
    at(i)->applyPBC();
    at(i)->globalIndex = i;
  }

  if (!fitInPBC())
  {
    cerr << "Atoms do not fit given PBC cell !" << endl << flush;
    throw Exception("Atoms do not fit given PBC cell !");
  }
}

void
AtomsContainer::setAttributesByElementID()
{
  for(size_t i = 0; i < size(); i++)
    at(i)->setAttributesByElementID();
}

AtomsContainer::AtomsContainer()
  :std::vector<Atom*>(),
   createdAtoms()
{
}

AtomsContainer::AtomsContainer(const AtomsContainer &c)
  :std::vector<Atom*>(),
   createdAtoms()
{
  for(size_t i = 0; i < c.size(); i++)
  {
    Atom& a = *(c[i]);
    push_back(createAtom(a));
  }
}

AtomsContainer&
AtomsContainer::operator =(const AtomsContainer &c)
{
  if (this == &c) return *this;

  clear();

  for(size_t i = 0; i < c.size(); i++)
  {
    Atom& a = *(c[i]);
    push_back(createAtom(a));
  }

  return *this;
}

void
AtomsContainer::addAtoms(const AtomsContainer &ac)
{
  for(size_t i = 0; i < ac.size(); i++)
  {
    Atom& a = *(ac[i]);
    push_back(createAtom(a));
  }
}

AtomsContainer::~AtomsContainer()
{
  for(size_t i = 0; i < createdAtoms.size(); i++)
    delete createdAtoms[i];
}

void
AtomsContainer::saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
{
//      YAATK_FSTREAM_WRITE(os,PBC,smode);
}

void
AtomsContainer::loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
{
//      YAATK_FSTREAM_READ(is,PBC,smode);
}

void
AtomsContainer::normalize()
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

}


