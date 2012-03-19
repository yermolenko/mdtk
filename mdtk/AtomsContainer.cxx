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
AtomsArray::applyPBC()
{
  for(size_t i = 0; i < size(); i++)
    at(i).applyPBC();
}

void
AtomsArray::unfoldPBC()
{
  for (size_t i = 0; i < size(); i++)
    at(i).unfoldPBC();
}

void
AtomsArray::PBC(Vector3D newPBC)
{
  arrayPBC = newPBC;
  for (size_t i = 0; i < size(); i++)
  {
    Atom& a = at(i);

    a.unfoldPBC();
    a.PBC = arrayPBC;
    a.applyPBC();
  }
}

Vector3D
AtomsArray::PBC() const
{
  return arrayPBC;
}

bool
AtomsArray::PBCEnabled() const
{
  return PBC() != NO_PBC;
}

bool
AtomsArray::checkMIC(Float Rc) const
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
AtomsArray::fitInPBC() const
{
  for(size_t i = 0; i < size(); i++)
  {
    Vector3D PBC(at(i).PBC);
    Vector3D aci = at(i).coords;
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
AtomsArray::prepareForSimulatation()
{
  REQUIRE(size() > 0);
  for(size_t i = 0; i < size(); i++)
  {
    at(i).PBC = arrayPBC;
    at(i).applyPBC();
    at(i).globalIndex = i;
  }

  if (!fitInPBC())
  {
    cerr << "Atoms do not fit given PBC cell !" << endl << flush;
    throw Exception("Atoms do not fit given PBC cell !");
  }
}

void
AtomsArray::setAttributesByElementID()
{
  for(size_t i = 0; i < size(); i++)
    at(i).setAttributesByElementID();
}

AtomsArray::AtomsArray(size_t size)
  :std::vector<Atom>(size),
   arrayPBC(NO_PBC)
{
}

AtomsArray::AtomsArray(const AtomsArray &c)
  :std::vector<Atom>(c),
   arrayPBC(c.arrayPBC)
{
}

AtomsArray&
AtomsArray::operator =(const AtomsArray &c)
{
  if (this == &c) return *this;

  std::vector<Atom>::operator =(c);
  arrayPBC = c.arrayPBC;

  return *this;
}

void
AtomsArray::addAtoms(const AtomsArray &ac)
{
  for(size_t i = 0; i < ac.size(); i++)
    push_back(ac[i]);

  PBC(arrayPBC);
}

AtomsArray::~AtomsArray()
{
}

void
AtomsArray::saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
{
  int i,atoms_count = size();
  YAATK_FSTREAM_WRITE(os,atoms_count,smode);
  for(i = 0; i < atoms_count; i++)
    YAATK_FSTREAM_WRITE(os,operator[](i),smode);

  YAATK_FSTREAM_WRITE(os,arrayPBC,smode);
}

void
AtomsArray::loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
{
  int i,atoms_count;
  YAATK_FSTREAM_READ(is,atoms_count,smode);
  cout << "Reading info about " << atoms_count << " atoms..." << endl;
  resize(atoms_count);
  for(i = 0; i < atoms_count; i++)
    YAATK_FSTREAM_READ(is,operator[](i),smode);
  cout << " done." << endl;

  YAATK_FSTREAM_READ(is,arrayPBC,smode);
}

Float
AtomsArray::mass() const
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < size(); ai++)
  {
    const mdtk::Atom& atom = at(ai);
    moleculeMass += atom.M;
  }
  return moleculeMass;
}

Vector3D
AtomsArray::velocity() const
{
  REQUIRE(size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < size(); ai++)
  {
    const mdtk::Atom& atom = at(ai);
    if (atom.isFixed()) continue;
    sumOfM += atom.M;
    sumOfP += atom.V*atom.M;
  };
  return sumOfP/sumOfM;
}

Vector3D
AtomsArray::massCenter() const
{
  REQUIRE(size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < size(); ai++)
  {
    const mdtk::Atom& atom = at(ai);
    sumOfM += atom.M;
    sumOfP += atom.coords*atom.M;
  };
  return sumOfP/sumOfM;
}

Vector3D
AtomsArray::geomCenter() const
{
  Float clusterRadius = 0.0;
  Vector3D clusterCenter(0,0,0);

  for(size_t i = 0; i < size(); i++)
    clusterCenter += at(i).coords;
  clusterCenter /= size();

  return clusterCenter;
}

Float
AtomsArray::maxDistanceFrom(Vector3D point) const
{
  Float clusterRadius = 0.0;

  for(size_t i = 0; i < size(); i++)
  {
    Float currentDist = (at(i).coords-point).module();
    clusterRadius = (currentDist>clusterRadius)?currentDist:clusterRadius;
  }

  return clusterRadius;
}

Float
AtomsArray::radius() const
{
  return maxDistanceFrom(geomCenter());
}

void
AtomsArray::removeMomentum()
{
  Vector3D v = velocity();
  for(size_t ai = 0; ai < size(); ai++)
  {
    mdtk::Atom& atom = at(ai);
    if (atom.isFixed()) continue;
    atom.V -= v;
  };
}

void
AtomsArray::addTranslationalEnergy(Float energy, Vector3D direction)
{
  direction.normalize();
  Vector3D v = velocity();
  for(size_t ai = 0; ai < size(); ai++)
  {
    mdtk::Atom& a = at(ai);
    if (a.isFixed()) continue;
    a.V += sqrt(2.0*energy/(mass()))*direction;
  };
}

void
AtomsArray::shiftToOrigin()
{
  if (size() == 0) return;

  Vector3D clusterCenter = massCenter();

  for(size_t i = 0; i < size(); i++)
    at(i).coords -= clusterCenter;
}

void
AtomsArray::shiftToPosition(Vector3D v)
{
  shiftToOrigin();
  for(size_t i = 0; i < size(); i++)
    at(i).coords += v;
}

AtomsArray::Dimensions
AtomsArray::dimensions() const
{
  Dimensions d;

  const Atom& a = at(0);

  Float x_max = a.coords.x;
  Float x_min = a.coords.x;
  Float y_max = a.coords.y;
  Float y_min = a.coords.y;
  Float z_max = a.coords.z;
  Float z_min = a.coords.z;

  for(size_t i = 0; i < size(); i++)
  {
    const Atom& a = at(i);

    if (a.coords.x > d.x_max)
      d.x_max = a.coords.x;
    if (a.coords.x < d.x_min)
      d.x_min = a.coords.x;

    if (a.coords.y > d.y_max)
      d.y_max = a.coords.y;
    if (a.coords.y < d.y_min)
      d.y_min = a.coords.y;

    if (a.coords.z > d.z_max)
      d.z_max = a.coords.z;
    if (a.coords.z < d.z_min)
      d.z_min = a.coords.z;
  }

  d.x_len = d.x_max - d.x_min;
  d.y_len = d.y_max - d.y_min;
  d.z_len = d.z_max - d.z_min;

  return d;
}

std::vector<size_t>
AtomsArray::fixNotFixedAtoms(const size_t begin, const size_t end)
{
  std::vector<size_t> fixated;
  for(size_t i = 0; i < end; i++)
    if (!at(i).isFixed())
    {
      at(i).fix();
      fixated.push_back(i);
    }
  return fixated;
}

std::vector<size_t>
AtomsArray::unfixFixedAtoms(const size_t begin, const size_t end)
{
  std::vector<size_t> unfixated;
  for(size_t i = 0; i < end; i++)
    if (at(i).isFixed())
    {
      at(i).unfix();
      unfixated.push_back(i);
    }
  return unfixated;
}

std::vector<size_t>
AtomsArray::fixUnfixedCHAtoms(const size_t begin, const size_t end)
{
  std::vector<size_t> fixated;
  for(size_t i = 0; i < end; i++)
    if (!at(i).isFixed())
      if (at(i).ID == C_EL || at(i).ID == H_EL)
      {
        at(i).fix();
        fixated.push_back(i);
      }
  return fixated;
}

void
AtomsArray::unfixAtoms(const std::vector<size_t> fixedAtoms)
{
  for(size_t i = 0; i < fixedAtoms.size(); i++)
    at(fixedAtoms[i]).unfix();
}

void
AtomsArray::fixAtoms(const std::vector<size_t> atomsToFix)
{
  for(size_t i = 0; i < atomsToFix.size(); i++)
    at(atomsToFix[i]).fix();
}

bool
AtomsArray::hasTag(unsigned int tagMask) const
{
  for(size_t i = 0; i < size(); i++)
    if (at(i).hasTag(tagMask))
      return true;
  return false;
}

void
AtomsArray::tag(unsigned int tagMask)
{
  for(size_t i = 0; i < size(); i++)
    at(i).tag(tagMask);
}

void
AtomsArray::untag(unsigned int tagMask)
{
  for(size_t i = 0; i < size(); i++)
    at(i).untag(tagMask);
}

void
AtomsArray::clearTags()
{
  for(size_t i = 0; i < size(); i++)
    at(i).clearTags();
}

AtomRefsContainer::AtomRefsContainer()
  :std::vector<Atom*>()
{
}

AtomRefsContainer::AtomRefsContainer(const AtomRefsContainer &c)
  :std::vector<Atom*>(c)
{
}

AtomRefsContainer::AtomRefsContainer(AtomsArray& c)
  :std::vector<Atom*>()
{
  for(size_t i = 0; i < c.size(); i++)
    push_back(&(c[i]));
}

AtomsArray
AtomRefsContainer::genAtomsArray()
{
  AtomsArray ar;
  for(size_t i = 0; i < size(); i++)
    ar.push_back(*at(i));
  return ar;
}

AtomRefsContainer&
AtomRefsContainer::operator =(const AtomRefsContainer &c)
{
  if (this == &c) return *this;

  std::vector<Atom*>::operator =(c);

  return *this;
}

void
AtomRefsContainer::addAtoms(const AtomRefsContainer &ac)
{
  for(size_t i = 0; i < ac.size(); i++)
    push_back(ac[i]);
}

AtomRefsContainer::~AtomRefsContainer()
{
}

void
AtomRefsContainer::saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
{
}

void
AtomRefsContainer::loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
{
}

}

