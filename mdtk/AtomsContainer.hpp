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

class AtomsArray:public std::vector<Atom>
{
  Vector3D arrayPBC;
public:
  void applyPBC();
  void unfoldPBC();

  void PBC(Vector3D newPBC);
  Vector3D PBC() const;
  bool PBCEnabled() const;

  bool checkMIC(Float Rc) const; // check Minimum Image Criteria
  bool fitInPBC() const;

  void prepareForSimulatation();
  void setAttributesByElementID();

  AtomsArray(size_t size = 0);
  AtomsArray(const AtomsArray &c);
  AtomsArray& operator =(const AtomsArray &c);
  void addAtoms(const AtomsArray &ac);
  virtual ~AtomsArray();
//  AtomRefsContainer refs() {return AtomRefsContainer(*this);};

  void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode);
  void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode);

  Float mass() const;
  Vector3D velocity() const;
  Vector3D massCenter() const;
  Vector3D geomCenter() const;
  Float radius() const;
  Float maxDistanceFrom(Vector3D point) const;

  void removeMomentum();
  void addTranslationalEnergy(Float energy, Vector3D direction);
  void shiftToOrigin();
  void shiftToPosition(Vector3D v);

  struct Dimensions
  {
    Float x_max;
    Float x_min;
    Float x_len;

    Float y_max;
    Float y_min;
    Float y_len;

    Float z_max;
    Float z_min;
    Float z_len;
  };

  Dimensions dimensions() const;

  std::vector<size_t> fixNotFixedAtoms(const size_t begin, const size_t end);
  std::vector<size_t> unfixFixedAtoms(const size_t begin, const size_t end);
  std::vector<size_t> fixUnfixedCHAtoms(const size_t begin, const size_t end);
  void unfixAtoms(const std::vector<size_t> fixedAtoms);
  void fixAtoms(const std::vector<size_t> atomsToFix);

  bool hasTag(unsigned int tagMask) const;
  void tag(unsigned int tagMask);
  void untag(unsigned int tagMask);
  void clearTags();
};

class AtomRefsContainer:public std::vector<Atom*>
{
public:
  AtomRefsContainer();
  AtomRefsContainer(const AtomRefsContainer &c);
  AtomRefsContainer(AtomsArray &c);
  AtomsArray genAtomsArray();
  AtomRefsContainer& operator =(const AtomRefsContainer &c);
  void addAtoms(const AtomRefsContainer &ac);
  virtual ~AtomRefsContainer();

  void saveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode);
  void loadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode);
};

}  // namespace mdtk


#endif


