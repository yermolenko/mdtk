#ifndef mdpp_AtomGroup_hpp
#define mdpp_AtomGroup_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>
#include <vector>

namespace mdepp
{
  using namespace mdtk;

  class Molecule;

class AtomGroup
{
public:
  std::vector<mdtk::Atom> atoms;

  AtomGroup();
  virtual ~AtomGroup();
  AtomGroup(const AtomGroup &c);

  AtomGroup& operator =(const AtomGroup &c);

  virtual void get(std::istream& is);
  virtual void put(std::ostream& os) const;

  void addAtom(const mdtk::Atom& a);
  bool hasAtom(const mdtk::Atom&) const;

  size_t size() const {return atoms.size();};

  void build(const mdtk::SimLoop& ml, 
	     const double SPOTTED_DISTANCE = +1000.0*mdtk::Ao);
  void update(const mdtk::SimLoop& ml);

  Molecule molecule(size_t atomIndex) const;
  Molecule molecule(const mdtk::Atom& a) const;
  Molecule maxMolecule() const;
  bool isMolecule() const;
  mdtk::Vector3D massCenter() const;

  Float          mass() const;
  mdtk::Vector3D velocity() const;

  friend std::istream&  operator>> (std::istream& is, AtomGroup& vec);
  friend std::ostream&  operator<< (std::ostream& os, const AtomGroup& vec);
};

std::istream& operator>> (std::istream& is, AtomGroup& vec);
std::ostream& operator<< (std::ostream& os, const AtomGroup& vec);

}

#endif
