#ifndef mdpp_AtomGroup_hpp
#define mdpp_AtomGroup_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>
#include <vector>

namespace mdepp
{
  using namespace mdtk;

  Float          getMassInAMU(const std::vector<mdtk::Atom>& atoms);
  mdtk::Vector3D getVelocity(const std::vector<mdtk::Atom>& atoms);
  void printGlobalIndexes(const std::vector<mdtk::Atom>& atoms,
			  std::ostream& fo);

//#define REF_POT_OF(ML_FPOT) ML_FPOT.potentials[0]

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

  Molecule molecule(size_t atomIndex);
  Molecule molecule(const mdtk::Atom& a);
  Molecule maxMolecule();

  friend std::istream&  operator>> (std::istream& is, AtomGroup& vec);
  friend std::ostream&  operator<< (std::ostream& os, const AtomGroup& vec);
};

std::istream& operator>> (std::istream& is, AtomGroup& vec);
std::ostream& operator<< (std::ostream& os, const AtomGroup& vec);

}

#endif
