#include "AtomGroup.hpp"
#include "Molecule.hpp"

namespace mdepp
{

AtomGroup::AtomGroup()
 :atoms()
{
}

AtomGroup::~AtomGroup()
{
}

AtomGroup::AtomGroup(const AtomGroup &c)
 :atoms()
{
  atoms = c.atoms;
}

AtomGroup& 
AtomGroup::operator =(const AtomGroup &c) 
{
  if (this == &c) return *this;
  atoms = c.atoms;
  return *this;
}

void
AtomGroup::put(std::ostream& os) const
{
  os << atoms.size() << "\n";
  for(size_t i = 0; i < atoms.size(); i++)
    os << atoms[i] << "\n";
}

void
AtomGroup::get(std::istream& is)
{
  size_t sz;
  is >> sz;
  atoms.resize(sz);
  for(size_t i = 0; i < atoms.size(); i++)
    is >> atoms[i];
}

std::ostream&
operator<<(std::ostream& os, const AtomGroup& v)
{
  v.put(os);
  return os;
}

std::istream&
operator>>(std::istream& is, AtomGroup& v)
{
  v.get(is);
  return is;
}

bool
AtomGroup::hasAtom(const mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
      found = true;
  return found;
}  

void
AtomGroup::addAtom(const mdtk::Atom& a)
{
  atoms.push_back(a);
//  atoms.back().container = NULL;
}

void
AtomGroup::update(const mdtk::SimLoop& ml)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    mdtk::Atom& a = atoms[i];
    REQUIRE(ac.size() > a.globalIndex);
    const mdtk::Atom& aUpd = ac[a.globalIndex];
    a.V = aUpd.V;
    a.coords = aUpd.coords;
  }
}  

void
AtomGroup::build(const mdtk::SimLoop& ml, 
		 const double SPOTTED_DISTANCE)
{
  const AtomsArray &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = ac[i];
    if (a.coords.z < SPOTTED_DISTANCE)
      addAtom(a);
  }
}  

Molecule
AtomGroup::molecule(size_t atomIndex) const
{
  REQUIRE(atomIndex < atoms.size());
  Molecule m;
  m.buildFromAtom(atoms[atomIndex],*this);
  return m;
}

Molecule
AtomGroup::molecule(const mdtk::Atom& a) const
{
  REQUIRE(hasAtom(a));
  return molecule(a.globalIndex);
}

bool
AtomGroup::isMolecule() const
{
  return maxMolecule().size() == atoms.size();
}

Molecule
AtomGroup::maxMolecule() const
{
  Molecule mMax;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    Molecule m = molecule(i);
    if (m.size() > mMax.size())
      mMax = m;
  }
  return mMax;
}

Float
AtomGroup::mass() const
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    moleculeMass += atom.M;
  } 
  return moleculeMass;
}

mdtk::Vector3D
AtomGroup::velocity() const
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.V*atom.M;
  };
  return sumOfP/sumOfM;    
}  

mdtk::Vector3D
AtomGroup::massCenter() const
{
  REQUIRE(atoms.size() > 0);
  mdtk::Vector3D sumOfP = 0.0;
  Float sumOfM = 0.0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    sumOfM += atom.M;
    sumOfP += atom.coords*atom.M;
  };
  return sumOfP/sumOfM;    
}  

}
