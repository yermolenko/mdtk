#include "AtomGroup.hpp"

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
  atoms.back().container = NULL;
}

void
AtomGroup::update(const mdtk::SimLoop& ml)
{
  const AtomsContainer &ac = ml.atoms;
  for(size_t i = 0; i < atoms.size(); i++)
  {
    mdtk::Atom& a = atoms[i];
    REQUIRE(ac.size() > a.globalIndex);
    const mdtk::Atom& aUpd = *ac[a.globalIndex];
    a.V = aUpd.V;
    a.coords = aUpd.coords;
  }
}  

Float getMassInAMU(const std::vector<mdtk::Atom>& atoms)
{
  Float moleculeMass = 0;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    moleculeMass += atom.M;
  } 
  return mdtk::academic_round(moleculeMass/mdtk::amu);
}

mdtk::Vector3D getVelocity(const std::vector<mdtk::Atom>& atoms)
{
  using mdtk::Exception;    
    
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

void printGlobalIndexes(const std::vector<mdtk::Atom>& atoms,
			std::ostream& fo)
{
  fo << atoms.size() << std::endl;
  for(size_t ai = 0; ai < atoms.size(); ai++)
  {
    const mdtk::Atom& atom = atoms[ai];
    fo << atom.globalIndex << std::endl;
  } 
} 


}
