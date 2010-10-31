#include "Fullerene.hpp"

namespace mdepp
{

Fullerene::Fullerene()
 :AtomGroup()
{
}

Fullerene::~Fullerene()
{
}

Fullerene::Fullerene(const Fullerene &c)
 :AtomGroup(c)
{
}

Fullerene& 
Fullerene::operator =(const Fullerene &c) 
{
  if (this == &c) return *this;
  AtomGroup::operator=(c);
  return *this;
}

void
Fullerene::put(std::ostream& os) const
{
  AtomGroup::put(os);
}

void
Fullerene::get(std::istream& is)
{
  AtomGroup::get(is);
}

void
Fullerene::build(const mdtk::SimLoop& ml)
{
  const AtomsContainer &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = *ac[i];
    if (a.ID == C_EL && a.coords.z < -3.615*Ao*1.5)
      addAtom(a);
  }
  REQUIRE(atoms.size() == 60);
}  

}
