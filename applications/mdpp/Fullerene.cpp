#include "Fullerene.hpp"

namespace mdepp
{

Fullerene::Fullerene()
  :AtomGroup(),cluster()
{
}

Fullerene::~Fullerene()
{
}

Fullerene::Fullerene(const Fullerene &c)
  :AtomGroup(c),cluster(c.cluster)
{
}

Fullerene& 
Fullerene::operator =(const Fullerene &c) 
{
  if (this == &c) return *this;
  AtomGroup::operator=(c);
  cluster = c.cluster;
  return *this;
}

void
Fullerene::put(std::ostream& os) const
{
  AtomGroup::put(os);
  cluster.put(os);
}

void
Fullerene::get(std::istream& is)
{
  AtomGroup::get(is);
  cluster.get(is);
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
  cluster.build(ml);
}  

bool
Fullerene::isEndoFullerene()
{
  return cluster.atoms.size() != 0;
}

}
