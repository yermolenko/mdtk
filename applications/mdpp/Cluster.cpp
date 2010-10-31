#include "Cluster.hpp"

namespace mdepp
{

Cluster::Cluster()
 :AtomGroup()
{
}

Cluster::~Cluster()
{
}

Cluster::Cluster(const Cluster &c)
 :AtomGroup(c)
{
}

Cluster& 
Cluster::operator =(const Cluster &c) 
{
  if (this == &c) return *this;
  AtomGroup::operator=(c);
  return *this;
}

void
Cluster::put(std::ostream& os) const
{
  AtomGroup::put(os);
}

void
Cluster::get(std::istream& is)
{
  AtomGroup::get(is);
}

void
Cluster::build(const mdtk::SimLoop& ml)
{
  const AtomsContainer &ac = ml.atoms;
  for(size_t i = 0; i < ac.size(); i++)
  {
    const mdtk::Atom a = *ac[i];
    if (a.ID == Cu_EL && a.coords.z < -3.615*Ao*1.5)
      addAtom(a);
  }
  REQUIRE(atoms.size() == 60);
}  

}
