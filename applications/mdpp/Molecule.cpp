#include "Molecule.hpp"

namespace mdepp
{

Molecule::Molecule()
 :AtomGroup()
{
  initParams();
}

Molecule::~Molecule()
{
}

Molecule::Molecule(const Molecule &c)
 :AtomGroup(c)
{
  initParams();    
}

Molecule& 
Molecule::operator =(const Molecule &c) 
{
  if (this == &c) return *this;
  AtomGroup::operator=(c);
  return *this;
}

void
Molecule::put(std::ostream& os) const
{
  AtomGroup::put(os);
}

void
Molecule::get(std::istream& is)
{
  AtomGroup::get(is);
}

void
Molecule::buildFromAtom(const mdtk::Atom& a, const AtomGroup& ag)
{
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);

  for(size_t i = 0; i < ag.atoms.size(); i++)
  {
    const mdtk::Atom& nb_a = ag.atoms[i];
    if (!isHandled(nb_a)) continue;
    if (Rc(a,nb_a) >= depos(a,nb_a).module())
    {
      buildFromAtom(nb_a,ag);
//      if (atoms.size() == 0) break;
    }
  }  
}  

}
