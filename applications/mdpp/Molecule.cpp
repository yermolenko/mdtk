/*
   The Molecule class.

   Copyright (C) 2010, 2012 Oleksandr Yermolenko
   <oleksandr.yermolenko@gmail.com>

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
