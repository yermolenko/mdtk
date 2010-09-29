#ifndef MDBUILDER_H2_HPP
#define MDBUILDER_H2_HPP

#include "../common.hpp"

namespace mdbuilder
{

inline
void
place_H2(mdtk::SimLoop& sl)
{
  sl.atoms.push_back(new mdtk::Atom());
  sl.atoms.push_back(new mdtk::Atom());
  sl.atoms[1]->coords.x = 1.0*mdtk::Ao;
}

}

#endif
