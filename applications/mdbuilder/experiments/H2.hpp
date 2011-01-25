#ifndef MDBUILDER_H2_HPP
#define MDBUILDER_H2_HPP

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

inline
void
place_H2_simple(mdtk::SimLoop& sl)
{
  place(H_EL,sl,Vector3D(0,0,0));
  place(H_EL,sl,Vector3D(1.0*mdtk::Ao,0,0));
}

inline
void
place_H2(mdtk::SimLoop& sl)
{
  TRACE("Ok");
  glPushMatrix();
  place(H_EL,sl);
  glTranslated(1.0*mdtk::Ao,0,0);
  place(H_EL,sl);
  glPopMatrix();
}

}

#endif
