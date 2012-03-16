/*
   Some physical constants (header file).

   Copyright (C) 2004, 2005, 2009, 2012 Oleksandr Yermolenko
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

#ifndef mdtk_consts_hpp
#define mdtk_consts_hpp

#include <mdtk/config.hpp>
#include <mdtk/Exception.hpp>
#include <mdtk/Vector3D.hpp>
#include <string>

namespace mdtk
{

extern const Float amu;
extern const Float e;
extern const Float Ao;
extern const Float eV;
extern const Float kb;
extern const Float K;

extern const Float fs;
extern const Float ps;

extern const Float Deg;

enum ElementID
{
  H_EL = 1,
  C_EL = 12,
  Cu_EL = 64,
  Ag_EL = 108,
  Au_EL = 197,
  Ar_EL = 40,
  Xe_EL = 131,
  DUMMY_EL = -1
};

#define EL_ID_size 200
//extern const int EL_ID_size;

inline
std::string
ElementIDtoString(ElementID id)
{
  std::string str;
  switch (id)
  {
  case H_EL  : str = "H"; break;
  case C_EL  : str = "C"; break;
  case Cu_EL : str = "Cu"; break;
  case Ag_EL : str = "Ag"; break;
  case Au_EL : str = "Au"; break;
  case Ar_EL : str = "Ar"; break;
  case Xe_EL : str = "Xe"; break;
  case DUMMY_EL :
  default : Exception("Unknown element");
  }
  return str;
}

extern const Vector3D NO_PBC;

}


#endif
