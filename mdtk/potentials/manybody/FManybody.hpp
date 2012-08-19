/*
   The generalized manybody interatomic potential class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2012 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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

#ifndef mdtk_FManybody_hpp
#define mdtk_FManybody_hpp

#include <mdtk/potentials/FGeneral.hpp>

namespace mdtk
{

#define FMANYBODY_PAIRS_RESERVE_ADD 5

class FManybody : public FGeneral
{
public:
  FManybody();
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::SaveToStream(os,smode);
  }
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    FGeneral::LoadFromStream(is,smode);
  }
};

} // namespace apme

#endif

