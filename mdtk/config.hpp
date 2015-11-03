/*
   Some configuration constants for MDTK (header file).

   Copyright (C) 2004, 2005, 2009, 2013, 2015 Oleksandr Yermolenko
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

#ifndef mdtk_config_hpp
#define mdtk_config_hpp

#include <cmath>
#include <iostream>

#include <yaatk/yaatk.hpp>
#include <mdtk/release_info.hpp>

#ifndef M_PI
#define M_PI        3.14159265358979323846
#endif

#ifdef MDTK_PARALLEL
#define YAATK_PARALLEL
#endif

namespace mdtk
{

using yaatk::verboseTrace;
using yaatk::Exception;
using yaatk::MPI_Exception;

typedef double Float;

extern const int FLOAT_PRECISION;

extern std::string buildID;

}

#endif

