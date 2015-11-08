#  GSL detection module for CMake
#
#  Copyright (C) 2012, 2015 Oleksandr Yermolenko
#  <oleksandr.yermolenko@gmail.com>
#
#  This file is part of MDTK, the Molecular Dynamics Toolkit.
#
#  MDTK is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MDTK is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MDTK.  If not, see <http://www.gnu.org/licenses/>.
#

if (NOT GSL_INCLUDE_DIRS)
  find_path(GSL_INCLUDE_DIRS gsl/gsl_rng.h
    /usr/openwin/share/include
    )
endif (NOT GSL_INCLUDE_DIRS)

if (NOT GSL_LIBRARY_CORE)
  find_library(GSL_LIBRARY_CORE gsl
    /usr/openwin/lib
    )
endif (NOT GSL_LIBRARY_CORE)

if (NOT GSL_LIBRARY_CBLAS)
  find_library(GSL_LIBRARY_CBLAS NAMES gslcblas cblas
    /usr/openwin/lib
    )
endif (NOT GSL_LIBRARY_CBLAS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GSL DEFAULT_MSG GSL_LIBRARY_CORE GSL_LIBRARY_CBLAS GSL_INCLUDE_DIRS)

set(GSL_LIBRARIES "${GSL_LIBRARY_CORE}" "${GSL_LIBRARY_CBLAS}")

mark_as_advanced(GSL_INCLUDE_DIRS GSL_LIBRARY_CORE GSL_LIBRARY_CBLAS GSL_LIBRARIES)
