#  OSMesa detection module for CMake
#
#  Copyright (C) 2012 Oleksandr Yermolenko
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

if (NOT OSMESA_INCLUDE_DIR)
  find_path(OSMESA_INCLUDE_DIR GL/osmesa.h
    /usr/openwin/share/include
    /opt/graphics/OpenGL/include
    )
endif (NOT OSMESA_INCLUDE_DIR)

if (NOT OSMESA_LIBRARIES)
  find_library(OSMESA_LIBRARIES OSMesa
    /opt/graphics/OpenGL/lib
    /usr/openwin/lib
    )
endif (NOT OSMESA_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OSMesa DEFAULT_MSG OSMESA_LIBRARIES OSMESA_INCLUDE_DIR)

mark_as_advanced(OSMESA_INCLUDE_DIR OSMESA_LIBRARIES)
