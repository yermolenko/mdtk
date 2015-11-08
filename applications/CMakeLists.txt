#  CMakeLists.txt file for MDTK applications and tools
#
#  Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2015 Oleksandr
#  Yermolenko <oleksandr.yermolenko@gmail.com>
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

add_subdirectory (mdtrajsim)

IF(WIN32)
  message(STATUS "Molecular Dynamics Postprocessor will not be compiled on Win32. Use *nix to build it.")
ELSE(WIN32)
  add_subdirectory (mdpp)
ENDIF(WIN32)

IF(MSVC)
  SET(FLTK_LIBRARIES fltk fltkgl)
  SET(FLTK_FOUND true)
ELSE(MSVC)
  SET(FLTK_SKIP_FLUID true)
  FIND_PACKAGE(FLTK)
ENDIF(MSVC)

FIND_PACKAGE(OpenGL)

IF(OPENGL_FOUND AND FLTK_FOUND)
  add_subdirectory (mdtrajview)
  add_subdirectory (mdbuilder)
ELSE(OPENGL_FOUND AND FLTK_FOUND)
  message(STATUS "FLTK or OpenGL not found! Install it in order to visualize MD trajectories and prepare experiments using mdbuilder.")
ENDIF(OPENGL_FOUND AND FLTK_FOUND)
