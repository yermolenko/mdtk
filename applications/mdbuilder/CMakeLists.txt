#  CMakeLists.txt file for the Molecular Dynamics Experiments Builder
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

SET(MDTK_CurrentTarget "mdbuilder${MDTK_BINARY_SUFFIX}")

IF(OSMESA_FOUND)
  include_directories (${MDTK_SOURCE_DIR} ${GSL_INCLUDE_DIRS} ${OSMESA_INCLUDE_DIR})
ELSE(OSMESA_FOUND)
  include_directories (${MDTK_SOURCE_DIR} ${GSL_INCLUDE_DIRS} ${FLTK_INCLUDE_DIR})
ENDIF(OSMESA_FOUND)

link_directories (${MDTK_BINARY_DIR})

add_executable (${MDTK_CurrentTarget}
  main.cpp
  common.cxx
  experiments/H2.cpp
  experiments/FCC.cpp
  experiments/Clusters.cpp
  experiments/Polyethylene.cpp
  experiments/Graphite.cpp
  experiments/Fullerite.cpp
  experiments/Fulleride.cpp
  )

target_link_libraries (${MDTK_CurrentTarget} mdtk yaatk ${YAATK_COMPRESSION_LIBRARIES} ${GSL_LIBRARIES})

IF(WIN32)
  target_link_libraries (${MDTK_CurrentTarget} ole32 uuid comctl32 wsock32 gdi32)
ENDIF(WIN32)

IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(WIN32)
    SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES LINK_FLAGS "-s -static")
  ELSE(WIN32)
    SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES LINK_FLAGS "-s")
  ENDIF(WIN32)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(OSMESA_FOUND)
  MESSAGE(STATUS "OSMesa found. Using OSMesa in mdbuilder.")
  target_link_libraries (${MDTK_CurrentTarget} ${OSMESA_LIBRARIES} GLU)
  SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES COMPILE_FLAGS "-DMDBUILDER_OSMESA")
ELSE(OSMESA_FOUND)
  MESSAGE(STATUS "OSMesa not found. MDBuilder will require a graphics environment.")
  target_link_libraries (${MDTK_CurrentTarget} ${FLTK_LIBRARIES} ${OPENGL_LIBRARIES})
ENDIF(OSMESA_FOUND)

install(TARGETS ${MDTK_CurrentTarget}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

