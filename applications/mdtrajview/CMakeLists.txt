#  CMakeLists.txt file for the Molecular Dynamics Trajectory Viewer
#
#  Copyright (C) 2007, 2008, 2009, 2010, 2012 Oleksandr Yermolenko
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

SET(MDTK_CurrentTarget "mdtrajview${MDTK_BINARY_SUFFIX}")

include_directories (${MDTK_SOURCE_DIR} ${ZLIB_INCLUDE_DIR})

link_directories (${MDTK_BINARY_DIR})

add_executable (${MDTK_CurrentTarget} main.cpp MainWindow.cpp VisBox.cpp bmpImage.cxx MDTrajectory.cxx)

IF(WIN32)
  target_link_libraries (${MDTK_CurrentTarget} mdtk ${FLTK_LIBRARIES} ${OPENGL_LIBRARIES} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES} ole32 uuid comctl32 wsock32 gdi32)
ELSE(WIN32)
  target_link_libraries (${MDTK_CurrentTarget} mdtk ${FLTK_LIBRARIES} ${OPENGL_LIBRARIES} ${ZLIB_LIBRARIES} ${GSL_LIBRARIES})
ENDIF(WIN32)

IF(CMAKE_COMPILER_IS_GNUCXX)
  IF(WIN32)
    SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES LINK_FLAGS "-s -static")
  ELSE(WIN32)
    SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES LINK_FLAGS "-s")
  ENDIF(WIN32)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(PNG_FOUND)
  MESSAGE("libpng found. Using libpng in mdtrajview.")
  target_link_libraries (${MDTK_CurrentTarget} ${PNG_LIBRARIES})
  SET_TARGET_PROPERTIES(${MDTK_CurrentTarget} PROPERTIES COMPILE_FLAGS "-DMDTRAJVIEW_PNG")
ELSE(PNG_FOUND)
  MESSAGE("libpng not found. mdtrajview will write BMP files instead of PNG.")
ENDIF(PNG_FOUND)

install(TARGETS ${MDTK_CurrentTarget}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

