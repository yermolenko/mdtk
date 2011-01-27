#  Top-level CMakeLists.txt file for MDTK
#
#  Copyright (C) 2007, 2008, 2009, 2010 Oleksandr Yermolenko
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

project (MDTK)

cmake_minimum_required(VERSION 2.4.5)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

SET(MDTK_BINARY_SUFFIX "")

IF(WIN32)
  ADD_DEFINITIONS(-D__WIN32__)
ENDIF(WIN32)

IF(MDTK_HEAVY_OPTIMZE)

  IF(BUILD_INTELCC)
    SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -march=core2")
    SET(MDTK_BINARY_SUFFIX "_intel")
    MESSAGE("Using Intel C++ Compiler with optimizations!")
  ENDIF(BUILD_INTELCC)

  IF(CMAKE_COMPILER_IS_GNUCXX)
    SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-pedantic -W -Wall -Wno-long-long -fexpensive-optimizations -O3 -march=nocona -mfpmath=sse -funroll-loops -static")
    SET(CMAKE_C_FLAGS_DISTRIBUTION "${CMAKE_CXX_FLAGS_DISTRIBUTION}")
    MESSAGE("Using GNU C++ Compiler with optimizations!")
  ENDIF(CMAKE_COMPILER_IS_GNUCXX)

  SET(CMAKE_BUILD_TYPE distribution)

ELSEIF(MDTK_HEAVY_OPTIMZE)

  SET(CMAKE_BUILD_TYPE Release)

ENDIF(MDTK_HEAVY_OPTIMZE)

#SET(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_BUILD_TYPE Release)

#FIND_PACKAGE(GSL REQUIRED) # no such macro in cmake yet
IF(MSVC)
  SET(GSL_LIBRARIES gsl cblas)
  SET(GSL_FOUND true)
ELSE(MSVC)
  SET(GSL_LIBRARIES gsl gslcblas)
  SET(GSL_FOUND true)
ENDIF(MSVC)

IF(WIN32 AND CMAKE_COMPILER_IS_GNUCXX) # for that combination cmake fails to detect zlib :(
  SET(ZLIB_LIBRARIES z)
  SET(ZLIB_FOUND true)
ELSE(WIN32 AND CMAKE_COMPILER_IS_GNUCXX)
  FIND_PACKAGE(ZLIB REQUIRED)
ENDIF(WIN32 AND CMAKE_COMPILER_IS_GNUCXX)

add_subdirectory (mdtk)
add_subdirectory (applications)
add_subdirectory (examples)

# may be useful for some old versions of Intel compilers
#  SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O3 -D \"__sync_fetch_and_add(ptr,addend)=_InterlockedExchangeAdd(const_cast<void*>(reinterpret_cast<volatile void*>(ptr)), addend)\"")

