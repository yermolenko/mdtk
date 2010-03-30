#  CMakeLists.txt file for MDTK programming library
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

include_directories (${MDTK_SOURCE_DIR})
add_library (mdtk
  Atom.cxx
  config.cxx
  consts.cxx
  Exception.cxx
  release_info.cxx
  SimLoop.cxx
  SplineAux.cxx
  Spline.cxx
  Spline5n.cxx
  Spline2D.cxx
  Spline3D.cxx
  tools.cxx
  Vector3D.cxx
  potentials/NeighbourList.cxx
  potentials/FProxy.cxx
  potentials/FGeneral.cxx
  potentials/pairwise/FPairwise.cxx
  potentials/pairwise/FLJ.cxx
  potentials/pairwise/FBZL.cxx
  potentials/pairwise/FBM.cxx
  potentials/manybody/FManybody.cxx
  potentials/manybody/TightBinding/TightBinding.cxx
  potentials/manybody/Ackland/Ackland.cxx
  potentials/manybody/Brenner/cubici.cxx
  potentials/manybody/Brenner/Brenner.cxx
  potentials/manybody/AIREBO/cubici.cxx
  potentials/manybody/AIREBO/REBO.cxx
  potentials/manybody/AIREBO/AIREBO_LJ.cxx
  potentials/manybody/AIREBO/AIREBO_ETors.cxx
  yaatk.cxx
  procmon.cxx
)

install(TARGETS mdtk
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

install(DIRECTORY ./ DESTINATION include/mdtk FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ./ DESTINATION include/mdtk FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ./ DESTINATION include/mdtk FILES_MATCHING PATTERN "*.inc")

