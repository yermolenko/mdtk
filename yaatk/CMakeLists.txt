#  CMakeLists.txt file for YAATK programming library
#
#  Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2015
#  Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>
#
#  This file is part of YAATK, Yet another auxiliary toolkit.
#
#  YAATK is free software: you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  YAATK is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with YAATK.  If not, see <http://www.gnu.org/licenses/>.
#

include_directories (../ ${YAATK_COMPRESSION_INCLUDE_DIRS})
add_library (yaatk
  config.cxx
  Exception.cxx
  yaatk.cxx
)

install(TARGETS yaatk
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)

install(DIRECTORY ./ DESTINATION include/yaatk FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ./ DESTINATION include/yaatk FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ./ DESTINATION include/yaatk FILES_MATCHING PATTERN "*.inc")

