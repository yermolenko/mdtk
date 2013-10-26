/*
   Some configuration constants for YAATK (header file).

   Copyright (C) 2004, 2005, 2009, 2013 Oleksandr Yermolenko
   <oleksandr.yermolenko@gmail.com>

   This file is part of YAATK, Yet another auxiliary toolkit.

   YAATK is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   YAATK is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with YAATK.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef yaatk_config_hpp
#define yaatk_config_hpp

#include <iostream>

namespace yaatk
{

extern bool verboseTrace;

struct VerboseOutput
{
  bool prevVerboseOutputState;
  VerboseOutput(bool newState)
  :prevVerboseOutputState(verboseTrace)
    {
      verboseTrace = newState;
    }
  ~VerboseOutput()
    {
      verboseTrace = prevVerboseOutputState;
    }
};

}

#endif

