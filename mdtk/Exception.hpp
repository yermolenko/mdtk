/*
   Custom exception classes header file.

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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

#ifndef	mdtk_Exception_hpp
#define	mdtk_Exception_hpp

#include <string>

namespace mdtk
{

class Exception
{
  std::string _msg;
public:
  Exception() : _msg("Unknown MDTK exception") { }
  Exception(const char* msg) : _msg(msg) { }
  Exception(std::string msg) : _msg(msg) { }
  const char* what() const
  {
    return _msg.c_str();
  }
}; 

}

class MPI_Exception
{
  std::string _msg;
public:
  MPI_Exception() : _msg("Unknown internal MPI exception") { }
  MPI_Exception(const char* msg) : _msg(msg) { }
  MPI_Exception(std::string msg) : _msg(msg) { }
  MPI_Exception(int errorCode) : _msg("Error Code ") { _msg += errorCode; }
  const char* what() const
  {
    return _msg.c_str();
  }
}; 

#endif

