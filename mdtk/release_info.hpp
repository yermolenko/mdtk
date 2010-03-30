/*
   The ReleaseInfo class header file.

   Copyright (C) 2004, 2009, 2010 Oleksandr Yermolenko
   <oleksandr.yermolenko@gmail.com>

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

#ifndef	_MDTK_RELEASE_INFO_H_
#define	_MDTK_RELEASE_INFO_H_

#include <string>
#include <iostream>

namespace mdtk
{


struct ReleaseInfo
{
  const std::string PRODUCT_NAME;
  const std::string PRODUCT_VERSION;
  const std::string PRODUCT_BUILDINFO;
  const std::string PRODUCT_COPYRIGHT;
  const std::string PRODUCT_NOTES;
  ReleaseInfo(std::string name,
              std::string version,
              std::string buildinfo,
              std::string copyright,
              std::string notes):
  PRODUCT_NAME(name),
  PRODUCT_VERSION(std::string("")+version),
  PRODUCT_BUILDINFO(buildinfo+std::string("Built ")+std::string(__DATE__)+" at "+__TIME__),
  PRODUCT_COPYRIGHT("Copyright (C) "+copyright),
  PRODUCT_NOTES(notes)
  {
  }
  void print()
  {
    std::cout << getInfo();
  }
  void print_stderr()
  {
    std::cerr << getInfo();
  }
  std::string getInfo()
  {
    return /*PRODUCT_NAME + " " +*/ PRODUCT_VERSION
      /*               + " ("
		       + PRODUCT_BUILDINFO + ")\n" */ + "\n"
               + PRODUCT_COPYRIGHT + "\n" + PRODUCT_NOTES+"\n";
  } 
};

extern ReleaseInfo release_info;

}

#endif

