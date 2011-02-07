/* 
   mdtrajsim (molecular dynamics trajectory simulator)

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011
   Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>

#include "../common.h"

using namespace mdtk;

class CustomSimLoop : public SimLoop
{
public:
  CustomSimLoop();
  void doBeforeIteration();
  void doAfterIteration();
};

CustomSimLoop::CustomSimLoop()
  :SimLoop()
{
}

void
CustomSimLoop::doBeforeIteration()
{
  if (simTime < 4.0*ps) simTimeSaveTrajInterval = 0.05*ps;
  else simTimeSaveTrajInterval = 0.2*ps;
}

void
CustomSimLoop::doAfterIteration()
{
}

bool
isAlreadyFinished();

int
main(int argc , char *argv[])
{
  if (argc > 1 && !strcmp(argv[1],"--version"))
  {
    std::cout << "mdtrajsim (Molecular dynamics trajectory simulator) ";
    mdtk::release_info.print();
    return 0;
  }

  if (argc > 1 && (!std::strcmp(argv[1],"--help") || !std::strcmp(argv[1],"-h")))
  {
    std::cout << "\
Usage: mdtrajsim [OPTION]... [FILE]\n\
Simulates molecular dynamics trajectory described in the FILE\
 (in.mde.gz file in the current directory by default)\n\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
    return 0;
  }

  if (isAlreadyFinished()) return 0;

try
{
  CustomSimLoop mdloop;

  setupPotentials(mdloop);
  if (yaatk::exists("simloop.conf") || yaatk::exists("simloop.conf.bak")) // have to continue interrupted simulation ?
  {
    mdloop.loadstate();
  }
  else
  {
    std::string inputFile = "in.mde";
    if (argc > 1) inputFile = argv[1];

    yaatk::text_ifstream fi(inputFile.c_str());
    mdloop.loadFromMDE(fi);
//    mdloop.loadFromMDE_OLD(fi);
    fi.close();
  yaatk::text_ofstream fo1("mde""_init");
    mdloop.saveToStream(fo1);
  fo1.close();
  }

  mdloop.execute();

  if (mdloop.simTime >= mdloop.simTimeFinal) // is simulation really finished ?
  {
    yaatk::text_ofstream fo2("mde""_final");
    mdloop.saveToStream(fo2);
    fo2.close();
    mdloop.writetrajXVA();
  }
}  
catch(mdtk::Exception& e)
{ 
  std::cerr << "MDTK exception in the main thread: "
            << e.what() << std::endl;
  std::cerr << "Probably wrong input file format or no input file (default is in.mde)" << std::endl;
  return 1;
}
catch(...)
{
  std::cerr << "Unknown exception in the main thread." << std::endl; 
  return 2;
}  

  return 0;
}

bool
isAlreadyFinished()
{
  {
    if (yaatk::exists("mde_final")) return true;
  }
  {
    std::ifstream ifinal("in.mde.after_crash");
    if (ifinal)
      return true;
    else
      ifinal.close();
  }
  return false;
}

