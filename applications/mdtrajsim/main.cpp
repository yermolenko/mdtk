/* 
   mdtrajsim (molecular dynamics trajectory simulator)

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010 Oleksandr
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>

#include "../common.h"

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

#ifdef MDE_PARALLEL
  MPI_TEST_SUCCESS(MPI_Init(&argc,&argv));
#endif

try
{
  mdtk::SimLoop mdloop;

  setupPotentials(mdloop);
  std::ifstream simloop_conf("simloop.conf.xz");
  std::ifstream simloop_conf_bak("simloop.conf.bak.xz");
  if (simloop_conf || simloop_conf_bak) // have to continue interrupted simulation ?
  {
#ifdef MDE_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
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
#ifdef MDE_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
    if (mdtk::comm_rank == 0) {
#endif
  yaatk::text_ofstream fo1("mde""_init");
    mdloop.saveToStream(fo1);
  fo1.close();
#ifdef MDE_PARALLEL
    } // if (rank == 0)
#endif
  }

  mdloop.execute();

#ifdef MDE_PARALLEL
    if (mdtk::comm_rank == 0) {
#endif

  if (mdloop.simTime >= mdloop.simTimeFinal) // is simulation really finished ?
  {
    yaatk::text_ofstream fo2("mde""_final");
    mdloop.saveToStream(fo2);
    fo2.close();
    mdloop.writetrajXVA();
  }

#ifdef MDE_PARALLEL
    } // if (rank == 0)
#endif

#ifdef MDE_PARALLEL
  MPI_TEST_SUCCESS(MPI_Finalize());
#endif
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
    std::ifstream ifinal("mde_final.xz");
    if (ifinal)
      return true;
    else
      ifinal.close();
  }
  {
    std::ifstream ifinal("mde_final");
    if (ifinal)
      return true;
    else
      ifinal.close();
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

