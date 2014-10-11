/* 
   mdpp (molecular dynamics postprocessor)

   Copyright (C) 2010, 2011 Oleksandr Yermolenko
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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>

#include "StatPostProcess.hpp"
#include "BatchPostProcess.hpp"

int
main(int argc , char *argv[])
{
  if (argc > 1 && !strcmp(argv[1],"--version"))
  {
    std::cout << "mdpp (Molecular dynamics postprocessor) ";
    mdtk::release_info.print();
    return 0;
  }

  if (argc > 1 && (!std::strcmp(argv[1],"--help") || !std::strcmp(argv[1],"-h")))
  {
    std::cout << "\
Usage: mdpp [OPTION]... [FILE]\n\
Postprocess molecular dynamics experiment data.\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
    return 0;
  }

try
{
  bool fullpp = !yaatk::exists("pp.state.orig");

  if (!fullpp)
  {
    yaatk::text_ifstream fi("pp.state.orig");
    mdepp::BatchPostProcess pp;
    pp.loadFromStream(fi);
    fi.close();

    pp.printResults();

    {
      yaatk::text_ofstream fo("pp.state.after");
      pp.saveToStream(fo);
      fo.close();
    }
  }
  else
  {
    mdepp::BatchPostProcess pp("../mdepp.in");
    pp.execute();

    pp.printResults();

    {
      yaatk::text_ofstream fo("pp.state.orig");
      pp.saveToStream(fo);
      fo.close();
    }

    {
      yaatk::text_ofstream fo("pp.state.after");
      pp.saveToStream(fo);
      fo.close();
    }
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

