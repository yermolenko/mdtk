/*
   mdbuilder (molecular dynamics experiments preparation tool)

   Copyright (C) 2010 Oleksandr Yermolenko
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

#include "../mdtrajview/VisBox.hpp"
#include <FL/Fl.H>
#include <cstring>

//#include "pentacene.hpp"

#include <iostream>
#include <fstream>

#include "experiments/H2.hpp"

int main(int argc, char *argv[])
{
  if (argc > 1 && !strcmp(argv[1],"--version"))
  {
    std::cout << "mdbuilder (molecular dynamics experiments preparation tool) ";
    mdtk::release_info.print();
    return 0;
  }

  if (argc > 1 && (!std::strcmp(argv[1],"--help") || !std::strcmp(argv[1],"-h")))
  {
    std::cout << "\
Usage: mdbuilder [OPTION]... \n\
Prepares molecular dynamics experiments.\n\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
    return 0;
  }

  xmde::VisBox visualizer(15,35,500,500,"","");

  mdbuilder::place_H2(*visualizer.ml_);

  {
    std::ofstream fomde("two_atoms.mde");
    visualizer.ml_->saveToMDE(fomde);
    fomde.close();
    YAATK_ZIP_FILE("two_atoms.mde");
  }

  visualizer.updateData();

  visualizer.show();

  Fl::run();

  return 0;
}
