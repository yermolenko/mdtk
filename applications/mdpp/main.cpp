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
  std::vector<mdepp::BatchPostProcess> bpps;
  bpps.resize(4);

  for(size_t haloIndex = 0; haloIndex < bpps.size(); haloIndex++)
  {
    std::stringstream resultsDir;
    if (haloIndex == 0)
      resultsDir << "area" << haloIndex;
    else
      resultsDir << "halo" << haloIndex;
    yaatk::ChDir cd(resultsDir.str());

    bool fullpp = !yaatk::exists("pp.state.orig");

    if (!fullpp)
    {
      yaatk::text_ifstream fi("pp.state.orig");
      bpps[haloIndex].loadFromStream(fi);
      fi.close();

      bpps[haloIndex].printResults();

      {
        yaatk::text_ofstream fo("pp.state.after");
        bpps[haloIndex].saveToStream(fo);
        fo.close();
      }
    }
    else
    {
      bpps[haloIndex] = mdepp::BatchPostProcess("../../mdepp.in", haloIndex);
      bpps[haloIndex].execute();

      bpps[haloIndex].printResults();

      {
        yaatk::text_ofstream fo("pp.state.orig");
        bpps[haloIndex].saveToStream(fo);
        fo.close();
      }

      {
        yaatk::text_ofstream fo("pp.state.after");
        bpps[haloIndex].saveToStream(fo);
        fo.close();
      }
    }
  }

  {
    yaatk::ChDir cd("area0-reproducibility-test");
    mdepp::BatchPostProcess bppsArea1;
    for(size_t haloIndex = 0; haloIndex <= 0; haloIndex++)
      bppsArea1 = bpps[haloIndex];

    bppsArea1.printResults();

    {
      yaatk::text_ofstream fo("pp.state.after");
      bppsArea1.saveToStream(fo);
      fo.close();
    }
  }
/*
  {
    yaatk::ChDir cd("area1");
    mdepp::BatchPostProcess bppsArea1;
    for(size_t haloIndex = 0; haloIndex <= 1; haloIndex++)
      bppsArea1.addHalo(bpps[haloIndex]);
  }
  {
    yaatk::ChDir cd("area2");
    mdepp::BatchPostProcess bppsArea2;
    for(size_t haloIndex = 0; haloIndex <= 2; haloIndex++)
      bppsArea2.addHalo(bpps[haloIndex]);
  }
  {
    yaatk::ChDir cd("area3");
    mdepp::BatchPostProcess bppsArea3;
    for(size_t haloIndex = 0; haloIndex <= 3; haloIndex++)
      bppsArea3.addHalo(bpps[haloIndex]);
  }
*/
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

