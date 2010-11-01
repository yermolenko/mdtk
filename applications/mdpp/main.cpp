/* 
   mdpp (molecular dynamics postprocessor)

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <mdtk/SimLoop.hpp>

#include "StatPostProcess.hpp"

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
  mdepp::FProcessTrajectory fpt = mdepp::trajProcess_Cu13_at_C60_Only;

  std::vector<std::string> trajDirNames;

  mdepp::findTrajDirs("../mdepp.in","../../trajset"DIR_DELIMIT_STR,trajDirNames,fpt);

  mdepp::StatPostProcess* pp = NULL;
  bool fullpp = !yaatk::exists("pp.state.orig");

  if (!fullpp)
  {
    yaatk::text_ifstream fi("pp.state.orig");
    pp = new mdepp::StatPostProcess();
    pp->loadFromStream(fi);
    fi.close();
  }
  else
  {
    pp = new mdepp::StatPostProcess(trajDirNames);
    pp->execute();
  }

  yaatk::mkdir("results");
  yaatk::chdir("results");

  for(size_t i = 0; i < pp->trajData.size(); i++)
  {
    TRACE(pp->trajData[i].trajDir);
//    TRACE(pp->getAboveSpottedHeight(i,mdepp::StatPostProcess::ProcessAll));

    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessCluster));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessSubstrate));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessClusterAndSubstrate));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessProjectile));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessAll));
  }

//  TRACE(pp->getAboveSpottedHeightTotal());
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessCluster));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessSubstrate));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessClusterAndSubstrate));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessProjectile));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessAll));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessCluster));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessSubstrate));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessClusterAndSubstrate));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessProjectile));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessAll));
  TRACE(pp->getAverageYieldProgress(mdepp::StatPostProcess::ProcessCluster));
  TRACE(pp->getAverageYieldProgress(mdepp::StatPostProcess::ProcessSubstrate));
  TRACE(pp->getAverageYieldProgress(mdepp::StatPostProcess::ProcessClusterAndSubstrate));
  TRACE(pp->getAverageYieldProgress(mdepp::StatPostProcess::ProcessProjectile));
  TRACE(pp->getAverageYieldProgress(mdepp::StatPostProcess::ProcessAll));

  TRACE(pp->getAverageEnergyOfSputtered(mdepp::StatPostProcess::ProcessCluster)/mdtk::eV);
  TRACE(pp->getAverageEnergyOfSputtered(mdepp::StatPostProcess::ProcessSubstrate)/mdtk::eV);
  TRACE(pp->getAverageEnergyOfSputtered(mdepp::StatPostProcess::ProcessClusterAndSubstrate)/mdtk::eV);
  TRACE(pp->getAverageEnergyOfSputtered(mdepp::StatPostProcess::ProcessProjectile)/mdtk::eV);
  TRACE(pp->getAverageEnergyOfSputtered(mdepp::StatPostProcess::ProcessAll)/mdtk::eV);

  pp->printMoleculesTotal();

  pp->printFullereneInfo();

  pp->printClusterDynamicsTotal();

  yaatk::chdir("..");

  if (fullpp)
  {
    yaatk::text_ofstream fo("pp.state.orig");
    pp->saveToStream(fo);
    fo.close();
  }

  yaatk::text_ofstream fo("pp.state.after");
  pp->saveToStream(fo);
  fo.close();

  delete pp;
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

