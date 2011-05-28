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
  mdepp::FProcessTrajectory fpt = mdepp::trajProcess_Custom1;
//  mdepp::FProcessTrajectory fpt = mdepp::trajProcess_Cu13_at_C60_Only;

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

    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessFullerene));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessCluster));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessProjectile));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessSubstrate));
    TRACE(pp->getYield(i,mdepp::StatPostProcess::ProcessAll));
  }

  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessFullerene));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessCluster));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessProjectile));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessSubstrate));
  TRACE(pp->getYieldSum(mdepp::StatPostProcess::ProcessAll));

  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessFullerene));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessCluster));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessProjectile));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessSubstrate));
  TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessAll));

  pp->printClassicMoleculesTotal();

  pp->printFullereneInfo();

  using namespace mdtk;

  Float integralThreshold = 0.1*Ao;
  while (integralThreshold < 21.0*Ao)
  {
    pp->plotFullereneImplantDepth(false,"010",integralThreshold);
    pp->plotFullereneImplantDepth(false,"001",integralThreshold);
    pp->plotFullereneImplantDepth(false,"011",integralThreshold);
    pp->plotFullereneImplantDepth(true,"010",integralThreshold);
    pp->plotFullereneImplantDepth(true,"001",integralThreshold);
    pp->plotFullereneImplantDepth(true,"011",integralThreshold);

    if (integralThreshold < 3.0*Ao)
      integralThreshold += 0.1*Ao;
    else if (integralThreshold < 10.0*Ao)
      integralThreshold += 0.5*Ao;
    else
      integralThreshold += 5.0*Ao;
  }
/*
  pp->plotFullereneIntegrityEvolutions(10.0*ps);
  pp->plotFullereneIntegrityEvolutions(10.0*ps,16.0*Ao);
  pp->plotFullereneIntegrityEvolutions(5.0*ps);
  pp->plotFullereneIntegrityEvolutions(5.0*ps,16.0*Ao);
*/
  pp->plotFullereneIntegrity(false,"010");
  pp->plotFullereneIntegrity(false,"001");
  pp->plotFullereneIntegrity(false,"011");
  pp->plotFullereneIntegrity(true,"010");
  pp->plotFullereneIntegrity(true,"001");
  pp->plotFullereneIntegrity(true,"011");

  pp->plotFullereneIntegrityHistogram(false,"010",false);
  pp->plotFullereneIntegrityHistogram(false,"001",false);
  pp->plotFullereneIntegrityHistogram(false,"011",false);
  pp->plotFullereneIntegrityHistogram(true,"010",false);
  pp->plotFullereneIntegrityHistogram(true,"001",false);
  pp->plotFullereneIntegrityHistogram(true,"011",false);

  pp->plotFullereneIntegrityHistogram(false,"010",true);
  pp->plotFullereneIntegrityHistogram(false,"001",true);
  pp->plotFullereneIntegrityHistogram(false,"011",true);
  pp->plotFullereneIntegrityHistogram(true,"010",true);
  pp->plotFullereneIntegrityHistogram(true,"001",true);
  pp->plotFullereneIntegrityHistogram(true,"011",true);

//  pp->printClusterDynamicsTotal();

//  pp->spottedTotalByMass();

//  pp->spottedByDepth();

  pp->buildMassSpectrum();

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

