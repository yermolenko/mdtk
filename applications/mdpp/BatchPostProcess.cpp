/* 
   Molecular dynamics postprocessor, BatchPostProcess  class

   Copyright (C) 2007, 2008, 2009, 2010, 2011 Oleksandr Yermolenko
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

#include "BatchPostProcess.hpp"

#include <algorithm>

#include <fstream>
#include <mdtk/tools.hpp>

namespace mdepp
{

using mdtk::Exception;

BatchPostProcess::BatchPostProcess(std::string mdeppinPath)
  :pps()
{
  char trajsetDir[10000];

  REQUIRE(yaatk::exists(mdeppinPath));

  std::ifstream fi(mdeppinPath.c_str());

  while(fi.getline(trajsetDir, 1000-1, '\n'))
  {
    if (strlen(trajsetDir) > 0)
    {
      pps.push_back(new mdepp::StatPostProcess(trajsetDir));
    }
  };

  fi.close();
}

BatchPostProcess::BatchPostProcess()
  :pps()
{
}

void
BatchPostProcess::saveToStream(std::ostream& os) const
{
  os << pps.size() << "\n";
  for(size_t i = 0; i < pps.size(); i++)
    pps[i]->saveToStream(os);
}

void
BatchPostProcess::loadFromStream(std::istream& is)
{
  size_t sz;
  is >> sz;
  pps.resize(sz);
  for(size_t i = 0; i < pps.size(); i++)
  {
    pps[i] = new mdepp::StatPostProcess();
    pps[i]->loadFromStream(is);
  }
}

void
BatchPostProcess::execute()
{
  for(size_t i = 0; i < pps.size(); ++i)
    pps[i]->execute();
}

void
BatchPostProcess::printResults()
{
  yaatk::mkdir("results");
  yaatk::chdir("results");

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess* pp = pps[i];

    yaatk::mkdir(pp->id.str.c_str());
    yaatk::chdir(pp->id.str.c_str());

    yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

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

    using namespace mdtk;

//  pp->printClusterDynamicsTotal();

//  pp->spottedTotalByMass();

//  pp->spottedByDepth();

    pp->buildMassSpectrum();

    yaatk::chdir("..");
  }

  {
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessCluster,"yields-Cluster");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessProjectile,"yields-Projectile");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessSubstrate,"yields-Substrate");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessAll,"yields-All");
  }

  yaatk::chdir("..");

  {
    yaatk::StreamToFileRedirect cout_redir(std::cout,"yieldsout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"yieldserr.txt");

    for(size_t i = 0; i < pps.size(); ++i)
    {
      mdepp::StatPostProcess* pp = pps[i];

      TRACE(pp->id.str);

      TRACE(pp->trajData.size());

      TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessCluster));
      TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessProjectile));
      TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessSubstrate));
      TRACE(pp->getAverageYield(mdepp::StatPostProcess::ProcessAll));

      TRACE("--------------------");
    }
  }
}

void
BatchPostProcess::plotYieldsAgainstIonEnergy(StatPostProcess::FProcessClassicMolecule fpm,
                                             std::string idStr) const
{
  std::stringstream fnb;
  fnb << idStr;

  ofstream fplt((fnb.str()+".plt").c_str());
  fplt << "\
reset\n\
set xrange [0:500]\n\
set yrange [0:*]\n\
set xtics (100,200,400)\n\
set pointsize 1.5\n\
#set grid ytics\n\
#set key left top\n\
#set key right top\n\
set key spacing 1.5\n\
set xlabel \"Енергія іона, еВ\"\n\
set ylabel \"Коефіцієнт зворотного розсіювання, атом/іон\"\n\
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 8cm, 8cm \"Arial,18\" enhanced\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess* pp = pps[i];

    data << pp->id.ionEnergy << " " << pp->getAverageYield(fpm) << "\n";

    if (i%3 == 2)
    {
      std::ostringstream cmd;
      cmd << "'-' title \"{/Italic "
          << ElementIDtoString(pp->id.ionElement) << " "
          << ElementIDtoString(pp->id.clusterElement)
          << "_{" << pp->id.clusterSize << "}"
          << "}\" "
          << "with linespoints";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

  REQUIRE(plotCmds.size() > 0);
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fplt << plotCmds[i] << ",\\\n";
    else
      fplt << plotCmds[i] << "\n";
  }

  fplt << data.str();

  fplt.close();
}

}
