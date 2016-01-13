/* 
   Molecular dynamics postprocessor, BatchPostProcess  class

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2014, 2015, 2016
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

#include "BatchPostProcess.hpp"

#include <algorithm>

#include <fstream>
#include <mdtk/tools.hpp>
#include <mdtk/SnapshotList.hpp>

namespace mdepp
{

using mdtk::Exception;

size_t
BatchPostProcess::getHaloIndex(std::string trajsetDir)
{
  std::string s = yaatk::extractItemFromEnd(trajsetDir,0);

  bool haloSizeRecognized = false;
  Float haloSizeAo = 0.0;

  {
    size_t istart = s.find("halo-");
    REQUIRE(istart == 0);
    istart += 5;

    size_t iend = s.find("Ao");
    REQUIRE(iend == s.size()-2);

    std::string haloSizeString = s.substr(istart,iend-istart);
    {
      istringstream is(haloSizeString);
      is >> haloSizeAo;
    }

    haloSizeRecognized = true;
  }

  REQUIRE(haloSizeRecognized);

  size_t haloIndex = 3;
  if (haloSizeAo < 11)
    haloIndex = 2;
  if (haloSizeAo < 6)
    haloIndex = 1;
  if (haloSizeAo < 2)
    haloIndex = 0;

  return haloIndex;
}

BatchPostProcess::BatchPostProcess(const std::string mdeppinPath, const size_t haloIndex)
  :pps()
{
  char trajsetDir[10000];

  REQUIRE(yaatk::exists(mdeppinPath));

  std::ifstream fi(mdeppinPath.c_str());

  while(fi.getline(trajsetDir, 1000-1, '\n'))
  {
    if (strlen(trajsetDir) > 0 && trajsetDir[0] != '#' &&
        haloIndex == getHaloIndex(trajsetDir))
    {
      std::vector<std::string> pp;
      pp.push_back(trajsetDir);
      pps.push_back(pp);
    }
  }

  fi.close();
}

BatchPostProcess::BatchPostProcess()
  :pps()
{
}

void
BatchPostProcess::addHalo(const std::string mdeppinPath, const size_t haloIndex)
{
  REQUIRE(pps.size() != 0);
  for(size_t i = 0; i < pps.size(); ++i)
  {
    REQUIRE(pps[i].size() > 0);
    REQUIRE(getHaloIndex(pps[i][0]) != haloIndex);

    char trajsetDir[10000];

    REQUIRE(yaatk::exists(mdeppinPath));

    std::ifstream fi(mdeppinPath.c_str());

    size_t haloMatches = 0;

    while(fi.getline(trajsetDir, 1000-1, '\n'))
    {
      if (strlen(trajsetDir) > 0 && trajsetDir[0] != '#' &&
          haloIndex == getHaloIndex(trajsetDir))
      {
        TRACE(pps[i][0]);
        TRACE(trajsetDir);
        mdepp::StatPostProcess::Id id0(yaatk::extractItemFromEnd(pps[i][0],1));
        mdepp::StatPostProcess::Id idi(yaatk::extractItemFromEnd(trajsetDir,1));

        if (id0.clusterElement == idi.clusterElement &&
            id0.clusterSize == idi.clusterSize &&
            id0.ionElement == idi.ionElement &&
            id0.ionEnergy == idi.ionEnergy)
        {
          pps[i].push_back(trajsetDir);
          ++haloMatches;
          TRACE(haloMatches);
        }
      }
    }

    REQUIRE(haloMatches <= 1);

    fi.close();
  }
}

void
BatchPostProcess::printResults() const
{
  yaatk::mkdir("results");
  yaatk::chdir("results");

  for(size_t i = 0; i < pps.size(); ++i)
  {
    const mdepp::StatPostProcess pp(pps[i]);

    yaatk::mkdir(pp.id.str.c_str());
    yaatk::chdir(pp.id.str.c_str());

    yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

    for(size_t i = 0; i < pp.trajData.size(); i++)
    {
      TRACE(pp.trajData[i]->trajDir);

      TRACE(pp.getYield(i,mdepp::StatPostProcess::ProcessFullerene));
      TRACE(pp.getYield(i,mdepp::StatPostProcess::ProcessCluster));
      TRACE(pp.getYield(i,mdepp::StatPostProcess::ProcessProjectile));
      TRACE(pp.getYield(i,mdepp::StatPostProcess::ProcessSubstrate));
      TRACE(pp.getYield(i,mdepp::StatPostProcess::ProcessAll));
    }

    TRACE(pp.getYieldSum(mdepp::StatPostProcess::ProcessFullerene));
    TRACE(pp.getYieldSum(mdepp::StatPostProcess::ProcessCluster));
    TRACE(pp.getYieldSum(mdepp::StatPostProcess::ProcessProjectile));
    TRACE(pp.getYieldSum(mdepp::StatPostProcess::ProcessSubstrate));
    TRACE(pp.getYieldSum(mdepp::StatPostProcess::ProcessAll));

    TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessFullerene));
    TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessCluster));
    TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessProjectile));
    TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessSubstrate));
    TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessAll));

    pp.printClassicMoleculesTotal();

    using namespace mdtk;

//  pp.printCoefficients();

//  pp.printClusterDynamicsTotal();

//  pp.spottedTotalByMass();

//    pp.spottedByDepth();

    // pp.buildMassSpectrum();

    yaatk::chdir("..");
  }

  std::vector<size_t> clusterSizes;
  clusterSizes.push_back(0);
  clusterSizes.push_back(13);
  clusterSizes.push_back(27);
  clusterSizes.push_back(39);
  clusterSizes.push_back(75);
//  clusterSizes.push_back(195);
  for(size_t i = 0; i < clusterSizes.size(); i++)
  {
    std::vector<ElementID> elements;
    elements.push_back(DUMMY_EL);
    elements.push_back(Ar_EL);
    elements.push_back(Xe_EL);
    for(size_t j = 0; j < elements.size(); j++)
    {
      std::vector<Float> energies;
      energies.push_back(100*eV);
      energies.push_back(200*eV);
      energies.push_back(400*eV);
      for(size_t k = 0; k < energies.size(); k++)
      {
        std::vector<ElementID> clusterElements;
        clusterElements.push_back(DUMMY_EL);
        clusterElements.push_back(Cu_EL);
        clusterElements.push_back(Au_EL);
        for(size_t l = 0; l < clusterElements.size(); l++)
        {
          // plotEnergyLoss(elements[j], clusterSizes[i], energies[k], clusterElements[l]);

          std::vector<StatPostProcess::FProcessClassicMolecule>
            moleculeFilters;
          moleculeFilters.push_back(mdepp::StatPostProcess::ProcessCluster);
          moleculeFilters.push_back(mdepp::StatPostProcess::ProcessProjectile);
          moleculeFilters.push_back(mdepp::StatPostProcess::ProcessSubstrate);
          moleculeFilters.push_back(mdepp::StatPostProcess::ProcessAll);

          std::vector<StatPostProcess::FProcessClassicMolecule>::const_iterator
            moleculeFilter;
          for(moleculeFilter = moleculeFilters.begin();
              moleculeFilter != moleculeFilters.end();
              ++moleculeFilter)
          {
            std::vector<StatPostProcess::FMoleculeAttribute> moleculeAttributes;
            // moleculeAttributes.push_back(StatPostProcess::moleculeAtomsCount);
            moleculeAttributes.push_back(StatPostProcess::moleculeCount);
            moleculeAttributes.push_back(StatPostProcess::moleculeEnergy);

            std::vector<StatPostProcess::FMoleculeAttribute>::const_iterator
              moleculeAttribute;
            for(moleculeAttribute = moleculeAttributes.begin();
                moleculeAttribute != moleculeAttributes.end();
                ++moleculeAttribute)
            {
              plotAngular(true, *moleculeFilter, *moleculeAttribute,
                          elements[j], clusterSizes[i], clusterElements[l]);
              plotAngular(false, *moleculeFilter, *moleculeAttribute,
                          elements[j], clusterSizes[i], clusterElements[l]);
            }

            plotYieldsAgainstIonEnergy(
              *moleculeFilter,
              elements[j], clusterSizes[i], clusterElements[l]);

            plotMassSpectrum(
              *moleculeFilter,
              elements[j], clusterSizes[i], clusterElements[l]);

            ClassicMolecule mol;
            mol.atoms.push_back(Atom(H_EL));
            Float binWidth = mol.getMass();
            if (*moleculeFilter == StatPostProcess::ProcessCluster ||
                *moleculeFilter == StatPostProcess::ProcessAll)
            {
              mol.atoms.clear();
              mol.atoms.push_back(Atom(clusterElements[l]));
              binWidth = mol.getMass();
            }
            if (*moleculeFilter == StatPostProcess::ProcessProjectile)
            {
              mol.atoms.clear();
              mol.atoms.push_back(Atom(elements[j]));
              binWidth = mol.getMass();
            }
            if (*moleculeFilter == StatPostProcess::ProcessSubstrate)
            {
              mol.atoms.clear();
              mol.atoms.push_back(Atom(C_EL));
              binWidth = mol.getMass();
            }
            plotMassSpectrumHistogram(
              binWidth,
              *moleculeFilter,
              StatPostProcess::moleculeCount,
              elements[j], clusterSizes[i], clusterElements[l]);
            plotMassSpectrumHistogram(
              binWidth,
              *moleculeFilter,
              StatPostProcess::moleculeMass,
              elements[j], clusterSizes[i], clusterElements[l]);
            plotMassSpectrumHistogram(
              binWidth,
              *moleculeFilter,
              StatPostProcess::moleculeEnergy,
              elements[j], clusterSizes[i], clusterElements[l]);

            plotEnergySpectrum(
              *moleculeFilter,
              StatPostProcess::moleculeCount,
              elements[j], clusterSizes[i], clusterElements[l]);
            plotEnergySpectrum(
              *moleculeFilter,
              StatPostProcess::moleculeEnergy,
              elements[j], clusterSizes[i], clusterElements[l]);
          }
        }
      }
    }
  }

  yaatk::chdir("..");

  {
    yaatk::StreamToFileRedirect cout_redir(std::cout,"yieldsout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"yieldserr.txt");

    for(size_t i = 0; i < pps.size(); ++i)
    {
      const mdepp::StatPostProcess pp(pps[i]);

      TRACE(pp.id.str);

      TRACE(pp.trajData.size());

      TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessCluster));
      TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessProjectile));
      TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessSubstrate));
      TRACE(pp.getYieldAverage(mdepp::StatPostProcess::ProcessAll));

      TRACE("--------------------");
    }
  }

  {
    std::ofstream trajdb("trajdb.txt");

    for(size_t i = 0; i < pps.size(); ++i)
    {
      const mdepp::StatPostProcess pp(pps[i]);

      for(size_t trajIndex = 0; trajIndex < pp.trajData.size(); trajIndex++)
      {
        const mdepp::StatPostProcess::TrajData& td = *pp.trajData[trajIndex];
        trajdb << td.trajDir << "\t"
               << pp.getYield(trajIndex, mdepp::StatPostProcess::ProcessCluster) << "\t"
               << pp.getYieldNormalizedByClusterSize(trajIndex, mdepp::StatPostProcess::ProcessCluster) << "\t"
               << pp.getYieldFragmentsCount(trajIndex, mdepp::StatPostProcess::ProcessCluster) << "\t"
               << pp.getYieldFragmentsCountNormalizedByClusterSize(trajIndex, mdepp::StatPostProcess::ProcessCluster) << "\n";
      }
    }
  }

  plot_Ekin_t(mdepp::StatPostProcess::ProcessCluster, 100);
  plot_Ekin_t(mdepp::StatPostProcess::ProcessProjectile, 100);

  plot_Ekin_t(mdepp::StatPostProcess::ProcessCluster, 200);
  plot_Ekin_t(mdepp::StatPostProcess::ProcessProjectile, 200);

  plot_Ekin_t(mdepp::StatPostProcess::ProcessCluster, 400);
  plot_Ekin_t(mdepp::StatPostProcess::ProcessProjectile, 400);
}

void
BatchPostProcess::plot_Ekin_t(mdepp::StatPostProcess::FProcessClassicMolecule fpm, Float ionEnergyFilter) const
{
  std::stringstream fnb;
  fnb << ionEnergyFilter << "eV-";
  fnb << "Ekin_t";

  bool pCluster = (fpm==mdepp::StatPostProcess::ProcessCluster);
  bool pProjectile = (fpm==mdepp::StatPostProcess::ProcessProjectile);

  REQUIRE(pCluster || pProjectile);

  string atomGroup = "";
  if (pCluster)
  {
    atomGroup += "cluster";
    fnb << "_cluster";
  }
  if (pProjectile)
  {
    atomGroup += "projectile";
    fnb << "_projectile";
  }

  ofstream fplt((fnb.str()+".plt").c_str());
  fplt << "\
reset\n\
set xrange [0:10]\n\
set yrange [0:*]\n\
#set xtics (100,200,400)\n\
set pointsize 1.5\n\
#set grid ytics\n\
" << (pProjectile?"#":"") << "set key left top\n\
set key right top\n\
set key spacing 1.5\n\
set xlabel \"time, ps\"\n\
set ylabel \"E_{kin} of " << atomGroup << ", eV\"\n\
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t ppi = 0; ppi < pps.size(); ++ppi)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[ppi][0],1));
    if (id.ionEnergy != ionEnergyFilter)
      continue;

    const mdepp::StatPostProcess pp(pps[ppi]);

    for(size_t trajIndex = 0; trajIndex < pp.trajData.size(); trajIndex++)
    {
      const mdepp::StatPostProcess::TrajData& td = *pp.trajData[trajIndex];

      REQUIRE(td.trajCluster.size() > 0);
      if (!td.trajCluster.rbegin()->second.isMolecule())
        continue;
      if (!isAmongSputtered(td.trajCluster.rbegin()->second,td.molecules))
        continue;

      std::map <Float,Float> singlePlot;
      if (fpm == mdepp::StatPostProcess::ProcessCluster)
        singlePlot = td.plot_Ekin_t(td.getTrajWithPartialSnapshots(td.trajCluster));
      else if (fpm == mdepp::StatPostProcess::ProcessProjectile)
        singlePlot = td.plot_Ekin_t(td.getTrajWithPartialSnapshots(td.trajProjectile));
      else
        throw Exception("plot_Ekin_t : cannot process requested atom groups");

      std::map<Float,Float>::const_iterator i;
      for(i = singlePlot.begin(); i != singlePlot.end(); ++i)
        data << i->first/mdtk::ps << " " << i->second/mdtk::eV << "\n";

      {
        std::ostringstream cmd;
        cmd << "'-' title \"{/Italic "
            << pp.id.ionEnergy << "eV "
            << ElementIDtoString(pp.id.ionElement) << " -> "
            << ElementIDtoString(pp.id.clusterElement)
            << "_{" << pp.id.clusterSize << "}"
            << " " << yaatk::extractItemFromEnd(td.trajDir,1) << "/"
            << yaatk::extractLastItem(td.trajDir)
            << "}\" "
            << " with linespoints pt -1 ";
        plotCmds.push_back(cmd.str());

        data << "e\n";
      }
    }
  }

  if (!(plotCmds.size() > 0))
  {
    TRACE(atomGroup);
  }

//  REQUIRE(plotCmds.size() > 0);
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fplt << plotCmds[i] << ",\\\n";
    else
      fplt << plotCmds[i] << "\n";
  }

  fplt << data.str();

  fplt << "\
set xrange [0:2]\n\
set output  \"" << fnb.str() << "-2ps" << ".eps\"\n\
replot\n\
";

  fplt << data.str();

  fplt.close();
}

void
BatchPostProcess::plotYieldsAgainstIonEnergy(StatPostProcess::FProcessClassicMolecule fpm,
                                             ElementID specIonElement,
                                             size_t specClusterSize,
                                             ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << "yields-" << StatPostProcess::FProcessClassicMoleculeToString(fpm);

  bool pCluster = (fpm==mdepp::StatPostProcess::ProcessCluster);
  bool pProjectile = (fpm==mdepp::StatPostProcess::ProcessProjectile);
  bool pSubstrate = (fpm==mdepp::StatPostProcess::ProcessSubstrate);
  bool pAll = (fpm==mdepp::StatPostProcess::ProcessAll);

  string yieldOfWhat = "розпилення";
  if (pCluster)
    yieldOfWhat = "розпилення кластера";
  if (pSubstrate)
    yieldOfWhat = "розпилення підкладинки";
  if (pProjectile)
    yieldOfWhat = "зворотного розсіювання";

  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
  {
    if (specClusterElement == DUMMY_EL)
      fnb << "_";
    fnb << specClusterSize;
  }

  ofstream fplt((fnb.str()+".plt").c_str());
  fplt << "\
reset\n\
set xrange [0:500]\n\
set yrange [0:*]\n\
set xtics (100,200,400)\n\
set pointsize 1.5\n\
#set grid ytics\n\
" << (pProjectile?"#":"") << "set key left top\n\
#set key right top\n\
set key spacing 1.5\n\
set xlabel \"Енергія іона, еВ\"\n\
set ylabel \"Коефіцієнт " << yieldOfWhat << ", атом/іон\"\n  \
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 8cm, 8cm \"Arial,18\" enhanced\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  std::map < std::string, std::map <Float, Float> > plot;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    const mdepp::StatPostProcess pp(pps[i]);

    std::stringstream expIdWithoutEnergy;
    expIdWithoutEnergy << ElementIDtoString(pp.id.ionElement) << " -> "
                       << ElementIDtoString(pp.id.clusterElement)
                       << "_{" << pp.id.clusterSize << "}";
    plot[expIdWithoutEnergy.str()][pp.id.ionEnergy] = pp.getYieldAverage(fpm);
  }

  std::map < std::string, std::map <Float, Float> >::const_iterator expId;
  for(expId = plot.begin(); expId != plot.end(); ++expId)
  {
    std::map <Float, Float>::const_iterator dataPoint;
    for(dataPoint = expId->second.begin(); dataPoint != expId->second.end(); ++dataPoint)
    {
      data << dataPoint->first << " " << dataPoint->second << "\n";
    }

    {
      std::ostringstream cmd;
      cmd << "'-' title \"{/Italic "
          << expId->first
          << "}\" "
          << "with linespoints";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

  if (!(plotCmds.size() > 0))
  {
    TRACE(yieldOfWhat);
    TRACE(StatPostProcess::FProcessClassicMoleculeToString(fpm));
    TRACE(specIonElement);
    TRACE(specClusterSize);
    TRACE(specClusterElement);
  }

//  REQUIRE(plotCmds.size() > 0);
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

void
BatchPostProcess::plotMassSpectrum(StatPostProcess::FProcessClassicMolecule fpm,
                                   ElementID specIonElement,
                                   size_t specClusterSize,
                                   ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << "MassSpectrum-" << StatPostProcess::FProcessClassicMoleculeToString(fpm);
/*
  bool pCluster = (fpm==mdepp::StatPostProcess::ProcessCluster);
  bool pProjectile = (fpm==mdepp::StatPostProcess::ProcessProjectile);
  bool pSubstrate = (fpm==mdepp::StatPostProcess::ProcessSubstrate);
  bool pAll = (fpm==mdepp::StatPostProcess::ProcessAll);

  string yieldOfWhat = "розпилення";
  if (pCluster)
  yieldOfWhat = "розпилення кластера";
  if (pSubstrate)
  yieldOfWhat = "розпилення підкладинки";
  if (pProjectile)
  yieldOfWhat = "зворотного розсіювання";
*/
  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
  {
    if (specClusterElement == DUMMY_EL)
      fnb << "_";
    fnb << specClusterSize;
  }

  std::ofstream fsameMass((fnb.str()+".sameMass").c_str());
  std::ofstream fo((fnb.str()+".txt").c_str());
  std::ofstream fplt((fnb.str()+".plt").c_str());
  std::ofstream fpltmulti((fnb.str()+"-multi.plt").c_str());

  fplt << "\
reset\n\
set yrange [0:*]\n\
set xrange [0:*]\n\
set border 1+2+4+8 lw 2\n\
\n\
set output \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
set xlabel \"Mass (amu)\"\n\
set ylabel ""\n\
\n\
#set label \"H\" at 2,279/1000.0+0.012\n\
#set label \"Cu_2\" at 63.5*2,0.43\n\
#set label \"Cu_1_3\" at 63.5*13,0.07\n\
#set label \"Cu_2_5\" at 63.5*25,0.05\n\
#set label \"Cu_3_5\" at 63.5*35,0.03\n\
#set label \"H_2\" at 2,83/1000.0+0.012\n\
#set label \"CH_3\" at 15,10/1000.0+0.012\n\
#set label \"C_2H_2\" at 26,27/1000.0+0.012 center\n\
#set label \"C_2H_4\" at 28,241/1000.0+0.012 center\n\
#set label \"C_3H_6\" at 42,16/1000.0+0.012 center\n\
#set label \"C_4H_7\" at 55,10/1000.0+0.012 center\n\
#set label \"C_6H_1_0\" at 82,2/1000.0+0.012 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,2/1000.0+0.03 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,(7+6+1)/500.0 center\n\
\n\
set xtics nomirror 200\n\
set tics scale -1\n\
\n\
plot \\\n\
";

  fpltmulti << "\
reset\n\
base_mass=64\n\
xrange_cut=2.5*base_mass\n\
set yrange [0:*]\n\
set border 1+2+4+8 lw 2\n\
\n\
set encoding koi8u\n\
set output \"" << fnb.str() << "-multi" << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
\n\
set multiplot\n\
set size 0.25,1\n\
set origin 0.0,0.0\n\
#set lmargin 10\n\
set rmargin 0\n\
set xrange [0:xrange_cut]\n\
set xlabel \" \"\n\
set ylabel \"Sputtering yield, specie/impact\"\n\
#set label \"Cu\" at 63.5+10,2.1\n\
#set label \"Cu_{2}\" at 63.5*2-20,0.40\n\
\n\
set key off\n\
set xtics nomirror 0,base_mass,2*base_mass\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    if (id.ionEnergy > 199.0 && id.ionEnergy < 201.0)
      continue;

    const mdepp::StatPostProcess pp(pps[i]);

    {
      std::ostringstream cmd;
      cmd << "'-' with points title \"{/Italic "
          << pp.id.ionEnergy << "eV "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\", "
          << "'-' with impulses notitle";
      plotCmds.push_back(cmd.str());
    }

    {
      fo << "#"
         << pp.id.ionEnergy << "eV "
         << ElementIDtoString(pp.id.ionElement) << " -> "
         << ElementIDtoString(pp.id.clusterElement)
         << "_{" << pp.id.clusterSize << "}\n";
      fsameMass << "#"
                << pp.id.ionEnergy << "eV "
                << ElementIDtoString(pp.id.ionElement) << " -> "
                << ElementIDtoString(pp.id.clusterElement)
                << "_{" << pp.id.clusterSize << "}\n";
    }

    std::map<ClassicMolecule, size_t> massSpectrum = pp.buildMassSpectrum(fpm);

    for(size_t c = 0; c < 2; ++c)
    {
      std::map<ClassicMolecule, size_t>::iterator i = massSpectrum.begin();
      size_t previousAMUMass = 0;
      while (i != massSpectrum.end())
      {
        size_t AMUMass = i->first.getAMUMass();

        fo << i->first.formula()
           << " (" << AMUMass << " amu) : "
           << i->second << "/" << pp.trajData.size()
           << " = " << double(i->second)/pp.trajData.size() << "\n";

        data << AMUMass << " " << double(i->second)/pp.trajData.size() << "\n";

        if (AMUMass == previousAMUMass)
        {
          fsameMass << i->first.formula()
                    << " (" << AMUMass << " amu) : "
                    << i->second << "/" << pp.trajData.size()
                    << " = " << double(i->second)/pp.trajData.size() << "\n";
        }

        previousAMUMass = AMUMass;

        i++;
      }
      data << "e\n";
    }
  }

  if (!(plotCmds.size() > 0))
  {
//    TRACE(yieldOfWhat);
    TRACE(StatPostProcess::FProcessClassicMoleculeToString(fpm));
    TRACE(specIonElement);
    TRACE(specClusterSize);
    TRACE(specClusterElement);
  }

//  REQUIRE(plotCmds.size() > 0);
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fplt << plotCmds[i] << ",\\\n";
    else
      fplt << plotCmds[i] << "\n";
  }

  fplt << data.str();

  fpltmulti << "plot \\\n";
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fpltmulti << plotCmds[i] << ",\\\n";
    else
      fpltmulti << plotCmds[i] << "\n";
  }

  fpltmulti << data.str();

  fpltmulti << "\
set size 0.7,1\n\
set origin 0.25,0.0\n\
set format y \"\"\n\
#set format y2 \"%7g\"\n\
set lmargin 0\n\
set rmargin 2\n\
#set nolog x\n\
set xrange [xrange_cut:*]\n\
#set xtic 0,10\n\
#set mxtic 5\n\
\n\
set xlabel \"Mass, a.m.u.\"\n\
set ylabel \"\"\n\
\n\
#set label \"H\" at 2,279/1000.0+0.012\n\
#set label \"Cu_{3}\" at 63.5*3+30,0.085\n\
#set label \"Cu_{11}\" at 63.5*11,0.040\n\
#set label \"Cu_{15}\" at 63.5*15,0.040\n\
#set label \"Cu_{21}\" at 63.5*21+25,0.076\n\
#set label \"Cu_{23}\" at 63.5*23-20,0.059\n\
#set label \"Cu_{25}\" at 63.5*25,0.059\n\
#set label \"Cu_{27}\" at 63.5*27,0.021\n\
#set label \"Cu_{35}\" at 63.5*35,0.055\n\
#set label \"H_2\" at 2,83/1000.0+0.012\n\
#set label \"CH_3\" at 15,10/1000.0+0.012\n\
#set label \"C_2H_2\" at 26,27/1000.0+0.012 center\n\
#set label \"C_2H_4\" at 28,241/1000.0+0.012 center\n\
#set label \"C_3H_6\" at 42,16/1000.0+0.012 center\n\
#set label \"C_4H_7\" at 55,10/1000.0+0.012 center\n\
#set label \"C_6H_1_0\" at 82,2/1000.0+0.012 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,2/1000.0+0.03 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,(7+6+1)/500.0 center\n\
\n\
set xtics nomirror 5*base_mass\n\
set xtics add (sprintf(\"%d\",3*base_mass) 3*base_mass)\n\
\n\
set y2tics\n\
\n\
set key samplen 1.0 spacing 1.3\n\
#set key left samplen 1.0 spacing 1.3\n\
";

  fpltmulti << "plot \\\n";
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fpltmulti << plotCmds[i] << ",\\\n";
    else
      fpltmulti << plotCmds[i] << "\n";
  }

  fpltmulti << data.str();

  fpltmulti << "set nomultiplot\n";

  fpltmulti.close();
  fplt.close();
  fo.close();
  fsameMass.close();
}

void
BatchPostProcess::plotMassSpectrumHistogram(
  Float binWidth,
  StatPostProcess::FProcessClassicMolecule fpm,
  StatPostProcess::FMoleculeAttribute fma,
  ElementID specIonElement,
  size_t specClusterSize,
  ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << "MassSpectrum-" << StatPostProcess::FProcessClassicMoleculeToString(fpm);

  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
  {
    if (specClusterElement == DUMMY_EL)
      fnb << "_";
    fnb << specClusterSize;
  }

  fnb << "-histogram";

  fnb << "-";

  fnb << StatPostProcess::FMoleculeAttributeToString(fma);

  std::ofstream fplt((fnb.str()+".plt").c_str());
  std::ofstream fpltmulti((fnb.str()+"-multi.plt").c_str());

  fplt << "\
reset\n\
set yrange [0:*]\n\
set xrange [0:*]\n\
set border 1+2+4+8 lw 2\n\
\n\
set output \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
set xlabel \"Mass (amu)\"\n\
set ylabel ""\n\
\n\
#set label \"H\" at 2,279/1000.0+0.012\n\
#set label \"Cu_2\" at 63.5*2,0.43\n\
#set label \"Cu_1_3\" at 63.5*13,0.07\n\
#set label \"Cu_2_5\" at 63.5*25,0.05\n\
#set label \"Cu_3_5\" at 63.5*35,0.03\n\
#set label \"H_2\" at 2,83/1000.0+0.012\n\
#set label \"CH_3\" at 15,10/1000.0+0.012\n\
#set label \"C_2H_2\" at 26,27/1000.0+0.012 center\n\
#set label \"C_2H_4\" at 28,241/1000.0+0.012 center\n\
#set label \"C_3H_6\" at 42,16/1000.0+0.012 center\n\
#set label \"C_4H_7\" at 55,10/1000.0+0.012 center\n\
#set label \"C_6H_1_0\" at 82,2/1000.0+0.012 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,2/1000.0+0.03 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,(7+6+1)/500.0 center\n\
\n\
set xtics nomirror 200\n\
set tics scale -1\n\
\n\
plot \\\n\
";

  fpltmulti << "\
reset\n\
base_mass=64\n\
xrange_cut=2.5*base_mass\n\
set yrange [0:*]\n\
set border 1+2+4+8 lw 2\n\
\n\
set encoding koi8u\n\
set output \"" << fnb.str() << "-multi" << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
\n\
set multiplot\n\
set size 0.25,1\n\
set origin 0.0,0.0\n\
#set lmargin 10\n\
set rmargin 0\n\
set xrange [0:xrange_cut]\n\
set xlabel \" \"\n\
set ylabel \"Sputtering yield, specie/impact\"\n\
#set label \"Cu\" at 63.5+10,2.1\n\
#set label \"Cu_{2}\" at 63.5*2-20,0.40\n\
\n\
set key off\n\
set xtics nomirror 0,base_mass,2*base_mass\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    if (id.ionEnergy > 199.0 && id.ionEnergy < 201.0)
      continue;

    const mdepp::StatPostProcess pp(pps[i]);

    {
      std::ostringstream cmd;
      cmd << "'-' with points title \"{/Italic "
          << pp.id.ionEnergy << "eV "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\", "
          << "'-' with boxes notitle";
      plotCmds.push_back(cmd.str());
    }

    Float maxMass = 0.0;
    {
      for(size_t trajIndex = 0; trajIndex < pp.trajData.size(); trajIndex++)
      {
        const StatPostProcess::TrajData& td = *(pp.trajData[trajIndex]);
        for(size_t mi = 0; mi < td.molecules.size(); mi++)
        {
          const ClassicMolecule& mol = td.molecules[mi];
          if (!fpm(mol)) continue;
          Float mass = mol.getMass();
          if (mass > maxMass)
            maxMass = mass;
        }
      }
    }

    std::map<Float, Float> histData;
    {
      const Float minHist = 0.0;
      const int n = (maxMass - minHist)/binWidth + 1;
      double range[n+1];
      range[0] = 0.0;
      for(size_t massIndex = 0; massIndex < n; ++massIndex)
        range[massIndex + 1] = range[massIndex] + binWidth;

      gsl_histogram * h = gsl_histogram_alloc (n);

      gsl_histogram_set_ranges(h, range, n+1);

      for(size_t trajIndex = 0; trajIndex < pp.trajData.size(); trajIndex++)
      {
        const StatPostProcess::TrajData& td = *(pp.trajData[trajIndex]);
        for(size_t mi = 0; mi < td.molecules.size(); mi++)
        {
          const ClassicMolecule& mol = td.molecules[mi];
          if (!fpm(mol)) continue;
          gsl_histogram_accumulate(
            h,
            mol.getMass(),
            fma(mol)/Float(pp.trajData.size()));
        }
      }

      for(size_t i = 0; i < gsl_histogram_bins(h); i++)
      {
        double lower, upper;
        gsl_histogram_get_range (h, i, &lower, &upper);
        histData[(lower+upper)/2.0] = gsl_histogram_get(h,i);
      }

      gsl_histogram_free (h);
    }

    for(size_t c = 0; c < 2; ++c)
    {
      std::map<Float, Float>::iterator i = histData.begin();
      while (i != histData.end())
      {
        data << i->first/mdtk::amu << " " << i->second << "\n";
        i++;
      }
      data << "e\n";
    }
  }

  if (!(plotCmds.size() > 0))
  {
    TRACE(StatPostProcess::FProcessClassicMoleculeToString(fpm));
    TRACE(specIonElement);
    TRACE(specClusterSize);
    TRACE(specClusterElement);
  }

//  REQUIRE(plotCmds.size() > 0);
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fplt << plotCmds[i] << ",\\\n";
    else
      fplt << plotCmds[i] << "\n";
  }

  fplt << data.str();

  fpltmulti << "plot \\\n";
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fpltmulti << plotCmds[i] << ",\\\n";
    else
      fpltmulti << plotCmds[i] << "\n";
  }

  fpltmulti << data.str();

  fpltmulti << "\
set size 0.7,1\n\
set origin 0.25,0.0\n\
set format y \"\"\n\
#set format y2 \"%7g\"\n\
set lmargin 0\n\
set rmargin 2\n\
#set nolog x\n\
set xrange [xrange_cut:*]\n\
#set xtic 0,10\n\
#set mxtic 5\n\
\n\
set xlabel \"Mass, a.m.u.\"\n\
set ylabel \"\"\n\
\n\
#set label \"H\" at 2,279/1000.0+0.012\n\
#set label \"Cu_{3}\" at 63.5*3+30,0.085\n\
#set label \"Cu_{11}\" at 63.5*11,0.040\n\
#set label \"Cu_{15}\" at 63.5*15,0.040\n\
#set label \"Cu_{21}\" at 63.5*21+25,0.076\n\
#set label \"Cu_{23}\" at 63.5*23-20,0.059\n\
#set label \"Cu_{25}\" at 63.5*25,0.059\n\
#set label \"Cu_{27}\" at 63.5*27,0.021\n\
#set label \"Cu_{35}\" at 63.5*35,0.055\n\
#set label \"H_2\" at 2,83/1000.0+0.012\n\
#set label \"CH_3\" at 15,10/1000.0+0.012\n\
#set label \"C_2H_2\" at 26,27/1000.0+0.012 center\n\
#set label \"C_2H_4\" at 28,241/1000.0+0.012 center\n\
#set label \"C_3H_6\" at 42,16/1000.0+0.012 center\n\
#set label \"C_4H_7\" at 55,10/1000.0+0.012 center\n\
#set label \"C_6H_1_0\" at 82,2/1000.0+0.012 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,2/1000.0+0.03 center\n\
#set label \"C_6H_1_0\\n(cyclohexene)\" at 82,(7+6+1)/500.0 center\n\
\n\
set xtics nomirror 5*base_mass\n\
set xtics add (sprintf(\"%d\",3*base_mass) 3*base_mass)\n\
\n\
set y2tics\n\
\n\
set key samplen 1.0 spacing 1.3\n\
#set key left samplen 1.0 spacing 1.3\n\
";

  fpltmulti << "plot \\\n";
  for(size_t i = 0; i < plotCmds.size(); ++i)
  {
    if (i != plotCmds.size()-1)
      fpltmulti << plotCmds[i] << ",\\\n";
    else
      fpltmulti << plotCmds[i] << "\n";
  }

  fpltmulti << data.str();

  fpltmulti << "set nomultiplot\n";

  fpltmulti.close();
  fplt.close();
}

void
BatchPostProcess::plotEnergySpectrum(
  StatPostProcess::FProcessClassicMolecule fpm,
  StatPostProcess::FMoleculeAttribute fma,
  ElementID specIonElement,
  size_t specClusterSize,
  ElementID specClusterElement) const
{
  std::stringstream fnb;

  fnb << "EnergySpectrum-";

  fnb << StatPostProcess::FProcessClassicMoleculeToString(fpm);

  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
  {
    if (specClusterElement == DUMMY_EL)
      fnb << "_";
    fnb << specClusterSize;
  }

  fnb << "-";

  fnb << StatPostProcess::FMoleculeAttributeToString(fma);

  ofstream fplt((fnb.str()+".plt").c_str());

  fplt << "\
reset\n\
set yrange [0:*]\n\
set xrange [0:*]\n\
set border 1+2+4+8 lw 2\n\
\n\
set output \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 16cm, 8cm \"Arial,18\" enhanced\n\
set xlabel \"Mass (amu)\"\n\
set ylabel ""\n\
\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    const mdepp::StatPostProcess pp(pps[i]);

    {
      std::map<Float, Float> histData =
        pp.distBy(
          StatPostProcess::moleculeEnergy,
          pp.suggestedBinWidth(
            StatPostProcess::moleculeEnergy, *fpm,
            specIonElement, specClusterSize, specClusterElement),
          0.0,
          pp.maxMoleculeAttribute(
            StatPostProcess::moleculeEnergy, *fpm),
          fma,
          fpm);

      std::map<Float, Float>::iterator i;
      for(i = histData.begin(); i != histData.end(); ++i)
        data << i->first << " " << histData[i->first] << "\n";
    }

    {
      std::ostringstream cmd;
      cmd << "'-' with boxes title \"{/Italic "
          << pp.id.ionEnergy << "eV "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\" "
          << "";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

//  REQUIRE(plotCmds.size() > 0);
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

void
BatchPostProcess::plotAngular(bool plotPolar,
                              StatPostProcess::FProcessClassicMolecule fpm,
                              StatPostProcess::FMoleculeAttribute fma,
                              ElementID specIonElement,
                              size_t specClusterSize,
                              ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << StatPostProcess::FProcessClassicMoleculeToString(fpm);

  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
  {
    if (specClusterElement == DUMMY_EL)
      fnb << "_";
    fnb << specClusterSize;
  }

  fnb << (plotPolar?"-polar":"-azimuthal");

  fnb << "-";

  fnb << StatPostProcess::FMoleculeAttributeToString(fma);

  ofstream fplt((fnb.str()+".plt").c_str());

  size_t n = plotPolar?3:12;
  char n_str[1024];
  sprintf(n_str,"%05lu",(long unsigned)n);

  if (plotPolar)
    fplt << "\
reset\n\
set xrange [-0.5:" << n-0.5 << "]\n\
set yrange [0:*]\n\
#set xtics 0,15,90\n\
set pointsize 1.5\n\
#set grid ytics\n\
set key left top\n\
#set key right top\n\
set key spacing 1.5\n\
set xlabel \"Полярний кут, градуси\"\n\
set ylabel \"{/Italic dY}_{x}, атом/іон\"\n\
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 8cm, 8cm \"Arial,18\" enhanced\n\
set style data histogram\n\
set boxwidth 0.8\n\
set style histogram cluster gap 1\n\
set style fill pattern border\n\
plot \\\n\
";
  else
    fplt << "\
reset\n\
set xrange [-0.5:" << n-0.5 << "]\n\
set yrange [0:*]\n\
#set xtics -180,45,180\n\
set pointsize 1.5\n\
#set grid ytics\n\
#set key left top\n\
#set key right top\n\
set key spacing 1.5\n\
set xlabel \"Азимутальний кут, градуси\"\n\
set ylabel \"{/Italic dY}_{x}, атом/іон\"\n\
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 8cm, 8cm \"Arial,18\" enhanced\n\
set style data histogram\n\
set boxwidth 0.8\n\
set style histogram cluster gap 1\n\
set style fill pattern border\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    const mdepp::StatPostProcess pp(pps[i]);

    {
      std::map<Float, Float> histData =
        pp.distByAngle(
          plotPolar ? StatPostProcess::ANGLE_POLAR : StatPostProcess::ANGLE_AZIMUTH,
          n,
          fma,
          fpm);

      std::map<Float, Float>::iterator i;
      for(i = histData.begin(); i != histData.end(); ++i)
        data << i->first << " " << histData[i->first] << "\n";
    }

    {
      std::ostringstream cmd;
      cmd << "'-' using 2:xticlabels(1) lt -1 title \"{/Italic "
          << pp.id.ionEnergy << "еВ "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\" "
          << "";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

//  REQUIRE(plotCmds.size() > 0);
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

Float Ek(ElementID id, SnapshotList::AtomSnapshot as)
{
  Atom a(id);
//  a.setAttributesByElementID();
  as.restoreToAtom(a);
  return a.M*SQR(a.V.module())/2.0;
}

void
BatchPostProcess::plotEnergyLoss(ElementID specIonElement,
                                 size_t specClusterSize,
                                 Float specIonEnergy,
                                 ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << "energy-loss";

  REQUIRE(specIonElement != DUMMY_EL);
  REQUIRE(specClusterSize > 0);
  REQUIRE(specIonEnergy > 0);
  REQUIRE(specClusterElement != DUMMY_EL);

  if (specIonElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specIonElement);

  if (specClusterElement != DUMMY_EL)
    fnb << "_" << ElementIDtoString(specClusterElement);

  if (specClusterSize > 0)
    fnb << "" << specClusterSize;

  if (specIonEnergy > 0)
    fnb << "_" << int(specIonEnergy/eV) << "eV";

  ofstream fplt((fnb.str()+".plt").c_str());

  std::vector<Float> bounds;

  bounds.push_back(-1000000.0*Ao);
  const Float c = 2.547*Ao;
  for(int i = 0; i < 11; ++i)
    bounds.push_back(c*i);
  bounds.push_back(+1000000.0*Ao);

  fplt << "# 1st bound = " << bounds[0]/Ao << " Ao" << "\n"
       << "# last depth = " << *(bounds.end()-1)/Ao << " Ao" << "\n"
       << "# number of bounds = " << bounds.size() << "\n"
       << "# number of bins = " << bounds.size()-1 << "\n";

  fplt << "\
reset\n\
set xrange [" << -c*1/Ao << ":" << c*9/Ao << "]\n\
#set yrange [0:*]\n\
set format x \"%.1f\"\n\
set xtics 0," << c*2/Ao << "\n\
#set grid xtics\n\
set pointsize 1.5\n\
#set grid ytics\n\
#set key left top\n\
#set key right top\n\
set key spacing 1.5\n\
set xlabel \"Глибина, ангстрем\"\n\
set ylabel \"Втрати енергії у шарі, еВ\"\n\
set encoding koi8u\n\
set output  \"" << fnb.str() << ".eps\"\n\
set terminal postscript eps size 8cm, 8cm \"Arial,18\" enhanced\n\
set style fill solid 0.5 border -1\n\
plot \\\n\
";

  std::vector<std::string> plotCmds;
  std::ostringstream data;

  for(size_t i = 0; i < pps.size(); ++i)
  {
    mdepp::StatPostProcess::Id id(yaatk::extractItemFromEnd(pps[i][0],1));

    if (specIonElement != DUMMY_EL)
    {
      if (id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (id.clusterSize != specClusterSize)
        continue;
    }

    if (specIonEnergy > 0)
    {
      if (fabs(id.ionEnergy*eV - specIonEnergy) > 0.05*eV)
        continue;
    }

    const mdepp::StatPostProcess pp(pps[i]);

    std::vector<Float> dEs;
    dEs.resize(bounds.size()-1);

    for(size_t ti = 0; ti < pp.trajData.size(); ti++)
    {
      const StatPostProcess::TrajData& td = *pp.trajData[ti];

      SnapshotList sn;
      {
        std::string cwd = yaatk::getcwd();
        yaatk::chdir("..");
        yaatk::chdir(td.trajDir.c_str());
        sn.loadstate();
        yaatk::chdir(cwd.c_str());
      }

      size_t dEindex = 0;
      Float prevEk;
      for(size_t shotIndex = 0; shotIndex < sn.snapshots.size(); ++shotIndex)
      {
        Float t = sn.snapshots[shotIndex].first;
        SnapshotList::SelectedAtomSnapshotList& asl
          = sn.snapshots[shotIndex].second;
        size_t ionIndex = asl.size()-1;

        Float depth = asl[ionIndex].pos[2];

        REQUIRE(asl.size() > 0);
        REQUIRE(shotIndex != 0 || depth < -3.0*mdtk::Ao);

        REQUIRE(dEindex < dEs.size());

        if (shotIndex == 0)
        {
          prevEk = Ek(pp.id.ionElement,asl[ionIndex]);
          TRACE(prevEk/eV);
        }

        size_t dEindex_next = dEindex;

        if (depth >= bounds[dEindex+1])
          dEindex_next = dEindex+1;
        if (depth < bounds[dEindex])
          dEindex_next = dEindex-1;

        bool transitionDetected = (dEindex != dEindex_next);

        if (transitionDetected || shotIndex == sn.snapshots.size()-1)
        {
          Float curEk = Ek(pp.id.ionElement,
                           asl[ionIndex]);

          dEs[dEindex] += curEk - prevEk;

          TRACE(shotIndex);
          TRACE(t/fs);
          TRACE(dEindex);
          TRACE((curEk - prevEk)/eV);
          TRACE(dEs[dEindex]/eV);

          prevEk = curEk;
          dEindex = dEindex_next;
          REQUIRE(/*dEindex >= 0 &&*/ dEindex < dEs.size());
        }
      }
      {
        Float dEsum = 0;
        for(int i = 0; i < dEs.size(); i++)
          dEsum += dEs[i];

        SnapshotList::SelectedAtomSnapshotList& asl_start
          = sn.snapshots[0].second;
        SnapshotList::SelectedAtomSnapshotList& asl_end
          = sn.snapshots[sn.snapshots.size()-1].second;

        size_t ionIndex = asl_start.size()-1;

        Float Ek_start = Ek(pp.id.ionElement,asl_start[ionIndex]);
        Float Ek_end   = Ek(pp.id.ionElement,asl_end[ionIndex]);

        TRACE(td.trajDir);

        TRACE(Ek_start/eV);
        TRACE(Ek_end/eV);
        TRACE((Ek_end-Ek_start)/eV);
        TRACE(dEsum/eV);
      }
    }

    for(int i = 0; i < dEs.size(); i++)
    {
      Float w = bounds[i+1]-bounds[i];
      Float x = bounds[i]+w/2;
      if (i == 0)
      {
        w = c;
        x = (bounds[i+1]-w/2);
      }
      if (i == dEs.size()-1)
      {
        w = c;
        x = (bounds[i]  +w/2);
      }
      data << x/Ao << " "
           << -dEs[i]/pp.trajData.size()/eV << " "
           << w/Ao << "\n";
    }

    {
      std::ostringstream cmd;
      cmd << "'-' title \"{/Italic "
          << pp.id.ionEnergy << "еВ "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\" "
          << "with boxes";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

//  REQUIRE(plotCmds.size() > 0);
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
