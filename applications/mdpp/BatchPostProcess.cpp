/* 
   Molecular dynamics postprocessor, BatchPostProcess  class

   Copyright (C) 2007, 2008, 2009, 2010, 2011, 2012, 2014 Oleksandr
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

BatchPostProcess::BatchPostProcess(std::string mdeppinPath, size_t haloIndex)
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
      pps.push_back(mdepp::StatPostProcess(trajsetDir));
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
    pps[i].saveToStream(os);
}

void
BatchPostProcess::loadFromStream(std::istream& is)
{
  size_t sz;
  is >> sz;
  for(size_t i = 0; i < sz; i++)
  {
    mdepp::StatPostProcess pp;
    pp.loadFromStream(is);
    pps.push_back(pp);
  }
}

void
BatchPostProcess::execute()
{
  for(size_t i = 0; i < pps.size(); ++i)
    pps[i].execute();
}

void
BatchPostProcess::addHalo(const mdepp::BatchPostProcess& bppHalo)
{
  REQUIRE(pps.size() != 0);
  for(size_t i = 0; i < pps.size(); ++i)
  {
    size_t haloMatches = 0;
    for(size_t j = 0; j < bppHalo.pps.size(); ++j)
      if (pps[i].id.str == bppHalo.pps[j].id.str)
      {
        ++haloMatches;
        pps[i].addHalo(bppHalo.pps[j]);
      }
    REQUIRE(haloMatches <= 1);
  }
}

void
BatchPostProcess::printResults() const
{
  yaatk::mkdir("results");
  yaatk::chdir("results");

  for(size_t i = 0; i < pps.size(); ++i)
  {
    const mdepp::StatPostProcess& pp = pps[i];

    yaatk::mkdir(pp.id.str.c_str());
    yaatk::chdir(pp.id.str.c_str());

    yaatk::StreamToFileRedirect cout_redir(std::cout,"stdout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"stderr.txt");

    for(size_t i = 0; i < pp.trajData.size(); i++)
    {
      TRACE(pp.trajData[i].trajDir);

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

    TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessFullerene));
    TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessCluster));
    TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessProjectile));
    TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessSubstrate));
    TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessAll));

    pp.printClassicMoleculesTotal();

    using namespace mdtk;

//  pp.printCoefficients();

//  pp.printClusterDynamicsTotal();

//  pp.spottedTotalByMass();

//    pp.spottedByDepth();

    pp.buildMassSpectrum();

#define MDPP_PROCESS_ONLY(FPM) \
    {\
      std::string s = "Process"#FPM;\
      yaatk::mkdir(s.c_str());\
      yaatk::chdir(s.c_str());\
      pp.buildMassSpectrum(mdepp::StatPostProcess::Process##FPM);\
      pp.buildAngular(mdepp::StatPostProcess::Process##FPM); \
      yaatk::chdir("..");\
    }

    MDPP_PROCESS_ONLY(Cluster);
    MDPP_PROCESS_ONLY(Projectile);
    MDPP_PROCESS_ONLY(Substrate);
    MDPP_PROCESS_ONLY(All);

    yaatk::chdir("..");
  }

  {
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessCluster,"yields-Cluster");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessProjectile,"yields-Projectile");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessSubstrate,"yields-Substrate");
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessAll,"yields-All");
  }


  std::vector<ElementID> elements;
  elements.push_back(Ar_EL);
  elements.push_back(Xe_EL);
  for(size_t i = 0; i < elements.size(); i++)
  {
//    plotEnergyLoss(elements[i]);

    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessCluster,
                               "yields-Cluster", elements[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessProjectile,
                               "yields-Projectile", elements[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessSubstrate,
                               "yields-Substrate", elements[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessAll,
                               "yields-All", elements[i]);

#define PLOT_ANGULAR_ION_TYPE(angleType)\
    plotAngular(angleType,mdepp::StatPostProcess::ProcessCluster,\
                "Cluster", elements[i]);\
    plotAngular(angleType,mdepp::StatPostProcess::ProcessProjectile,\
                "Projectile", elements[i]);\
    plotAngular(angleType,mdepp::StatPostProcess::ProcessSubstrate,\
                "Substrate", elements[i]);\
    plotAngular(angleType,mdepp::StatPostProcess::ProcessAll,\
                "All", elements[i]);

//    PLOT_ANGULAR_ION_TYPE(true);
//    PLOT_ANGULAR_ION_TYPE(false);
  }

  std::vector<size_t> clusterSizes;
  clusterSizes.push_back(13);
  clusterSizes.push_back(27);
  clusterSizes.push_back(39);
  clusterSizes.push_back(75);
//  clusterSizes.push_back(195);
  for(size_t i = 0; i < clusterSizes.size(); i++)
  {
//    plotEnergyLoss(DUMMY_EL, clusterSizes[i]);

    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessCluster,
                               "yields-Cluster", DUMMY_EL, clusterSizes[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessProjectile,
                               "yields-Projectile", DUMMY_EL, clusterSizes[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessSubstrate,
                               "yields-Substrate", DUMMY_EL, clusterSizes[i]);
    plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessAll,
                               "yields-All", DUMMY_EL, clusterSizes[i]);

#define PLOT_ANGULAR(angleType,elementType,clusterSize,clusterElement)  \
    plotAngular(angleType,mdepp::StatPostProcess::ProcessCluster,       \
                "Cluster", elementType, clusterSize, clusterElement);   \
    plotAngular(angleType,mdepp::StatPostProcess::ProcessProjectile,    \
                "Projectile", elementType, clusterSize, clusterElement); \
    plotAngular(angleType,mdepp::StatPostProcess::ProcessSubstrate,     \
                "Substrate", elementType, clusterSize, clusterElement); \
    plotAngular(angleType,mdepp::StatPostProcess::ProcessAll,           \
                "All", elementType, clusterSize, clusterElement);

//    PLOT_ANGULAR(true,DUMMY_EL,clusterSizes[i]);
//    PLOT_ANGULAR(false,DUMMY_EL,clusterSizes[i]);

    std::vector<ElementID> elements;
    elements.push_back(Ar_EL);
    elements.push_back(Xe_EL);
    for(size_t j = 0; j < elements.size(); j++)
    {
//      plotEnergyLoss(elements[j], clusterSizes[i]);

      std::vector<Float> energies;
      energies.push_back(100*eV);
      energies.push_back(200*eV);
      energies.push_back(400*eV);
      for(size_t k = 0; k < energies.size(); k++)
      {
        std::vector<ElementID> clusterElements;
        clusterElements.push_back(Cu_EL);
        clusterElements.push_back(Au_EL);
        for(size_t l = 0; l < clusterElements.size(); l++)
        {
          plotEnergyLoss(elements[j], clusterSizes[i], energies[k], clusterElements[l]);
          PLOT_ANGULAR(true,elements[j],clusterSizes[i], clusterElements[l]);
          PLOT_ANGULAR(false,elements[j],clusterSizes[i], clusterElements[l]);

          plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessCluster,
                                     "yields-Cluster", elements[j],0,clusterElements[l]);
          plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessProjectile,
                                     "yields-Projectile", elements[j],0,clusterElements[l]);
          plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessSubstrate,
                                     "yields-Substrate", elements[j],0,clusterElements[l]);
          plotYieldsAgainstIonEnergy(mdepp::StatPostProcess::ProcessAll,
                                     "yields-All", elements[j],0,clusterElements[l]);
        }
      }

      PLOT_ANGULAR(true,elements[j],clusterSizes[i],DUMMY_EL);
      PLOT_ANGULAR(false,elements[j],clusterSizes[i],DUMMY_EL);
    }
  }

  yaatk::chdir("..");

  {
    yaatk::StreamToFileRedirect cout_redir(std::cout,"yieldsout.txt");
    yaatk::StreamToFileRedirect cerr_redir(std::cerr,"yieldserr.txt");

    for(size_t i = 0; i < pps.size(); ++i)
    {
      const mdepp::StatPostProcess& pp = pps[i];

      TRACE(pp.id.str);

      TRACE(pp.trajData.size());

      TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessCluster));
      TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessProjectile));
      TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessSubstrate));
      TRACE(pp.getAverageYield(mdepp::StatPostProcess::ProcessAll));

      TRACE("--------------------");
    }
  }
}

void
BatchPostProcess::plotYieldsAgainstIonEnergy(StatPostProcess::FProcessClassicMolecule fpm,
                                             std::string idStr,
                                             ElementID specIonElement,
                                             size_t specClusterSize,
                                             ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << idStr;

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

  for(size_t i = 0; i < pps.size(); ++i)
  {
    const mdepp::StatPostProcess& pp = pps[i];

    if (specIonElement != DUMMY_EL)
    {
      if (pp.id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (pp.id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (pp.id.clusterSize != specClusterSize)
        continue;
    }

    data << pp.id.ionEnergy << " " << pp.getAverageYield(fpm) << "\n";

    if (i%3 == 2)
    {
      std::ostringstream cmd;
      cmd << "'-' title \"{/Italic "
          << ElementIDtoString(pp.id.ionElement) << " -> "
          << ElementIDtoString(pp.id.clusterElement)
          << "_{" << pp.id.clusterSize << "}"
          << "}\" "
          << "with linespoints";
      plotCmds.push_back(cmd.str());

      data << "e\n";
    }
  }

  if (!(plotCmds.size() > 0))
  {
    TRACE(yieldOfWhat);
    TRACE(idStr);
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
BatchPostProcess::plotAngular(bool plotPolar,
                              StatPostProcess::FProcessClassicMolecule fpm,
                              std::string idStr,
                              ElementID specIonElement,
                              size_t specClusterSize,
                              ElementID specClusterElement) const
{
  std::stringstream fnb;
  fnb << idStr;

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
    const mdepp::StatPostProcess& pp = pps[i];

    if (specIonElement != DUMMY_EL)
    {
      if (pp.id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (pp.id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (pp.id.clusterSize != specClusterSize)
        continue;
    }

    {
      ostringstream sfname;
      sfname << pp.id.str << "/" << "Process" << idStr << "/"
             << "_angular" << "/"
             << (plotPolar?"atomsCount_by_polar_":"atomsCount_by_azimuth_")
             << n_str
             << ".dat";
      TRACE(sfname.str());
      std::ifstream fiSingle(sfname.str().c_str());
      REQUIRE(fiSingle);
      for(size_t i = 0; i < n; ++i)
      {
        Float x,y;
        fiSingle >> x >> y;
        data << x << " " << y << "\n";
      }
      fiSingle.close();
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

  std::vector<Float> dEs;
  dEs.resize(bounds.size()-1);

  fplt << "# 1st bound = " << bounds[0]/Ao << " Ao" << "\n"
       << "# last depth = " << *(bounds.end()-1)/Ao << " Ao" << "\n"
       << "# number of bounds = " << bounds.size() << "\n"
       << "# number of bins = " << dEs.size() << "\n";

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
    const mdepp::StatPostProcess& pp = pps[i];

    if (specIonElement != DUMMY_EL)
    {
      if (pp.id.ionElement != specIonElement)
        continue;
    }

    if (specClusterElement != DUMMY_EL)
    {
      if (pp.id.clusterElement != specClusterElement)
        continue;
    }

    if (specClusterSize > 0)
    {
      if (pp.id.clusterSize != specClusterSize)
        continue;
    }

    if (specIonEnergy > 0)
    {
      if (fabs(pp.id.ionEnergy*eV - specIonEnergy) > 0.05*eV)
        continue;
    }

    for(size_t ti = 0; ti < pp.trajData.size(); ti++)
    {
      const StatPostProcess::TrajData& td = pp.trajData[ti];

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
