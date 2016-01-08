/* 
   Molecular dynamics postprocessor, BatchPostProcess classes, header

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

#ifndef mdtk_BatchPostProcess_hpp
#define mdtk_BatchPostProcess_hpp

#include "StatPostProcess.hpp"

namespace mdepp
{

using namespace mdtk;

class BatchPostProcess
{
  static size_t getHaloIndex(std::string);
public:
  std::vector< std::vector<std::string> > pps;

  BatchPostProcess(const std::string mdeppinPath, const size_t haloIndex = 0);
  BatchPostProcess();

  void addHalo(const std::string mdeppinPath, const size_t haloIndex);

  void printResults() const;

  void plot_Ekin_t(StatPostProcess::FProcessClassicMolecule fpm, Float ionEnergyFilter) const;

  void plotYieldsAgainstIonEnergy(StatPostProcess::FProcessClassicMolecule fpm,
                                  ElementID specIonElement = DUMMY_EL,
                                  size_t specClusterSize = 0,
                                  ElementID specClusterElement = DUMMY_EL) const;
  void plotMassSpectrum(StatPostProcess::FProcessClassicMolecule fpm,
                        ElementID specIonElement = DUMMY_EL,
                        size_t specClusterSize = 0,
                        ElementID specClusterElement = DUMMY_EL) const;
  void plotAngular(bool plotPolar,
                   StatPostProcess::FProcessClassicMolecule fpm,
                   ElementID specIonElement = DUMMY_EL,
                   size_t specClusterSize = 0,
                   ElementID specClusterElement = DUMMY_EL) const;
  void plotEnergyLoss(ElementID specIonElement = DUMMY_EL,
                      size_t specClusterSize = 0,
                      Float specIonEnergy = -100*eV,
                      ElementID specClusterElement = DUMMY_EL) const;
};

}


#endif
