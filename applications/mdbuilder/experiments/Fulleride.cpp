/*
   Building of Fullerides

   Copyright (C) 2011, 2012 Oleksandr Yermolenko
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

#include "Fulleride.hpp"
#include "FCC.hpp"
#include "Clusters.hpp"

namespace mdbuilder
{

using namespace mdtk;

SimLoop
build_Fulleride_C60(
  int a_num,
  int b_num,
  int c_num,
  ElementID id,
  size_t numberOfAtoms,
  bool fixBottomCellLayer,
  double a
  )
{
  double b = a;
  double c = a;

  mdtk::SimLoop sl;
  initialize_simloop(sl);

  mdtk::AtomsArray C60 = mdbuilder::C60();

  place_Generic_FCC_lattice(sl.atoms,C60,
                            a_num,b_num,c_num,
                            fixBottomCellLayer,0,
                            a,a,a);

  mdtk::AtomsArray intercal = mdbuilder::cluster(id, numberOfAtoms);

  // endo
  place_Generic_FCC_lattice(sl.atoms,intercal,
                            a_num,b_num,c_num,
                            fixBottomCellLayer,0,
                            a,a,a);

  // interstitial
  place_Generic_FCC_lattice(sl.atoms,intercal,
                            a_num,b_num,c_num,
                            fixBottomCellLayer,0,
                            a,a,a,
                            true);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBath.zMin = (c_num > 2)?(c*(c_num-2)-0.5*Ao):(0.0);
  sl.thermalBath.dBoundary = 3.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,1.0*K);

  sl.atoms.removeMomentum();

  return sl;
}

}
