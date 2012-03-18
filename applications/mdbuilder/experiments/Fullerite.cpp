/*
   Building of fullerite stuctures

   Copyright (C) 2007, 2008, 2009, 2011, 2012 Oleksandr Yermolenko
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

#include "Fullerite.hpp"
#include "FCC.hpp"
#include "Clusters.hpp"

namespace mdbuilder
{

using namespace mdtk;

SimLoop
build_Fullerite_C60(
  int a_num,
  int b_num,
  int c_num,
  bool fixBottomCellLayer,
  double a
  )
{
  double b = a;
  double c = a;

  SimLoopDump sl;
  initialize_simloop(sl);

  mdtk::AtomsArray C60 = mdbuilder::C60();

  C60.shiftToOrigin();

  place_Generic_FCC_lattice(sl.atoms,C60,
                            a_num,b_num,c_num,
                            fixBottomCellLayer,0,
                            a,a,a);

  sl.setPBC(Vector3D(a*a_num, b*b_num, c*c_num));

  std::vector<size_t> fixedAtoms =
    sl.atoms.unfixFixedAtoms(0,sl.atoms.size());

  {
    sl.enableDump();

    sl.dumpConst(0.95);
    relax_flush(sl,0.05*ps,"_tmp-W-relax095");

    sl.dumpConst(0.97);
    relax_flush(sl,0.05*ps,"_tmp-W-relax097");

    sl.dumpConst(0.99);
    relax_flush(sl,0.05*ps,"_tmp-W-relax099");

    sl.disableDump();

    relax_flush(sl,0.50*ps,"_tmp-W-relax100");
  }

  quench(sl,0.01*K,100000*ps,0.01*ps,"_tmp-X-quenchall");

  sl.atoms.fixAtoms(fixedAtoms);

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));

  {
    sl.enableDump();

    sl.dumpConst(0.95);
    relax_flush(sl,0.05*ps,"_tmp-X-relax095");

    sl.dumpConst(0.97);
    relax_flush(sl,0.05*ps,"_tmp-X-relax097");

    sl.dumpConst(0.99);
    relax_flush(sl,0.05*ps,"_tmp-X-relax099");

    sl.disableDump();

    relax_flush(sl,0.50*ps,"_tmp-X-relax100");
  }

  sl.setPBC(Vector3D(a*a_num, b*b_num, NO_PBC.z));
  sl.thermalBath.zMin = (c_num > 2)?(c*(c_num-2)-0.5*Ao):(0.0);
  sl.thermalBath.dBoundary = 3.0*Ao;

  relax(sl,0.01*ps);
  quench(sl,1.0*K);

  sl.atoms.removeMomentum();

  return sl;
}

}
