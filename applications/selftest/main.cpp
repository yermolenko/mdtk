/*
   selftest - program for the Molecular Dynamics Toolkit (MDTK) selftesting.

   Copyright (C) 2009, 2011, 2012 Oleksandr Yermolenko
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

#include <mdtk/SimLoop.hpp>
#include <mdtk/Atom.hpp>
#include <mdtk/potentials/pairwise/FLJ.hpp>
#include <mdtk/potentials/manybody/TightBinding/TightBinding.hpp>
#include <mdtk/potentials/manybody/Ackland/Ackland.hpp>
#include <fstream>

using namespace mdtk;
using namespace std;

void
pairtest()
{
  mdtk::verboseTrace = false;

  SimLoop mdloop;

  mdtk::FGeneral* pot = NULL;

  pot = new mdtk::TightBinding();
  mdloop.fpot.addPotential(pot);

  mdloop.atoms.push_back(Atom(Cu_EL,Vector3D(0.0*Ao,0.0*Ao,0.0*Ao)));
  mdloop.atoms.push_back(Atom(Cu_EL,Vector3D(0.5*Ao,0.0*Ao,0.0*Ao)));

  mdloop.executeDryRun();

  ofstream foe("e.dat");
  ofstream fof("f.dat");
  while (mdloop.atoms.back().coords.x < 6.0*Ao)
  {
    foe << mdloop.atoms.back().coords.x/Ao 
        << " " << (*pot)(mdloop.atoms)/eV
        << "\n";
    fof << mdloop.atoms.back().coords.x/Ao 
        << " " << (*pot).grad(mdloop.atoms.back(),mdloop.atoms).module() 
        << "\n";
    mdloop.atoms.back().coords.x += 0.001*Ao;
  }
  foe.close();
  fof.close();
}

int
main()
{
  pairtest();

  return 0;
}

