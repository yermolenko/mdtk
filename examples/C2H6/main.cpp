/*
   mdC2H6 - example program for the Molecular Dynamics Toolkit (MDTK).
   Molecular dynamics simulation of the ethane molecule.

   Copyright (C) 2009, 2012 Oleksandr Yermolenko
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
#include <mdtk/potentials/manybody/Brenner/Brenner.hpp>

using namespace mdtk;
using namespace std;

int
main()
{
  SimLoop mdloop;

  mdloop.fpot.addPotential(new Brenner());

  mdloop.atoms.push_back(Atom(C_EL,Vector3D(0.0*Ao,0.0*Ao,0.0*Ao)));
  mdloop.atoms.push_back(Atom(C_EL,Vector3D(1.5*Ao,0.0*Ao,0.0*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(1.9*Ao,0.0*Ao,1.0*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(1.9*Ao,-0.9*Ao,-0.5*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(1.9*Ao,0.9*Ao,-0.5*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(-0.4*Ao,0.0*Ao,1.0*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(-0.4*Ao,-0.9*Ao,-0.2*Ao)));
  mdloop.atoms.push_back(Atom(H_EL,Vector3D(-0.4*Ao,0.8*Ao,-0.2*Ao)));

  mdloop.simTimeSaveTrajInterval = 0.001*ps;
  mdloop.simTimeFinal = 0.5*ps;

  yaatk::text_ofstream fomde("input.mde");
  mdloop.saveToMDE(fomde);
  fomde.close();

  mdloop.execute();

  puts("Simulation finished.");
  puts("To visualize MD trajectory run \"mdtrajview\" from current directory.");
  puts("Run \"mdtrajview -a\" to run animation instantly.");

  return 0;
}

