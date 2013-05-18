/*
   Building of various clusters (header file)

   Copyright (C) 2007, 2008, 2010, 2011, 2012, 2013 Oleksandr
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

#ifndef MDBUILDER_Clusters_HPP
#define MDBUILDER_Clusters_HPP

#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>

#include "../common.hpp"

namespace mdbuilder
{

using namespace mdtk;

void
place_C60(AtomsArray& sl);

AtomsArray
C60();

AtomsArray
cluster(ElementID id, int clusterSize);

AtomsArray
clusterFromCrystal(const AtomsArray& atoms, int clusterSize, Vector3D c = NO_PBC);

AtomsArray
clusterFromFCCCrystal(ElementID id, int clusterSize);

AtomsArray
embed(AtomsArray cluster, AtomsArray shell);

void
add_rotational_motion(
  AtomsArray& atoms,
  Float totalRotEnergy = 100.0*eV,
  Vector3D rotAxis = Vector3D(0.0,1.0,0.0)
  );

SimLoop
build_target_by_cluster_bombardment(
  const SimLoop& sl_target,
  AtomsArray sl_cluster,
  Float clusterEnergy = 100*eV,
  Float interactionDistance = 3.3*Ao
  );

void
prepare_Cu_by_Cu_at_C60_bobardment();

void
prepare_Graphite_by_Cu_at_C60_bombardment();

SimLoop
build_Cluster_Landed_on_Substrate(
  const mdtk::SimLoop sl_Substrate,
  AtomsArray cluster,
  bool applyPBCtoCluster = false
  );

void
bomb_Cluster_with_Ions(
  std::string dirname,
  const SimLoop& target,
  std::vector<size_t> clusterAtomIndices,
  ElementID ionElement,
  Float ionEnergy,
  std::set<Float> halos,
  size_t numberOfImpacts = 1024
  );

void
bomb_landedCluster_with_Ions(
  const SimLoop& target,
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  std::set<Float> halos,
  size_t numberOfImpacts = 1024
);

void
bomb_orthorhombic_with_clusters(
  std::string dirname,
  mdtk::SimLoop cluster,
  const mdtk::SimLoop target,
  int a_num,
  int b_num,
  double a,
  double b,
  size_t numberOfImpacts = 1024
  );

void
build_FCC_metal_bombardment_with_ions(
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  size_t numberOfImpacts = 1024,
  int a_num = 7,
  int b_num = 7,
  int c_num = 7,
  double a = 3.615*Ao,
  double b = 3.615*Ao,
  double c = 3.615*Ao,
  ElementID metalElement = Cu_EL
  );

void
build_FCC_metal_bombardment_with_C60(
  std::vector<Float> fullereneEnergies,
  size_t numberOfImpacts = 1024,
  int a_num = 7,
  int b_num = 7,
  int c_num = 7,
  double a = 3.615*Ao,
  double b = 3.615*Ao,
  double c = 3.615*Ao,
  ElementID metalElement = Cu_EL
  );

void
build_fullerite_bombardment_with_ions(
  std::vector<ElementID> ionElements,
  std::vector<Float> ionEnergies,
  size_t numberOfImpacts = 1024,
  int a_num = 8,
  int b_num = 8,
  int c_num = 10,
  double a = 14.17*Ao
  );

void
build_metal_C60_mixing(
  std::vector<Float> impactEnergies,
  ElementID metalElement = Cu_EL
  );

}

#endif
