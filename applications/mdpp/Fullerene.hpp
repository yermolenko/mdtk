#ifndef mdpp_Fullerene_hpp
#define mdpp_Fullerene_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>
#include <vector>
#include "AtomGroup.hpp"
#include "Cluster.hpp"
#include "Molecule.hpp"

namespace mdepp
{
  using namespace mdtk;

class Fullerene : public AtomGroup
{
public:
  Cluster cluster;
  Fullerene();
  ~Fullerene();
  Fullerene(const Fullerene &c);

  Fullerene& operator =(const Fullerene &c);

  virtual void get(std::istream& is);
  virtual void put(std::ostream& os) const;

  void build(const mdtk::SimLoop& ml);

  bool isEndoFullerene() const;
  bool isUnparted() const {return maxMolecule().size() == atoms.size();}
  Float maxDistanceFromMassCenter() const;
  Float minDistanceFromMassCenter() const;
  bool isIntegral(Float threshold = 3.0*Ao) const;
};

}

#endif
