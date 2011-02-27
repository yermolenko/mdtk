#ifndef mdpp_Cluster_hpp
#define mdpp_Cluster_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>
#include <vector>
#include "AtomGroup.hpp"

namespace mdepp
{
  using namespace mdtk;

class Cluster : public AtomGroup
{
public:
  Cluster();
  ~Cluster();
  Cluster(const Cluster &c);

  Cluster& operator =(const Cluster &c);

  virtual void get(std::istream& is);
  virtual void put(std::ostream& os) const;

  void build(const mdtk::SimLoop& ml);
};

}

#endif
