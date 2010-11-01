#include "Molecule.hpp"

namespace mdepp
{

bool
Molecule::hasAtom(mdtk::Atom& a) const
{
  bool found = false;
  for(size_t i = 0; i < atoms.size(); i++)
    if (atoms[i].globalIndex == a.globalIndex)
      found = true;
  return found;
}  

void
Molecule::buildFromAtom(mdtk::Atom& a, NeighbourList& nl,double SPOTTED_DISTANCE)
{
//TRACE("adding atom");
  if (hasAtom(a) || !isHandled(a)) return;
  atoms.push_back(a);atoms[atoms.size()-1].container = &dummy_ac;

//  TRACE(a.fixed);
//  TRACE(atoms[atoms.size()-1].fixed);
  

//TRACE(a.globalIndex);  
//TRACE(a.ID);  
//TRACE("1");  
//TRACE(a.nl(&(ml.fpot)).size());
  for(size_t i = 0; i < nl.nl[a.globalIndex].size(); i++)
  {
//TRACE("2");
    mdtk::Atom& nb_a = *nl.nl[a.globalIndex][i];
//    nb_a.container = &dummy_ac;
//       a.container = &dummy_ac;
    if (!isHandled(nb_a)) continue;
    Float distance = depos(a,nb_a).module();//sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
    if (Rc(a,nb_a) >= distance)
    {
      if (nb_a.coords.z < SPOTTED_DISTANCE)
      {
//TRACE("Next dive");        
        buildFromAtom(nb_a,nl,SPOTTED_DISTANCE);
      }
      else
      {
//TRACE("atoms.clear");        
        atoms.clear();
      }  
      if (atoms.size() == 0) break;
    }  
  }  
}  

}
