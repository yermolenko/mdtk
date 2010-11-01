#ifndef mdpp_ClassicMolecule_hpp
#define mdpp_ClassicMolecule_hpp

#include <iostream>
#include <mdtk/SimLoop.hpp>
#include <algorithm>

#include "NeighbourList.hpp"

namespace mdepp
{

using namespace mdtk;

#define ATOMTAG_FIXED 1<<0
#define ATOMTAG_SUBSTRATE 1<<1
#define ATOMTAG_CLUSTER   1<<2
#define ATOMTAG_NOTAG 0

//#define REF_POT_OF(ML_FPOT) ML_FPOT.potentials[0]


inline
void
setTags(mdtk::SimLoop* ml)
{
  for(size_t i = 0; i < ml->atoms_.size(); i++)
  {
    mdtk::Atom& atom = *(ml->atoms_[i]);
    atom.tag = 0;
    if (atom.M > 1000.0*mdtk::amu) atom.tag |= ATOMTAG_FIXED;
    if (atom.ID == mdtk::C_EL) atom.tag |= ATOMTAG_CLUSTER;
//    if (atom.ID == mdtk::Cu_EL) atom.tag |= ATOMTAG_CLUSTER;
//    if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) atom.tag |= ATOMTAG_SUBSTRATE;
//    if (atom.ID == mdtk::Ar_EL) atom.tag |= ATOMTAG_PROJECTILE;
  }
}

inline
bool  isProjectileAtom(const mdtk::Atom& atom)// const
{
  if (atom.ID == mdtk::Ar_EL) return true;
//  if (atom.tag & ATOMTAG_CLUSTER) return true;
  return false;
}  

inline
bool  isClusterAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::Cu_EL) return true;
  if (atom.tag & ATOMTAG_CLUSTER) return true;
  return false;
}  

inline
bool  isSubstrateAtom(const mdtk::Atom& atom)// const
{
//  if (atom.ID == mdtk::C_EL || atom.ID == mdtk::H_EL) return true;
  if (atom.tag & ATOMTAG_SUBSTRATE) return true;
  return false;
}  

class ClassicMolecule
{
private:
  enum {ECOUNT = 4};
  enum {C = 0};
  enum {H = 1};
  enum {Ar = 2};
  enum {Cu = 3};
  Float Rc_[ECOUNT][ECOUNT];
  size_t e2i(const mdtk::Atom &atom) const
  {
    using namespace mdtk;
    switch (atom.ID)
    {
      case H_EL : return H; break;
      case C_EL : return C; break;
      case Ar_EL : return Ar; break;
      case Cu_EL : return Cu; break;
      default : throw Exception("ClassicMolecule::e2i() : unknown element");
    };  
  }  
public:
  Float Rc(const mdtk::Atom &atom1,const mdtk::Atom &atom2) const
  {
    return Rc_[e2i(atom1)][e2i(atom2)];
  }  
  std::set<mdtk::ElementID> handledElements;
  bool
  isHandled(mdtk::Atom& atom) const
  {
    if (handledElements.find(atom.ID) != handledElements.end())
      return true;
    else
      return false;
  }  
public:
  Float formationTime;
  Float escapeTime;
  std::vector<mdtk::Atom> atoms;
  std::vector<mdtk::Atom> atoms_init;
//  AtomsContainer atoms;
//  AtomsContainer atoms_init;
  int  trajectory_dummy;
/*
ClassicMolecule( const ClassicMolecule &C )
 :atoms(),atoms_init()
{
  initParams();    

  formationTime = C.formationTime;
  escapeTime = C.escapeTime;
  trajectory = C.trajectory;
  atoms = C.atoms;
  atoms_init = C.atoms_init;
  TRACE(atoms.size());
  for(size_t i = 0; i < atoms.size(); i++)
  {
    TRACE(i);
    atoms[i].container = &dummy_ac;
    TRACE("1");
    atoms_init[i].container = &dummy_ac;
    TRACE("2");
  }
}

inline
ClassicMolecule&
operator =(const ClassicMolecule &C) 
{
  if (this == &C) return *this;

  formationTime = C.formationTime;
  escapeTime = C.escapeTime;
  trajectory = C.trajectory;
  atoms = C.atoms;
  atoms_init = C.atoms_init;
//  TRACE(atoms.size());
  for(size_t i = 0; i < atoms.size(); i++)
  {
    atoms[i].container = &dummy_ac;
    atoms_init[i].container = &dummy_ac;
  }

  return *this;
}
*/
  void initParams()
  {
    using namespace mdtk;
    handledElements.insert(H_EL);
    handledElements.insert(C_EL);
    handledElements.insert(Ar_EL);
    handledElements.insert(Cu_EL);

    for(size_t e1 = 0; e1 < ECOUNT; e1++)
      for(size_t e2 = 0; e2 < ECOUNT; e2++)
        Rc_[e1][e2] = 0.001*Ao;

    Rc_[C][C] = 2.0*Ao;
    Rc_[H][H] = 1.7*Ao;
    Rc_[C][H] = 1.8*Ao;
      Rc_[H][C] = Rc_[C][H];

    Rc_[Cu][Cu] = 3.0*Ao;


    Rc_[Cu][H] = 2.75*Ao;
      Rc_[H][Cu] = Rc_[Cu][H];

    Rc_[Cu][C] = 2.75*Ao;
      Rc_[C][Cu] = Rc_[Cu][C];
      
//    Rc_[Ar][Ar] = 0.001*Ao;
//      Rc_[Ar][ H] = Rc_[ H][Ar] = Rc_[ C][Ar] = Rc_[Ar][ C] = Rc_[Ar][Ar];
  }

  void addAtom(mdtk::Atom& a) {atoms.push_back(a);}
  void buildFromAtom(mdtk::Atom&, NeighbourList& nl,double SPOTTED_DISTANCE);
  bool hasAtom(mdtk::Atom&) const;
  ClassicMolecule():handledElements(),formationTime(-1),escapeTime(-1),
    atoms(/*0*/),atoms_init(),trajectory_dummy(-1)
  {
    initParams();    
  }
  ~ClassicMolecule(){;}
/*  
  bool  hasBackScatteredAtoms() const
  {
    return hasProjectileAtoms();
  }  
*/  

  bool  hasProjectileAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isProjectileAtom(atoms[ai])) return true;
    return false;
  }  
  bool  hasClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isClusterAtom(atoms[ai])) return true;
    return false;
  }  
  bool  hasSubstrateAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (isSubstrateAtom(atoms[ai])) return true;
    return false;
  }

  bool  hasOnlyProjectileAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isProjectileAtom(atoms[ai])) return false;
    return true;
  }  
  bool  hasOnlyClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isClusterAtom(atoms[ai])) return false;
    return true;
  }  
  bool  hasOnlySubstrateAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!isSubstrateAtom(atoms[ai])) return false;
    return true;
  }
  bool  hasOnlySubstrateOrClusterAtoms() const
  {
    for(size_t ai = 0; ai < atoms.size(); ai++)
      if (!(isSubstrateAtom(atoms[ai]) || isClusterAtom(atoms[ai]))) return false;
    return true;
  }
  Float getAMUMass() const
  {
    Float moleculeMass = 0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      moleculeMass += atom.M;
    } 
    return mdtk::academic_round(moleculeMass/mdtk::amu);
  }
  mdtk::Vector3D getVelocity() const
  {
    using mdtk::Exception;    
    
    REQUIRE(atoms.size() > 0);
    mdtk::Vector3D sumOfP = 0.0;
    Float sumOfM = 0.0;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      sumOfM += atom.M;
      sumOfP += atom.V*atom.M;
    };
    return sumOfP/sumOfM;    
  }  
  void printGlobalIndexes(std::ostream& fo) const
  {
    fo << atoms.size() << std::endl;
    for(size_t ai = 0; ai < atoms.size(); ai++)
    {
      const mdtk::Atom& atom = atoms[ai];
      fo << atom.globalIndex << std::endl;
    } 
  } 
  void saveToStream(std::ostream& os) const
  {
    os << formationTime << "\n";
    os << escapeTime << "\n";
    os << atoms.size() << "\n";
    for(size_t i = 0; i < atoms.size(); i++)
      os << atoms[i] << "\n";
    os << atoms_init.size() << "\n";
    for(size_t i = 0; i < atoms_init.size(); i++)
      os << atoms_init[i] << "\n";
    os << trajectory_dummy << "\n";
  }  
  void loadFromStream(std::istream& is)
  {
//    ERRTRACE("Loading ClassicMolecule....");
    is >> formationTime;
    is >> escapeTime;
    size_t sz;
    is >> sz;
    atoms.resize(sz);
    for(size_t i = 0; i < atoms.size(); i++)
    {
//      ERRTRACE(i);
      is >> atoms[i];
    }
    is >> sz;
    atoms_init.resize(sz);
    for(size_t i = 0; i < atoms_init.size(); i++)
      is >> atoms_init[i];
    is >> trajectory_dummy;
  }  

  friend int operator<(const ClassicMolecule& left, const ClassicMolecule& right);
  friend int operator<(ClassicMolecule& left, ClassicMolecule& right);
};  


inline
int operator<(const ClassicMolecule& left, const ClassicMolecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}

inline
int operator<(ClassicMolecule& left, ClassicMolecule& right)
{
  return left.getAMUMass() < right.getAMUMass();
}


}

#endif
