/*
   Implementation of the many-body interatomic potential for
   hydrocarbons.
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009, 2012 Oleksandr
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

#include "Brenner.hpp"
#include <algorithm>

namespace mdtk
{

Float
Brenner::operator()(AtomsArray& gl)
{
  Float Ei = 0;
  for(size_t ii = 0; ii < gl.size(); ii++)
  {
    Atom &atom_i = gl[ii];
    if (isHandled(atom_i))
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        if (isHandled(atom_j))
          if (&atom_i != &atom_j)
          {
            if (!probablyAreNeighbours(atom_i,atom_j)) continue;
            AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));

            Float VAvar = VA(ij);
            Ei += VR(ij,1.0);
            if (VAvar != 0.0)
            {
              Float BaverVal = Baver(ij,-VAvar);
              Ei += -BaverVal*VA(ij,-BaverVal);
            }
          }
      }
  }

  return Ei;
}

inline
Float
Brenner::VR(AtomsPair& ij, const Float V)
{
  return VR_Exp(ij,V);
}

inline
Float
Brenner::VR_Exp(AtomsPair& ij, const Float V)
{
  Float fvar = ij.f();

#ifdef BRENNER_OPTIMIZED
  if (fvar == 0.0) return 0.0;
#endif

  Float r = ij.r();

  Float De=   this->De(ij);
  Float S=    this->S(ij);
  Float beta= this->beta(ij);
  Float Re=   this->Re(ij);

  Float a = De/(S-1.0);
  Float b = sqrt(2.0*S)*beta;
  Float c = Re;

  Float Val = a*exp(-b*(r-c));

  if (V != 0)
  {
    Float Der = Val*(-b);
    ij.r(Der*fvar*V);
    ij.f(Val*V);
  }

  return fvar*Val;
}

inline
Float
Brenner::VA(AtomsPair& ij, const Float V)
{
  return VA_Exp(ij,V);
}

inline
Float
Brenner::VA_Exp(AtomsPair& ij, const Float V)
{
  Float fvar = ij.f();

#ifdef BRENNER_OPTIMIZED
  if (fvar == 0.0) return 0.0;
#endif

  Float r = ij.r();

  Float De=   this->De(ij);
  Float S=    this->S(ij);
  Float beta= this->beta(ij);
  Float Re=   this->Re(ij);

  Float a = De*S/(S-1.0);
  Float b = sqrt(2.0/S)*beta;
  Float c = Re;

  Float Val = a*exp(-b*(r-c));

  if (V != 0)
  {
    Float Der = Val*(-b);
    ij.r(Der*fvar*V);
    ij.f(Val*V);
  }

  return fvar*Val;
}

inline
Float
Brenner::Baver(AtomsPair& ij, const Float V)
{
  AtomsPair ji(-ij);

   return (B(ij,V/2.0)+B(ji,V/2.0))/2.0
         +
         (
           FN(ij,V/2.0)/2.0
         );
}

inline
Float
Brenner::ExpTerm(AtomsPair& ij, AtomsPair& ik, const Float V)
{
#ifdef BRENNER_OPTIMIZED
  if (alpha(ij,ik) == 0.0) return 1.0;
#endif
  Float r_ij = ij.r();
  Float r_ik = ik.r();

  Float Val = exp(alpha(ij,ik)*((r_ij-Re(ij))-(r_ik-Re(ik))));

  Float Dertmp = Val*alpha(ij,ik);

  if (V != 0.0)
  {
    ij.r(Dertmp*V);
    ik.r(-Dertmp*V);
/*
    ij.atom1.grad += ij.dr(ij.atom1)*Dertmp*V;
    ij.atom2.grad += ij.dr(ij.atom2)*Dertmp*V;
    ik.atom2.grad += ij.dr(ik.atom2)*Dertmp*V;
*/
  }

  return Val;
}

Float
Brenner::D(AtomsPair& ij, const Float V)
{
  Float Dij = 1.0;
  for(size_t k = 0; k < NL(ij.atom1).size(); k++)
  {
    Atom& atom_k = *(NL(ij.atom1)[k]);
    if (isHandled(atom_k))
    if (&atom_k != &ij.atom2/* && &atom_k != &atom1*/)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));

      Float ExpTermvar = ExpTerm(ij,ik);
      if (ExpTermvar != 0.0)
      {
        Float Gvar = G(ij,ik);
        Float fvar = ik.f();
        if (V != 0.0)
        {
          G(ij,ik,fvar*ExpTermvar*V);
          ik.f(Gvar*ExpTermvar*V);
          ExpTerm(ij,ik,Gvar*fvar*V);
        }
        Dij += Gvar*fvar*ExpTermvar;
      }
    }
  }
  Dij += HN(ij,V);
  return Dij;
}

inline
Float
Brenner::B(AtomsPair& ij, const Float V)
{
  Float DVal = D(ij);
  if (V != 0)
  {
    Float Der = -delta(ij.atom1)*pow(DVal,-delta(ij.atom1)-1.0);
    D(ij,Der*V);
  };
  return pow(DVal,-delta(ij.atom1));
}

inline
Float
Brenner::G(AtomsPair& ij, AtomsPair& ik, const Float V)
{
  if (ij.atom1.ID == C_EL)
  {
    Float CosT = CosTheta(ij,ik);
    if (V != 0)
    {
      CosTheta(ij,ik,dGdCT(CosT)*V);
    }
    return ao_*(1.0+SQR(co_)/SQR(do_)-SQR(co_)/(SQR(do_)+SQR(1.0+CosT)));
  }
  else if (ij.atom1.ID == H_EL)
  {
    return G_HH_;
  }
  else
  {
    throw Exception("Brenner::G() : unknown element");
  }
}

inline
Float
Brenner::dGdCT(Float CosT) const
{
  return 2.0*ao_*SQR(co_)*(1.0+CosT)/SQR(SQR(do_)+SQR(1.0+CosT));
}

Float
Brenner::Nt(AtomsPair& ij, const Float V)
{
  REQUIREM(ij.atom1.ID == C_EL,"Nt is only for C !!!");

  Float Nt_i = 0;
  for(size_t k = 0; k < NL(ij.atom1).size(); k++)
  {
    Atom &atom_k = *NL(ij.atom1)[k];
    if (isHandled(atom_k))
    if (&atom_k != &ij.atom2)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      Nt_i += ik.f(V);
    }
  }
  return Nt_i;
}

Float
Brenner::NH(AtomsPair& ij, const Float V)
{
  REQUIREM(ij.atom1.ID == C_EL,"NH is only for C !!!");
  Float NH_i = 0;
  for(size_t k = 0; k < NL(ij.atom1).size(); k++)
  {
    Atom &atom_k = *NL(ij.atom1)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == H_EL && &atom_k != &ij.atom2)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      NH_i += ik.f(V);
    }
  }
  return NH_i;
}

Float
Brenner::NC(AtomsPair& ij, const Float V)
{
  REQUIREM(ij.atom1.ID == C_EL,"NC is only for C !!!");
  Float NC_i = 0;
  for(size_t k = 0; k < NL(ij.atom1).size(); k++)
  {
    Atom &atom_k = *NL(ij.atom1)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == C_EL && &atom_k != &ij.atom2)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      NC_i += ik.f(V);
    }
  }
  return NC_i;
}

inline
Float
Brenner::F(Float x) const
{
  if (x<=2.0)
  {
    return 1.0;
  }
  else if (x>=3.0)
  {
    return 0.0;
  }
  else
  {
    return (1.0+cos(M_PI*(x-2.0)))/2.0;
  }
}

inline
Float
Brenner::dF(Float x) const
{
  if (x<=2.0)
  {
    return 0.0;
  }
  else if (x>=3.0)
  {
    return 0.0;
  }
  else
  {
    return (-M_PI*sin(M_PI*(x-2.0)))/2.0;
  }
}

Float
Brenner::Nconj(AtomsPair& ij, const Float V)
{
  REQUIREM(ij.atom1.ID == C_EL && ij.atom2.ID == C_EL,"Nconj is only for C-C !!!");

  Float Nconj_ij = 1.0;

  for(size_t k = 0; k < NL(ij.atom1).size(); k++)
  {
    Atom& atom_k = *(NL(ij.atom1)[k]);
//    if (isHandled(atom_k))
    if (&atom_k != &ij.atom2 && atom_k.ID == C_EL)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
//      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      AtomsPair ik(atom_k,ij.atom1,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));

      Float f_ik = ik.f();

      Float x_ik = Nt(ik);

      Float F_x_ik = F(x_ik);

      if (V != 0.0)
      {
        ik.f(F_x_ik*V);
        Nt(ik,f_ik*dF(x_ik)*V);
      }

      Nconj_ij += f_ik*F_x_ik;
    }
  }

  for(size_t l = 0; l < NL(ij.atom2).size(); l++)
  {
    Atom& atom_l = *(NL(ij.atom2)[l]);
//    if (isHandled(atom_l))
    if (&atom_l != &ij.atom1 && atom_l.ID == C_EL)
    {
      if (!probablyAreNeighbours(ij.atom2,atom_l)) continue;
//      AtomsPair jl(ij.atom2,atom_l,R(0,ij.atom2,atom_l),R(1,ij.atom2,atom_l));
      AtomsPair jl(atom_l,ij.atom2,R(0,ij.atom2,atom_l),R(1,ij.atom2,atom_l));

      Float f_jl = jl.f();

      Float x_jl = Nt(jl);

      Float F_x_jl = F(x_jl);

      if (V != 0.0)
      {
        jl.f(F_x_jl*V);
        Nt(jl,f_jl*dF(x_jl)*V);
      }

      Nconj_ij += f_jl*F_x_jl;
    }
  }

  return Nconj_ij;
}

inline
Float
Brenner::HN(AtomsPair& ij, const Float V)
{
RETURN_BRENNER_0;

  if (ij.atom1.ID == C_EL && ij.atom2.ID == C_EL)
  {
    Float NHVar = NH(ij);
    Float NCVar = NC(ij);
    Float Val = HN_CC_num(NHVar,NCVar);
    if (V != 0.0)
    {
      NH(ij,HN_CC_dH_num(NHVar,NCVar)*V);
      NC(ij,HN_CC_dC_num(NHVar,NCVar)*V);
    }
    return Val;
  }
  else if ( (ij.atom1.ID == C_EL && ij.atom2.ID == H_EL)/* || (atom_i.ID == H_EL && atom_j.ID == C_EL)*/)
  {
    Float NHVar = NH(ij);
    Float NCVar = NC(ij);
    Float Val = HN_CH_num(NHVar,NCVar);
    if (V != 0.0)
    {
      NH(ij,HN_CH_dH_num(NHVar,NCVar)*V);
      NC(ij,HN_CH_dC_num(NHVar,NCVar)*V);
    }
    return Val;
  }
  else return 0.0;
}

inline
Float
Brenner::FN(AtomsPair& ij, const Float V)
{
  RETURN_BRENNER_0;

  AtomsPair ji(-ij);

  if (ij.atom1.ID != C_EL || ij.atom2.ID != C_EL)
    return 0.0;
  else
  {
    Float Nti = Nt(ij);
    Float Ntj = Nt(ji);
    Float Nconj_ij = Nconj(ij);
    Float Val = FN_num(Nti,Ntj,Nconj_ij);
    if (V != 0.0)
    {
      Nt(ij,FN_dNt_i_num(Nti,Ntj,Nconj_ij)*V);
      Nt(ji,FN_dNt_j_num(Nti,Ntj,Nconj_ij)*V);
      Nconj(ij,FN_dNconj_num(Nti,Ntj,Nconj_ij)*V);
    }
    return Val;
  }
}

Brenner::Brenner(ParamSet parSet):
  FManybody(),
  paramSet(parSet),
  funcH_CC((paramSet==POTENTIAL1)?0:1),
  funcH_CH((paramSet==POTENTIAL1)?0:1),
  funcF((paramSet==POTENTIAL1)?0:1)
{
  handledElements.insert(H_EL);
  handledElements.insert(C_EL);
  handledElementPairs.insert(std::make_pair(H_EL,C_EL));
  handledElementPairs.insert(std::make_pair(C_EL,H_EL));
  switch (paramSet)
  {
    case POTENTIAL1:  setupPotential1(); break;
    case POTENTIAL2:  setupPotential2(); break;
    default : throw Exception("Unknown ParamSet for Brenner");
  };

  nl.Rcutoff = getRcutoff();
}

void
Brenner::setupPotential1()
{
  PTRACE("Setup Brenner1");
  
  Re_[C][C] = 1.315*Ao;
  Re_[H][H] = 0.74144*Ao;
  Re_[C][H] = 1.1199*Ao;
    Re_[H][C] = Re_[C][H]; 

  De_[C][C] = 6.325*eV;
  De_[H][H] = 4.7509*eV;
  De_[C][H] = 3.6422*eV;
    De_[H][C] = De_[C][H];

  beta_[C][C] = 1.5/Ao;
  beta_[H][H] = 1.9436/Ao;
  beta_[C][H] = 1.9583/Ao;
    beta_[H][C] = beta_[C][H];

  S_[C][C] = 1.29;
  S_[H][H] = 2.3432;
  S_[C][H] = 1.7386;
    S_[H][C] = S_[C][H];

  delta_[C][C] = 0.80469;
  delta_[H][H] = 0.80469;
  delta_[C][H] = 0.0; // undefined and unneeded
    delta_[H][C] = delta_[C][H];

  R_[C][C][0] = 1.7*Ao;
  R_[H][H][0] = 1.1*Ao;
  R_[C][H][0] = 1.3*Ao;
    R_[H][C][0] = R_[C][H][0];

// max R[*][*][1] == Rcutoff // see f()

  R_[C][C][1] = 2.0*Ao;
  R_[H][H][1] = 1.7*Ao;
  R_[C][H][1] = 1.8*Ao;
    R_[H][C][1] = R_[C][H][1];

  alpha_[H][H][H] = 3.0/Ao;
  alpha_[H][H][C] = 3.0/Ao;
  alpha_[H][C][H] = 3.0/Ao;
  alpha_[H][C][C] = 3.0/Ao;
  alpha_[C][H][H] = 3.0/Ao;
  alpha_[C][H][C] = 0.0/Ao;
  alpha_[C][C][H] = 0.0/Ao;
  alpha_[C][C][C] = 0.0/Ao;
    
  ao_ = 0.011304;
  co_ = 19.0;
  do_ = 2.5;

  G_HH_ = 4.0;
}  

void
Brenner::setupPotential2()
{
  PTRACE("Setup Brenner2");

  Re_[C][C] = 1.39*Ao;
  Re_[H][H] = 0.74144*Ao;
  Re_[C][H] = 1.1199*Ao;
    Re_[H][C] = Re_[C][H]; 

  De_[C][C] = 6.0*eV;
  De_[H][H] = 4.7509*eV;
  De_[C][H] = 3.6422*eV;
    De_[H][C] = De_[C][H];

  beta_[C][C] = 2.1/Ao;
  beta_[H][H] = 1.9436/Ao;
  beta_[C][H] = 1.9583/Ao;
    beta_[H][C] = beta_[C][H];

  S_[C][C] = 1.22;
  S_[H][H] = 2.3432;
  S_[C][H] = 1.69077;
    S_[H][C] = S_[C][H];

  delta_[C][C] = 0.5;
  delta_[H][H] = 0.5;
  delta_[C][H] = 0.0; // undefined and unneeded
    delta_[H][C] = delta_[C][H];

  R_[C][C][0] = 1.7*Ao;
  R_[H][H][0] = 1.1*Ao;
  R_[C][H][0] = 1.3*Ao;
    R_[H][C][0] = R_[C][H][0];

// max R[*][*][1] == Rcutoff // see f()

  R_[C][C][1] = 2.0*Ao;
  R_[H][H][1] = 1.7*Ao;
  R_[C][H][1] = 1.8*Ao;
    R_[H][C][1] = R_[C][H][1];

  alpha_[H][H][H] = 4.0/Ao;
  alpha_[H][H][C] = 4.0/Ao;
  alpha_[H][C][H] = 4.0/Ao;
  alpha_[H][C][C] = 4.0/Ao;
  alpha_[C][H][H] = 4.0/Ao;
  alpha_[C][H][C] = 0.0/Ao;
  alpha_[C][C][H] = 0.0/Ao;
  alpha_[C][C][C] = 0.0/Ao;
    
  ao_ = 0.00020813;
  co_ = 330.0;
  do_ = 3.5;

  G_HH_ = 12.33;
}  

inline
Float
Brenner::HN_CC_num(Float a1, Float a2) const
{
  return funcH_CC(a1,a2);
} 

inline
Float
Brenner::HN_CC_dH_num(Float a1, Float a2) const
{
  return funcH_CC.dH(a1,a2);
} 

inline
Float
Brenner::HN_CC_dC_num(Float a1, Float a2) const
{
  return funcH_CC.dC(a1,a2);
} 

inline
Float
Brenner::HN_CH_num(Float a1, Float a2) const
{
  return funcH_CH(a1,a2);
} 

inline
Float
Brenner::HN_CH_dH_num(Float a1, Float a2) const
{
  return funcH_CH.dH(a1,a2);
} 

inline
Float
Brenner::HN_CH_dC_num(Float a1, Float a2) const
{
  return funcH_CH.dC(a1,a2);
} 

inline
Float
Brenner::FN_num(Float a1, Float a2, Float a3) const 
{
  return funcF(a1,a2,a3);
} 

inline
Float
Brenner::FN_dNconj_num(Float a1, Float a2, Float a3) const 
{
  return funcF.dk(a1,a2,a3);
} 

inline
Float
Brenner::FN_dNt_i_num(Float a1, Float a2, Float a3) const 
{
  return funcF.di(a1,a2,a3);
} 

inline
Float
Brenner::FN_dNt_j_num(Float a1, Float a2, Float a3) const 
{
  return funcF.dj(a1,a2,a3);
} 

}


