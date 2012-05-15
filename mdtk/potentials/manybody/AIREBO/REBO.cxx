/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. The REBO potential part.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2011, 2012 Oleksandr
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

#include "REBO.hpp"
#include <algorithm>

namespace mdtk
{

Float
REBO::operator()(AtomsArray& gl)
{
  countNeighbours(gl);
  Float Ei = 0;
  if (gl.size() != pairs.size()) pairs.resize(gl.size());
  size_t ii;
  for(ii = 0; ii < gl.size(); ii++)
  {
    size_t prevSize = pairs[ii].size();
    pairs[ii].clear();
    pairs[ii].reserve(prevSize+FMANYBODY_PAIRS_RESERVE_ADD);
  }
  for(ii = 0; ii < gl.size(); ii++)
  {
    Atom &atom_i = gl[ii];
    if (isHandled(atom_i))
    {
      AtomRefsContainer& nli = NL(atom_i);
      for(size_t jj = 0; jj < nli.size(); jj++)
      {
        Atom &atom_j = *(nli[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        if (&atom_i != &atom_j)
        {
          if (!probablyAreNeighbours(atom_i,atom_j)) continue;
          AtomsPair ij(atom_i,atom_j,R(0,atom_i,atom_j),R(1,atom_i,atom_j));

          Float VAvar = VA(ij);
          Ei += VR(ij,1.0);
          if (VAvar != 0.0)
          {
            Float BaverVal = Baver(ij,VAvar);
            Ei += BaverVal*VA(ij,BaverVal);
          }
        }
      }
    }
  }

  return Ei;
}

inline
Float
REBO::Sprime(Float arg, Float arg_min, Float arg_max) const
{
  REQUIRE(arg_min < arg_max);

  if (arg<arg_min)
  {
    return 1.0;
  }
  else if (arg>arg_max)
  {
    return 0.0;
  }
  else
  {
    return (1.0+cos(M_PI*(arg-arg_min)
            /(arg_max-arg_min)))/2.0;
  }
}

inline
Float
REBO::dSprime(Float arg, Float arg_min, Float arg_max) const
{
  REQUIRE(arg_min < arg_max);

  if (arg<arg_min)
  {
    return 0.0;
  }
  else if (arg>arg_max)
  {
    return 0.0;
  }
  else
  {
    return (-(M_PI/(arg_max-arg_min))*sin(M_PI*(arg-arg_min)/(arg_max-arg_min)))/2.0;
  }
}

inline
Float
REBO::fprime(AtomsPair& ij, const Float V)
{
  Float R1=Rprime(0,ij);
  Float R2=Rprime(1,ij);

  Float r = ij.r();

  Float Val = Sprime(r,R1,R2);

  if (r >= R1 && r <= R2)
  {
    ij.r(dSprime(r,R1,R2)*V);
  }

  return Val;
}

inline
Float
REBO::VR(AtomsPair& ij, const Float V)
{
  return VR_Exp(ij,V);
}

inline
Float
REBO::VR_Exp(AtomsPair& ij, const Float V)
{
  Float f = ij.f();

#ifdef REBO_OPTIMIZED
  if (f == 0.0) return 0.0;
#endif

  Float r = ij.r();

  Float Q=     this->Q(ij);
  Float A=     this->A(ij);
  Float alpha= this->alpha(ij);

  Float t1 = A*exp(-alpha*r);

  Float Val = (1+Q/r)*t1;

  if (V != 0)
  {
    Float Der = -t1*(Q+alpha*r*r+alpha*r*Q)/(r*r);
    ij.r(Der*f*V);
    ij.f(Val*V);
  }

  return f*Val;
}

inline
Float
REBO::VA(AtomsPair& ij, const Float V)
{
  Float f = ij.f();

#ifdef REBO_OPTIMIZED
  if (f == 0.0) return 0.0;
#endif

  Float r = ij.r();

  Float B1=   this->B(1,ij);
  Float B2=   this->B(2,ij);
  Float B3=   this->B(3,ij);

  Float beta1=   this->beta(1,ij);
  Float beta2=   this->beta(2,ij);
  Float beta3=   this->beta(3,ij);

  Float t1 = B1*exp(-beta1*r);
  Float t2 = B2*exp(-beta2*r);
  Float t3 = B3*exp(-beta3*r);

  Float Val = t1+t2+t3;

  if (V != 0)
  {
    Float Der = (-beta1)*t1+(-beta2)*t2+(-beta3)*t3;
    ij.r(-Der*f*V);
    ij.f(-Val*V);
  }

  return -f*Val;
}

Float
REBO::Baver(AtomsPair& ij, const Float V)
{
  AtomsPair ji(-ij);

  return (p_sigma_pi(ij,V/2.0)+p_sigma_pi(ji,V/2.0))/2.0
    +pi_rc(ij,V)
#ifdef REBO_DIHEDRAL
    +pi_dh(ij,V)
#endif
    ;
}

inline
Float
REBO::ExpTerm(AtomsPair& ij, AtomsPair& ik, const Float V)
{
  if (ij.atom1.ID != H_EL) return 1.0;

  Float r_ij = ij.r_alter();
  Float r_ik = ik.r();

  Float Val = exp((4.0/Ao)*((rho(ik) - r_ik)-(rho(ij) - r_ij)));

  Float Dertmp = Val*(4.0/Ao);

  if (V != 0.0)
  {
    ij.r_alter(Dertmp*V);
    ik.r(-Dertmp*V);
  }

  return Val;
}

Float
REBO::D(AtomsPair& ij, const Float V)
{
  Float Dij = 1.0;
  AtomRefsContainer& nli = NL(ij.atom1);
  for(size_t k = 0; k < nli.size(); k++)
  {
    Atom& atom_k = *(nli[k]);
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
  Dij += P(ij,V);
  return Dij;
}

inline
Float
REBO::p_sigma_pi(AtomsPair& ij, const Float V)
{
  Float DVal = D(ij);
  if (V != 0)
  {
    Float Der = -0.5*pow(DVal,-0.5-1.0);
    D(ij,Der*V);
  };
  return pow(DVal,-0.5);
}

inline
Float
REBO::G(AtomsPair& ij, AtomsPair& ik, const Float V)
{
  Float CosT = CosTheta(ij,ik);

  if (ij.atom1.ID == C_EL)
  {
    Float Nt1=Ntot_[0];
    Float Nt2=Ntot_[1];

    Float N = Nt(ij);
    Float Sp = Sprime(N,Nt1,Nt2);

    Float funcGC1 = funcG_C1(CosT);
    Float funcGC2 = funcG_C2(CosT);

    if (V != 0.0)
    {
      Float SpDer = dSprime(N,Nt1,Nt2);
      Float funcGC1Der = funcG_C1.dCosT(CosT);
      Float funcGC2Der = funcG_C2.dCosT(CosT);

      if (funcGC2Der != 0.0 || funcGC1Der != 0.0)
        CosTheta(ij,ik,funcGC2Der*V + funcGC1Der*Sp*V + (-funcGC2Der*Sp*V));

      if ((funcGC2 != 0.0 || funcGC1 != 0.0) && SpDer != 0.0)
        Nt_donly(ij,funcGC1*SpDer*V + (-funcGC2*SpDer*V));
    }

    return funcGC2+Sp*(funcGC1 - funcGC2);
  }
  else if (ij.atom1.ID == H_EL)
  {
    if (V != 0.0)
    {
      Float funcGHDer = funcG_H.dCosT(CosT);
      if (funcGHDer != 0.0)
        CosTheta(ij,ik,funcGHDer*V);
    }
    return funcG_H(CosT);
  }
  else
  {
    throw Exception("REBO::G() : unknown element");
  }
}

void
REBO::countNeighbours(AtomsArray& gl)
{
  Nts.clear();
  Nts.resize(gl.size());
  Nts2touch.clear();
  Nts2touch.resize(gl.size());

  NCs.clear();
  NCs.resize(gl.size());
  NCs2touch.resize(gl.size());
  NCs2touch.clear();

  NHs.clear();
  NHs.resize(gl.size());
  NHs2touch.clear();
  NHs2touch.resize(gl.size());

  for(size_t i = 0; i < gl.size(); i++)
  {
    Atom &atom_i = gl[i];
    if (!isHandled(atom_i)) continue;
    AtomRefsContainer& nli = NL(atom_i);
    for(size_t k = 0; k < nli.size(); k++)
    {
      Atom &atom_k = *nli[k];

      Float fval = 0.0;
      bool  touchNeeded = false;

      Float r = depos(atom_i,atom_k).module();
      Float R2 = R(1,atom_i,atom_k);
      if (r>R2)
        continue;
      else
      {
        Float R1=R(0,atom_i,atom_k);
        if (r<R1)
          fval = 1.0;
        else
        {
          touchNeeded = true;
          fval = (1.0+cos(M_PI*(r-R1)/(R2-R1)))/2.0;
        }
      }

      {
        Nts[i] += fval;
        if (touchNeeded)
          Nts2touch[i].push_back(&atom_k);
      }
      if (atom_k.ID == H_EL)
      {
        NHs[i] += fval;
        if (touchNeeded)
          NHs2touch[i].push_back(&atom_k);
      }
      if (atom_k.ID == C_EL)
      {
        NCs[i] += fval;
        if (touchNeeded)
          NCs2touch[i].push_back(&atom_k);
      }
    }
  }
}

inline
bool
REBO::areNeighbours(Atom &atom_i, Atom &atom_j)
{
  AtomRefsContainer& nli = NL(atom_i);
  for(size_t k = 0; k < nli.size(); k++)
  {
    Atom &atom_k = *nli[k];
    if (&atom_k == &atom_j)
      return true;
  }
  return false;
}

inline
void
REBO::dfThem(const std::vector<std::vector<Atom*> >& patoms, AtomsPair& ij, const Float V)
{
  if (V == 0.0) return;
  const std::vector<Atom*>& they = patoms[ij.atom1.globalIndex];
  for(size_t k = 0; k < they.size(); k++)
    if (they[k] != &ij.atom2)
    {
      AtomsPair ik(ij.atom1,*they[k],R(0,ij.atom1,*they[k]),R(1,ij.atom1,*they[k]));
      ik.f(V);
    }
}

Float
REBO::Nt(AtomsPair& ij, const Float V)
{
//  REQUIRE(areNeighbours(atom_i,atom_j));
  if (V != 0.0)
  {
    dfThem(Nts2touch,ij,V);
  }
  Float N = Nts[ij.atom1.globalIndex];
  if (ij.r() < R(1,ij))
    N -= ij.f();
  return N;
}

void
REBO::Nt_donly(AtomsPair& ij, const Float V)
{
//  REQUIRE(areNeighbours(atom_i,atom_j));
  if (V != 0.0)
  {
    dfThem(Nts2touch,ij,V);
  }
}

Float
REBO::NH(AtomsPair& ij, const Float V)
{
  if (V != 0.0)
  {
    dfThem(NHs2touch,ij,V);
  }
  Float N = NHs[ij.atom1.globalIndex];
  if (ij.atom2.ID == H_EL && ij.r() < R(1,ij))
  {
    N -= ij.f();
  }
  return N;
}

void
REBO::NH_donly(AtomsPair& ij, const Float V)
{
  if (V != 0.0)
  {
    dfThem(NHs2touch,ij,V);
  }
}

Float
REBO::NC(AtomsPair& ij, const Float V)
{
  if (V != 0.0)
  {
    dfThem(NCs2touch,ij,V);
  }
  Float N = NCs[ij.atom1.globalIndex];
  if (ij.atom2.ID == C_EL && ij.r() < R(1,ij))
  {
    N -= ij.f();
  }
  return N;
}

void
REBO::NC_donly(AtomsPair& ij, const Float V)
{
  if (V != 0.0)
  {
    dfThem(NCs2touch,ij,V);
  }
}

inline
Float
REBO::F(Float x) const
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
REBO::dF(Float x) const
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
REBO::NconjSum1(AtomsPair& ij, const Float V)
{
  Float sum1 = 0.0;
  AtomRefsContainer& nli = NL(ij.atom1);
  for(Index k = 0; k < nli.size(); k++)
  {
    Atom& atom_k = *(nli[k]);
    if (&atom_k != &ij.atom2 && atom_k.ID == C_EL)
    {
      if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
//      AtomsPair ik(ij.atom1,atom_k,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      AtomsPair ik(atom_k,ij.atom1,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
      Float f_ik = ik.f();

      Float N   = Nt(ik);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];

      Float Sp = Sprime(N,Nt1,Nt2);

      if (V != 0.0)
      {
        ik.f(Sp*V);
        Nt_donly(ik,f_ik*dSprime(N,Nt1,Nt2)*V);
      }

      sum1 += f_ik*Sp;
    }
  }

  return sum1;
}

Float
REBO::NconjSum2(AtomsPair& ij, const Float V)
{
  Float sum2 = 0.0;
  AtomRefsContainer& nlj = NL(ij.atom2);
  for(Index l = 0; l < nlj.size(); l++)
  {
    Atom& atom_l = *(nlj[l]);
    if (&atom_l != &ij.atom1 && atom_l.ID == C_EL)
    {
      if (!probablyAreNeighbours(ij.atom2,atom_l)) continue;
//      AtomsPair jl(ij.atom2,atom_l,R(0,ij.atom2,atom_l),R(1,ij.atom2,atom_l));
      AtomsPair jl(atom_l,ij.atom2,R(0,ij.atom2,atom_l),R(1,ij.atom2,atom_l));
      Float f_jl = jl.f();

      Float N   = Nt(jl);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];

      Float Sp = Sprime(N,Nt1,Nt2);

      if (V != 0.0)
      {
        jl.f(Sp*V);
        Nt_donly(jl,f_jl*dSprime(N,Nt1,Nt2)*V);
      }

      sum2 += f_jl*Sp;
    }
  }

  return sum2;
}

Float
REBO::Nconj(AtomsPair& ij, const Float V)
{
  Float Nconj_ij = 1.0;

  Float sum1 = NconjSum1(ij);
  NconjSum1(ij,2.0*sum1*V);
  Nconj_ij += SQR(sum1);

  Float sum2 = NconjSum2(ij);
  NconjSum2(ij,2.0*sum2*V);
  Nconj_ij += SQR(sum2);

  return Nconj_ij;
}

inline
Float
REBO::P(AtomsPair& ij, const Float V)
{
RETURN_REBO_0;

  if (ij.atom1.ID == C_EL && ij.atom2.ID == C_EL)
  {
    Float NHVar = NH(ij);
    Float NCVar = NC(ij);
    Float Val = P_CC_num(NHVar,NCVar);
    if (V != 0.0)
    {
      NH_donly(ij,P_CC_dH_num(NHVar,NCVar)*V);
      NC_donly(ij,P_CC_dC_num(NHVar,NCVar)*V);
    }
    return Val;
  }
  else if (ij.atom1.ID == C_EL && ij.atom2.ID == H_EL)
  {
    Float NHVar = NH(ij);
    Float NCVar = NC(ij);
    Float Val = P_CH_num(NHVar,NCVar);
    if (V != 0.0)
    {
      NH_donly(ij,P_CH_dH_num(NHVar,NCVar)*V);
      NC_donly(ij,P_CH_dC_num(NHVar,NCVar)*V);
    }
    return Val;
  }
  else return 0.0;
}

inline
Float
REBO::pi_rc(AtomsPair& ij, const Float V)
{
  RETURN_REBO_0;

  AtomsPair ji(-ij);

  if (ij.atom1.ID == C_EL && ij.atom2.ID == C_EL)
  {
    Float Nti = Nt(ij);
    Float Ntj = Nt(ji);
    Float Nconj_ij = Nconj(ij);
    Float Val = pi_rc_CC_num(Nti,Ntj,Nconj_ij);
    if (V != 0.0)
    {
      Nt_donly(ij,pi_rc_CC_dNt_i_num(Nti,Ntj,Nconj_ij)*V);
      Nt_donly(ji,pi_rc_CC_dNt_j_num(Nti,Ntj,Nconj_ij)*V);
      Nconj(ij,pi_rc_CC_dNconj_num(Nti,Ntj,Nconj_ij)*V);
    }
    return Val;
  }
  if (ij.atom1.ID == C_EL && ij.atom2.ID == H_EL)
  {
    Float Nti = Nt(ij);
    Float Ntj = Nt(ji);
    Float Nconj_ij = Nconj(ij);
    Float Val = pi_rc_CH_num(Nti,Ntj,Nconj_ij);
    if (V != 0.0)
    {
      Nt_donly(ij,pi_rc_CH_dNt_i_num(Nti,Ntj,Nconj_ij)*V);
      Nt_donly(ji,pi_rc_CH_dNt_j_num(Nti,Ntj,Nconj_ij)*V);
      Nconj(ij,pi_rc_CH_dNconj_num(Nti,Ntj,Nconj_ij)*V);
    }
    return Val;
  }
  if (ij.atom1.ID == H_EL && ij.atom2.ID == C_EL)
  {
    return pi_rc(ji,V);
  }
  if (ij.atom1.ID == H_EL && ij.atom2.ID == H_EL)
  {
    Float Nti = Nt(ij);
    Float Ntj = Nt(ji);
    Float Nconj_ij = Nconj(ij);
    Float Val = pi_rc_HH_num(Nti,Ntj,Nconj_ij);
    if (V != 0.0)
    {
      Nt_donly(ij,pi_rc_HH_dNt_i_num(Nti,Ntj,Nconj_ij)*V);
      Nt_donly(ji,pi_rc_HH_dNt_j_num(Nti,Ntj,Nconj_ij)*V);
      Nconj(ij,pi_rc_HH_dNconj_num(Nti,Ntj,Nconj_ij)*V);
    }
    return Val;
  }
  return 0.0;
}

Float
REBO::pi_dh_sum(AtomsPair& ij, const Float V)
{
RETURN_REBO_0;

  Float  temp_sum = 0.0;
  AtomRefsContainer& nli = NL(ij.atom1);
  for(Index k = 0; k < nli.size(); k++)
  {
    Atom& atom_k = *(nli[k]);
    if (&atom_k == &ij.atom2 /* && atom_k.ID == C_EL*/) continue;
    if (!probablyAreNeighbours(ij.atom1,atom_k)) continue;
    AtomRefsContainer& nlj = NL(ij.atom2);
    for(Index l = 0; l < nlj.size(); l++)
    {
      Atom& atom_l = *(nlj[l]);
      if (&atom_l == &ij.atom1 /* && atom_l.ID == C_EL*/) continue;
      if (!probablyAreNeighbours(ij.atom2,atom_l)) continue;
      if (&atom_k != &atom_l) // otherwise cos=1 -> temp_sum=0
      {
        AtomsPair ik(atom_k,ij.atom1,R(0,ij.atom1,atom_k),R(1,ij.atom1,atom_k));
        if (fabs(SinTheta(ij,ik))<0.1) continue;
        AtomsPair ki(-ik);
        Float f_ik = fprime(ik);
        AtomsPair jl(atom_l,ij.atom2,R(0,ij.atom2,atom_l),R(1,ij.atom2,atom_l));
        if (fabs(SinTheta(ij,jl))<0.1) continue;
        AtomsPair lj(-jl);
        Float f_jl = fprime(jl);
        Float CosDh = - CosDihedral(ij,ki,lj);

        if (V != 0)
        {
          - CosDihedral(ij,ki,lj,+1.0*2*CosDh*f_ik*f_jl*V);
          fprime(ik,(1.0-SQR(CosDh))*f_jl*V);
          fprime(jl,(1.0-SQR(CosDh))*f_ik*V);
        }

        temp_sum += (1.0-SQR(CosDh))*f_ik*f_jl;
      }
    }
  }

  return temp_sum;
}

Float
REBO::pi_dh_Tij(AtomsPair& ij, const Float V)
{
RETURN_REBO_0;

  AtomsPair ji(-ij);

  if (ij.atom1.ID == H_EL || ij.atom2.ID == H_EL) return 0.0;

  Float Nti = Nt(ij);
  Float Ntj = Nt(ji);
  Float Nconj_ij = Nconj(ij);
  Float Tij_num_val = Tij_num(Nti,Ntj,Nconj_ij);
  if (V != 0.0)
  {
    Float Tij_dNt_i_num_val = Tij_dNt_i_num (Nti,Ntj,Nconj_ij);
    Float Tij_dNt_j_num_val = Tij_dNt_j_num (Nti,Ntj,Nconj_ij);
    Float Tij_dNconj_num_val = Tij_dNconj_num(Nti,Ntj,Nconj_ij);

    if (
      Tij_num_val == 0.0 &&
      Tij_dNt_i_num_val == 0.0 &&
      Tij_dNt_j_num_val == 0.0 &&
      Tij_dNconj_num_val == 0.0
      ) return 0.0;

    Nt_donly(ij,Tij_dNt_i_num_val*V);
    Nt_donly(ji,Tij_dNt_j_num_val*V);
    Nconj(ij,Tij_dNconj_num_val*V);
  }

  return Tij_num_val;
}

Float
REBO::pi_dh(AtomsPair& ij, const Float V)
{
  Float  TijVal = pi_dh_Tij(ij);
  if (TijVal == 0.0) return 0.0;

  Float temp_sum = pi_dh_sum(ij);

  if (V != 0.0)
  {
    pi_dh_Tij(ij,temp_sum*V);
    pi_dh_sum(ij,TijVal*V);
  }

  return TijVal*temp_sum;
}

REBO::REBO(ParamSet /*parSet*/):
  FManybody(),
  funcP_CC(0),
  funcP_CH(0),
  func_pi_rc_CC(),
  func_pi_rc_CH(),
  func_pi_rc_HH(),
  funcG_C1(),
  funcG_C2(),
  funcG_H(),
  func_Tij(),
  Nts(),
  NCs(),
  NHs()
{

  func_pi_rc_CC.init(0);
  func_pi_rc_CH.init(0);
  func_pi_rc_HH.init(0);

  funcG_C1.init();
  funcG_C2.init();
  funcG_H.init();
  
  func_Tij.init();

  handledElements.insert(H_EL);
  handledElements.insert(C_EL);
  handledElementPairs.insert(std::make_pair(H_EL,C_EL));
  handledElementPairs.insert(std::make_pair(C_EL,H_EL));
  setupPotential1();

  nl.Rcutoff = getRcutoff();
}

void
REBO::setupPotential1()
{
  PTRACE("Setup REBO");

  Q_[C][C] = 0.313460*Ao;
  Q_[C][H] = 0.340776*Ao;
    Q_[H][C] = Q_[C][H];
  Q_[H][H] = 0.370471*Ao;
  
  alpha_[C][C] = 4.7465391/Ao;
  alpha_[C][H] = 4.1025498/Ao;
    alpha_[H][C] = alpha_[C][H];
  alpha_[H][H] = 3.5362986/Ao;
  
  A_[C][C] = 10953.544*eV;
  A_[C][H] = 149.94099*eV;
    A_[H][C] = A_[C][H];
  A_[H][H] = 32.817356*eV;
  
  B_[0][C][C] = 12388.792*eV;
  B_[0][C][H] = 32.355187*eV;
    B_[0][H][C] = B_[0][C][H];
  B_[0][H][H] = 29.632593*eV;
  
  B_[1][C][C] = 17.567065*eV;
  B_[1][C][H] = 0.0*eV; // AIREBO undefined
    B_[1][H][C] = B_[1][C][H]; // AIREBO undefined
  B_[1][H][H] = 0.0*eV; // AIREBO undefined
  
  B_[2][C][C] = 30.714932*eV;
  B_[2][C][H] = 0.0*eV; // AIREBO undefined
    B_[2][H][C] = B_[2][C][H]; // AIREBO undefined
  B_[2][H][H] = 0.0*eV; // AIREBO undefined
  
  beta_[0][C][C] = 4.7204523/Ao;
  beta_[0][C][H] = 1.4344581/Ao;
    beta_[0][H][C] = beta_[0][C][H];
  beta_[0][H][H] = 1.7158922/Ao;
  
  beta_[1][C][C] = 1.4332132/Ao;
  beta_[1][C][H] = 0.0/Ao; // AIREBO undefined
    beta_[1][H][C] = beta_[1][C][H]; // AIREBO undefined
  beta_[1][H][H] = 0.0/Ao; // AIREBO undefined
  
  beta_[2][C][C] = 1.3826913/Ao;
  beta_[2][C][H] = 0.0/Ao; // AIREBO undefined
    beta_[2][H][C] = beta_[2][C][H]; // AIREBO undefined
  beta_[2][H][H] = 0.0/Ao; // AIREBO undefined
  
  rho_[C][C] = 0.0; // AIREBO undefined
  rho_[C][H] = 1.09*Ao;
    rho_[H][C] = rho_[C][H];
  rho_[H][H] = 0.7415887*Ao;
  
  R_[C][C][0] = 1.7*Ao;
  R_[H][H][0] = 1.1*Ao;
  R_[C][H][0] = 1.3*Ao;
    R_[H][C][0] = R_[C][H][0];

  R_[C][C][1] = 2.0*Ao;
  R_[H][H][1] = 1.7*Ao;
//  R_[H][H][1] = 1.5*Ao;
  R_[C][H][1] = 1.8*Ao;
    R_[H][C][1] = R_[C][H][1];

  Rprime_[C][C][0] = 1.7*Ao;
  Rprime_[H][H][0] = 1.1*Ao;
  Rprime_[C][H][0] = 1.3*Ao;
    Rprime_[H][C][0] = Rprime_[C][H][0];

  Rprime_[C][C][1] = 2.0*Ao;  REQUIRE(Rprime_[C][C][1] <= R_[C][C][1]);
  Rprime_[H][H][1] = 1.7*Ao;  REQUIRE(Rprime_[H][H][1] <= R_[H][H][1]);
//  Rprime_[H][H][1] = 1.5*Ao;  REQUIRE(Rprime_[H][H][1] <= R_[H][H][1]);
  Rprime_[C][H][1] = 1.6*Ao;  REQUIRE(Rprime_[C][H][1] <= R_[C][H][1]);
    Rprime_[H][C][1] = Rprime_[C][H][1];  REQUIRE(Rprime_[H][C][1] <= R_[H][C][1]);

  Ntot_[0] = 2;
  Ntot_[1] = 3;
    
  ao_ = 0.011304;
  co_ = 19.0;
  do_ = 2.5;

  G_HH_ = 4.0;
}  

inline
Float
REBO::P_CC_num(Float a1, Float a2) const
{
  return funcP_CC(a1,a2);
} 

inline
Float
REBO::P_CC_dH_num(Float a1, Float a2) const
{
  return funcP_CC.dH(a1,a2);
} 

inline
Float
REBO::P_CC_dC_num(Float a1, Float a2) const
{
  return funcP_CC.dC(a1,a2);
} 

inline
Float
REBO::P_CH_num(Float a1, Float a2) const
{
  return funcP_CH(a1,a2);
} 

inline
Float
REBO::P_CH_dH_num(Float a1, Float a2) const
{
  return funcP_CH.dH(a1,a2);
} 

inline
Float
REBO::P_CH_dC_num(Float a1, Float a2) const
{
  return funcP_CH.dC(a1,a2);
} 

inline
Float
REBO::pi_rc_CC_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CC(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CC_dNconj_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CC.dk(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CC_dNt_i_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CC.di(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CC_dNt_j_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CC.dj(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CH_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CH(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CH_dNconj_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CH.dk(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CH_dNt_i_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CH.di(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_CH_dNt_j_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_CH.dj(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_HH_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_HH(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_HH_dNconj_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_HH.dk(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_HH_dNt_i_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_HH.di(a1,a2,a3);
} 

inline
Float
REBO::pi_rc_HH_dNt_j_num(Float a1, Float a2, Float a3) const 
{
  return func_pi_rc_HH.dj(a1,a2,a3);
} 

inline
Float
REBO::Tij_num(Float a1, Float a2, Float a3) const 
{
  return func_Tij(a1,a2,a3);
} 

inline
Float
REBO::Tij_dNconj_num(Float a1, Float a2, Float a3) const 
{
  return func_Tij.dk(a1,a2,a3);
} 

inline
Float
REBO::Tij_dNt_i_num(Float a1, Float a2, Float a3) const 
{
  return func_Tij.di(a1,a2,a3);
} 

inline
Float
REBO::Tij_dNt_j_num(Float a1, Float a2, Float a3) const 
{
  return func_Tij.dj(a1,a2,a3);
} 

}


