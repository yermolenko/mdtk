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

  Float  REBO::buildPairs(AtomsArray& gl)
  {
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
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        std::pair<int,int> sample_pair(atom_i.globalIndex,atom_j.globalIndex);
        if (&atom_i != &atom_j)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
        {

    currentPairPtr = &sample_pair;
    ontouch_enabled = true;

#ifdef REBO_OPTIMIZED  
        Float VAvar = VA(atom_i,atom_j);
        if (VAvar != 0.0)
          Ei += VR(atom_i,atom_j) + Baver(atom_i,atom_j)*VAvar;
        else
          Ei += VR(atom_i,atom_j);
#else
        Ei += VR(atom_i,atom_j) + Baver(atom_i,atom_j)*VA(atom_i,atom_j);
#endif      

    currentPairPtr = NULL;
    ontouch_enabled = false;
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
REBO::fprime(Atom &atom1,Atom &atom2)
{
  Float r = r_vec_module(atom1,atom2);

  Float R1=Rprime(0,atom1,atom2);
  Float R2=Rprime(1,atom1,atom2);

  return Sprime(r,R1,R2);
}

inline
Vector3D
REBO::dfprime(Atom &atom1,Atom &atom2, Atom &datom)
{
  if (&datom != &atom1 && &datom != &atom2) return 0.0;

  Float r = r_vec_module(atom1,atom2);

  Float R1=Rprime(0,atom1,atom2);
  Float R2=Rprime(1,atom1,atom2);

#ifdef FGENERAL_OPTIMIZED  
  Vector3D dvar = dr_vec_module(atom1,atom2,datom);
  if (dvar != 0.0)
    return dSprime(r,R1,R2)*dvar;
  else
    return 0.0;
#else
  return dSprime(r,R1,R2)*dr_vec_module(atom1,atom2,datom);  
#endif
}

inline
Float
REBO::VR(Atom &atom1,Atom &atom2)
{
  return VR_Exp(atom1, atom2);
}  

inline
Vector3D
REBO::dVR(Atom &atom1,Atom &atom2, Atom &datom)
{
  return dVR_Exp(atom1, atom2, datom);
}
  
inline
Float
REBO::VR_Exp(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef REBO_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float Q=     this->Q(atom1,atom2);
  Float A=     this->A(atom1,atom2);
  Float alpha= this->alpha(atom1,atom2);
  
  return fvar*(1+Q/r)*A*exp(-alpha*r);
}

inline
Vector3D
REBO::dVR_Exp(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef REBO_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float Q=     this->Q(atom1,atom2);
  Float A=     this->A(atom1,atom2);
  Float alpha= this->alpha(atom1,atom2);

  Float t1 = A*exp(-alpha*r);

  return (
           df(atom1,atom2,datom)
           *(1+Q/r)*t1
           +
           f(atom1,atom2)       
           *(-t1*(Q+alpha*r*r+alpha*r*Q)/(r*r))*dr_vec_module(atom1,atom2,datom)
         );
}

inline
Float
REBO::VA(Atom &atom1,Atom &atom2)
{
  Float fvar = f(atom1,atom2);

#ifdef REBO_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float B1=   this->B(1,atom1,atom2);
  Float B2=   this->B(2,atom1,atom2);
  Float B3=   this->B(3,atom1,atom2);

  Float beta1=   this->beta(1,atom1,atom2);
  Float beta2=   this->beta(2,atom1,atom2);
  Float beta3=   this->beta(3,atom1,atom2);
  
  return -fvar*(
          B1*exp(-beta1*r)+
          B2*exp(-beta2*r)+
          B3*exp(-beta3*r)
          );
}

inline
Vector3D
REBO::dVA(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef REBO_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float B1=   this->B(1,atom1,atom2);
  Float B2=   this->B(2,atom1,atom2);
  Float B3=   this->B(3,atom1,atom2);

  Float beta1=   this->beta(1,atom1,atom2);
  Float beta2=   this->beta(2,atom1,atom2);
  Float beta3=   this->beta(3,atom1,atom2);

  Float t1 = B1*exp(-beta1*r);
  Float t2 = B2*exp(-beta2*r);
  Float t3 = B3*exp(-beta3*r);
  
  return -dfvar*(
          t1+
          t2+
          t3
          )
         -f(atom1,atom2)*(
          (-beta1)*t1+
          (-beta2)*t2+
          (-beta3)*t3
         )*drmodvar;
}

Float
REBO::operator()(AtomsArray& gl)
{
  return buildPairs(gl);
}  

 
 
Vector3D
REBO::grad(Atom &atom,AtomsArray &gl)
{
  Vector3D dEi(0.0,0.0,0.0);

  if (isHandled(atom))
  {
    Index i;

    std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];

    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = gl[acnt[i].first];
      {
        Atom &atom_j = gl[acnt[i].second];

        REQUIREM(&atom_j != &atom_i,"must be (&atom_j != &atom_i)");
#if defined(REBO_OPTIMIZED_EVEN_BETTER)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
#endif
        {
          Vector3D dEij;
#ifdef REBO_OPTIMIZED  
          dEij = dVR(atom_i,atom_j,atom);
          Float VAvar = VA(atom_i,atom_j);
          Vector3D dVAvar = dVA(atom_i,atom_j,atom);

          if (VAvar != 0.0) dEij += dBaver(atom_i,atom_j,atom)* VAvar;
          if (dVAvar != 0.0) dEij += Baver(atom_i,atom_j)*dVAvar;
#else
          dEij = 
                 (dVR(atom_i,atom_j,atom) + 
                (dBaver(atom_i,atom_j,atom)* VA(atom_i,atom_j) + 
                  Baver(atom_i,atom_j)*dVA(atom_i,atom_j,atom)));
#endif
          dEi += dEij;

        }  
      }
    }    
  }  

  return  dEi;
}  

Float
REBO::Baver(Atom &atom1,Atom &atom2)
{
  return (p_sigma_pi(atom1,atom2)+p_sigma_pi(atom2,atom1))/2.0
         +pi_rc(atom1,atom2)
#ifdef REBO_DIHEDRAL
         +pi_dh(atom1,atom2)
#endif
         ;   
}

Vector3D
REBO::dBaver(Atom &atom1,Atom &atom2, Atom &datom)
{
  return (dp_sigma_pi(atom1,atom2,datom)+dp_sigma_pi(atom2,atom1,datom))/2.0
         +dpi_rc(atom1,atom2,datom)
#ifdef REBO_DIHEDRAL
         +dpi_dh(atom1,atom2,datom)
#endif
         ; 
}

inline
Float
REBO::ExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k)
{
  if (atom_i.ID != H_EL) return 1.0;

  Float r_ij = r_vec_module(atom_i,atom_j);
  Float r_ik = r_vec_module(atom_i,atom_k);

  return exp(
             (4.0/Ao)*(
             (rho(atom_i, atom_k) - r_ik)
            -(rho(atom_i, atom_j) - r_ij)
             )
            );
}

inline
Vector3D
REBO::dExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom)
{
  if (atom_i.ID != H_EL) return 0.0;

  return  ExpTerm(atom_i,atom_j,atom_k)
          *(
            (4.0/Ao)
            *(
              dr_vec_module(atom_i,atom_j,datom)-dr_vec_module(atom_i,atom_k,datom)
             )
           );
}

Float
REBO::D(Atom &atom1,Atom &atom2)
{
  Index k;
  Float Dij = 1.0;
  for(k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
    if (&atom_k != &atom2/* && &atom_k != &atom1*/)
    {
#ifdef REBO_OPTIMIZED  
      Float fvar = f(atom1,atom_k);
      if (fvar != 0.0)
      {
        Float ExpTermvar = ExpTerm(atom1,atom2,atom_k);
        if (ExpTermvar != 0.0)
          Dij += G(atom1,atom2,atom_k)*fvar*ExpTermvar;
      }    
#else
      Dij += G(atom1,atom2,atom_k)*f(atom1,atom_k)*ExpTerm(atom1,atom2,atom_k);
#endif
    }  
  }  
  Dij += P(atom1,atom2);  
  return Dij;
}

Vector3D
REBO::dD(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D DerD = 0.0;
  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
    if (&atom_k != &atom2/* && &atom_k != &atom1*/)
#ifdef REBO_OPTIMIZED_EVEN_BETTER
    if (r_vec_module(atom1,atom_k) < R(1,atom1,atom_k))
#endif
    {
#ifdef REBO_OPTIMIZED  
      Float ExpTermvar = ExpTerm(atom1,atom2,atom_k);
      Vector3D dExpTermvar = dExpTerm(atom1,atom2,atom_k,datom);

      if (ExpTermvar != 0.0 || dExpTermvar != 0.0)
      {
      Float fvar = f(atom1,atom_k);
      Float Gvar = G(atom1,atom2,atom_k);
      Vector3D dfvar = df(atom1,atom_k,datom);
      Vector3D dGvar = dG(atom1,atom2,atom_k,datom);

      DerD += dGvar* fvar* ExpTermvar;
      DerD +=  Gvar*dfvar* ExpTermvar;
      DerD +=  Gvar* fvar*dExpTermvar;
      }  
#else
      DerD += dG(atom1,atom2,atom_k,datom)* f(atom1,atom_k)* ExpTerm(atom1,atom2,atom_k);
      DerD +=  G(atom1,atom2,atom_k)*df(atom1,atom_k,datom)* ExpTerm(atom1,atom2,atom_k);
      DerD +=  G(atom1,atom2,atom_k)* f(atom1,atom_k)*dExpTerm(atom1,atom2,atom_k,datom);
#endif
    }  
  }  
  DerD += dP(atom1,atom2,datom); 

  return DerD;
}

inline
Float
REBO::p_sigma_pi(Atom &atom1,Atom &atom2)
{
  return pow(D(atom1,atom2),-0.5);
}

inline
Vector3D
REBO::dp_sigma_pi(Atom &atom1,Atom &atom2, Atom &datom)
{
  return -0.5*pow(D(atom1,atom2),-0.5-1.0)*dD(atom1,atom2,datom);
}

inline
Float
REBO::G(Atom &atom_i,Atom &atom_j,Atom &atom_k)
{
  Float CosT = CosThetaJIK(atom_i,atom_j,atom_k);

  if (atom_i.ID == C_EL)
  {
    Float Nt1=Ntot_[0];
    Float Nt2=Ntot_[1];

    Float N = Nt(atom_i,atom_j);

    return funcG_C2(CosT)+Sprime(N,Nt1,Nt2)*(funcG_C1(CosT) - funcG_C2(CosT));    
  }
  else if (atom_i.ID == H_EL)
  {
    return funcG_H(CosT);
  }  
  else
  {
    throw Exception("REBO::G() : unknown element");
  }  
}


inline
Vector3D
REBO::dG(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom)
{
  Float CosT = CosThetaJIK(atom_i,atom_j,atom_k);
  Vector3D dCosT = dCosThetaJIK(atom_i,atom_j,atom_k,datom);

  if (atom_i.ID == C_EL)
  {
    Float Nt1=Ntot_[0];
    Float Nt2=Ntot_[1];

    Float N = Nt(atom_i,atom_j);
    Vector3D dN = dNt(atom_i,atom_j,datom);

    return funcG_C2.dCosT(CosT)*dCosT
           + 
           (
           (dN!=0.0)?
           (
           dSprime(N,Nt1,Nt2)*dN*(funcG_C1(CosT) - funcG_C2(CosT))
           ):(0.0)
           )
           +
           Sprime(N,Nt1,Nt2)*(funcG_C1.dCosT(CosT) - funcG_C2.dCosT(CosT))*dCosT
           ;    
  }
  else if (atom_i.ID == H_EL)
  {
    return funcG_H.dCosT(CosT)*dCosT;
  }  
  else
  {
    throw Exception("REBO::dG() : unknown element");
  }  
}

Float
REBO::Nt(Atom &atom_i, Atom &atom_j)
{
  Float Nt_i = 0;

  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (&atom_k != &atom_j)
    {
      Nt_i += f(atom_i,atom_k);
    }  
  }     

  return Nt_i;
} 

Vector3D
REBO::dNt(Atom &atom_i, Atom &atom_j, Atom &datom)
{
  Vector3D Nt_i = 0;

  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (&atom_k != &atom_j)
    {
      Nt_i += df(atom_i,atom_k,datom);
    }  
  }     

  return Nt_i;
}

Float
REBO::NH(Atom &atom_i, Atom &atom_j)
{ 
  Float NH_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (atom_k.ID == H_EL && &atom_k != &atom_j)
    {
      NH_i += f(atom_i,atom_k);
    }  
  }    

  return NH_i;
}

Vector3D
REBO::dNH(Atom &atom_i, Atom &atom_j, Atom &datom)
{ 
  Vector3D NH_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (atom_k.ID == H_EL && &atom_k != &atom_j)
    {
      NH_i += df(atom_i,atom_k,datom);
    }  
  }    

  return NH_i;
}

Float
REBO::NC(Atom &atom_i, Atom &atom_j)
{
  Float NC_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (atom_k.ID == C_EL && &atom_k != &atom_j)
    {
      NC_i += f(atom_i,atom_k);
    }  
  }    

  return NC_i;
}

Vector3D
REBO::dNC(Atom &atom_i, Atom &atom_j, Atom &datom)
{
  Vector3D NC_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (atom_k.ID == C_EL && &atom_k != &atom_j)
    {
      NC_i += df(atom_i,atom_k,datom);
    }  
  }    

  return NC_i;
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
REBO::Nconj(Atom &atom1,Atom &atom2)
{
  Float Nconj_ij = 1.0;

  Float sum1 = 0.0;
  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
    if (&atom_k != &atom2 && 
        atom_k.ID == C_EL)
    {
      Float f_ik = f(atom_k,atom1);
#ifdef REBO_OPTIMIZED_EVEN_BETTER
      if (f_ik == 0.0) continue;
#endif

      Float N   = Nt(atom_k,atom1);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];
      sum1 += f_ik*Sprime(N,Nt1,Nt2);
    }  
  }    
  Nconj_ij += SQR(sum1);

  Float sum2 = 0.0;
  for(Index l = 0; l < NL(atom2).size(); l++)
  {
    Atom& atom_l = *(NL(atom2)[l]);
    if (&atom_l != &atom1 &&
        atom_l.ID == C_EL)
    {
      Float f_jl = f(atom_l,atom2);
#ifdef REBO_OPTIMIZED_EVEN_BETTER
      if (f_jl == 0.0) continue;
#endif

      Float N   = Nt(atom_l,atom2);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];
      sum2 += f_jl*Sprime(N,Nt1,Nt2);
    }  
  }    
  Nconj_ij += SQR(sum2);

  return Nconj_ij; 
}


Vector3D
REBO::dNconj(Atom &atom1,Atom &atom2, Atom &datom)
{
  Vector3D dNconj = 0.0;

  Float    sum1 = 0.0;
  Vector3D dsum1 = 0.0;
  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);

    if (&atom_k != &atom2 && 
        atom_k.ID == C_EL)
    {
      Float f_ik = f(atom_k,atom1);
#ifdef REBO_OPTIMIZED_EVEN_BETTER
      if (f_ik == 0.0) continue;
#endif
      Vector3D df_ik = df(atom_k,atom1,datom);
      
      Float N   = Nt(atom_k,atom1);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];
      Vector3D dN  = dNt(atom_k,atom1,datom);

       sum1 += f_ik*Sprime(N,Nt1,Nt2);

      dsum1 += df_ik* Sprime(N,Nt1,Nt2)+
               (
               (dN!=0.0)?(f_ik*dSprime(N,Nt1,Nt2)*dN):(0.0)
               );
    }  
  }    
  dNconj += 2.0*sum1*dsum1;

  Float     sum2 = 0.0;
  Vector3D dsum2 = 0.0;
  for(Index l = 0; l < NL(atom2).size(); l++)
  {
    Atom& atom_l = *(NL(atom2)[l]);

    if (&atom_l != &atom1 &&
        atom_l.ID == C_EL)
    {
      Float f_jl = f(atom_l,atom2);
#ifdef REBO_OPTIMIZED_EVEN_BETTER
      if (f_jl == 0.0) continue;
#endif
      Vector3D df_jl = df(atom_l,atom2,datom);

      Float N   = Nt(atom_l,atom2);
      Float Nt1 = Ntot_[0];
      Float Nt2 = Ntot_[1];
      Vector3D dN  = dNt(atom_l,atom2,datom);

       sum2 +=  f_jl* Sprime(N,Nt1,Nt2);
      dsum2 += df_jl* Sprime(N,Nt1,Nt2)+
               (
               (dN!=0.0)?(f_jl*dSprime(N,Nt1,Nt2)*dN):(0.0)
               );
    }
  }
  dNconj += 2.0*sum2*dsum2;

  return dNconj;
}

inline
Float
REBO::P(Atom &atom_i, Atom &atom_j)
{
RETURN_REBO_0;

  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
    return P_CC_num(NH(atom_i,atom_j),NC(atom_i,atom_j));
  else if (atom_i.ID == C_EL && atom_j.ID == H_EL)
    return P_CH_num(NH(atom_i,atom_j),NC(atom_i,atom_j));  
  else return 0.0;  
}

inline
Vector3D
REBO::dP(Atom &atom_i, Atom &atom_j, Atom &datom)
{
RETURN_REBO_0;

#ifdef REBO_OPTIMIZED  
  Vector3D RESULT = 0.0;
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
  {
    Float NHvar = NH(atom_i,atom_j);
    Float NCvar = NC(atom_i,atom_j);

    Float p1 = P_CC_dH_num (NHvar,NCvar);
    Float p2 = P_CC_dC_num (NHvar,NCvar);

    if (p1 != 0.0) 
      RESULT += p1*dNH(atom_i,atom_j,datom);
    if (p2 != 0.0) 
      RESULT += p2*dNC(atom_i,atom_j,datom);
  }  
  else if (atom_i.ID == C_EL && atom_j.ID == H_EL)
  {
    Float NHvar = NH(atom_i,atom_j);
    Float NCvar = NC(atom_i,atom_j);

    Float p1 = P_CH_dH_num (NHvar,NCvar);
    Float p2 = P_CH_dC_num (NHvar,NCvar);

    if (p1 != 0.0) 
      RESULT += p1*dNH(atom_i,atom_j,datom);
    if (p2 != 0.0) 
      RESULT += p2*dNC(atom_i,atom_j,datom);
  }
  else RESULT = 0.0;
  return RESULT;
#else
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
    return P_CC_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNH(atom_i,atom_j,datom)+
           P_CC_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNC(atom_i,atom_j,datom);
  else if (atom_i.ID == C_EL && atom_j.ID == H_EL)         
    return P_CH_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNH(atom_i,atom_j,datom)+
           P_CH_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNC(atom_i,atom_j,datom);
  else return 0.0;
#endif  
}

inline
Float
REBO::pi_rc(Atom &atom_i, Atom &atom_j)
{
RETURN_REBO_0;
  
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
    return pi_rc_CC_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j));
  if (atom_i.ID == C_EL && atom_j.ID == H_EL)
    return pi_rc_CH_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j));
  if (atom_i.ID == H_EL && atom_j.ID == C_EL)
  {
    Atom& atom_i_new = atom_j;
    Atom& atom_j_new = atom_i;
    return pi_rc(atom_i_new,atom_j_new);
  };
  if (atom_i.ID == H_EL && atom_j.ID == H_EL)
    return pi_rc_HH_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j));
  return 0.0;
}

inline
Vector3D
REBO::dpi_rc(Atom &atom_i, Atom &atom_j, Atom &datom)
{
RETURN_REBO_0;
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
#ifdef REBO_OPTIMIZED  
{
  Vector3D RESULT = 0.0;
  {  
    Float Ntijvar = Nt(atom_i,atom_j);
    Float Ntjivar = Nt(atom_j,atom_i);
    Float Nconjvar = Nconj(atom_i,atom_j);

    Float p1 = pi_rc_CC_dNt_i_num (Ntijvar,Ntjivar,Nconjvar);
    Float p2 = pi_rc_CC_dNt_j_num (Ntijvar,Ntjivar,Nconjvar);
    Float p3 = pi_rc_CC_dNconj_num(Ntijvar,Ntjivar,Nconjvar);

    if (p1 != 0.0) 
      RESULT += p1*dNt   (atom_i,atom_j,datom);
    if (p2 != 0.0) 
      RESULT += p2*dNt   (atom_j,atom_i,datom);
    if (p3 != 0.0)
      RESULT += p3*dNconj(atom_i,atom_j,datom);
  }  
  return RESULT;
}  
#else
    return pi_rc_CC_dNt_i_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_i,atom_j,datom)+
           pi_rc_CC_dNt_j_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_j,atom_i,datom)+
           pi_rc_CC_dNconj_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNconj(atom_i,atom_j,datom);
#endif
  if (atom_i.ID == C_EL && atom_j.ID == H_EL)
#ifdef REBO_OPTIMIZED  
{
  Vector3D RESULT = 0.0;
  {  
    Float Ntijvar = Nt(atom_i,atom_j);
    Float Ntjivar = Nt(atom_j,atom_i);
    Float Nconjvar = Nconj(atom_i,atom_j);

    Float p1 = pi_rc_CH_dNt_i_num (Ntijvar,Ntjivar,Nconjvar);
    Float p2 = pi_rc_CH_dNt_j_num (Ntijvar,Ntjivar,Nconjvar);
    Float p3 = pi_rc_CH_dNconj_num(Ntijvar,Ntjivar,Nconjvar);

    if (p1 != 0.0) 
      RESULT += p1*dNt   (atom_i,atom_j,datom);
    if (p2 != 0.0) 
      RESULT += p2*dNt   (atom_j,atom_i,datom);
    if (p3 != 0.0)
      RESULT += p3*dNconj(atom_i,atom_j,datom);
  }  
  return RESULT;
}  
#else
    return pi_rc_CH_dNt_i_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_i,atom_j,datom)+
           pi_rc_CH_dNt_j_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_j,atom_i,datom)+
           pi_rc_CH_dNconj_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNconj(atom_i,atom_j,datom);
#endif
  if (atom_i.ID == H_EL && atom_j.ID == C_EL)
  {
    Atom& atom_i_new = atom_j;
    Atom& atom_j_new = atom_i;
    return dpi_rc(atom_i_new,atom_j_new,datom);
  }; 
  if (atom_i.ID == H_EL && atom_j.ID == H_EL)
#ifdef REBO_OPTIMIZED  
{
  Vector3D RESULT = 0.0;
  {  
    Float Ntijvar = Nt(atom_i,atom_j);
    Float Ntjivar = Nt(atom_j,atom_i);
    Float Nconjvar = Nconj(atom_i,atom_j);

    Float p1 = pi_rc_HH_dNt_i_num (Ntijvar,Ntjivar,Nconjvar);
    Float p2 = pi_rc_HH_dNt_j_num (Ntijvar,Ntjivar,Nconjvar);
    Float p3 = pi_rc_HH_dNconj_num(Ntijvar,Ntjivar,Nconjvar);

    if (p1 != 0.0) 
      RESULT += p1*dNt   (atom_i,atom_j,datom);
    if (p2 != 0.0) 
      RESULT += p2*dNt   (atom_j,atom_i,datom);
    if (p3 != 0.0)
      RESULT += p3*dNconj(atom_i,atom_j,datom);
  }  
  return RESULT;
}  
#else
    return pi_rc_HH_dNt_i_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_i,atom_j,datom)+
           pi_rc_HH_dNt_j_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_j,atom_i,datom)+
           pi_rc_HH_dNconj_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNconj(atom_i,atom_j,datom);
#endif

  return 0.0;
}

Float
REBO::pi_dh(Atom &atom_i, Atom &atom_j)
{
RETURN_REBO_0;

  if (atom_i.ID == H_EL || atom_j.ID == H_EL) return 0.0;

  Float  pi_dh_val = 1.0;
  pi_dh_val *= Tij_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j));
  if (pi_dh_val == 0.0) return 0.0;

  Float  temp_sum = 0.0;
  for(Index k = 0; k < NL(atom_i).size(); k++)
  for(Index l = 0; l < NL(atom_j).size(); l++)
  {
    Atom& atom_k = *(NL(atom_i)[k]);
    Atom& atom_l = *(NL(atom_j)[l]);
    if (&atom_k != &atom_l) // otherwise cos=1 -> temp_sum=0
    if (&atom_k != &atom_j /* && atom_k.ID == C_EL*/)
    if (&atom_l != &atom_i /* && atom_l.ID == C_EL*/)
#ifdef REBO_OPTIMIZED_EVEN_BETTER
    if (probablyAreNeighbours(atom_i,atom_k))
#endif
#ifdef REBO_OPTIMIZED_EVEN_BETTER
    if (probablyAreNeighbours(atom_j,atom_l))
#endif
    {
      if (fabs(SinTheta(atom_i,atom_j,atom_k))<0.1) continue;
      if (fabs(SinTheta(atom_i,atom_j,atom_l))<0.1) continue;
      Float f_ik = fprime(atom_k,atom_i);
      Float f_jl = fprime(atom_l,atom_j);
      
      Float CosDh = - CosDihedral(atom_i,atom_j,atom_k,atom_l);

      temp_sum += (1.0-SQR(CosDh))
                  *f_ik 
                  *f_jl 
/*                  *Heaviside(sin_jik-smin)
                  *Heaviside(sin_ijl-smin)*/;
    }  
  }    

  pi_dh_val *= temp_sum;

  return pi_dh_val;
}

Vector3D
REBO::dpi_dh(Atom &atom_i, Atom &atom_j, Atom &datom)
{
RETURN_REBO_0;

  if (atom_i.ID == H_EL || atom_j.ID == H_EL) return 0.0;

  Vector3D  dpi_dh_val = 0.0;

    Float Nt_i = Nt(atom_i,atom_j);
    Float Nt_j = Nt(atom_j,atom_i);
    Float Nconj_ij = Nconj(atom_i,atom_j);

    Float Tij_num_val = Tij_num(Nt_i,Nt_j,Nconj_ij);

    Float Tij_dNt_i_num_val = Tij_dNt_i_num (Nt_i,Nt_j,Nconj_ij);
    Float Tij_dNt_j_num_val = Tij_dNt_j_num (Nt_i,Nt_j,Nconj_ij);
    Float Tij_dNconj_num_val = Tij_dNconj_num(Nt_i,Nt_j,Nconj_ij);

    if (
        Tij_num_val == 0.0 &&
        Tij_dNt_i_num_val == 0.0 &&
        Tij_dNt_j_num_val == 0.0 &&
        Tij_dNconj_num_val == 0.0
       ) return 0.0;

  Float  temp_sum = 0.0;
  Vector3D  dtemp_sum = 0.0;
  for(Index k = 0; k < NL(atom_i).size(); k++)
  for(Index l = 0; l < NL(atom_j).size(); l++)
  {
    Atom& atom_k = *(NL(atom_i)[k]);
    Atom& atom_l = *(NL(atom_j)[l]);
    if (&atom_k != &atom_l) // otherwise cos=1 -> temp_sum=0
    if (&atom_k != &atom_j /* && atom_k.ID == C_EL*/)
    if (&atom_l != &atom_i /* && atom_l.ID == C_EL*/)
#ifdef REBO_OPTIMIZED_EVEN_BETTER
    if (probablyAreNeighbours(atom_i,atom_k))
#endif
#ifdef REBO_OPTIMIZED_EVEN_BETTER
    if (probablyAreNeighbours(atom_j,atom_l))
#endif
    {
      if (fabs(SinTheta(atom_i,atom_j,atom_k))<0.1) continue;
      if (fabs(SinTheta(atom_i,atom_j,atom_l))<0.1) continue;
      Float f_ik = fprime(atom_k,atom_i);
      Float f_jl = fprime(atom_l,atom_j);

      Float CosDh = - CosDihedral(atom_i,atom_j,atom_k,atom_l);

      temp_sum += (1.0-SQR(CosDh))  *f_ik   *f_jl 
/*                  *Heaviside(sin_jik-smin)
                  *Heaviside(sin_ijl-smin)*/;
                  
      Vector3D df_ik = dfprime(atom_k,atom_i,datom);
      Vector3D df_jl = dfprime(atom_l,atom_j,datom);
      
      Vector3D dCosDh = - dCosDihedral(atom_i,atom_j,atom_k,atom_l,datom);

      dtemp_sum += (-2.0*CosDh*dCosDh) * f_ik   * f_jl 
                  +(1.0-SQR(CosDh)   ) *df_ik   * f_jl 
                  +(1.0-SQR(CosDh)   ) * f_ik   *df_jl 
/*                  *Heaviside(sin_jik-smin)
                  *Heaviside(sin_ijl-smin)*/;
    }  
  }    

    dpi_dh_val += 
            Tij_num_val
            *dtemp_sum;
  if (temp_sum != 0.0)
    dpi_dh_val += 
            (
             Tij_dNt_i_num_val*dNt   (atom_i,atom_j,datom)+
             Tij_dNt_j_num_val*dNt   (atom_j,atom_i,datom)+
             Tij_dNconj_num_val*dNconj(atom_i,atom_j,datom)
            )
            *temp_sum;

  return dpi_dh_val;
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
  func_Tij()
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


