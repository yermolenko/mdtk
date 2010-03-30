/*
   Implementation of the many-body interatomic potential for
   hydrocarbons.
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009 Oleksandr Yermolenko
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

#include "Brenner.hpp"
#include <algorithm>

namespace mdtk
{

  Float  Brenner::buildPairs(AtomsContainer& gl)
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
      Atom &atom_i = *(gl[ii]);
      if (isHandled(atom_i))
      for(size_t jj = 0; jj < NL(atom_i).size(); jj++)
      {
        Atom &atom_j = *(NL(atom_i)[jj]);
        if (atom_i.globalIndex > atom_j.globalIndex) continue;
        std::pair<int,int> sample_pair(atom_i.globalIndex,atom_j.globalIndex);
        if (isHandled(atom_j))
        if (&atom_i != &atom_j)
        if (r_vec_module(atom_i,atom_j) < R(1,atom_i,atom_j))
        {
    currentPairPtr = &sample_pair;
    ontouch_enabled = true;

#ifdef BRENNER_OPTIMIZED  
        Float VAvar = VA(atom_i,atom_j);
        if (VAvar != 0.0)
          Ei += VR(atom_i,atom_j) - Baver(atom_i,atom_j)*VAvar;
        else
          Ei += VR(atom_i,atom_j);
#else
        Ei += VR(atom_i,atom_j) - Baver(atom_i,atom_j)*VA(atom_i,atom_j);
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
Brenner::VR(Atom &atom1,Atom &atom2) 
{
  return VR_Exp(atom1, atom2);
}  

inline
Vector3D
Brenner::dVR(Atom &atom1,Atom &atom2, Atom &datom) 
{
  return dVR_Exp(atom1, atom2, datom);
}
  
inline
Float
Brenner::VR_Exp(Atom &atom1,Atom &atom2) 
{
  Float fvar = f(atom1,atom2);

#ifdef BRENNER_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float De=   this->De(atom1,atom2);
  Float S=    this->S(atom1,atom2);
  Float beta= this->beta(atom1,atom2);
  Float Re=   this->Re(atom1,atom2);
  
  Float a = De/(S-1.0);
  Float b = sqrt(2.0*S)*beta;
  Float c = Re;

  return fvar*a*exp(-b*(r-c));
}

inline
Vector3D
Brenner::dVR_Exp(Atom &atom1,Atom &atom2, Atom &datom) 
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef BRENNER_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float De=   this->De(atom1,atom2);
  Float S=    this->S(atom1,atom2);
  Float beta= this->beta(atom1,atom2);
  Float Re=   this->Re(atom1,atom2);
  
  Float a = De/(S-1.0);
  Float b = sqrt(2.0*S)*beta;
  Float c = Re;

  return a*exp(-b*(r-c))*(dfvar-b*drmodvar*f(atom1,atom2));
}

inline
Float
Brenner::VA(Atom &atom1,Atom &atom2) 
{
  return VA_Exp(atom1,atom2);
}  

inline
Float
Brenner::VA_Exp(Atom &atom1,Atom &atom2) 
{
  Float fvar = f(atom1,atom2);

#ifdef BRENNER_OPTIMIZED  
  if (fvar == 0.0) return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float De=   this->De(atom1,atom2);
  Float S=    this->S(atom1,atom2);
  Float beta= this->beta(atom1,atom2);
  Float Re=   this->Re(atom1,atom2);
  
  Float a = De*S/(S-1.0);
  Float b = sqrt(2.0/S)*beta;
  Float c = Re;

  return fvar*a*exp(-b*(r-c));
}

inline
Vector3D
Brenner::dVA(Atom &atom1,Atom &atom2, Atom &datom) 
{
  return dVA_Exp(atom1,atom2, datom);
}
  
inline
Vector3D
Brenner::dVA_Exp(Atom &atom1,Atom &atom2, Atom &datom) 
{
  Vector3D dfvar = df(atom1,atom2,datom);
  Vector3D drmodvar = dr_vec_module(atom1,atom2,datom);

#ifdef BRENNER_OPTIMIZED  
  if (dfvar == 0.0 && drmodvar == 0.0)  return 0.0;
#endif

  Float r = r_vec_module(atom1,atom2);

  Float De=   this->De(atom1,atom2);
  Float S=    this->S(atom1,atom2);
  Float beta= this->beta(atom1,atom2);
  Float Re=   this->Re(atom1,atom2);
  
  Float a = De*S/(S-1.0);
  Float b = sqrt(2.0/S)*beta;
  Float c = Re;

  return a*exp(-b*(r-c))*(dfvar-b*drmodvar*f(atom1,atom2));
}

Float
Brenner::operator()(AtomsContainer& gl) 
{
  return buildPairs(gl);
}  

static size_t Bij_count = 0;
 
Vector3D
Brenner::grad(Atom &atom,AtomsContainer&gl) 
{
  Bij_count = 0;
  
  Index i;
  
  Vector3D dEi(0.0,0.0,0.0);

  if (isHandled(atom))
  {

  std::vector<std::pair<int,int> >& acnt = pairs[atom.globalIndex];

    for(i = 0; i < acnt.size(); i++)
    {
      Atom &atom_i = *(gl[acnt[i].first]);
      if (isHandled(atom_i))
      {
        Atom &atom_j = *(gl[acnt[i].second]);

        REQUIREM(&atom_j != &atom_i,"must be (&atom_j != &atom_i)");
        if (isHandled(atom_j))
#if defined(BRENNER_OPTIMIZED_EVEN_BETTER)
        if (r_vec_module_no_touch(atom_i,atom_j) < R(1,atom_i,atom_j))
#endif
        {
          Vector3D dEij;
#ifdef BRENNER_OPTIMIZED  
          dEij = dVR(atom_i,atom_j,atom);
          Float VAvar = VA(atom_i,atom_j);
          Vector3D dVAvar = dVA(atom_i,atom_j,atom);
          if (VAvar != 0.0) dEij -= dBaver(atom_i,atom_j,atom)* VAvar;
          if (dVAvar != 0.0) dEij -= Baver(atom_i,atom_j)*dVAvar;
          ++Bij_count;
#else
          dEij = 
                 (dVR(atom_i,atom_j,atom) - 
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

inline
Float
Brenner::Baver(Atom &atom1,Atom &atom2) 
{
  return (B(atom1,atom2)+B(atom2,atom1))/2.0
         + 
         (
          FN(atom1,atom2)/2.0
         );   

}

Vector3D
Brenner::dBaver(Atom &atom1,Atom &atom2, Atom &datom) 
{
  return (dB(atom1,atom2,datom)+dB(atom2,atom1,datom))/2.0
         + 
         (
          dFN(atom1,atom2,datom)/2.0
         ); 
}

inline
Float
Brenner::ExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k) 
{
#ifdef BRENNER_OPTIMIZED  
  if (alpha(atom_i,atom_j,atom_k) == 0.0) return 1.0;
#endif
  Float r_ij = r_vec_module(atom_i,atom_j);
  Float r_ik = r_vec_module(atom_i,atom_k);

  return exp(alpha(atom_i,atom_j,atom_k)*
      ((r_ij-Re(atom_i,atom_j))-(r_ik-Re(atom_i,atom_k)))); 
}

inline
Vector3D
Brenner::dExpTerm(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) 
{
#ifdef BRENNER_OPTIMIZED  
  if (alpha(atom_i,atom_j,atom_k) == 0.0) return 0;
#endif
  return  ExpTerm(atom_i,atom_j,atom_k)
          *(
            alpha(atom_i,atom_j,atom_k)
            *(
              dr_vec_module(atom_i,atom_j,datom)-dr_vec_module(atom_i,atom_k,datom)
             )
           );
}

Float
Brenner::D(Atom &atom1,Atom &atom2) 
{
  Index k;
  Float Dij = 1.0;
  for(k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
    if (isHandled(atom_k))
    if (&atom_k != &atom2/* && &atom_k != &atom1*/)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom1,atom_k) < R(1,atom1,atom_k))
#endif
    {
#ifdef BRENNER_OPTIMIZED  
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
  Dij += HN(atom1,atom2);  
  return Dij;
}

Vector3D
Brenner::dD(Atom &atom1,Atom &atom2, Atom &datom) 
{
  Vector3D DerD = 0.0;
  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
    if (isHandled(atom_k))
    if (&atom_k != &atom2/* && &atom_k != &atom1*/)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom1,atom_k) < R(1,atom1,atom_k))
#endif
    {
#ifdef BRENNER_OPTIMIZED  
      Float fvar = f(atom1,atom_k);
      Float Gvar = G(atom1,atom2,atom_k);
      Float ExpTermvar = ExpTerm(atom1,atom2,atom_k);
      Vector3D dfvar = df(atom1,atom_k,datom);
      Vector3D dGvar = dG(atom1,atom2,atom_k,datom);
      Vector3D dExpTermvar = dExpTerm(atom1,atom2,atom_k,datom);

      DerD += dGvar* fvar* ExpTermvar;
      DerD +=  Gvar*dfvar* ExpTermvar;
      DerD +=  Gvar* fvar*dExpTermvar;
#else
      DerD += dG(atom1,atom2,atom_k,datom)* f(atom1,atom_k)* ExpTerm(atom1,atom2,atom_k);
      DerD +=  G(atom1,atom2,atom_k)*df(atom1,atom_k,datom)* ExpTerm(atom1,atom2,atom_k);
      DerD +=  G(atom1,atom2,atom_k)* f(atom1,atom_k)*dExpTerm(atom1,atom2,atom_k,datom);
#endif
    }  
  }  

  DerD += dHN(atom1,atom2,datom); 
  return DerD;
}

inline
Float
Brenner::B(Atom &atom1,Atom &atom2) 
{
  return pow(D(atom1,atom2),-delta(atom1));
}

inline
Vector3D
Brenner::dB(Atom &atom1,Atom &atom2, Atom &datom) 
{
  return -delta(atom1)*pow(D(atom1,atom2),-delta(atom1)-1.0)*dD(atom1,atom2,datom);
}

inline
Float
Brenner::G(Atom &atom_i,Atom &atom_j,Atom &atom_k) 
{
  if (atom_i.ID == C_EL)
  {
    Float CosT = CosTheta(atom_i,atom_j,atom_k);
    return ao_*(1.0+SQR(co_)/SQR(do_)-SQR(co_)/(SQR(do_)+SQR(1.0+CosT)));    
  }
  else if (atom_i.ID == H_EL)
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

inline
Vector3D
Brenner::dG(Atom &atom_i,Atom &atom_j,Atom &atom_k, Atom &datom) 
{
  if (atom_i.ID == C_EL)
  {
    Float CosT = CosTheta(atom_i,atom_j,atom_k);
    Vector3D dCosT = dCosTheta(atom_i,atom_j,atom_k,datom);

    return dCosT*dGdCT(CosT);    
  }
  else if (atom_i.ID == H_EL)
  {
    return 0.0;
  }  
  else
  {
    throw Exception("Brenner::dG() : unknown element");
  }  
}

Float
Brenner::Nt(Atom &atom_i, Atom &atom_j) 
{
  REQUIREM(atom_i.ID == C_EL,"Nt is only for C !!!");

  Float Nt_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (isHandled(atom_k))
    if (&atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      Nt_i += f(atom_i,atom_k);
    }  
  }    
  return Nt_i;

} 

Vector3D
Brenner::dNt(Atom &atom_i, Atom &atom_j, Atom &datom) 
{
  REQUIREM(atom_i.ID == C_EL,"dNt is only for C !!!");

  Vector3D Nt_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
    if (isHandled(atom_k))
    if (&atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      Nt_i += df(atom_i,atom_k,datom);
    }  
  }    
  return Nt_i;

}

Float
Brenner::NH(Atom &atom_i, Atom &atom_j) 
{ 
  REQUIREM(atom_i.ID == C_EL,"NH is only for C !!!");
  Float NH_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == H_EL && &atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      NH_i += f(atom_i,atom_k);
    }  
  }    
  return NH_i;
}

Vector3D
Brenner::dNH(Atom &atom_i, Atom &atom_j, Atom &datom) 
{ 
  REQUIREM(atom_i.ID == C_EL,"dNH is only for C !!!");
  Vector3D NH_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == H_EL && &atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      NH_i += df(atom_i,atom_k,datom);
    }  
  }    
  return NH_i;
}

Float
Brenner::NC(Atom &atom_i, Atom &atom_j) 
{
  REQUIREM(atom_i.ID == C_EL,"NC is only for C !!!");
  Float NC_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == C_EL && &atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      NC_i += f(atom_i,atom_k);
    }  
  }    
  return NC_i;
}

Vector3D
Brenner::dNC(Atom &atom_i, Atom &atom_j, Atom &datom) 
{
  REQUIREM(atom_i.ID == C_EL,"dNC is only for C !!!");
  Vector3D NC_i = 0;
  for(Index k = 0; k < NL(atom_i).size(); k++)  
  {
    Atom &atom_k = *NL(atom_i)[k];
//    if (isHandled(atom_k))
    if (atom_k.ID == C_EL && &atom_k != &atom_j)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom_i,atom_k) < R(1,atom_i,atom_k))
#endif
    {
      NC_i += df(atom_i,atom_k,datom);
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
Brenner::Nconj(Atom &atom1,Atom &atom2) 
{
  REQUIREM(atom1.ID == C_EL && atom2.ID == C_EL,"Nconj is only for C-C !!!");

  Float Nconj_ij = 1.0;

  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);
//    if (isHandled(atom_k))
    if (&atom_k != &atom2 && 
        atom_k.ID == C_EL)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom1,atom_k) < R(1,atom1,atom_k))
#endif
    {
//      Float f_ik = f(atom1,atom_k);
      Float f_ik = f(atom_k,atom1);

//      Float x_ik = Nt(atom1,atom_k);
      Float x_ik = Nt(atom_k,atom1);

      Nconj_ij += f_ik*F(x_ik);
    }  
  }    

  for(Index l = 0; l < NL(atom2).size(); l++)
  {
    Atom& atom_l = *(NL(atom2)[l]);
//    if (isHandled(atom_l))
    if (&atom_l != &atom1 &&
        atom_l.ID == C_EL)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom2,atom_l) < R(1,atom2,atom_l))
#endif
    {
//      Float f_jl = f(atom2,atom_l);
      Float f_jl = f(atom_l,atom2);

//      Float x_jl = Nt(atom2,atom_l);
      Float x_jl = Nt(atom_l,atom2);

      Nconj_ij += f_jl*F(x_jl);
    }  
  }    

  return Nconj_ij; 
}

Vector3D
Brenner::dNconj(Atom &atom1,Atom &atom2, Atom &datom) 
{
  REQUIREM(atom1.ID == C_EL && atom2.ID == C_EL,"dNconj is only for C-C !!!");

  Vector3D dNconj = 0.0;

  for(Index k = 0; k < NL(atom1).size(); k++)
  {
    Atom& atom_k = *(NL(atom1)[k]);

//    if (isHandled(atom_k))
    if (&atom_k != &atom2 && 
        atom_k.ID == C_EL)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom1,atom_k) < R(1,atom1,atom_k))
#endif
    {
//      Float f_ik = f(atom1,atom_k);
      Float f_ik = f(atom_k,atom1);

//      Float x_ik = Nt(atom1,atom_k);
      Float x_ik = Nt(atom_k,atom1);

//      Vector3D df_ik = df(atom1,atom_k,datom);
      Vector3D df_ik = df(atom_k,atom1,datom);

//      Vector3D dF_ik = dF(x_ik)*(dNt(atom1,atom_k,datom));
      Vector3D dF_ik = dF(x_ik)*(dNt(atom_k,atom1,datom));

      dNconj += df_ik*F(x_ik)+f_ik*dF_ik;
    }  
  }    

  for(Index l = 0; l < NL(atom2).size(); l++)
  {
    Atom& atom_l = *(NL(atom2)[l]);

//    if (isHandled(atom_l))
    if (&atom_l != &atom1 &&
        atom_l.ID == C_EL)
#ifdef BRENNER_OPTIMIZED_EVEN_BETTER
    if (r_vec_module_no_touch(atom2,atom_l) < R(1,atom2,atom_l))
#endif
    {
//      Float f_jl = f(atom2,atom_l);
      Float f_jl = f(atom_l,atom2);

//      Float x_jl = Nt(atom2,atom_l);
      Float x_jl = Nt(atom_l,atom2);

//      Vector3D df_jl = df(atom2,atom_l,datom);
      Vector3D df_jl = df(atom_l,atom2,datom);

//      Vector3D dF_jl = dF(x_jl)*(dNt(atom2,atom_l,datom));
      Vector3D dF_jl = dF(x_jl)*(dNt(atom_l,atom2,datom));

      dNconj += df_jl*F(x_jl)+f_jl*dF_jl;
    }  
  }    

  return dNconj;
}

inline
Float
Brenner::HN(Atom &atom_i, Atom &atom_j) 
{
RETURN_BRENNER_0;

  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
    return HN_CC_num(NH(atom_i,atom_j),NC(atom_i,atom_j));
  else if ( (atom_i.ID == C_EL && atom_j.ID == H_EL)/* || (atom_i.ID == H_EL && atom_j.ID == C_EL)*/)
    return HN_CH_num(NH(atom_i,atom_j),NC(atom_i,atom_j));  
  else return 0.0;  
}

inline
Vector3D
Brenner::dHN(Atom &atom_i, Atom &atom_j, Atom &datom) 
{
RETURN_BRENNER_0;

#ifdef BRENNER_OPTIMIZED  
  Vector3D RESULT = 0.0;
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
  {
    Vector3D dNHvar = dNH(atom_i,atom_j,datom);
    Vector3D dNCvar = dNC(atom_i,atom_j,datom);
    if (dNHvar != 0.0) RESULT += HN_CC_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNHvar;
    if (dNCvar != 0.0) RESULT += HN_CC_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNCvar;
  }  
  else if ( (atom_i.ID == C_EL && atom_j.ID == H_EL)/* || (atom_i.ID == H_EL && atom_j.ID == C_EL)*/)         
  {
    Vector3D dNHvar = dNH(atom_i,atom_j,datom);
    Vector3D dNCvar = dNC(atom_i,atom_j,datom);
    if (dNHvar != 0.0) RESULT += HN_CH_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNHvar;
    if (dNCvar != 0.0) RESULT += HN_CH_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNCvar;
  }  
  else RESULT = 0.0;
  return RESULT;
#else
  if (atom_i.ID == C_EL && atom_j.ID == C_EL)
    return HN_CC_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNH(atom_i,atom_j,datom)+
           HN_CC_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNC(atom_i,atom_j,datom);
  else if (atom_i.ID == C_EL && atom_j.ID == H_EL)         
    return HN_CH_dH_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNH(atom_i,atom_j,datom)+
           HN_CH_dC_num(NH(atom_i,atom_j),NC(atom_i,atom_j))*dNC(atom_i,atom_j,datom);
  else return 0.0;
#endif  
}

inline
Float
Brenner::FN(Atom &atom_i, Atom &atom_j) 
{
RETURN_BRENNER_0;
  
  if (atom_i.ID != C_EL || atom_j.ID != C_EL)
    return 0.0;
  else  
    return FN_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j));
}

inline
Vector3D
Brenner::dFN(Atom &atom_i, Atom &atom_j, Atom &datom) 
{
RETURN_BRENNER_0;

#ifdef BRENNER_OPTIMIZED  
  Vector3D RESULT = 0.0;
  if (atom_i.ID != C_EL || atom_j.ID != C_EL)
    RESULT = 0.0;
  else
  {  
    Vector3D dNtijvar = dNt   (atom_i,atom_j,datom);
    Vector3D dNtjivar = dNt   (atom_j,atom_i,datom);
    Vector3D dNconjvar = dNconj(atom_i,atom_j,datom);
    if (dNtijvar != 0.0) 
      RESULT += FN_dNt_i_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNtijvar;
    if (dNtjivar != 0.0) 
      RESULT += FN_dNt_j_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNtjivar;
    if (dNconjvar != 0.0)
      RESULT += FN_dNconj_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNconjvar;
  }  
  return RESULT;
#else
  if (atom_i.ID != C_EL || atom_j.ID != C_EL)
    return 0.0;
  else  
    return FN_dNt_i_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_i,atom_j,datom)+
           FN_dNt_j_num (Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNt   (atom_j,atom_i,datom)+
           FN_dNconj_num(Nt(atom_i,atom_j),Nt(atom_j,atom_i),Nconj(atom_i,atom_j))*dNconj(atom_i,atom_j,datom);
#endif
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


