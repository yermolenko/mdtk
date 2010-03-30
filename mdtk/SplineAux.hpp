/*
   Auxillary N-cubic spline classes and functions (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009 Oleksandr
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

#ifndef mdtk_SplineAux_hpp
#define mdtk_SplineAux_hpp

#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>

namespace mdtk
{

class SplineTerm
{
private:
public:
  std::vector<int> xpowers;
  Float c;
  Float multiplier;
  SplineTerm(const size_t xCount=1, const Float cval = 0.0);
  Float operator()(const std::vector<Float>& x) const;
  Float woC(const std::vector<Float>& x) const;
  void differentiate(size_t xIndex);
  virtual ~SplineTerm() {}
  void print() const;
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    int i, xpowers_count = xpowers.size();
    YAATK_FSTREAM_WRITE(os,xpowers_count,smode);
    for(i = 0; i < xpowers_count; i++)
    {
      YAATK_FSTREAM_WRITE(os,xpowers[i],smode);
    }

    YAATK_FSTREAM_WRITE(os,c,smode);
    YAATK_FSTREAM_WRITE(os,multiplier,smode);
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    int i,xpowers_count;
    YAATK_FSTREAM_READ(is,xpowers_count,smode);
    xpowers.resize(xpowers_count);
    for(i = 0; i < xpowers_count; i++)
    {
      YAATK_FSTREAM_READ(is,xpowers[i],smode);
    }

    YAATK_FSTREAM_READ(is,c,smode);
    YAATK_FSTREAM_READ(is,multiplier,smode);
  }  

};

inline
SplineTerm::SplineTerm(const size_t xCount, const Float cval)
:xpowers(),c(cval),multiplier(1.0)
{
  xpowers.resize(xCount);
}

inline
Float 
SplineTerm::operator()(const std::vector<Float>& x) const
{
//  REQUIRE(x.size() == xpowers.size());
  if (multiplier == 0.0 || c == 0.0) return 0.0;
  Float val = multiplier*c;
  size_t i;
  for(i = 0; i < x.size(); i++)
  {
//    if (xpowers[i] != 0.0)
    val *= intpow(x[i],xpowers[i]);
  }
  return val;
}

inline
Float 
SplineTerm::woC(const std::vector<Float>& x) const
{
  REQUIRE(x.size() == xpowers.size());
  if (multiplier == 0) return 0.0;
  Float val = multiplier;
  size_t i;
  for(i = 0; i < x.size(); i++)
  {
//    if (xpowers[i] != 0)
    val *= intpow(x[i],xpowers[i]);
  }
  return val;
}

inline
void
SplineTerm::print() const
{
  if (multiplier == 0) return;
  std::cout <<  "+" << multiplier << "*" << "c" << "*";
  size_t i;
  for(i = 0; i < xpowers.size(); i++)
  {
    std::cout << "x[" << i << "]^" << xpowers[i];
  }
}

inline
void
SplineTerm::differentiate(size_t xIndex)
{
  multiplier *= xpowers[xIndex];
  xpowers[xIndex] -= 1;
}

class SplineMultiD
{
private:
public:
  std::vector<SplineTerm> terms;
//  std::vector<std::vector<Float> > bounds;
  SplineMultiD(size_t dimNum = 1, int splineOrder=3);
  virtual ~SplineMultiD();
  void differentiate(size_t xIndex);
  Float operator()(const std::vector<Float>& x) const;
  void print() const;
  void computeCoeff();
  void addEquations(std::vector<std::vector<Float> >& xvs, std::vector<Float>& fvs);
  void SaveToStream(std::ostream& os, YAATK_FSTREAM_MODE smode)
  {
    int i, terms_count = terms.size();
    YAATK_FSTREAM_WRITE(os,terms_count,smode);
    for(i = 0; i < terms_count; i++)
    {
      terms[i].SaveToStream(os, smode);
    }
  }  
  void LoadFromStream(std::istream& is, YAATK_FSTREAM_MODE smode)
  {
    int i,terms_count;
    YAATK_FSTREAM_READ(is,terms_count,smode);
    terms.resize(terms_count);
    for(i = 0; i < terms_count; i++)
    {
      terms[i].LoadFromStream(is, smode);
    }

  }  
};


void
addEquation(SplineMultiD& s, std::vector<Float>& xv, Float& fv,
   std::vector<std::vector<Float> >& xvs, std::vector<Float>& fvs);

inline
SplineMultiD::~SplineMultiD()
{
}

inline
SplineMultiD::SplineMultiD(size_t dimNum, int splineOrder)
:terms()/*,bounds()*/
{
/*
  bounds.resize(dimNum);
  for(size_t i = 0; i < dimNum; i++)
  {
    bounds[i].resize(2);
  }
*/
  std::vector<int> xpowers;
  xpowers.resize(dimNum);
  int c = 0;
  while(!c)
  {
    terms.push_back(SplineTerm(dimNum));
    terms[terms.size()-1].xpowers = xpowers;

    size_t i;
    xpowers[0]++;
    for(i = 0; i < dimNum; i++)
    {
      int prevc = c;
      c = (prevc+xpowers[i])/(splineOrder+1);
      xpowers[i] = (prevc+xpowers[i])%(splineOrder+1);
    }
  }
}

inline
void
SplineMultiD::print() const
{
  TRACE(terms.size());
  size_t i;
  for(i = 0; i < terms.size(); i++)
    terms[i].print();
  std::cout << "\n";
}

inline
void
SplineMultiD::differentiate(size_t xIndex)
{
  size_t i;
  for(i = 0; i < terms.size(); i++)
  {
    terms[i].differentiate(xIndex);
  }
}

inline
Float 
SplineMultiD::operator()(const std::vector<Float>& x) const
{
  Float val = 0;
  size_t i;
  for(i = 0; i < terms.size(); i++)
  {
    val += terms[i](x);
  }
  return val;
}

} 

#endif

