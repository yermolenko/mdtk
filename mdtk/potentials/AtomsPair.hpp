/*
   The AtomsPair class (header file).

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2011, 2012
   Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#ifndef mdtk_AtomsPair_hpp
#define mdtk_AtomsPair_hpp

#include <mdtk/config.hpp>
#include <mdtk/consts.hpp>
#include <mdtk/Atom.hpp>

#include <cmath>

namespace mdtk
{

struct AtomsPair
{
  Atom& atom1;
  Atom& atom2;
  Vector3D rv;
  Float r_;
  Float r_squared_;
  Vector3D dr_vec_module_template;
public:
  Float r(const Float V = 0.0)
    {
      if (V != 0.0)
      {
        atom1.grad += dr(atom1)*V;
        atom2.grad += dr(atom2)*V;
      }
      return r_;
    }
  Float r_alter(const Float V = 0.0)
    {
      if (!alter_depos_enabled)
        return r(V);
      if (V != 0.0)
      {
        atom1.grad += dr_alter(atom1)*V;
        atom2.grad += dr_alter(atom2)*V;
      }
      return r_alter_;
    }
  Vector3D dr(const Atom& datom) const
    {
      Vector3D v;
      if (&atom1 == &datom)
        v = dr_vec_module_template;
      else if (&atom2 == &datom)
        v = -dr_vec_module_template;
      return v;
    }
  Vector3D dr_alter(const Atom& datom) const
    {
      if (!alter_depos_enabled)
        return dr(datom);
      Vector3D temp = rv.normalized()*r_alter_;
      Vector3D v;
      if (&atom1 == &datom)
        v = temp/temp.module();
      else if (&atom2 == &datom)
        v = -temp/temp.module();
      return v;
    }
  Float r_squared(const Float V = 0.0)
    {
      if (V != 0.0)
        r(2.0*V*r_);
      return r_squared_;
    }
private:
  Float f_;
public:
  Float f(const Float V = 0.0)
    {
      if (df_template != 0.0)
        r(V*df_template);
      return f_;
    }
private:
  Float df_template;
public:
  Vector3D df(const Atom& datom) const
    {
      return df_template*dr(datom);
    }
  Float r_alter_;
  bool alter_depos_enabled;
  AtomsPair operator-() const
    {
      AtomsPair p(atom2,atom1,false);
      p.rv = -rv;
      p.r_ = r_;
      p.r_squared_ = r_squared_;
      p.dr_vec_module_template = -dr_vec_module_template;
      p.f_ = f_;
      p.df_template = df_template;
      p.r_alter_ = r_alter_;
      p.alter_depos_enabled = alter_depos_enabled;
      return p;
    }
  AtomsPair(Atom& ai, Atom& aj, const bool init)
    :atom1(ai),atom2(aj)
    {
      if (init)
      {
        rv = depos(ai,aj);
        r_ = rv.module();
        r_squared_ = (/*rv.module_squared()*/r_*r_);
        dr_vec_module_template = (rv/r_);
        f_ = 1.0;
        df_template = 0.0;
        r_alter_ = 0.0;
        alter_depos_enabled = false;
      }
    }
  AtomsPair(Atom& ai, Atom& aj, const Float R1, const Float R2, const Float alter_depos = 0.0)
    :atom1(ai),atom2(aj),
     rv(depos(ai,aj)),
     r_(rv.module()),
     r_squared_(/*rv.module_squared()*/r_*r_),
     dr_vec_module_template(rv/r_),
     f_(),
     df_template(),
     r_alter_(alter_depos),
     alter_depos_enabled(alter_depos != 0.0)
    {
      if (r_<R1)
        f_ = 1.0;
      else
      {
        if (r_<R2)
        {
          f_ = (1.0+cos(M_PI*(r_-R1)/(R2-R1)))/2.0;
          df_template = (-(M_PI/(R2-R1))*sin(M_PI*(r_-R1)/(R2-R1)))/2.0;
        }
      }
    }
};

} // namespace mdtk

#endif
