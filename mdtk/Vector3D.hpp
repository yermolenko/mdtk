/*
   The Vector3D class header file.

   Copyright (C) 2004, 2005, 2009 Oleksandr Yermolenko
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

#ifndef mdtk_Vector3D_hpp
#define mdtk_Vector3D_hpp

#include <iostream>
#include <iomanip>
#include <mdtk/config.hpp>
#include <mdtk/tools.hpp>

namespace mdtk
{

class Vector3D
{
public:
  Float x,y,z;
  Float& X(size_t i)
  {
        switch (i)
        { 
          case 0: return x;break;
          case 1: return y;break;
          case 2: return z;break;
          default: throw;
        };
  }
  const Float& X(size_t i) const
  {
        switch (i)
        { 
          case 0: return x;break;
          case 1: return y;break;
          case 2: return z;break;
          default: throw;
        };
  }
  Vector3D(Float xc/* = 0.0*/, Float yc/* = 0.0*/, Float zc/* = 0.0*/);
  Vector3D();
  Vector3D& operator=(const Float&);
//  Vector3D& operator=(const Vector3D&);
  Vector3D(const Float&);

  Float    module() const;
  void     normalize();
  Vector3D normalized() const;
  friend Vector3D operator+(const Vector3D& a,const Vector3D& b);
  friend bool     operator==(const Vector3D& v1, const Vector3D& v2);
  friend bool     operator!=(const Vector3D& v1, const Vector3D& v2);
  friend void     operator+=(Vector3D& v, const Vector3D& a);
  friend void     operator-=(Vector3D& v, const Vector3D& a);
  friend Vector3D operator-(const Vector3D& v,const Vector3D& b);
  friend Vector3D operator*(const Vector3D& v,const Float& m);
  friend Vector3D operator/(const Vector3D& v,const Float& m);
  friend Vector3D operator*(const Float& m,const Vector3D& v);
  friend void     operator*=(Vector3D& v, const Float& m);
  friend void     operator/=(Vector3D& v, const Float& m);
  Vector3D operator-() const;

  friend Float    distance(const Vector3D &v1,const Vector3D &v2);
  friend Vector3D vectormul(const Vector3D& a,const Vector3D& b);
  friend Float    scalarmul(const Vector3D& a,const Vector3D& b);
  friend Float    relangle(const Vector3D& a, const Vector3D& b);

  friend std::istream&  operator>> (std::istream& is, Vector3D& vec);
  friend std::ostream&  operator<< (std::ostream& os, const Vector3D& vec);
};

std::istream&
operator>> (std::istream& is, Vector3D& vec);

std::ostream&
operator<< (std::ostream& os, const Vector3D& vec);

bool
operator==(const Vector3D& v1, const Vector3D& v2);

bool
operator!=(const Vector3D& v1, const Vector3D& v2);

Vector3D
operator+(const Vector3D& a, const Vector3D& b);

void
operator+=(Vector3D& v, const Vector3D& a);

Vector3D
operator-(const Vector3D& v, const Vector3D& b);

Vector3D operator*(const Vector3D& v, const Float& m);

Vector3D
operator/(const Vector3D& v, const Float& m);

Vector3D
operator*(const Float& m,const Vector3D& v);

void operator*=(Vector3D& v, const Float& m);

Float
distance(Vector3D &v1,Vector3D &v2);

Vector3D
vectormul(const Vector3D& a, const Vector3D& b);

Float
scalarmul(const Vector3D& a, const Vector3D& b);

Float
relangle(const Vector3D& a, const Vector3D& b);

using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::sqrt;
using std::acos;
using std::FILE;

inline
Vector3D::Vector3D(Float xc, Float yc, Float zc):
  x(xc),y(yc),z(zc)
{
}

inline
Vector3D::Vector3D():
  x(0.0),y(0.0),z(0.0)
{
}

inline
Vector3D::Vector3D(const Float& obj):
  x(obj),y(obj),z(obj)
{
}

/*
inline
Vector3D&
Vector3D::operator=(const Vector3D& obj)
{
  if (this == &obj)
    return (*this);
  x = obj.x;
  y = obj.y;
  z = obj.z;
  return (*this); 
}
*/

inline
Vector3D&
Vector3D::operator=(const Float& obj)
{
  x = obj;
  y = obj;
  z = obj;
  return (*this); 
}

inline
bool
operator==(const Vector3D& a, const Vector3D& b)
{
  return ((a.x == b.x) && (a.y == b.y) && (a.z == b.z));
}

inline
bool
operator!=(const Vector3D& a, const Vector3D& b)
{
  return ((a.x != b.x) || (a.y != b.y) || (a.z != b.z));
}

inline
Vector3D
operator+(const Vector3D& a, const Vector3D& b)
{
  return Vector3D(a.x+b.x, a.y+b.y, a.z+b.z);
}

inline
void
operator+=(Vector3D& v, const Vector3D& a)
{
  v.x+=a.x; v.y+=a.y; v.z+=a.z;
}

inline
void
operator-=(Vector3D& v, const Vector3D& a)
{
  v.x-=a.x; v.y-=a.y; v.z-=a.z;
}


inline
Vector3D
operator-(const Vector3D& v, const Vector3D& b)
{
  return Vector3D(v.x-b.x, v.y-b.y ,v.z-b.z);
}


inline
Vector3D
operator*(const Vector3D& v, const Float& m)
{
  return Vector3D(v.x*m, v.y*m, v.z*m);
}

inline
Vector3D
operator/(const Vector3D& v, const Float& m)
{
  return Vector3D(v.x/m, v.y/m, v.z/m);
}

inline
Vector3D operator*(const Float& m,const Vector3D& v)
{
  return Vector3D(v.x*m, v.y*m, v.z*m);
}

inline
void
operator*=(Vector3D& v, const Float& m)
{
  v.x*=m; v.y*=m;; v.z*=m;
}

inline
void
operator/=(Vector3D& v, const Float& m)
{
  v.x/=m; v.y/=m;; v.z/=m;
}

inline
Float
distance(const Vector3D &v1, const Vector3D &v2)
{
  return sqrt(SQR(v1.x-v2.x)+SQR(v1.y-v2.y)+SQR(v1.z-v2.z));
}

inline
Vector3D
vectormul(const Vector3D& a, const Vector3D& b)
{
  return Vector3D(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

inline
Float
scalarmul(const Vector3D& a, const Vector3D& b)
{
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

inline
Float
relangle(const Vector3D& a, const Vector3D& b)
{
  return acos(scalarmul(a,b)/(a.module()*b.module()));
}

inline
Float
Vector3D::module() const
{
  return sqrt(x*x+y*y+z*z);
}

inline
void
Vector3D::normalize()
{
  Float L=module();
  x/=L; y/=L; z/=L;
}

inline
Vector3D
Vector3D::normalized() const
{
  Float L=module();
  return Vector3D(x/L, y/L, z/L);
}

inline
Vector3D
Vector3D::operator-() const
{
  return Vector3D(-x,-y,-z);
}

inline
istream&
operator>>(istream& is,  Vector3D& vec)
{
  Float u = 0, v = 0, w = 0;
  is >> u >> v >> w;
  if (is)
    vec = Vector3D(u,v,w);
  else
    cerr << " Error in reading Vector3D " << endl << flush;
  return is;
}

inline
ostream&
operator<<(ostream& os, const Vector3D& vec)
{
  os << vec.x << " "
     << vec.y << " "
     << vec.z << " ";
  return os;
}

} // namespace mdtk

#endif


