/*
   Implementation of the many-body interatomic potential for
   hydrocarbons. Splines.
   See [D.W. Brenner, Phys. Rev. B 42, 9458 (1990)]

   Copyright (C) 2004, 2005, 2006, 2007, 2009, 2010 Oleksandr
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

#include "cubici.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace mdtk
{

using namespace std;

FuncH_CC::FuncH_CC(int paramSet)
{
  Float HCC[5][5];

  int i,j;
  for(i = 0; i < 5; i++)
    for(j = 0; j < 5; j++)
      HCC[i][j] = 0.0;

  switch (paramSet)
  {
    case 0:
  HCC[1][1] = -0.0175;
  HCC[2][0] = -0.0070;
  HCC[3][0] =  0.0119;
  HCC[1][2] =  0.0115;
  HCC[2][1] =  0.0118;
    break;
    case 1:
  HCC[1][1] = -0.0226;
  HCC[2][0] = -0.0061;
  HCC[3][0] =  0.0173;
  HCC[1][2] =  0.0149;
  HCC[2][1] =  0.0160;
    break;
    default : throw Exception("Unknown ParamSet for FuncH_CC");
  };


  Float       x[2];
  Float       y[2];
  Float    z[2][2];
  Float dzdx[2][2];
  Float dzdy[2][2];

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4;j++)
    {
      x[0] = i;
      y[0] = j;
      x[1] = i+1;
      y[1] = j+1;

      z[0][0] = HCC[i][j];
      z[0][1] = HCC[i][j+1];
      z[1][0] = HCC[i+1][j];
      z[1][1] = HCC[i+1][j+1];

      dzdx[0][0] = 0;
      dzdx[0][1] = 0;
      dzdx[1][0] = 0;
      dzdx[1][1] = 0;

      dzdy[0][0] = 0;
      dzdy[0][1] = 0;
      dzdy[1][0] = 0;
      dzdy[1][1] = 0;

      spline[i][j] = Spline2D(x,y,z,dzdx,dzdy);
    }  
}  

FuncH_CH::FuncH_CH(int paramSet)
{
  Float HCH[5][5];
  Float DHCHDH[5][5];
  Float DHCHDC[5][5];

  int i,j;
  for(i = 0; i < 5; i++)
    for(j = 0; j < 5; j++)
    {
      HCH[i][j] = 0.0;
      DHCHDH[i][j] = 0.0;
      DHCHDC[i][j] = 0.0;
    }  

  switch (paramSet)
  {
    case 0:
  HCH[1][0] = -0.0760;
  HCH[2][0] = -0.2163;
  HCH[3][0] = -0.3375;
  HCH[0][1] = -0.1792;
  HCH[0][2] = -0.2407;
  HCH[1][1] = -0.2477;
  HCH[2][1] = -0.3320;
  HCH[0][3] = -0.3323;
  HCH[1][2] = -0.3321;

  DHCHDH[2][0] = -0.13075; // -0.130745;
  DHCHDH[1][1] = -0.0764;

  DHCHDC[0][2] = -0.07655;
  DHCHDC[1][1] = -0.12805;
    break;
    case 1:
  HCH[1][0] = -0.0984;
  HCH[2][0] = -0.2878;
  HCH[3][0] = -0.4507;
  HCH[0][1] = -0.2479;
  HCH[0][2] = -0.3221;
  HCH[1][1] = -0.3344;
  HCH[2][1] = -0.4438;
  HCH[0][3] = -0.4460;
  HCH[1][2] = -0.4449;

  DHCHDH[2][0] = -0.17615;
  DHCHDH[1][1] = -0.09795;

  DHCHDC[0][2] = -0.09905;
  DHCHDC[1][1] = -0.17325;
    break;
    default : throw Exception("Unknown ParamSet for FuncH_CH");
  };

  Float       x[2];
  Float       y[2];
  Float    z[2][2];
  Float dzdx[2][2];
  Float dzdy[2][2];

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4;j++)
    {
      x[0] = i;
      y[0] = j;
      x[1] = i+1;
      y[1] = j+1;

      z[0][0] = HCH[i][j];
      z[0][1] = HCH[i][j+1];
      z[1][0] = HCH[i+1][j];
      z[1][1] = HCH[i+1][j+1];

      dzdx[0][0] = DHCHDH[i][j];
      dzdx[0][1] = DHCHDH[i][j+1];
      dzdx[1][0] = DHCHDH[i+1][j];
      dzdx[1][1] = DHCHDH[i+1][j+1];

      dzdy[0][0] = DHCHDC[i][j];
      dzdy[0][1] = DHCHDC[i][j+1];
      dzdy[1][0] = DHCHDC[i+1][j];
      dzdy[1][1] = DHCHDC[i+1][j+1];

      spline[i][j] = Spline2D(x,y,z,dzdx,dzdy);
    }  
}  

  
FuncF::FuncF(int paramSet)
{
  Float FCC[5][5][6];
  Float DFCCDI[5][5][6];
  Float DFCCDJ[5][5][6];

  int i,j,k;

  for(i = 0; i < 5; i++)
    for(j = 0; j < 5; j++)
      for(k = 0; k < 6; k++)
      {
        FCC[i][j][k] = 0.0;
        DFCCDI[i][j][k] = 0.0;
        DFCCDJ[i][j][k] = 0.0;
      }  

  switch (paramSet)
  {
    case 0:
  FCC[2][3][1] = FCC[3][2][1] = -0.0465;
  FCC[2][3][2] = FCC[3][2][2] = -0.0465; // -0.01860
  FCC[1][2][2] = FCC[2][1][2] = -0.0355; 
  FCC[1][2][1] = FCC[2][1][1] =  0.0126;
  FCC[1][3][1] = FCC[3][1][1] = -0.1130;
  FCC[1][3][2] = FCC[3][1][2] = -0.1130;
  FCC[0][3][1] = FCC[3][0][1] = -0.1220;
  FCC[0][3][2] = FCC[3][0][2] = -0.1220;
  FCC[0][2][2] = FCC[2][0][2] = -0.0445;
  FCC[0][2][1] = FCC[2][0][1] =  0.0320;
  FCC[0][1][1] = FCC[1][0][1] =  0.1100;

  FCC[1][1][2] = FCC[1][1][2] =  0.0074;
  FCC[1][1][1] = FCC[1][1][1] =  0.1511;
  FCC[2][2][1] = FCC[2][2][1] =  0.0750;

  DFCCDI[2][0][1] = -0.1160;
  DFCCDI[2][1][1] = -0.13205;
  DFCCDI[2][0][2] = -0.0610;
  DFCCDI[1][2][2] =  0.02225;
//  DFCCDI[1][3][2] =  -0.03775; // 0.05170 GK;
  DFCCDI[1][3][2] =  0.03775; // 0.05170 GK;
  DFCCDI[2][3][2] =  0.0565;
  DFCCDI[2][3][1] =  0.0565;
//  DFCCDI[2][1][2] = -0.1065;
  DFCCDI[2][1][2] = -0.0602;

  DFCCDJ[0][2][1] = -0.1160;
  DFCCDJ[1][2][1] = -0.13205;
  DFCCDJ[0][2][2] = -0.0610;
  DFCCDJ[2][1][2] =  0.02225;
//  DFCCDJ[3][1][2] =  -0.03775; // 0.05170 GK;
  DFCCDJ[3][1][2] =  0.03775; // 0.05170 GK;
  DFCCDJ[3][2][2] =  0.0565;
  DFCCDJ[3][2][1] =  0.0565;
//  DFCCDJ[1][2][2] = -0.1065;
  DFCCDJ[1][2][2] = -0.0602;
    break;
    case 1:
  FCC[2][3][1] = FCC[3][2][1] = -0.0363;
  FCC[2][3][2] = FCC[3][2][2] = -0.0363;
  FCC[1][2][2] = FCC[2][1][2] = -0.0243; 
  FCC[1][2][1] = FCC[2][1][1] =  0.0120;
  FCC[1][3][1] = FCC[3][1][1] = -0.0903;
  FCC[1][3][2] = FCC[3][1][2] = -0.0903;
  FCC[0][3][1] = FCC[3][0][1] = -0.0904;
  FCC[0][3][2] = FCC[3][0][2] = -0.0904;
  FCC[0][2][2] = FCC[2][0][2] = -0.0269;
  FCC[0][2][1] = FCC[2][0][1] =  0.0427;
  FCC[0][1][1] = FCC[1][0][1] =  0.0996;

  FCC[1][1][2] = FCC[1][1][2] =  0.0108;
  FCC[1][1][1] = FCC[1][1][1] =  0.1264;
  FCC[2][2][1] = FCC[2][2][1] =  0.0605;

  DFCCDI[2][0][1] = -0.0950;
  DFCCDI[2][1][1] = -0.10835;
  DFCCDI[2][0][2] = -0.0452;
  DFCCDI[1][2][2] =  0.01345;
  DFCCDI[1][3][2] =  0.02705;
  DFCCDI[2][3][2] =  0.04515;
  DFCCDI[2][3][1] =  0.04515;
  DFCCDI[2][1][2] = -0.05055;

  DFCCDJ[0][2][1] = -0.0950;
  DFCCDJ[1][2][1] = -0.10835;
  DFCCDJ[0][2][2] = -0.0452;
  DFCCDJ[2][1][2] =  0.01345;
  DFCCDJ[3][1][2] =  0.02705;
  DFCCDJ[3][2][2] =  0.04515;
  DFCCDJ[3][2][1] =  0.04515;
  DFCCDJ[1][2][2] = -0.05055;
    break;
    default : throw Exception("Unknown ParamSet for FuncC");
  };

  for(i = 0; i < 5; i++)
    for(j = 0; j < 5; j++)
      for(k = 3; k < 6; k++)
      {
        FCC[i][j][k] = FCC[i][j][2];
      }  

  Float          x[2];
  Float          y[2];
  Float          z[2];
  Float    v[2][2][2];
  Float dvdx[2][2][2];
  Float dvdy[2][2][2];
  Float dvdz[2][2][2];

  for(i = 0; i < 4; i++)
   for(j = 0; j < 4;j++)
    for(k = 0; k < 5;k++)
    {
      x[0] = i;
      y[0] = j;
      z[0] = k;
      x[1] = i+1;
      y[1] = j+1;
      z[1] = k+1;

      v[0][0][0] = FCC[i]  [j]  [k];
      v[0][0][1] = FCC[i]  [j]  [k+1];
      v[0][1][0] = FCC[i]  [j+1][k];
      v[0][1][1] = FCC[i]  [j+1][k+1];
      v[1][0][0] = FCC[i+1][j]  [k];
      v[1][0][1] = FCC[i+1][j]  [k+1];
      v[1][1][0] = FCC[i+1][j+1][k];
      v[1][1][1] = FCC[i+1][j+1][k+1];

      dvdx[0][0][0] = DFCCDI[i]  [j]  [k];
      dvdx[0][0][1] = DFCCDI[i]  [j]  [k+1];
      dvdx[0][1][0] = DFCCDI[i]  [j+1][k];
      dvdx[0][1][1] = DFCCDI[i]  [j+1][k+1];
      dvdx[1][0][0] = DFCCDI[i+1][j]  [k];
      dvdx[1][0][1] = DFCCDI[i+1][j]  [k+1];
      dvdx[1][1][0] = DFCCDI[i+1][j+1][k];
      dvdx[1][1][1] = DFCCDI[i+1][j+1][k+1];

      dvdy[0][0][0] = DFCCDJ[i]  [j]  [k];
      dvdy[0][0][1] = DFCCDJ[i]  [j]  [k+1];
      dvdy[0][1][0] = DFCCDJ[i]  [j+1][k];
      dvdy[0][1][1] = DFCCDJ[i]  [j+1][k+1];
      dvdy[1][0][0] = DFCCDJ[i+1][j]  [k];
      dvdy[1][0][1] = DFCCDJ[i+1][j]  [k+1];
      dvdy[1][1][0] = DFCCDJ[i+1][j+1][k];
      dvdy[1][1][1] = DFCCDJ[i+1][j+1][k+1];

      dvdz[0][0][0] = 0.0;
      dvdz[0][0][1] = 0.0;
      dvdz[0][1][0] = 0.0;
      dvdz[0][1][1] = 0.0;
      dvdz[1][0][0] = 0.0;
      dvdz[1][0][1] = 0.0;
      dvdz[1][1][0] = 0.0;
      dvdz[1][1][1] = 0.0;

      spline[i][j][k] = Spline3D(x,y,z,v,dvdx,dvdy,dvdz);
    }  
}    

}



