/*
   Implementation of the AIREBO interatomic potential for
   hydrocarbons. Splines.
   See [S.J. Stuart, A.B. Tutein and J.A. Harrison,
   J. Chem. Phys. 112, 6472 (2000)]

   Copyright (C) 2005, 2006, 2007, 2008, 2009 Oleksandr Yermolenko
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

#include "cubici.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

namespace mdtk
{

using namespace std;

#define SET_FUNC_G(i,x_val,y_val,dy_val,d2y_val) \
{ \
  x[i] = x_val; y[i] = y_val; dy[i] = dy_val; d2y[i] = d2y_val; \
}  

void
FuncG_C1::init()
{
#define DATA_POINTS_G_C1  5
  splineCount = DATA_POINTS_G_C1-1;
  Float x[DATA_POINTS_G_C1];
  Float y[DATA_POINTS_G_C1];
  Float dy[DATA_POINTS_G_C1];
  Float d2y[DATA_POINTS_G_C1];
  
  SET_FUNC_G(0,    -1.0,  -0.010000,  0.104000,  0.000000);
  SET_FUNC_G(1,-2.0/3.0,   0.028207,  0.131443,  0.140229);
  SET_FUNC_G(2,-1.0/2.0,   0.052804,  0.170000,  0.370000);
  SET_FUNC_G(3,-1.0/3.0,   0.097321,  0.400000,  1.980000);
  SET_FUNC_G(4,     1.0,   1.000000,  2.834570, 10.264700);

  for(int i = 0; i < splineCount; i++)
    spline[i] = Spline5n(x[i  ],y[i  ],dy[i  ],d2y[i  ],
                         x[i+1],y[i+1],dy[i+1],d2y[i+1]);
}  

void
FuncG_C2::init()
{
#define DATA_POINTS_G_C2  5
  splineCount = DATA_POINTS_G_C2-1;
  Float x[DATA_POINTS_G_C2];
  Float y[DATA_POINTS_G_C2];
  Float dy[DATA_POINTS_G_C2];
  Float d2y[DATA_POINTS_G_C2];
  
  SET_FUNC_G(0,    -1.0,  -0.010000,  0.104000,  0.000000);
  SET_FUNC_G(1,-2.0/3.0,   0.028207,  0.131443,  0.140229);
  SET_FUNC_G(2,-1.0/2.0,   0.052804,  0.170000,  0.370000);
  SET_FUNC_G(3,-1.0/3.0,   0.097321,  0.400000,  1.980000);
  SET_FUNC_G(4,     1.0,   8.000000, 20.243600,43.9336000);
  
  for(int i = 0; i < splineCount; i++)
    spline[i] = Spline5n(x[i  ],y[i  ],dy[i  ],d2y[i  ],
                         x[i+1],y[i+1],dy[i+1],d2y[i+1]);
}  

void
FuncG_H::init()
{
#define DATA_POINTS_G_H  4
  splineCount = DATA_POINTS_G_H-1;
  Float x[DATA_POINTS_G_H];
  Float y[DATA_POINTS_G_H];
  Float dy[DATA_POINTS_G_H];
  Float d2y[DATA_POINTS_G_H];
  
  SET_FUNC_G(0,    -1.0,  11.235700,  0.000000,115.115000);
  SET_FUNC_G(1,-5.0/6.0,  12.595300, 13.854300, 32.361800);
  SET_FUNC_G(2,-1.0/2.0,  16.811100,  8.641230,-25.061700);
  SET_FUNC_G(3,     1.0,  19.991800,  0.333013, -0.474189);
  
  for(int i = 0; i < splineCount; i++)
    spline[i] = Spline5n(x[i  ],y[i  ],dy[i  ],d2y[i  ],
                         x[i+1],y[i+1],dy[i+1],d2y[i+1]);
}  

FuncP_CC::FuncP_CC(int paramSet)
{
  Float HCC[5][5];

  int i,j;

  for(i = 0; i < 5; i++)
    for(j = 0; j < 5; j++)
      HCC[i][j] = 0.0;

  switch (paramSet)
  {
    case 0:
  HCC[2][0] = -0.000500;
  HCC[3][0] =  0.016125;
  HCC[1][1] = -0.010960;
  HCC[2][1] =  0.006326;
  HCC[0][2] = -0.027603;
  HCC[1][2] =  0.003180;
    break;
/*
    case 1:
  HCC[1][1] = -0.0226;
  HCC[2][0] = -0.0061;
  HCC[3][0] =  0.0173;
  HCC[1][2] =  0.0149;
  HCC[2][1] =  0.0160;
    break;
*/
    default : throw Exception("Unknown ParamSet for FuncP_CC");
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

FuncP_CH::FuncP_CH(int paramSet)
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
  HCH[1][0] =  0.209337;
  HCH[2][0] = -0.064450;
  HCH[3][0] = -0.303928;
  HCH[0][1] =  0.010000;
  HCH[1][1] = -0.125123;
  HCH[2][1] = -0.298905;
  HCH[0][2] = -0.122042;
  HCH[1][2] = -0.300529;
  HCH[0][3] = -0.307585;

/*
  DHCHDH[2][0] = -0.13075; // -0.130745;
  DHCHDH[1][1] = -0.0764;

  DHCHDC[0][2] = -0.07655;
  DHCHDC[1][1] = -0.12805;
*/
    break;
/*
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
*/
    default : throw Exception("Unknown ParamSet for FuncP_CH");
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

#define FCCfillUp(iBase, jBase, kBase, i_size, j_size, k_size) \
{ \
  for(int i = 0; i < i_size; i++) \
    for(int j = 0; j < j_size; j++) \
      if ((i==iBase && j==jBase) || (i==jBase && j==iBase)) \
        for(int k = kBase+1; k < k_size; k++)  FCC[i][j][k] = FCC[i][j][kBase]; \
}  

#define FCCfillSym(i, j, k, value) \
  FCC[i][j][k] = FCC[j][i][k] =  value;

#define DFCCDIfillSym(i, j, k, value) \
  DFCCDI[i][j][k] = DFCCDI[j][i][k] =  value;

#define DFCCDJfillSym(i, j, k, value) \
  DFCCDJ[i][j][k] = DFCCDJ[j][i][k] =  value;


#define DFCCDKfillSym(i, j, k, value) \
  DFCCDK[i][j][k] = DFCCDK[j][i][k] =  value;

#define pi_rc_c0_COMMON_INC \
\
  int i,j,k;\
\
  Float FCC[i_size_pi_rc][j_size_pi_rc][k_size_pi_rc];\
  Float DFCCDI[i_size_pi_rc][j_size_pi_rc][k_size_pi_rc];\
  Float DFCCDJ[i_size_pi_rc][j_size_pi_rc][k_size_pi_rc];\
  Float DFCCDK[i_size_pi_rc][j_size_pi_rc][k_size_pi_rc];\
\
  for(i = 0; i < i_size_pi_rc; i++)\
    for(j = 0; j < j_size_pi_rc; j++)\
      for(k = 0; k < k_size_pi_rc; k++)\
      {\
        FCC[i][j][k] = 0.0;\
        DFCCDI[i][j][k] = 0.0;\
        DFCCDJ[i][j][k] = 0.0;\
        DFCCDK[i][j][k] = 0.0;\
      }  \\


#define pi_rc_c1_COMMON_INC \
\
  Float          x[2];\
  Float          y[2];\
  Float          z[2];\
  Float    v[2][2][2];\
  Float dvdx[2][2][2];\
  Float dvdy[2][2][2];\
  Float dvdz[2][2][2];\
\
  for(i = 0; i < i_size_pi_rc-1; i++)\
   for(j = 0; j < j_size_pi_rc-1;j++)\
    for(k = 0; k < k_size_pi_rc-1;k++)\
    {\
      x[0] = i;\
      y[0] = j;\
      z[0] = k;\
      x[1] = i+1;\
      y[1] = j+1;\
      z[1] = k+1;\
\
      v[0][0][0] = FCC[i]  [j]  [k];\
      v[0][0][1] = FCC[i]  [j]  [k+1];\
      v[0][1][0] = FCC[i]  [j+1][k];\
      v[0][1][1] = FCC[i]  [j+1][k+1];\
      v[1][0][0] = FCC[i+1][j]  [k];\
      v[1][0][1] = FCC[i+1][j]  [k+1];\
      v[1][1][0] = FCC[i+1][j+1][k];\
      v[1][1][1] = FCC[i+1][j+1][k+1];\
\
      dvdx[0][0][0] = DFCCDI[i]  [j]  [k];\
      dvdx[0][0][1] = DFCCDI[i]  [j]  [k+1];\
      dvdx[0][1][0] = DFCCDI[i]  [j+1][k];\
      dvdx[0][1][1] = DFCCDI[i]  [j+1][k+1];\
      dvdx[1][0][0] = DFCCDI[i+1][j]  [k];\
      dvdx[1][0][1] = DFCCDI[i+1][j]  [k+1];\
      dvdx[1][1][0] = DFCCDI[i+1][j+1][k];\
      dvdx[1][1][1] = DFCCDI[i+1][j+1][k+1];\
\
      dvdy[0][0][0] = DFCCDJ[i]  [j]  [k];\
      dvdy[0][0][1] = DFCCDJ[i]  [j]  [k+1];\
      dvdy[0][1][0] = DFCCDJ[i]  [j+1][k];\
      dvdy[0][1][1] = DFCCDJ[i]  [j+1][k+1];\
      dvdy[1][0][0] = DFCCDJ[i+1][j]  [k];\
      dvdy[1][0][1] = DFCCDJ[i+1][j]  [k+1];\
      dvdy[1][1][0] = DFCCDJ[i+1][j+1][k];\
      dvdy[1][1][1] = DFCCDJ[i+1][j+1][k+1];\
\
      dvdz[0][0][0] = DFCCDK[i]  [j]  [k];\
      dvdz[0][0][1] = DFCCDK[i]  [j]  [k+1];\
      dvdz[0][1][0] = DFCCDK[i]  [j+1][k];\
      dvdz[0][1][1] = DFCCDK[i]  [j+1][k+1];\
      dvdz[1][0][0] = DFCCDK[i+1][j]  [k];\
      dvdz[1][0][1] = DFCCDK[i+1][j]  [k+1];\
      dvdz[1][1][0] = DFCCDK[i+1][j+1][k];\
      dvdz[1][1][1] = DFCCDK[i+1][j+1][k+1];\
\
      spline[i][j][k] = Spline3D(x,y,z,v,dvdx,dvdy,dvdz);\
    }  \\



void  
Func_pi_rc_CC::init(int paramSet)
{
pi_rc_c0_COMMON_INC;

  switch (paramSet)
  {
    case 0:
    break;
    default : throw Exception("Unknown ParamSet for FuncC");
  };
  FCCfillSym(0,0,3, 0.004959);
   FCCfillUp(0,0,3, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(1,0,1, 0.021694);
  FCCfillSym(1,0,2, 0.004959);
   FCCfillUp(1,0,2, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(1,1,1, 0.052500);
  FCCfillSym(1,1,2,-0.002089);
  FCCfillSym(1,1,3,-0.008043);
   FCCfillUp(1,1,3, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(2,0,1, 0.024699);
  FCCfillSym(2,0,2,-0.005971);
  FCCfillSym(2,0,3, 0.004959);
   FCCfillUp(2,0,3, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(2,1,1, 0.004825);  DFCCDIfillSym(2,1,1,-0.026250);
  FCCfillSym(2,1,2, 0.015000);
  FCCfillSym(2,1,3,-0.010000);
  FCCfillSym(2,1,4,-0.011689);  DFCCDKfillSym(2,1,4,-0.010022);
  FCCfillSym(2,1,5,-0.013378);  DFCCDIfillSym(2,1,5,-0.027188);  DFCCDKfillSym(2,1,5,-0.010022);
  FCCfillSym(2,1,6,-0.015067);  DFCCDIfillSym(2,1,6,-0.027188);
  FCCfillSym(2,1,7,-0.015067);  DFCCDIfillSym(2,1,7,-0.027188);
   FCCfillUp(2,1,7, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(2,2,1, 0.047225);
  FCCfillSym(2,2,2, 0.011000);
  FCCfillSym(2,2,3, 0.019853);
  FCCfillSym(2,2,4, 0.016544);  DFCCDKfillSym(2,2,4,-0.003309);
  FCCfillSym(2,2,5, 0.013235);  DFCCDKfillSym(2,2,5,-0.003309);
  FCCfillSym(2,2,6, 0.009926);  DFCCDKfillSym(2,2,6,-0.003309);
  FCCfillSym(2,2,7, 0.006618);  DFCCDKfillSym(2,2,7,-0.003309);
  FCCfillSym(2,2,8, 0.003309);  DFCCDKfillSym(2,2,8,-0.003309);
  FCCfillSym(3,0,1,-0.099899);
  FCCfillSym(3,0,2,-0.099899);
  FCCfillSym(3,0,3, 0.004959);
   FCCfillUp(3,0,3, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(3,1,2,-0.062418);  DFCCDJfillSym(3,1,2, 0.037545);
  FCCfillSym(3,1,3,-0.062418);
   FCCfillUp(3,1,3, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(3,2,1,-0.022355);
  FCCfillSym(3,2,2,-0.022355);  DFCCDJfillSym(3,2,2, 0.062418);
   FCCfillUp(3,2,2, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);

pi_rc_c1_COMMON_INC;

}    
  

void  
Func_pi_rc_CH::init(int paramSet)
{
pi_rc_c0_COMMON_INC;

  switch (paramSet)
  {
    case 0:
    break;
    default : throw Exception("Unknown ParamSet for FuncC");
  };

  FCCfillSym(1,1,1,-0.050000);
  FCCfillSym(1,1,2,-0.050000);
  FCCfillSym(1,1,3,-0.300000);
  FCCfillSym(1,1,4,-0.050000);
   FCCfillUp(1,1,4, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(2,0,5,-0.004524);
   FCCfillUp(2,0,5, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);
  FCCfillSym(2,1,2,-0.250000);
  FCCfillSym(2,1,3,-0.250000);
  FCCfillSym(3,1,1,-0.100000);
  FCCfillSym(3,1,2,-0.125000);
  FCCfillSym(3,1,3,-0.125000);
  FCCfillSym(3,1,4,-0.100000);
   FCCfillUp(3,1,4, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);

pi_rc_c1_COMMON_INC;

}    
  
void  
Func_pi_rc_HH::init(int paramSet)
{
pi_rc_c0_COMMON_INC;

  switch (paramSet)
  {
    case 0:
    break;
    default : throw Exception("Unknown ParamSet for FuncC");
  };

  FCCfillSym(1,1,1, 0.124916);

pi_rc_c1_COMMON_INC;

}    
  

void  
Func_Tij::init(int paramSet)
{
pi_rc_c0_COMMON_INC;

  switch (paramSet)
  {
    case 0:
    break;
    default : throw Exception("Unknown ParamSet for FuncTij");
  };
  FCCfillSym(2,2,1,-0.035140);
  FCCfillSym(2,2,2,-0.004048);
   FCCfillUp(2,2,2, i_size_pi_rc, j_size_pi_rc, k_size_pi_rc);

pi_rc_c1_COMMON_INC

}    

}



