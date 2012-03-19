/*
   The VisBox class for the molecular dynamics trajectory viewer

   Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
   2012 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include <FL/Fl.H>

#include "VisBox.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include <algorithm>

#include "../common.h"
namespace xmde
{

using namespace mdtk;

void
VisBox::loadNewSnapshot(std::string base_state_filename,std::string file)
{
  using mdtk::Exception;

  TRACE(file);

  std::vector<std::string>::iterator xvai = 
    std::find(xvaList.begin(),xvaList.end(),file);

  REQUIRE(xvai != xvaList.end());

  if (xvai == xvaList.begin())
  {
    
    TRACE("*********LOADING INITIAL STATE *******");
    
    TRACE(base_state_filename);

    if (base_state_filename.find("simloop.conf") != std::string::npos) 
    {
      ml_->loadstate();
    }
    else
    {
      yaatk::text_ifstream fi(base_state_filename.c_str()); 

      ml_->initNLafterLoading = false;

      if (base_state_filename.find("mde_init") != std::string::npos)
	ml_->loadFromStream(fi);
      else
      {
	ml_->loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
	ml_->allowPartialLoading = true; // hack, disables essential checks
	ml_->atoms.prepareForSimulatation();
      }
      fi.close(); 
    }
    if (file != "shot")
    {
      yaatk::text_ifstream fixva(file.c_str()); 
      ml_->loadFromStreamXVA(fixva);
      fixva.close(); 
      /*
	yaatk::binary_ifstream fixva(file.c_str()); 
	ml_->loadFromStreamXVA_bin(fixva);
	fixva.close(); 
      */
    }
  }
  else
  {
    TRACE("********* UPDATING FROM MDT ***********");
    MDTrajectory::const_iterator t = mdt.begin();
    int xvaCount = 0;
    TRACE(xvai-xvaList.begin());
    while (xvaCount < xvai-xvaList.begin())
    {
      ++t;
      ++xvaCount;
    }
    TRACE(xvaCount);
    const std::vector<Atom>& atoms = t->second;
    TRACE(atoms.size());
    TRACE(ml_->atoms.size());
    REQUIRE(atoms.size() == ml_->atoms.size());
    for(size_t i = 0; i < ml_->atoms.size(); ++i)
    {
      ml_->atoms[i] = atoms[i];
    }
    ml_->simTime = t->first;
  }

  Float sc = scale;
  Float msc = maxScale;

  setData(*ml_);

  if (maxScale < msc) maxScale = msc;
  scale = sc;
  redraw();
}  

VisBox::VisBox(int x,int y,int w,int h,std::string base_state_filename,
	       const std::vector<std::string>& xvas)
  : Fl_Gl_Window(x,y,w,h,"MDTK Trajectory Viewer - 3D View"),
    allowRescale(true),
    vertexColor(combineRGB(255,255,255)),
    edgeColor(combineRGB(255,255,255)),
    bgColor(combineRGB(255,255,255)),
    showAxes(false),
    showCTree(false),
    showCTreeConnected(true),
    showCTreeAtoms(true),
    showCTreeAllTimes(false),
    tinyAtoms(false),
    downscaleCTree(3.0),
    energyThresholdCTree(5.0),
    showAtoms(true),
    showBath(false),
    showBathSketch(false),
    unfoldPBC(false),
    showCustom1(false),
    showCustom2(false),
    showCustom3(false),
    showSelected(false),
    showBarrier(false),
    nativeVertexColors(true),
    atomsQuality(14),
    atomsQualityInHQMode(50),
    hqMode(false),
    selectedAtomIndex(0),
    nRange(50),
    vertexRadius(1.0), axesRadius(1.0), scale(1.0), maxScale(1.0),
    R(),Ro(),
    ml_(NULL),
    mdt(),
    xvaList(xvas),
    ctree(NULL),
    zbar(0.0),
    lstBall(0),
    lstBallHQ(0),
    lstStick(0),
    lstStickHQ(0),
    old_rot_x(0.0), old_rot_y(0.0), old_rot_z(0.0),
    tiledMode(false),tileCount(8)
{
  mode(FL_RGB | FL_DOUBLE | FL_ACCUM | FL_ALPHA | FL_DEPTH | FL_MULTISAMPLE);
  end();

  light0_dir[0] = 1.0;
  light0_dir[1] = 1.0;
  light0_dir[2] = -1.0;
  light0_dir[3] = 0.0;

  using mdtk::Exception;

  ml_ = new mdtk::SimLoop();

  setupPotentials(*ml_);
  if (base_state_filename != "")
  {
  if (base_state_filename.find("simloop.conf") != std::string::npos) 
  {
    ml_->loadstate();
  }
  else
  {
    ml_->initNLafterLoading = false;
    yaatk::text_ifstream fi(base_state_filename.c_str()); 
    if (base_state_filename.find("mde_init") != std::string::npos)
      ml_->loadFromStream(fi);
    else
      ml_->loadFromMDE(fi);
//	  ml_->loadFromMDE_OLD(fi);
    fi.close(); 
  }
  setData(*ml_);
  if (xvas.size() > 0 && xvas[0] != "shot")
    MDTrajectory_read(mdt,base_state_filename,xvas);
  else
  {
    if (yaatk::exists("snapshots.conf"))
      MDTrajectory_read_from_SnapshotList(mdt,base_state_filename);
  }
//    ctree = new CollisionTree(*(ml_->atoms.back()),mdt.begin(),mdt);
  }
  size_range(100, 100, 5000, 5000, 3*4, 3*4, 1);


  callback(window_cb);
}

void
VisBox::reArrange(double xmin, double xmax,
		  double ymin, double ymax,
		  double zmin, double zmax)
{
  double xmin_ = XMin+(XMax-XMin)*xmin;
  double xmax_ = XMin+(XMax-XMin)*xmax;
  double ymin_ = YMin+(YMax-YMin)*ymin;
  double ymax_ = YMin+(YMax-YMin)*ymax;
  double zmin_ = ZMin+(ZMax-ZMin)*zmin;
  double zmax_ = ZMin+(ZMax-ZMin)*zmax;

  int i;
  int VC = Ro.size();
	
  R.resize(0);
	
  for(i=0;i<VC;i++)
  {
    if (Ro[i].coords.x >= xmin_ && Ro[i].coords.x <= xmax_ &&
	Ro[i].coords.y >= ymin_ && Ro[i].coords.y <= ymax_ &&
	Ro[i].coords.z >= zmin_ && Ro[i].coords.z <= zmax_)
    {
      R.push_back(Ro[i]);
    }
  };	

  GLdouble /*DistMin,*/DistMax,XMin,XMax,YMin,YMax,ZMin,ZMax;
	
  VC = R.size();
	
  XMin=R[0].coords.x;
  XMax=R[0].coords.x;
  YMin=R[0].coords.y;
  YMax=R[0].coords.y;
  ZMin=R[0].coords.z;
  ZMax=R[0].coords.z;
  for(int i=0;i<VC;i++)
  {
    if (R[i].coords.x<XMin)
    {
      XMin=R[i].coords.x;
    }
    else
      if (R[i].coords.x>XMax)
      {
	XMax=R[i].coords.x;
      };
    if (R[i].coords.y<YMin)
    {
      YMin=R[i].coords.y;
    }
    else
      if (R[i].coords.y>YMax)
      {
	YMax=R[i].coords.y;
      };
    if (R[i].coords.z<ZMin)
    {
      ZMin=R[i].coords.z;
    }
    else
      if (R[i].coords.z>ZMax)
      {
	ZMax=R[i].coords.z;
      };
  };

  if (allowRescale)
  {
    XCenter=(XMin+XMax)/2;
    YCenter=(YMin+YMax)/2;
    ZCenter=(ZMin+ZMax)/2;
  }	
  /*DistMin=*/DistMax=sqrt(SQR(XMax-XMin)+SQR(YMax-YMin)+SQR(ZMax-ZMin));
  //	VertexRadius=DistMin/4;
  vertexRadius = 2.57*mdtk::Ao/2.0/3.0/2.0;
  axesRadius = vertexRadius*3;
  if (allowRescale)
  {
    scale=(2.0*nRange)/(DistMax+2.0*vertexRadius);
  }  
  maxScale=2.0*(scale);
  redraw();
}  

void
VisBox::setData(mdtk::SimLoop &ml)
{
  if (unfoldPBC)
    ml.atoms.unfoldPBC();

  Ro = ml.atoms;

  int VC = Ro.size();
	
  XMin=Ro[0].coords.x;
  XMax=Ro[0].coords.x;
  YMin=Ro[0].coords.y;
  YMax=Ro[0].coords.y;
  ZMin=Ro[0].coords.z;
  ZMax=Ro[0].coords.z;
  for(int i=0;i<VC;i++)
  {
    if (Ro[i].coords.x<XMin)
    {
      XMin=Ro[i].coords.x;
    }
    else
      if (Ro[i].coords.x>XMax)
      {
	XMax=Ro[i].coords.x;
      };
    if (Ro[i].coords.y<YMin)
    {
      YMin=Ro[i].coords.y;
    }
    else
      if (Ro[i].coords.y>YMax)
      {
	YMax=Ro[i].coords.y;
      };
    if (Ro[i].coords.z<ZMin)
    {
      ZMin=Ro[i].coords.z;
    }
    else
      if (Ro[i].coords.z>ZMax)
      {
	ZMax=Ro[i].coords.z;
      };
  };

  reArrange(-1,101,-1,101,-1,101);
}

void
VisBox::drawObjects()
{
  glInitNames();

  glPushMatrix();

  glLightfv(GL_LIGHT0, GL_POSITION, light0_dir);

  glScaled(scale,scale,scale);
  glTranslated(-XCenter, -YCenter, -ZCenter);
  glEnable(GL_LIGHTING);

  if (showAtoms)
    listVertexes();

  if (showCTree)
    listCTree();
  if (showAxes)
  {
//    glDisable(GL_LIGHTING);
    listAxes();
//    glEnable(GL_LIGHTING);
  }
  if (showBarrier)
  {
    glDisable(GL_LIGHTING);
    listBarrier();
    glEnable(GL_LIGHTING);
  }
  if (showBath)
  {
    glDisable(GL_LIGHTING);
    listThermalBath();
    glEnable(GL_LIGHTING);
  }
  if (showBathSketch)
  {
    glDisable(GL_LIGHTING);
    listThermalBathSketch();
    glEnable(GL_LIGHTING);
  }
  if (showCustom1)
  {
    glDisable(GL_LIGHTING);
    listCustom1();
    glEnable(GL_LIGHTING);
  }
  if (showCustom2)
  {
    glDisable(GL_LIGHTING);
    listCustom2();
    glEnable(GL_LIGHTING);
  }
  if (showCustom3)
  {
//    glDisable(GL_LIGHTING);
    listCustom3();
//    glEnable(GL_LIGHTING);
  }
  glPopMatrix();
}

void
VisBox::listBarrier()
{
  glPushMatrix();
  glColor4ub(127,127,127,127);
  glBegin(GL_QUADS);
  glVertex3d(XMax,YMax,zbar);
  glVertex3d(XMin,YMax,zbar);
  glVertex3d(XMin,YMin,zbar);
  glVertex3d(XMax,YMin,zbar);
  glEnd();
  glBegin(GL_QUADS);
  glVertex3d(XMax,YMin,zbar);
  glVertex3d(XMin,YMin,zbar);
  glVertex3d(XMin,YMax,zbar);
  glVertex3d(XMax,YMax,zbar);
  glEnd();
  glPopMatrix();  
}

void
VisBox::listThermalBath()
{
  Float tb[3][2];

  tb[0][0] = ml_->thermalBath.dBoundary;
  tb[0][1] = -ml_->thermalBath.dBoundary+ml_->atoms.PBC().x;
  tb[1][0] = ml_->thermalBath.dBoundary;
  tb[1][1] = -ml_->thermalBath.dBoundary+ml_->atoms.PBC().y;
  tb[2][0] = ml_->thermalBath.zMinOfFreeZone;
  tb[2][1] = ml_->thermalBath.zMin;

  Float cb[3][2] = {{0,0},{0,0},{0,0}};
  cb[0][1] = ml_->atoms.PBC().x;
  cb[1][1] = ml_->atoms.PBC().y;
  cb[2][0] = tb[2][0];
  cb[2][1] = tb[2][1]+30*Ao;

  GLubyte tbc[4]={0,227,127,127};

  glColor4ubv(tbc);

  glPushMatrix();
  glBegin(GL_QUADS);
  glVertex3d(tb[0][0],tb[1][0],tb[2][1]);
  glVertex3d(tb[0][0],tb[1][1],tb[2][1]);
  glVertex3d(tb[0][1],tb[1][1],tb[2][1]);
  glVertex3d(tb[0][1],tb[1][0],tb[2][1]);
  glEnd();

  if (ml_->thermalBath.dBoundary != 0.0)
  {
    int vi = 0, vj;
    do
    {
      int i[2],j[2];
      if (vi == 0) vj = 1; else
	if (vi == 1) vj = 3; else
	  if (vi == 3) vj = 2; else
	    if (vi == 2) vj = 0;
      i[0] = vi/2; i[1] = vi%2;
      j[0] = vj/2; j[1] = vj%2;

      glBegin(GL_QUADS);
      glVertex3d(tb[0][i[0]],tb[1][i[1]],tb[2][1]);
      glVertex3d(tb[0][i[0]],tb[1][i[1]],tb[2][0]);
      glVertex3d(tb[0][j[0]],tb[1][j[1]],tb[2][0]);
      glVertex3d(tb[0][j[0]],tb[1][j[1]],tb[2][1]);
      glEnd();

      ++vi;
    }while(vi != 4);
  }

  glPopMatrix();  
}

void
VisBox::listThermalBathSketch()
{
  Float tb[3][2];

  tb[0][0] = ml_->thermalBath.dBoundary;
  tb[0][1] = -ml_->thermalBath.dBoundary+ml_->atoms.PBC().x;
  tb[1][0] = ml_->thermalBath.dBoundary;
  tb[1][1] = -ml_->thermalBath.dBoundary+ml_->atoms.PBC().y;
  tb[2][0] = ml_->thermalBath.zMinOfFreeZone;
  tb[2][1] = ml_->thermalBath.zMin;

  Float cb[3][2] = {{0,0},{0,0},{0,0}};
  cb[0][1] = ml_->atoms.PBC().x;
  cb[1][1] = ml_->atoms.PBC().y;
  cb[2][0] = tb[2][0];
  cb[2][1] = tb[2][1]+30*Ao;

  GLubyte tbc[4]={0,227,127,127};

  glColor4ubv(tbc);

  glPushMatrix();
  glBegin(GL_QUADS);
  glVertex3d(cb[0][0],cb[1][1]/2.0+cb[1][0],cb[2][1]);
  glVertex3d(cb[0][0],cb[1][1]/2.0+cb[1][0],tb[2][1]);
  glVertex3d(cb[0][1],cb[1][1]/2.0+cb[1][0],tb[2][1]);
  glVertex3d(cb[0][1],cb[1][1]/2.0+cb[1][0],cb[2][1]);
  glEnd();

  glBegin(GL_QUADS);
  glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][0],cb[2][1]);
  glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][0],tb[2][1]);
  glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][1],tb[2][1]);
  glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][1],cb[2][1]);
  glEnd();

  if (ml_->thermalBath.dBoundary != 0.0)
  {
    glBegin(GL_QUADS);
    glVertex3d(cb[0][0],cb[1][1]/2.0+cb[1][0],cb[2][0]);
    glVertex3d(cb[0][0],cb[1][1]/2.0+cb[1][0],tb[2][1]);
    glVertex3d(tb[0][0],cb[1][1]/2.0+cb[1][0],tb[2][1]);
    glVertex3d(tb[0][0],cb[1][1]/2.0+cb[1][0],cb[2][0]);
    glEnd();
    glBegin(GL_QUADS);
    glVertex3d(cb[0][1],cb[1][1]/2.0+cb[1][0],cb[2][0]);
    glVertex3d(cb[0][1],cb[1][1]/2.0+cb[1][0],tb[2][1]);
    glVertex3d(tb[0][1],cb[1][1]/2.0+cb[1][0],tb[2][1]);
    glVertex3d(tb[0][1],cb[1][1]/2.0+cb[1][0],cb[2][0]);
    glEnd();

    glBegin(GL_QUADS);
    glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][0],cb[2][0]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][0],tb[2][1]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],tb[1][0],tb[2][1]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],tb[1][0],cb[2][0]);
    glEnd();
    glBegin(GL_QUADS);
    glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][1],cb[2][0]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],cb[1][1],tb[2][1]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],tb[1][1],tb[2][1]);
    glVertex3d(cb[0][1]/2.0+cb[0][0],tb[1][1],cb[2][0]);
    glEnd();
  }

  glPopMatrix();  
}

void
VisBox::listAxes()
{
  Vector3D xyz000(0,0,0);
  Vector3D x(XMax,0,0);
  Vector3D y(0,YMax,0);
  Vector3D z(0,0,ZMax);

  drawArrow(x,xyz000,0x0000FF,
            axesRadius,1.0/20);
  drawArrow(y,xyz000,0x00FF00,
            axesRadius,1.0/20);
  drawArrow(z,xyz000,0xFF0000,
            axesRadius,1.0/20);
}

void
VisBox::listVertexes()
{
  size_t i;
  Color c;
  for(i=0;i<R.size();i++)
  {
    glPushName(i);
    glPushMatrix();
    switch (nativeVertexColors)
    {
    case false: 
      c = vertexColor;
      break;
    case true:  
      switch (R[i].ID)
      {
      case H_EL:  c = (0x00FF00); break; 
      case C_EL:  c = (0x0000FF); break; 
      case Ar_EL: c = (0xFF00FF); break; 
      case Xe_EL: c = (0xFF0000); break; 
      case Cu_EL: c = (0x00FFFF); break; 
      case Ag_EL: c = (0x00FFFF); break; 
      case Au_EL: c = (0x00FFFF); break; 
      default: c = R[i].tag;
      };

      if (showSelected)
      {
	if (i!=selectedAtomIndex)
	{
	  if ((R[i].coords-R[selectedAtomIndex].coords).module() < 5.5*mdtk::Ao)
	  {
	    c = (0x777777);
	    if ((R[i].coords-R[selectedAtomIndex].coords).module() < 2.5*mdtk::Ao)
	      c = (0xFFFF00);
	  }  
	}  
	else  
	  c = (0xFFFFFF);
      }
    }
    myglColor(c);
    if (R[i].M > 1000.0*amu) myglColor(0x0);

    if (atomsQuality > 2)
    {
      glTranslated(R[i].coords.x,R[i].coords.y,R[i].coords.z);
      Atom a = R[i]; a.setAttributesByElementID();
      Float scale = 1.0*vertexRadius*pow(a.M/mdtk::amu,1.0/3.0);
      if (tinyAtoms) scale /= 5;
      glScaled(scale,scale,scale);
      glCallList(hqMode?lstBallHQ:lstBall);
    }
    else
    {
      glPointSize(vertexRadius*scale);
      glBegin(GL_POINTS);
      glVertex3d(R[i].coords.x,R[i].coords.y,R[i].coords.z);
      glEnd();
    }
    glPopMatrix();
    glPopName();
  }
}

Vector3D
_vectorMul(const Vector3D& a, const Vector3D& b)
{
  return Vector3D(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

Float
_scalarMul(const Vector3D& a, const Vector3D& b)
{
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

Float
_relAngle(const Vector3D& a, const Vector3D& b)
{
  return acos(_scalarMul(a,b)/(a.module()*b.module()));
}

void
VisBox::drawEdge(const Vector3D& vi, const Vector3D& vj, 
		 unsigned int color, double radius)
{
  Vector3D TempRotVector;
  double    TempRotAngle;

  glPushMatrix();
  myglColor(color);
  glTranslated(vi.x,vi.y,vi.z);
  TempRotVector=_vectorMul(Vector3D(0,0,1.0L),vj-vi);
  TempRotAngle=(_relAngle(Vector3D(0,0,1.0L),vj-vi)/M_PI)*180.0L;
  glRotated(TempRotAngle,TempRotVector.x,TempRotVector.y,TempRotVector.z);
  glScaled(radius,radius,(vi-vj).module());
  glCallList(hqMode?lstStickHQ:lstStick);
  glPopMatrix();
}

void
VisBox::drawArrow(const Vector3D& vi, const Vector3D& vj, 
		  unsigned int color, double radius, Float arrowPart)
{
  Vector3D TempRotVector;
  double    TempRotAngle;

  GLUquadricObj *quadObj;

  glPushMatrix();
  quadObj = gluNewQuadric ();
  gluQuadricDrawStyle (quadObj, GLU_FILL);
  myglColor(color);
  Vector3D vish = vi+(vj-vi)*arrowPart;
  glTranslated(vish.x,vish.y,vish.z);
  TempRotVector=_vectorMul(Vector3D(0,0,1.0L),vj-vi);
  TempRotAngle=(_relAngle(Vector3D(0,0,1.0L),vj-vi)/M_PI)*180.0L;
  glRotated(TempRotAngle,TempRotVector.x,TempRotVector.y,TempRotVector.z);
  gluCylinder (quadObj,
	       radius/3.0,
	       radius/3.0,
	       (vi-vj).module()*(1.0-arrowPart), 
	       hqMode?atomsQualityInHQMode:atomsQuality, 
	       hqMode?atomsQualityInHQMode:atomsQuality);
  gluDeleteQuadric(quadObj);
  glPopMatrix();

  glPushMatrix();
  quadObj = gluNewQuadric ();
  gluQuadricDrawStyle (quadObj, GLU_FILL);
  myglColor(color);
  glTranslated(vi.x,vi.y,vi.z);
  TempRotVector=_vectorMul(Vector3D(0,0,1.0L),vj-vi);
  TempRotAngle=(_relAngle(Vector3D(0,0,1.0L),vj-vi)/M_PI)*180.0L;
  glRotated(TempRotAngle,TempRotVector.x,TempRotVector.y,TempRotVector.z);
  gluCylinder (quadObj,
	       0,
	       radius,
	       (vi-vj).module()*arrowPart, 
	       hqMode?atomsQualityInHQMode:atomsQuality, 
	       hqMode?atomsQualityInHQMode:atomsQuality);
  gluDeleteQuadric(quadObj);
  glPopMatrix();
}

void
VisBox::drawCTree(CollisionTree* ct)
{
  Atom& a = ct->a;
  if (ct->t1)
  {
    Atom& a1 = ct->t1->a;
    drawEdge(a.coords,a1.coords,
	     0xFF0000,vertexRadius*pow(a.M/mdtk::amu,1.0/3.0));
    TRACE(a.globalIndex);
    TRACE(a1.globalIndex);
    drawCTree(ct->t1);
  }
  if (ct->t2)
  {
    Atom& a2 = ct->t2->a;
    drawEdge(a.coords,a2.coords,
	     0xFF0000,vertexRadius*pow(a.M/mdtk::amu,1.0/3.0));
    TRACE(a.globalIndex);
    TRACE(a2.globalIndex);
    drawCTree(ct->t2);
  }
}

void
VisBox::listCTree()
{
//  CTree_List(ctree);
/*
  size_t i = 0;
  size_t j = 10695;
  Draw_Edge(R[i].coords,R[j].coords,0xFF0000);
*/
  std::vector<bool> ignore(mdt.begin()->second.size());
  std::vector<bool> hadEnteredCollision(mdt.begin()->second.size());

  MDTrajectory::const_iterator t = mdt.begin();
  while (t != mdt.end())
  {
    if (!showCTreeAllTimes && t->first > ml_->simTime)
      break;

    const std::vector<Atom>& atoms = t->second;

    MDTrajectory::const_iterator t_prev = t;
    if (t != mdt.begin()) --t_prev;
    const std::vector<Atom>& atoms_prev = t_prev->second;

    for(size_t i = 0; i < atoms.size(); ++i)
    {
      const Atom& a = atoms[i];
      if (a.ID == Cu_EL && showCustom3) continue;
      Float Ek = a.M*SQR(a.V.module())/2.0;
      if (Ek > energyThresholdCTree*eV && !ignore[i])
      {
	hadEnteredCollision[i] = true;
	Color c;
	switch (R[i].ID)
	{
	case H_EL:  c = (0x00FF00); break; 
	case C_EL:  c = (0x0000FF); break; 
	case Ar_EL: c = (0xFF00FF); break; 
	case Cu_EL: c = (0x00FFFF); break; 
	case Ag_EL: c = (0x00FFFF); break; 
	case Au_EL: c = (0x00FFFF); break; 
	default: c = R[i].tag;
	}

	if (showCTreeAtoms)
	{
	  myglColor(c);
	  glPushMatrix();
	  glTranslated(a.coords.x,a.coords.y,a.coords.z);
          Float scale = vertexRadius*pow(a.M/mdtk::amu,1.0/3.0)/downscaleCTree;
          glScaled(scale,scale,scale);
          glCallList(hqMode?lstBallHQ:lstBall);
	  glPopMatrix();
	}

	if (showCTreeConnected)
	{
	  const Atom& a_prev = atoms_prev[i];
	  if (t != mdt.begin())
	    drawEdge(a.coords,a_prev.coords,c,
		     vertexRadius*pow(a.M/mdtk::amu,1.0/3.0)/downscaleCTree);
	}
      }
      else
      {
	if (hadEnteredCollision[i])
	  ignore[i] = true;
      }
    }

    ++t;
  }
}

#include <gsl/gsl_qrng.h>

void
VisBox::listCustom1()
{
  std::vector<Atom> cluster;
  for(size_t i = 0; i < ml_->atoms.size(); i++)
  {
    Atom& a = ml_->atoms[i];
    if (a.ID == Cu_EL) cluster.push_back(a);
  }

  Float halo = 1.3*Ao;
  Float depth = 0.0;
  {
    halo = 5.0*Ao;
    glColor4ub(0,0,255,255);
    depth = -15.0*Ao;

    for(size_t i = 0; i < cluster.size(); i++)
    {
      glPushMatrix();
      glTranslated(cluster[i].coords.x,cluster[i].coords.y,depth);
      glScaled(halo,halo,halo);
      glCallList(hqMode?lstBallHQ:lstBall);
      glPopMatrix();  
    }

    glColor4ub(0,255,0,255);
    halo = 1.3*Ao;
    depth = -22.0*Ao;

    for(size_t i = 0; i < cluster.size(); i++)
    {
      glPushMatrix();
      glTranslated(cluster[i].coords.x,cluster[i].coords.y,depth);
      glScaled(halo,halo,halo);
      glCallList(hqMode?lstBallHQ:lstBall);
      glPopMatrix();  
    }

    glColor4ub(255,255,0,255);
    halo = 0.3*Ao;
    depth = -24.0*Ao;

    for(size_t i = 0; i < cluster.size(); i++)
    {
      glPushMatrix();
      glTranslated(cluster[i].coords.x,cluster[i].coords.y,depth);
      glScaled(halo,halo,halo);
      glCallList(hqMode?lstBallHQ:lstBall);
      glPopMatrix();  
    }
  }
}

void
VisBox::listCustom2()
{
  {
    const Atom& a = R.back();
    Vector3D c1 = a.coords;
    Vector3D c2 = c1; c1.z += 5.0*Ao;
    drawArrow(c1,c2,0xFF0000,
	      vertexRadius*pow(a.M/mdtk::amu,1.0/3.0));
  }
}

void
VisBox::listCustom3()
{
  MDTrajectory::const_iterator t = mdt.begin();
  while (t != mdt.end())
  {
    if (!showCTreeAllTimes && t->first > ml_->simTime)
      break;

    const std::vector<Atom>& atoms = t->second;

    MDTrajectory::const_iterator t_prev = t;
    if (t != mdt.begin()) --t_prev;
    const std::vector<Atom>& atoms_prev = t_prev->second;
    
    Vector3D clusterMassCenter;
    {
      mdtk::Vector3D sumOfC = 0.0;
      Float sumOfM = 0.0;
      for(size_t ai = 0; ai < atoms.size(); ai++)
      {
	const mdtk::Atom& atom = atoms[ai];
	if (atom.ID != Cu_EL) continue;
	sumOfM += atom.M;
	sumOfC += atom.coords*atom.M;
      };
      REQUIRE(sumOfM > 0.0);
      clusterMassCenter = sumOfC/sumOfM;
    }
/*
    Vector3D clusterMassCenterPrev;
    {
      mdtk::Vector3D sumOfC = 0.0;
      Float sumOfM = 0.0;
      for(size_t ai = 0; ai < atoms_prev.size(); ai++)
      {
	const mdtk::Atom& atom = atoms_prev[ai];
	if (atom.ID != Cu_EL) continue;
	sumOfM += atom.M;
	sumOfC += atom.coords*atom.M;
      };
      REQUIRE(sumOfM > 0.0);
      clusterMassCenterPrev = sumOfC/sumOfM;
    }
    
    Color c;
    c = (0x00FFFF);

    myglColor(c);
    
    if (t != mdt.begin())
      drawEdge(clusterMassCenter,clusterMassCenterPrev,c,
	       vertexRadius*pow(64,1.0/3.0));
*/
    Color c;
    c = (0x00FFFF);

    myglColor(c);

    glPushMatrix();
    glTranslated(clusterMassCenter.x,
		 clusterMassCenter.y,
		 clusterMassCenter.z);
    Float scale = vertexRadius*pow(64,1.0/3.0)/2;
    glScaled(scale,scale,scale);
    glCallList(hqMode?lstBallHQ:lstBall);
    glPopMatrix();
    
    ++t;
  }
}

void
VisBox::setupLighting()
{
  ////////// Material
  GLfloat MaterialAmbient[] = {0.5, 0.5, 0.5, 1.0};
  GLfloat MaterialDiffuse[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat MaterialSpecular[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat MaterialShininess[] = {100.0};

  glMaterialfv(GL_FRONT, GL_AMBIENT, MaterialAmbient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, MaterialDiffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, MaterialSpecular);
  glMaterialfv(GL_FRONT, GL_SHININESS, MaterialShininess);

  /////////// Global Lights

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
	
  GLfloat global_ambient[] = {0.2, 0.2, 0.2, 1.0};
	
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

  // Light0

  //	GLfloat	light0_pos[] = {100.0, 100.0, 100.0, 0.0}:
  //	GLfloat	light0_dir[] = {1.0, 1.0, 1.0, 0.0};
	
  GLfloat	light0_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  GLfloat	light0_ambient[] = {0.0, 0.0, 0.0, 1.0};
  GLfloat	light0_specular[] = {0.0, 0.0, 0.0, 1.0};
	
  glLightfv(GL_LIGHT0, GL_POSITION, light0_dir);
  glLightfv(GL_LIGHT0, GL_AMBIENT,  light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);

  // Final ...

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  //	glShadeModel(GL_SMOOTH);
	
  glEnable(GL_DEPTH_TEST);
  //	glEnable(GL_CULL_FACE);
  glEnable(GL_ALPHA_TEST);
  glEnable(GL_NORMALIZE);

  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_POINT_SMOOTH_HINT, GL_DONT_CARE);
  glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
}

void
VisBox::prepareBasicLists()
{
  glNewList(lstBall,GL_COMPILE);
  {
    GLUquadricObj *quadObj;
    quadObj = gluNewQuadric ();
    gluQuadricDrawStyle (quadObj, GLU_FILL);
    gluSphere (quadObj,
	       1.0,
	       atomsQuality,
	       atomsQuality);
    gluDeleteQuadric(quadObj);
  }
  glEndList();

  glNewList(lstBallHQ,GL_COMPILE);
  {
    GLUquadricObj *quadObj;
    quadObj = gluNewQuadric ();
    gluQuadricDrawStyle (quadObj, GLU_FILL);
    gluSphere (quadObj,
	       1.0,
	       atomsQualityInHQMode,
	       atomsQualityInHQMode);
    gluDeleteQuadric(quadObj);
  }
  glEndList();

  glNewList(lstStick,GL_COMPILE);
  {
    GLUquadricObj *quadObj;
    quadObj = gluNewQuadric ();
    gluQuadricDrawStyle (quadObj, GLU_FILL);
    gluCylinder (quadObj,
                 1.0,
                 1.0,
                 1.0, 
                 atomsQuality, 
                 atomsQuality);
    gluDeleteQuadric(quadObj);
  }
  glEndList();

  glNewList(lstStickHQ,GL_COMPILE);
  {
    GLUquadricObj *quadObj;
    quadObj = gluNewQuadric ();
    gluQuadricDrawStyle (quadObj, GLU_FILL);
    gluCylinder (quadObj,
                 1.0,
                 1.0,
                 1.0, 
                 atomsQualityInHQMode, 
                 atomsQualityInHQMode);
    gluDeleteQuadric(quadObj);
  }
  glEndList();
}

void
VisBox::onResizeGL()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,(w()>h())?h():w(),(w()>h())?h():w());
//  glOrtho (-nRange, nRange, -nRange, nRange, -nRange, nRange);
  GLfloat xr[3][2];

  if (!tiledMode)
  {
    xr[0][0] = -nRange;
    xr[0][1] = +nRange;
    xr[1][0] = -nRange;
    xr[1][1] = +nRange;
    xr[2][0] = -nRange;
    xr[2][1] = +nRange;
  }
  else
  {
    xr[0][0] = -nRange+2.0*nRange/tileCount*(tileIndex[0]);
    xr[0][1] = -nRange+2.0*nRange/tileCount*(tileIndex[0]+1);
    xr[1][0] = -nRange+2.0*nRange/tileCount*(tileIndex[1]);
    xr[1][1] = -nRange+2.0*nRange/tileCount*(tileIndex[1]+1);
    xr[2][0] = -nRange;
    xr[2][1] = +nRange;
  }
  glOrtho (xr[0][0], xr[0][1], xr[1][0], xr[1][1], xr[2][0], xr[2][1]);
}

void
VisBox::draw()
{
  if (!valid())
  {
    onResizeGL();
 
    setupLighting();
  
    if (!lstBall) lstBall = glGenLists(1);
    REQUIRE(lstBall);

    if (!lstBallHQ) lstBallHQ = glGenLists(1);
    REQUIRE(lstBallHQ);

    if (!lstStick) lstStick = glGenLists(1);
    REQUIRE(lstStick);

    if (!lstStickHQ) lstStickHQ = glGenLists(1);
    REQUIRE(lstStickHQ);

    prepareBasicLists();
  }

  setupLighting();

  glClearColor((bgColor%0x100)/255.0,
	       ((bgColor/0x100)%0x100)/255.0,
	       (bgColor/0x10000)/255.0,1.0f);
  glClearDepth(1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  drawObjects();
  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
}

void
VisBox::myglColor(Color color, GLubyte alpha)
{
  color = color % 0x1000000;
  glColor4ub(color%0x100,((color/0x100)%0x100),((color/0x10000)%0x100),
	     alpha);
  REQUIRE(((color/0x10000)%0x100) == color/0x10000);
//  glColor4ub(color%0x100,((color/0x100)%0x100),color/0x10000,0xFF);
}

void
VisBox::saveToMDE(char* filename)
{
  using namespace mdtk;  
  using namespace std;  
  yaatk::text_ofstream fo(filename);
  {
    ml_->saveToMDE(fo);
  }

  fo.close();
}  

void
VisBox::saveImageToFile(char* filename)
{
  hqMode = true;

  unsigned long width = w(); unsigned long height = h();
  
  unsigned char *d = new unsigned char[width*height*3];
  REQUIRE(d!=NULL);

  glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,d);
  for(unsigned long i = 0;i < width*height*3; i++)
  {
    if (i % 3 == 0)
    {
      unsigned char t;
      t = d[i];
      d[i] = d[i+2];
      d[i+2] = t;
    }	 
  }
  grctk::bmpImage* bmp = new grctk::bmpImage(width,height,d);

  bmp->SaveToFile(filename);


  delete [] d;
  delete bmp;

  hqMode = false;
}

void
VisBox::saveTiledImageToFile(char* filename)
{
  REQUIRE(w()==h());
  tiledMode = true;
  hqMode = true;

  unsigned long width = w()*tileCount; unsigned long height = h()*tileCount;
  
  unsigned char *d = new unsigned char[width*height*3];
  REQUIRE(d!=NULL);
  unsigned char *dtile = new unsigned char[w()*h()*3];
  REQUIRE(dtile!=NULL);

  grctk::bmpImage* bmp = new grctk::bmpImage(width,height,d);

  for(tileIndex[0] = 0; tileIndex[0] < tileCount; tileIndex[0]++)
  {
    for(tileIndex[1] = 0; tileIndex[1] < tileCount; tileIndex[1]++)
    {
      onResizeGL();
      redraw();
      while (!Fl::ready()) {};
      Fl::flush();
      while (!Fl::ready()) {};

      glReadPixels(0,0,w(),h(),GL_RGB,GL_UNSIGNED_BYTE,dtile);

      for(unsigned long i = 0;i < ((unsigned long) w())*h()*3; i++)
      {
	if (i % 3 == 0)
	{
	  unsigned char t;
	  t = dtile[i];
	  dtile[i] = dtile[i+2];
	  dtile[i+2] = t;
	}	 
      }

      grctk::bmpImage* tilebmp = new grctk::bmpImage(w(),h(),dtile);
//      tilebmp->SaveToFile("xxx.bmp");
      for(unsigned long i = 0;i < (unsigned long) h(); i++)
	for(unsigned long j = 0;j < (unsigned long) w(); j++)
	  bmp->setPixel(i+tileIndex[1]*w(),j+tileIndex[0]*h(),
			tilebmp->getPixel(i,j));
    }
  }

  bmp->SaveToFile(filename);

  delete [] dtile;
  delete [] d;
  delete bmp;

  hqMode = false;
  tiledMode = false;

  onResizeGL();
  redraw();
}

void
VisBox::rollAround(double angle,double x, double y,double z)
{
  glMatrixMode(GL_MODELVIEW);

  glRotated(angle,x,y,z);

  if (x == 1.0) old_rot_x += angle;
  if (y == 1.0) old_rot_y += angle;
  if (z == 1.0) old_rot_z += angle;

  redraw();
}

void
VisBox::window_cb(Fl_Widget* widget, void*)
{
//    if (fl_choice("Do you really want to exit?","No","Yes",NULL)==1)
  {
    ((Fl_Window*)widget)->hide();
    exit(0);
  }
}

int
VisBox::handle(int event)
{
  int x,y;
  int atom2select;
  switch(event)
  {
    case FL_PUSH:
      x = Fl::event_x();
      y = Fl::event_y();
      atom2select = pickAtom(x,y);
      if (atom2select >=0)
        MainWindow_GlobalPtr->setAtomViewIndex(atom2select);
      return 1;
    default:
      return Fl_Widget::handle(event);
  }
}

int
VisBox::pickAtom(int x, int y)
{
  int pickedAtom = -1;

  const int BUFSIZE = 4*R.size();
  GLuint* selectBuf = new GLuint[BUFSIZE];

  {
    glSelectBuffer(BUFSIZE,selectBuf);
    glRenderMode(GL_SELECT);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,viewport);
    gluPickMatrix(x,viewport[3]-y,
                  5,5,viewport);
    GLfloat xr[3][2];
    {
      xr[0][0] = -nRange;
      xr[0][1] = +nRange;
      xr[1][0] = -nRange;
      xr[1][1] = +nRange;
      xr[2][0] = -nRange;
      xr[2][1] = +nRange;
    }
    glOrtho (xr[0][0], xr[0][1], xr[1][0], xr[1][1], xr[2][0], xr[2][1]);
  }
  
  draw();

  {
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glFlush();
  }

  GLint hits;
  hits = glRenderMode(GL_RENDER);

  {
    GLint i, j, numberOfNames = 0;
    GLuint names, *ptr, minZ,*ptrNames;

    ptr = (GLuint *) selectBuf;
    minZ = 0xffffffff;
    for (i = 0; i < hits; i++)
    {
      names = *ptr;
      ptr++;
      if (*ptr < minZ)
      {
        numberOfNames = names;                                 
        minZ = *ptr;                                           
        ptrNames = ptr+2;                                      
      }                                                        
      
      ptr += names+2;                                          
    }
    if (numberOfNames > 0)
    {
      ptr = ptrNames;
      for (j = 0; j < numberOfNames; j++,ptr++)
        pickedAtom = *ptr;
    }
  }

  delete [] selectBuf;

  return pickedAtom;
}

}

