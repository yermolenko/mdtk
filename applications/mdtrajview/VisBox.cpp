/*
   The VisBox class for the molecular dynamics trajectory viewer

   Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010
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
	ml_->updateGlobalIndexes();
      }
      fi.close(); 
    }
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
    size_t xvaCount = 0;
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
      *(ml_->atoms[i]) = atoms[i];
    }
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
    showAxes(true),
    showBath(false),
    showSelected(false),
    showBarrier(false),
    nativeVertexColors(true),
    fixedLights(false),
    energyThreshold(10.0),
    nRange(50),
    vertexRadius(1.0), axesRadius(1.0), scale(1.0), maxScale(1.0),
    atomsQuality(14),
    selectedAtomIndex(0),
    R(),Ro(),
    ml_(NULL),
    mdt(),
    xvaList(xvas),
    ctree(NULL),
    zbar(0.0),
    old_rot_x(0.0), old_rot_y(0.0), old_rot_z(0.0), MM_orig(true)
{
  mode(FL_RGB | FL_DOUBLE | FL_ACCUM | FL_ALPHA | FL_DEPTH | FL_MULTISAMPLE);
  end();

  light0_dir[0] = 1.0;
  light0_dir[1] = 1.0;
  light0_dir[2] = 1.0;
  light0_dir[3] = 0.0;

  using mdtk::Exception;

  ml_ = new mdtk::SimLoop();

  setupPotentials(*ml_);
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
  size_range(100, 100, 5000, 5000, 3*4, 3*4, 1);

  MDTrajectory_read(mdt,base_state_filename,xvas);
//    ctree = new CollisionTree(*(ml_->atoms.back()),mdt.begin(),mdt);

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
    if (Ro[i]->coords.x >= xmin_ && Ro[i]->coords.x <= xmax_ &&
	Ro[i]->coords.y >= ymin_ && Ro[i]->coords.y <= ymax_ &&
	Ro[i]->coords.z >= zmin_ && Ro[i]->coords.z <= zmax_)
    {
      R.push_back(Ro[i]);
    }
  };	

  GLdouble /*DistMin,*/DistMax,XMin,XMax,YMin,YMax,ZMin,ZMax;
	
  VC = R.size();
	
  XMin=R[0]->coords.x;
  XMax=R[0]->coords.x;
  YMin=R[0]->coords.y;
  YMax=R[0]->coords.y;
  ZMin=R[0]->coords.z;
  ZMax=R[0]->coords.z;
  for(int i=0;i<VC;i++)
  {
    if (R[i]->coords.x<XMin)
    {
      XMin=R[i]->coords.x;
    }
    else
      if (R[i]->coords.x>XMax)
      {
	XMax=R[i]->coords.x;
      };
    if (R[i]->coords.y<YMin)
    {
      YMin=R[i]->coords.y;
    }
    else
      if (R[i]->coords.y>YMax)
      {
	YMax=R[i]->coords.y;
      };
    if (R[i]->coords.z<ZMin)
    {
      ZMin=R[i]->coords.z;
    }
    else
      if (R[i]->coords.z>ZMax)
      {
	ZMax=R[i]->coords.z;
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
  axesRadius = vertexRadius/(3*4);
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
  zbar = ml.barrier.z;
  Ro = ml.atoms_;

  int VC = Ro.size();
	
  XMin=Ro[0]->coords.x;
  XMax=Ro[0]->coords.x;
  YMin=Ro[0]->coords.y;
  YMax=Ro[0]->coords.y;
  ZMin=Ro[0]->coords.z;
  ZMax=Ro[0]->coords.z;
  for(int i=0;i<VC;i++)
  {
    if (Ro[i]->coords.x<XMin)
    {
      XMin=Ro[i]->coords.x;
    }
    else
      if (Ro[i]->coords.x>XMax)
      {
	XMax=Ro[i]->coords.x;
      };
    if (Ro[i]->coords.y<YMin)
    {
      YMin=Ro[i]->coords.y;
    }
    else
      if (Ro[i]->coords.y>YMax)
      {
	YMax=Ro[i]->coords.y;
      };
    if (Ro[i]->coords.z<ZMin)
    {
      ZMin=Ro[i]->coords.z;
    }
    else
      if (Ro[i]->coords.z>ZMax)
      {
	ZMax=Ro[i]->coords.z;
      };
  };

  reArrange(-1,101,-1,101,-1,101);
}

void
VisBox::drawObjects()
{
  glPushMatrix();

  glLightfv(GL_LIGHT0, GL_POSITION, light0_dir);

  glScaled(scale,scale,scale);
  glTranslated(-XCenter, -YCenter, -ZCenter);
  if (showBath)
  {
    glDisable(GL_LIGHTING);
    listThermalBath();
    glEnable(GL_LIGHTING);
  }
  glEnable(GL_LIGHTING);
  if (showAxes)
    listVertexes();
  listCTree();
  if (showAxes)
  {
    glDisable(GL_LIGHTING);
    listAxes();
    glEnable(GL_LIGHTING);
  }
  if (showBarrier)
  {
    glDisable(GL_LIGHTING);
    listBarrier();
    glEnable(GL_LIGHTING);
  }
  glPopMatrix();
}

void VisBox::listBarrier()
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
  Float tbXMin = ml_->thermalBath.dBoundary;
  Float tbXMax = -ml_->thermalBath.dBoundary+ml_->getPBC().x;
  Float tbYMin = ml_->thermalBath.dBoundary;
  Float tbYMax = -ml_->thermalBath.dBoundary+ml_->getPBC().y;
  Float tbZMin = ml_->thermalBath.zMinOfFreeZone;//ZMin;
  Float tbZMax = ml_->thermalBath.zMin;

  glPushMatrix();
  glColor4ub(0,227,127,127);
  glBegin(GL_QUADS);
  glVertex3d(tbXMax,tbYMax,tbZMax);
  glVertex3d(tbXMin,tbYMax,tbZMax);
  glVertex3d(tbXMin,tbYMin,tbZMax);
  glVertex3d(tbXMax,tbYMin,tbZMax);
  glEnd();

  if (ml_->thermalBath.dBoundary != 0.0)
  {
    glColor4ub(0,127,127,127);
    glBegin(GL_QUADS);
    glVertex3d(tbXMax,tbYMin,tbZMax);
    glVertex3d(tbXMin,tbYMin,tbZMax);
    glVertex3d(tbXMin,tbYMin,tbZMin);
    glVertex3d(tbXMax,tbYMin,tbZMin);
    glEnd();
    glColor4ub(0,127,127,127);
    glBegin(GL_QUADS);
    glVertex3d(tbXMax,tbYMax,tbZMax);
    glVertex3d(tbXMin,tbYMax,tbZMax);
    glVertex3d(tbXMin,tbYMax,tbZMin);
    glVertex3d(tbXMax,tbYMax,tbZMin);
    glEnd();
    glColor4ub(0,127,127,127);
    glBegin(GL_QUADS);
    glVertex3d(tbXMin,tbYMax,tbZMax);
    glVertex3d(tbXMin,tbYMin,tbZMax);
    glVertex3d(tbXMin,tbYMin,tbZMin);
    glVertex3d(tbXMin,tbYMax,tbZMin);
    glEnd();
    glColor4ub(0,127,127,127);
    glBegin(GL_QUADS);
    glVertex3d(tbXMax,tbYMax,tbZMax);
    glVertex3d(tbXMax,tbYMin,tbZMax);
    glVertex3d(tbXMax,tbYMin,tbZMin);
    glVertex3d(tbXMax,tbYMax,tbZMin);
    glEnd();
  }

  glPopMatrix();  
}

void
VisBox::listAxes()
{
  glBegin(GL_LINES);
  glColor3ub(255,0,0);
  glVertex3d(0,0,0);
  glVertex3d(100,0,0);
  glEnd();
  glBegin(GL_LINES);
  glColor3ub(0,255,0);
  glVertex3d(0,0,0);
  glVertex3d(0,100,0);
  glEnd();
  glBegin(GL_LINES);
  glColor3ub(0,0,255);
  glVertex3d(0,0,0);
  glVertex3d(0,0,100);
  glEnd();
}

void
VisBox::listVertexes()
{
  size_t i;
  GLUquadricObj *quadObj;
  Color c;
  for(i=0;i<R.size();i++)
  {
    glPushMatrix();
    switch (nativeVertexColors)
    {
    case false: 
      c = vertexColor;
      break;
    case true:  
      switch (R[i]->ID)
      {
      case H_EL:  c = (0x00FF00); break; 
      case C_EL:  c = (0x0000FF); break; 
      case Ar_EL: c = (0xFF00FF); break; 
      case Cu_EL: c = (0x00FFFF); break; 
      case Ag_EL: c = (0x00FFFF); break; 
      case Au_EL: c = (0x00FFFF); break; 
      default: c = R[i]->tag;
      };

      if (showSelected)
      {
	if (i!=selectedAtomIndex)
	{
	  if ((R[i]->coords-R[selectedAtomIndex]->coords).module() < 5.5*mdtk::Ao)
	  {
	    c = (0x777777);
	    if ((R[i]->coords-R[selectedAtomIndex]->coords).module() < 2.5*mdtk::Ao)
	      c = (0xFFFF00);
	  }  
	}  
	else  
	  c = (0xFFFFFF);
      }
    }
    myglColor(c);

    if (atomsQuality > 2)
    {
      quadObj = gluNewQuadric ();
      gluQuadricDrawStyle (quadObj, GLU_FILL);
      glTranslated(R[i]->coords.x,R[i]->coords.y,R[i]->coords.z);
      gluSphere (quadObj,
		 vertexRadius*pow(
		   ((R[i]->M < 1000.0*amu/*!=INFINITE_MASS*/)?
		    (R[i]->M):
		    (0.5*mdtk::amu))/mdtk::amu,1.0/3.0),
		 atomsQuality, atomsQuality);
      gluDeleteQuadric(quadObj);
    }
    else
    {
      glPointSize(vertexRadius*scale);
      glBegin(GL_POINTS);
      glVertex3d(R[i]->coords.x,R[i]->coords.y,R[i]->coords.z);
      glEnd();
    }
    glPopMatrix();
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
VisBox::drawEdge(const Vector3D& vi, const Vector3D& vj, unsigned int color)
{
  Vector3D TempRotVector;
  double    TempRotAngle;

  GLUquadricObj *quadObj;

  glPushMatrix();
  quadObj = gluNewQuadric ();
  gluQuadricDrawStyle (quadObj, GLU_FILL);
  myglColor(color);
  glTranslated(vi.x,vi.y,vi.z);
  TempRotVector=_vectorMul(Vector3D(0,0,1.0L),vj-vi);
  TempRotAngle=(_relAngle(Vector3D(0,0,1.0L),vj-vi)/M_PI)*180.0L;
  glRotated(TempRotAngle,TempRotVector.x,TempRotVector.y,TempRotVector.z);
  gluCylinder (quadObj,vertexRadius,vertexRadius,(vi-vj).module(), 6, 6);
  gluDeleteQuadric(quadObj);
  glPopMatrix();
}

void VisBox::drawCTree(CollisionTree* ct)
{
  Atom& a = ct->a;
  if (ct->t1)
  {
    Atom& a1 = ct->t1->a;
    drawEdge(a.coords,a1.coords,0xFF0000);
    TRACE(a.globalIndex);
    TRACE(a1.globalIndex);
    drawCTree(ct->t1);
  }
  if (ct->t2)
  {
    Atom& a2 = ct->t2->a;
    drawEdge(a.coords,a2.coords,0xFF0000);
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
  Draw_Edge(R[i]->coords,R[j]->coords,0xFF0000);
*/
  MDTrajectory::const_iterator t = mdt.begin();
  while (t != mdt.end())
  {
    const std::vector<Atom>& atoms = t->second;
    for(size_t i = 0; i < atoms.size(); ++i)
    {
      const Atom& a = atoms[i];
      Float Ek = a.M*SQR(a.V.module())/2.0;
      if (Ek > energyThreshold*eV)
      {
	Color c;
	switch (R[i]->ID)
	{
	case H_EL:  c = (0x00FF00); break; 
	case C_EL:  c = (0x0000FF); break; 
	case Ar_EL: c = (0xFF00FF); break; 
	case Cu_EL: c = (0x00FFFF); break; 
	case Ag_EL: c = (0x00FFFF); break; 
	case Au_EL: c = (0x00FFFF); break; 
	default: c = R[i]->tag;
	}
	myglColor(c);
	glPushMatrix();
	GLUquadricObj *quadObj;
	quadObj = gluNewQuadric ();
	gluQuadricDrawStyle (quadObj, GLU_FILL);
	glTranslated(a.coords.x,a.coords.y,a.coords.z);
	gluSphere (quadObj,
		   vertexRadius*pow(
		     ((a.M < 1000.0*amu)?
		      (a.M):
		      (0.5*mdtk::amu))/mdtk::amu,1.0/3.0),
		   4, 4);
	gluDeleteQuadric(quadObj);
	glPopMatrix();
      }
    }
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
VisBox::onResizeGL()
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glViewport(0,0,(w()>h())?h():w(),(w()>h())?h():w());
  glOrtho (-nRange, nRange, -nRange, nRange, -nRange, nRange);
/*
  glLoadIdentity();
  glViewport(0,0,w(),h());
  if (w() <= h())
  glOrtho (-nRange, nRange, -nRange*h()/w(), nRange*h()/w(), -nRange, nRange);
  else
  glOrtho (-nRange*w()/h(), nRange*w()/h(), -nRange, nRange, -nRange, nRange);
*/    
  setupLighting();
}

void
VisBox::resetProjection()
{
  glMatrixMode(GL_PROJECTION);
		
  if (MM_orig)
  {
    glLoadIdentity();
  }
  else
  {
    glLoadMatrixd(MM_orig_data);
  }
}

void
VisBox::draw()
{
  if (!valid())
  {
    onResizeGL();
    resetProjection();
  }

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
VisBox::myglColor(Color color)
{
  color = color % 0x1000000;
  glColor4ub(color%0x100,((color/0x100)%0x100),color/0x10000,0xFF);
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
}

void
VisBox::rollAround(double angle,double x, double y,double z)
{
  glMatrixMode(GL_PROJECTION);

  if (fixedLights)
  {
    glMatrixMode(GL_MODELVIEW);
  }
  else
  {
    glMatrixMode(GL_PROJECTION);
  }

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

}

