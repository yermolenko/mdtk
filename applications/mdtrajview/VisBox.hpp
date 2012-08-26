/*
   The VisBox class for the molecular dynamics trajectory viewer
   (header file)

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

#ifndef	mde_VisBox_h
#define	mde_VisBox_h

#include <FL/fl_ask.H>
#ifdef __WIN32__
  #include <windows.h>
#endif  
#include <GL/gl.h>
#include <GL/glu.h>
#include <FL/Fl_Gl_Window.H>
#include "bmpImage.hpp"
#include "mdtk/Atom.hpp"
#include "mdtk/config.hpp"
#include "mdtk/consts.hpp"
#include "mdtk/SimLoop.hpp"

#include "CollisionTree.hpp"
#include "MainWindow.hpp"

namespace xmde
{

typedef unsigned int Color;

class MainWindow;

class VisBox : public Fl_Gl_Window
{
public:
  bool allowRescale;
  unsigned long vertexColor,edgeColor,bgColor;
  bool showAxes;
  bool showCTree;
  bool showCTreeConnected;
  bool showCTreeAtoms;
  bool showCTreeAllTimes;
  bool tinyAtoms;
  Float downscaleCTree;
  Float energyThresholdCTree;
  bool showAtoms;
  bool showBath;
  bool showBathSketch;
  bool showBonds;
  bool unfoldPBC;
  bool showCustom1;
  bool showCustom2;
  bool showCustom3;
  bool showSelected;
  bool showBarrier;
  bool nativeVertexColors;
  int  atomsQuality;
  int  atomsQualityInHQMode;
  bool hqMode;
  size_t selectedAtomIndex;
  Atom* getSelectedAtomPtr() {return &(Ro[selectedAtomIndex]);}

  GLfloat nRange;
  Float vertexRadius, axesRadius, scale, maxScale,
    XMin,XMax,YMin,YMax,ZMin,ZMax, XCenter, YCenter, ZCenter;

  GLfloat  light0_dir[4];

public:
  mdtk::AtomsArray* getAtoms(){return &Ro;};
  size_t getAtomsCount(){return ml_->atoms.size();};
  void loadNewSnapshot(size_t);
private:
  mdtk::AtomsArray R,Ro;
  std::vector<bool> completeInfoPresent;
public:
  mdtk::SimLoop* ml_;
  MDTrajectory mdt;
private:
  CollisionTree *ctree;
private:
  mdtk::Float zbar;

  void draw();
  void drawcube(int wire);
  void drawObjects();
  void setupLighting();
  void myglColor(Color, GLubyte alpha = 0xFF);

  void onResizeGL();

  GLuint lstBall;
  GLuint lstBallHQ;
  GLuint lstStick;
  GLuint lstStickHQ;
  void prepareBasicLists();

  void listVertexes();
  void listAxes();
  void listBarrier();
  void listThermalBath();
  void listThermalBathSketch();
  void listCTree();
  void listBonds();
  void listCustom1();
  void listCustom2();
  void listCustom3();
  void drawCTree(CollisionTree* ct);
  void drawEdge(const Vector3D& vi, const Vector3D& vj, 
		unsigned int color, double radius, GLubyte alpha = 255);
  void drawArrow(const Vector3D& vi, const Vector3D& vj, 
		 unsigned int color, double radius, 
                 Float arrowPart = 0.2);

public:
  double old_rot_x;
  double old_rot_y;  
  double old_rot_z;  

  void  reArrange(double xmin, double xmax,
		  double ymin, double ymax,
		  double zmin, double zmax);  

  void  rollAround(double,double,double,double);

  double RX(int i) {return R[i].coords.x;}
  double RY(int i) {return R[i].coords.y;}
  double RZ(int i) {return R[i].coords.z;}
	
public:
  VisBox(int x,int y,int w,int h);
  void loadDataFromFiles(std::string base_state_filename,
                    const std::vector<std::string>& xvas);
  virtual ~VisBox(){delete ml_;};

  void setData(mdtk::SimLoop &);
  void updateData(){setData(*ml_);}

  void saveImageToFile(char* filename);
  void saveTiledImageToFile(char* filename);
  void saveToMDE(char* filename);
  void saveState(char* filename);

  static void window_cb(Fl_Widget *, void *);

  bool tiledMode;
  int  tileCount;
  int  tileIndex[2];

  int handle(int event);
  int pickAtom(int x, int y);
};


inline Color combineRGB(unsigned char r,unsigned char g,unsigned char b,unsigned char a = 0)
{
  return (a*0x1000000+b*0x10000+g*0x100+r);
}

inline Color combineRGBA(unsigned char r,unsigned char g,unsigned char b,unsigned char a = 0)
{
  return (a*0x1000000+b*0x10000+g*0x100+r);
}

inline void analyseRGB(Color c,unsigned char &r,unsigned char &g,unsigned char &b)
{
  r = c%0x100;
  g = (c/0x100)%0x100;
  b = (c/0x10000)%0x100;
}

inline void analyseRGBA(Color c,unsigned char &r,unsigned char &g,unsigned char &b,unsigned char &a)
{
  r = c%0x100;
  g = (c/0x100)%0x100;
  b = (c/0x10000)%0x100;
  a = (c/0x1000000)%0x100;
}

}

#endif

