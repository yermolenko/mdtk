/*
   The VisBox class for the molecular dynamics trajectory viewer
   (header file)

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

namespace xmde
{

  typedef unsigned int Color;

  class VisBox : public Fl_Gl_Window
  {
  public:
    bool alloweRescale;
    void SetAllowRescale(bool);
    bool GetAllowRescale();
  public:
    size_t getAtomsCount(){return ml_->atoms_.size();};
    void loadNewML(std::string base_state_filename,std::string);
  private:
    int atoms_quality;
    GLdouble DistMin,DistMax,XMin,XMax,YMin,YMax,ZMin,ZMax;

    GLfloat nRange;
    GLdouble		XCenter,YCenter,ZCenter,
      VertexRadius,AxesRadius,
      Scale,MaxScale;
	
    mdtk::AtomsContainer R,Ro;
    mdtk::SimLoop* ml_;

    MDTrajectory mdt;
    std::vector<std::string> xvaList;
    CollisionTree *ctree;

    unsigned long VertexColor,EdgeColor,BGColor;
    bool EnableAxes;
    bool EnableBath;
    bool ShowSelected;
    bool EnableBarrier;
    mdtk::Float zbar;
    bool NativeVertexColors;
    bool FixedLights;

    bool Fronted;
    int min_ind,max_ind;

    void draw();
    void drawcube(int wire);
    void draw_objects();
    void SetupLighting();
    void myglColor(Color);
    void ChangeAxesState(bool);
    void ChangeNativeVertexColorsState(bool);

    void OnResizeGL();
    void ResetProjection();

    void Vertexes_List();
    void Axes_List();
    void Barrier_List();
    void ThermalBath_List();
    void CoolEdges_List();
    void CTree_List(CollisionTree* ct);
    void Draw_Edge(const Vector3D& vi, const Vector3D& vj, unsigned int color);

    GLfloat	light0_dir[4];

    double old_rot_x;
    double old_rot_y;  
    double old_rot_z;  

  public:
    void  ReArrange(double xmin, double xmax,
		    double ymin, double ymax,
		    double zmin, double zmax);  

    int GetVertexCount() {return R.size();}
	
    double GetRX(int i) {return R[i]->coords.x;}
    double GetRY(int i) {return R[i]->coords.y;}
    double GetRZ(int i) {return R[i]->coords.z;}
	
    double	MM_orig_data[16];
    bool	MM_orig;

    VisBox(int x,int y,int w,int h,std::string base_state_filename,
	   const std::vector<std::string>& xvas);
    virtual ~VisBox(){delete ml_;};

    void SetData(mdtk::SimLoop &);

    Float EnergyThreshold;
    void SetEnergyThreshold(Float et);

    void SetFixedLightsState(bool);
    bool GetFixedLightsState();
    void SetAxesState(bool State);
    bool GetAxesState();
    void SetBathState(bool State);
    bool GetBathState();
    void SetShowSelectedState(bool State);
    bool GetShowSelectedState();
    void SetBarrierState(bool State);
    bool GetBarrierState();
    void SetNativeVertexColorsState(bool State);
    bool GetNativeVertexColorsState();
    void SetBGColor(Color);
    void SetVertexColor(Color);
    Color GetBGColor() {return BGColor;};
    Color GetVertexColor() {return VertexColor;};
    void	SetScale(double);
    double	GetMaxScale();
    void  	SetMaxScale(double);
    double	GetScale();

    void  	SetLightXDir(double val);
    void  	SetLightYDir(double val);
    void  	SetLightZDir(double val);

    double  Get_old_rot_x();
    double  Get_old_rot_y();
    double  Get_old_rot_z();

    void SetAtomsQuality(int);

    void	RollAround(double,double,double,double);

    void QSaveImageToFile(char* filename);
    void SaveToMDE(char* filename);
	
    size_t selectedIndex;
	
    void selectAtom(size_t);
	
    mdtk::AtomsContainer* getAtoms(){return &Ro;};

    static void window_cb(Fl_Widget *, void *);
  };

  inline Color CombineRGB(unsigned char r,unsigned char g,unsigned char b,unsigned char a = 0)
  {
    return (a*0x1000000+b*0x10000+g*0x100+r);
  }

  inline Color CombineRGBA(unsigned char r,unsigned char g,unsigned char b,unsigned char a = 0)
  {
    return (a*0x1000000+b*0x10000+g*0x100+r);
  }

  inline void AnalyseRGB(Color c,unsigned char &r,unsigned char &g,unsigned char &b)
  {
    r = c%0x100;
    g = (c/0x100)%0x100;
    b = (c/0x10000)%0x100;
  }

  inline void AnalyseRGBA(Color c,unsigned char &r,unsigned char &g,unsigned char &b,unsigned char &a)
  {
    r = c%0x100;
    g = (c/0x100)%0x100;
    b = (c/0x10000)%0x100;
    a = (c/0x1000000)%0x100;
  }


}

#endif

