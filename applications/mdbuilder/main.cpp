/*
   mdbuilder (molecular dynamics experiments preparation tool)

   Copyright (C) 2010, 2011 Oleksandr Yermolenko
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

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>

#include <cstring>

//#include "pentacene.hpp"

#include <iostream>
#include <fstream>

#include "applications/common.h"
#include "experiments/H2.hpp"
#include "experiments/FCC.hpp"
#include "experiments/Clusters.hpp"

class MDBuilderWindow : public Fl_Gl_Window
{
  void draw();
  void resize(int X,int Y,int W,int H)
    {
      Fl_Gl_Window::resize(X,Y,W,H);
      glLoadIdentity();
      glViewport(0,0,W,H);
      glOrtho(-W,W,-H,H,-1,1);
      redraw();
    }
public:
  MDBuilderWindow(int X,int Y,int W,int H,const char*L=0) : Fl_Gl_Window(X,Y,W,H,L)
    {
    }
};

void
MDBuilderWindow::draw()
{
  using namespace mdtk;
  if (!valid())
  {
    valid(1);
    glLoadIdentity();
    glViewport(0,0,w(),h());
    glOrtho(-w(),w(),-h(),h(),-1,1);
  }
  glClear(GL_COLOR_BUFFER_BIT);
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINE_STRIP); 
  glVertex2f(w(), h()); 
  glVertex2f(-w(),-h()); 
  glEnd();
  glBegin(GL_LINE_STRIP); 
  glVertex2f(w(),-h()); 
  glVertex2f(-w(), h());
  glEnd();

  glLoadIdentity();
  {
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_H2(sl);

      std::ofstream fomde("H2.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_FCC_cell(sl);

      std::ofstream fomde("CuCell.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_FCC_lattice(sl);

      std::ofstream fomde("Cu.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_C60(sl);

      std::ofstream fomde("C60.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      mdtk::SimLoop sl;

      mdbuilder::build_C60_optimized(sl);

      std::ofstream fomde("C60-optimized.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
//    if (0)
    {
      mdtk::SimLoop sl;

      mdbuilder::build_cluster(sl,Cu_EL,13);

      std::ofstream fomde("Cu13.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
  }
  exit(0);
}

int main(int argc, char *argv[])
{
  if (argc > 1 && !strcmp(argv[1],"--version"))
  {
    std::cout << "mdbuilder (molecular dynamics experiments preparation tool) ";
    mdtk::release_info.print();
    return 0;
  }

  if (argc > 1 && (!std::strcmp(argv[1],"--help") || !std::strcmp(argv[1],"-h")))
  {
    std::cout << "\
Usage: mdbuilder [OPTION]... \n\
Prepares molecular dynamics experiments.\n\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
    return 0;
  }

  Fl_Window win(500, 300, "MDBuilder");
  MDBuilderWindow mygl(10, 10, win.w()-20, win.h()-20);
  win.end();
  win.resizable(mygl);
  win.show();
  return(Fl::run());
}
