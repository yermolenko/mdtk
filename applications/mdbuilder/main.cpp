/*
   mdbuilder (molecular dynamics experiments preparation tool)

   Copyright (C) 2010, 2011, 2012, 2013 Oleksandr Yermolenko
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

#ifdef MDBUILDER_OSMESA
#include "GL/osmesa.h"
#include "GL/glu.h"
#else
#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#endif

#include <cstring>

//#include "pentacene.hpp"

#include <iostream>
#include <fstream>
#include <sstream>

#include "applications/common.h"
#include "experiments/H2.hpp"
#include "experiments/FCC.hpp"
#include "experiments/Graphite.hpp"
#include "experiments/Polyethylene.hpp"
#include "experiments/Clusters.hpp"
#include "experiments/Fullerite.hpp"
#include "experiments/Fulleride.hpp"

static bool generateClusters = false;
static int clusterSize = 0;

static std::string clusterConfigurationFile;
static size_t numberOfImpacts = 1024;

void
buildCommands()
{
  using namespace mdtk;
  glLoadIdentity();
  {
    if (generateClusters) // cluster generation was requested
    {
      PRINT("Generating clusters.\n");
      TRACE(clusterSize);

      std::vector<ElementID> ids;
      ids.push_back(Cu_EL);
      ids.push_back(Au_EL);
      ids.push_back(Ag_EL);

      std::vector<size_t> sizes;

      REQUIRE(1 < clusterSize && clusterSize < 2000);
      sizes.push_back(clusterSize);

      for(size_t idi = 0; idi < ids.size(); ++idi)
      for(size_t si = 0; si < sizes.size(); ++si)
      {
        ElementID id = ids[idi];
        int clusterSize = sizes[si];

        std::ostringstream dirname;
        dirname << ElementIDtoString(id);
        PRINT2STREAM_FW(dirname, clusterSize, '0', 6);
        yaatk::ChDir cd(dirname.str());

        AtomsArray cluster = mdbuilder::clusterFromFCCCrystal(id,clusterSize);
        SimLoop sl;
        mdbuilder::initialize_simloop(sl);
        sl.atoms = cluster;

        sl.iteration = 0;
        sl.simTime = 0.0*ps;
        sl.simTimeFinal = 10.0*ps;
        sl.simTimeSaveTrajInterval = 0.1*ps;

        yaatk::text_ofstream fomde("mde_init");
        sl.saveToStream(fomde);
        fomde.close();
      }
    }
    else // bombardment experiment preparation was requested
    {
      PRINT("Preparing MD experiment.\n");

      AtomsArray cluster;
      SimLoop clusterMDLoop;
      if (!yaatk::exists(clusterConfigurationFile.c_str()))
      {
        std::cerr << "Can not open cluster configuration file. Exiting." << std::endl;
        exit(1);
      }
      yaatk::text_ifstream fi(clusterConfigurationFile.c_str());
      clusterMDLoop.loadFromStream(fi);
      fi.close();
      cluster = clusterMDLoop.atoms;

      TRACE(cluster.size());
      REQUIRE(cluster.size() > 0);
      TRACE(ElementIDtoString(cluster[0].ID));
      TRACE(numberOfImpacts);

      glLoadIdentity();
      std::vector<ElementID> ionElements;
      ionElements.push_back(Ar_EL);
      ionElements.push_back(Xe_EL);
      std::vector<Float> ionEnergies;
      ionEnergies.push_back(100*eV);
      ionEnergies.push_back(200*eV);
      ionEnergies.push_back(300*eV);
      ionEnergies.push_back(400*eV);
      ionEnergies.push_back(500*eV);
      mdbuilder::bomb_MetalCluster_on_Polyethylene_with_Ions(8,12,17,
                                                             cluster,
                                                             ionElements,
                                                             ionEnergies,
                                                             numberOfImpacts);
    }
  }
  exit(0);
}

#ifdef MDBUILDER_OSMESA

int
buildWithOSMesa()
{
  int w = 500;
  int h = 500;

  OSMesaContext ctx;
  void *buffer;

  buffer = malloc(w*h*4*sizeof(GLubyte));
  if (!buffer)
  {
    printf("Alloc image buffer failed!\n");
    return -1;
  }

#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
  /* specify Z, stencil, accum sizes */
  ctx = OSMesaCreateContextExt(OSMESA_RGBA,16,0,0,NULL);
#else
  ctx = OSMesaCreateContext(OSMESA_RGBA,NULL);
#endif
  if (!ctx)
  {
    printf("OSMesaCreateContext failed!\n");
    return -1;
  }

  if (!OSMesaMakeCurrent(ctx,buffer,GL_UNSIGNED_BYTE,w,h))
  {
    printf("OSMesaMakeCurrent failed!\n");
    return -1;
  }

  {
    int z, s, a;
    glGetIntegerv(GL_DEPTH_BITS, &z);
    glGetIntegerv(GL_STENCIL_BITS, &s);
    glGetIntegerv(GL_ACCUM_RED_BITS, &a);
    if (z != 16 || s != 0 || a != 0)
    {
      fprintf(stderr,"Something is wrong with OSMesa setup\n");
      fprintf(stderr,"Depth=%d Stencil=%d Accum=%d\n", z, s, a);
      exit(1);
    }
  }


  using namespace mdtk;
  {
    glLoadIdentity();
    glViewport(0,0,w,h);
    glOrtho(-w,w,-h,h,-1,1);
  }
  glClear(GL_COLOR_BUFFER_BIT);
  glColor3f(1.0, 1.0, 1.0);
  glBegin(GL_LINE_STRIP);
  glVertex2f(w, h);
  glVertex2f(-w,-h);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glVertex2f(w,-h);
  glVertex2f(-w, h);
  glEnd();

  buildCommands();

  glFinish();

  OSMesaDestroyContext( ctx );

  return 0;
}

#else // MDBUILDER_OSMESA

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

  buildCommands();
}

#endif // MDBUILDER_OSMESA

int main(int argc, char *argv[])
{
  for(int argi = 0; argi < argc; ++argi)
  {
    if (!strcmp(argv[argi],"--cluster-size"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify cluster size, e.g. --cluster-size 13\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> clusterSize;
      if (!(clusterSize > 0 && clusterSize <= 200))
      {
        std::cerr << "Unsupported cluster size\n";
        return -1;
      }
      else
      {
        generateClusters = true;
      }
    }
    if (!strcmp(argv[argi],"--cluster"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify existing cluster configuration file, e.g. --cluster-size Cu13.mdloop\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> clusterConfigurationFile;
      generateClusters = false;
    }

    if (!strcmp(argv[argi],"--number-of-impacts"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify number of impacts, e.g. --number-of-impacts 300\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> numberOfImpacts;
      if (!(numberOfImpacts > 0 && numberOfImpacts <= 2048))
      {
        std::cerr << "Unsupported number of impacts\n";
        return -1;
      }
    }

    if (!strcmp(argv[argi],"--version"))
    {
      std::cout << "mdbuilder (molecular dynamics experiments preparation tool) ";
      mdtk::release_info.print();
      return 0;
    }

    if (!std::strcmp(argv[argi],"--help") || !std::strcmp(argv[argi],"-h"))
    {
      std::cout << "\
Usage: mdbuilder [OPTION]... \n\
Prepares molecular dynamics experiments.\n\
To prepare clusters of the specific size run:\n\
  mdbuilder --cluster-size  <cluster size>\n\
To prepare an experiment for the specific cluster run:\n\
  mdbuilder --cluster <mdloop file containing cluster configuration> --number-of-impacts  <number of impacts>\n\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
      --cluster-size  <cluster size>  generate Cu, Au and Ag clusters of the specific size\n\
      --cluster  <mdloop file containing cluster configuration>  prepare an experiment for the specific cluster\n\
      --number-of-impacts  <number of impacts>  generate specific number of impacts for each experiment (default : 1024)\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
      return 0;
    }
  }

  srand(12345);

#ifdef MDBUILDER_OSMESA
  PRINT("Using OSMesa-based backend.\n");
  buildWithOSMesa();
  return 0;
#else
  PRINT("Using FLTK OpenGL-based backend.\n");
  Fl_Window win(500, 300, "MDBuilder");
  MDBuilderWindow mygl(10, 10, win.w()-20, win.h()-20);
  win.end();
  win.resizable(mygl);
  win.show();
  return(Fl::run());
#endif
}
