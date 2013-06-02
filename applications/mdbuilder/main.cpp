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

static bool prepareClusters = false;
static int clusterSize = 0;

static bool prepareSubstrate = false;
// Polyethylene target size
// static int a_num = 8;
// static int b_num = 12;
// static int c_num = 17;

// static int a_num = 12;
// static int b_num = 18;
// static int c_num = 35;

// static int a_num = 12;
// static int b_num = 18;
// static int c_num = 17;

// static int a_num = 16;
// static int b_num = 24;
// static int c_num = 17;

// static int a_num = 16;
// static int b_num = 24;
// static int c_num = 20;

// static int a_num = 16;
// static int b_num = 24;
// static int c_num = 21;

static int a_num = 15;
static int b_num = 22;
static int c_num = 15+2;

static bool landCluster = false;
static std::string clusterDataId;
static std::string substrateDataId;

static bool prepareBombardment = false;
static std::string targetDataId;

void
buildCommands()
{
  using namespace mdtk;
  glLoadIdentity();
  {
    if (prepareClusters)
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

        yaatk::text_ofstream fomde("mdloop.opti.best");
        sl.saveToStream(fomde);
        fomde.close();

        SimLoopSaver mds(sl);
        mds.write("most-optimal");
        mds.removeAttributesButPosVel("most-optimal");
      }
    }

    if (prepareSubstrate)
    {
      PRINT("Preparing a substrate (polyethylene).\n");

      glLoadIdentity();
      VerboseOutput vo(true);

      SimLoop sl_Polyethylene =
        mdbuilder::build_Polyethylene_lattice_with_folds(a_num,b_num,c_num);
      sl_Polyethylene.atoms.tag(ATOMTAG_SUBSTRATE);

      {
        SimLoopSaver mds(sl_Polyethylene);
        mds.write("substrate-000-pbc-box");
      }

//      Float unfixedDepth = (c_num-2-2)*2.547*Ao;
//      Float dAboveSurface = 10.0*Ao;
//      mdbuilder::setupSpherical(sl_Polyethylene,dAboveSurface,unfixedDepth + dAboveSurface);

//      Float dAboveSurface = 25.0*Ao;
      Float dAboveSurface = 25.0*Ao;
      mdbuilder::setupSpherical(sl_Polyethylene,dAboveSurface);

      {
        SimLoopSaver mds(sl_Polyethylene);
        mds.write("substrate-010-spherical-unrelaxed");
      }

      {
        yaatk::ChDir cd("_build_PE_spherical_relax");
        mdbuilder::initialize_simloop(sl_Polyethylene);
        mdbuilder::relax(sl_Polyethylene,10.0*ps,"010-relax");
        mdbuilder::quench(sl_Polyethylene,0.01*K,200*ps,0.01*ps,"011-quench");
        sl_Polyethylene.atoms.removeMomentum();
      }

      {
        SimLoopSaver mds(sl_Polyethylene);
        mds.write("substrate");
      }
    }

    if (landCluster)
    {
      PRINT("Landing cluster on the surface.\n");

      SimLoop sl_Polyethylene;
      {
        SimLoopSaver mds(sl_Polyethylene);
        mds.load(substrateDataId);
      }

      AtomsArray cluster;
      {
        SimLoop clusterMDLoop;
        if (!yaatk::exists((clusterDataId+".r").c_str()))
        {
          std::cerr << "Can not open cluster configuration file. Exiting." << std::endl;
          exit(1);
        }
        SimLoopSaver mds(clusterMDLoop);
        mds.load(clusterDataId);
        cluster = clusterMDLoop.atoms;
      }

      TRACE(cluster.size());
      REQUIRE(cluster.size() > 0);
      TRACE(ElementIDtoString(cluster[0].ID));

      glLoadIdentity();

      VerboseOutput vo(true);

      SimLoop sl_Landed =
        mdbuilder::build_Cluster_Landed_on_Substrate(sl_Polyethylene, cluster);

      {
        SimLoopSaver mds(sl_Landed);
        mds.write("target");
      }
    }

    if (prepareBombardment)
    {
      PRINT("Preparing MD experiment.\n");

      std::vector<ElementID> ionElements;
      ionElements.push_back(Ar_EL);
      ionElements.push_back(Xe_EL);
      std::vector<Float> ionEnergies;
      ionEnergies.push_back(100*eV);
      ionEnergies.push_back(200*eV);
      ionEnergies.push_back(300*eV);
      ionEnergies.push_back(400*eV);
      ionEnergies.push_back(500*eV);

      std::set<Float> halos;
      halos.insert(5.5*Ao);
      halos.insert(1.3*Ao);

      SimLoop target;

      SimLoopSaver mds(target);
      mds.load("target");

      mdbuilder::bomb_landedCluster_with_Ions(target,
                                              ionElements,
                                              ionEnergies,
                                              halos,
                                              4096);
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
    if (!strcmp(argv[argi],"--prepare-clusters"))
    {
      prepareClusters = true;
    }
    if (!strcmp(argv[argi],"--cluster-size"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify cluster size, e.g. --generate-clusters --cluster-size 13\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> clusterSize;
      if (!(clusterSize > 0 && clusterSize <= 200))
      {
        std::cerr << "Unsupported cluster size\n";
        return -1;
      }
    }

    if (!strcmp(argv[argi],"--prepare-substrate"))
    {
      prepareSubstrate = true;
    }

    if (!strcmp(argv[argi],"--land-cluster"))
    {
      landCluster = true;
    }
    if (!strcmp(argv[argi],"--cluster"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify existing cluster configuration file id, e.g. --land-cluster --cluster Cu13 --substrate PE\nIn this example current directory should contain files Cu13.z and Cu13.r\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> clusterDataId;
    }
    if (!strcmp(argv[argi],"--substrate"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify existing substrate configuration file id, e.g. --land-cluster --cluster Cu13 --substrate PE\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> substrateDataId;
    }

    if (!strcmp(argv[argi],"--prepare-bombardment"))
    {
      prepareBombardment = true;
    }
    if (!strcmp(argv[argi],"--target"))
    {
      if (!(argi+1 < argc))
      {
        std::cerr << "You should specify existing target configuration file id, e.g. --prepare-bombardment --target Cu13_landed_on_PE\n";
        return -1;
      }
      std::istringstream iss(argv[argi+1]);
      iss >> targetDataId;
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
To prepare Cu, Ag and Au clusters of the specific size run:\n\
  mdbuilder --prepare-clusters --cluster-size <cluster size>\n\
To prepare polyethylene substrate run:\n\
  mdbuilder --prepare-substrate\n\
To land a cluster on substrate run:\n\
  mdbuilder --land-cluster --cluster <id of files containing the cluster configuration> --substrate <id of files containing the substrate configuration>\n\
To prepare bombardment of the specific target containing a cluster run:\n\
  mdbuilder --prepare-bombardment --target <id of files containing the target>\n\
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
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
