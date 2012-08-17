/*
   mdbuilder (molecular dynamics experiments preparation tool)

   Copyright (C) 2010, 2011, 2012 Oleksandr Yermolenko
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

static int clusterSize = 0;
static size_t numberOfImpacts = 1024;

void
buildCommands()
{
  using namespace mdtk;
  glLoadIdentity();
  {
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_H2(sl.atoms);

      yaatk::text_ofstream fomde("H2.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_FCC_cell(sl.atoms);

      yaatk::text_ofstream fomde("CuCell.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_FCC_lattice(sl.atoms);

      yaatk::text_ofstream fomde("Cu.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_C60(sl.atoms);

      yaatk::text_ofstream fomde("C60.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      mdtk::SimLoop sl;
      sl.atoms = mdbuilder::C60();

      yaatk::text_ofstream fomde("C60-optimized.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      mdtk::SimLoop sl;
      sl.atoms = mdbuilder::cluster(Cu_EL,13);

      yaatk::text_ofstream fomde("Cu13.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      {
        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::clusterFromFCCCrystal(Au_EL,13);

        yaatk::text_ofstream fomde("Au13.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::clusterFromFCCCrystal(Au_EL,27);

        yaatk::text_ofstream fomde("Au27.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::clusterFromFCCCrystal(Au_EL,39);

        yaatk::text_ofstream fomde("Au39.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::clusterFromFCCCrystal(Au_EL,75);

        yaatk::text_ofstream fomde("Au75.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::clusterFromFCCCrystal(Au_EL,195);

        yaatk::text_ofstream fomde("Au195.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
    }
    if (0)
    {
      mdtk::AtomsArray C60 = mdbuilder::C60();

      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,1);

        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::embed(cluster,C60);

        yaatk::text_ofstream fomde("Cu001_in_C60.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,5);

        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::embed(cluster,C60);

        yaatk::text_ofstream fomde("Cu005_in_C60.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,13);

        mdtk::SimLoop sl;
        sl.atoms = mdbuilder::embed(cluster,C60);

        yaatk::text_ofstream fomde("Cu013_in_C60.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
    }
    if (0)
    {
      mdtk::AtomsArray C60 = mdbuilder::C60();

      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,0);

        mdtk::AtomsArray atoms = mdbuilder::embed(cluster,C60);
        mdbuilder::add_rotational_motion(atoms,200*eV,Vector3D(0,0,1));
        mdtk::SimLoop sl;
        sl.atoms = atoms;

        yaatk::text_ofstream fomde("Cu000_in_C60-rot.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,1);

        mdtk::AtomsArray atoms = mdbuilder::embed(cluster,C60);
        mdbuilder::add_rotational_motion(atoms,200*eV,Vector3D(0,0,1));
        mdtk::SimLoop sl;
        sl.atoms = atoms;

        yaatk::text_ofstream fomde("Cu001_in_C60-rot.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,6);

        mdtk::AtomsArray atoms = mdbuilder::embed(cluster,C60);
        mdbuilder::add_rotational_motion(atoms,200*eV,Vector3D(0,0,1));
        mdtk::SimLoop sl;
        sl.atoms = atoms;

        yaatk::text_ofstream fomde("Cu006_in_C60-rot.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
      {
        mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,13);

        mdtk::AtomsArray atoms = mdbuilder::embed(cluster,C60);
        mdbuilder::add_rotational_motion(atoms,200*eV,Vector3D(0,0,1));
        mdtk::SimLoop sl;
        sl.atoms = atoms;

        yaatk::text_ofstream fomde("Cu013_in_C60-rot.mde");
        sl.saveToMDE(fomde);
        fomde.close();
      }
    }
    if (0)
    {
      mdtk::SimLoop sl_Cu = mdbuilder::build_FCC_lattice(14,14,7,Cu_EL);

      mdtk::AtomsArray C60 = mdbuilder::C60();

      mdtk::AtomsArray cluster = mdbuilder::cluster(Cu_EL,6);

      mdtk::AtomsArray endo = mdbuilder::embed(cluster,C60);
      mdbuilder::add_rotational_motion(endo,50*eV,Vector3D(0,0,1));

      mdtk::SimLoop sl = mdbuilder::build_target_by_cluster_bombardment(sl_Cu,endo,200*eV);

      TRACE(sl.energyKin()/eV);

      sl.iteration = 0;
      sl.simTime = 0.0*ps;
      sl.simTimeFinal = 6.0*ps;
      sl.simTimeSaveTrajInterval = 0.05*ps;

      yaatk::text_ofstream fomde("Cu100_by_Cu06@C60_n200eV_z050eV.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      mdbuilder::prepare_Cu_by_Cu_at_C60_bobardment();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl;
      mdbuilder::place_Graphite_lattice(sl.atoms);

      yaatk::text_ofstream fomde("Graphite.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl = mdbuilder::build_Graphite_lattice();

      yaatk::text_ofstream fomde("Graphite_PBC.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      mdbuilder::prepare_Graphite_by_Cu_at_C60_bombardment();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl = mdbuilder::build_Polyethylene_lattice_without_folds();

      yaatk::text_ofstream fomde("Polyethylene_without_folds.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl = mdbuilder::build_Polyethylene_lattice_with_folds(4,6,10);

      yaatk::text_ofstream fomde("Polyethylene_with_folds.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
//    if (0)
    {
      TRACE(clusterSize);
      glLoadIdentity();
      std::vector<int> clusterSizes;
      clusterSizes.push_back(clusterSize);
/*
      clusterSizes.push_back(1);
      clusterSizes.push_back(13);
      clusterSizes.push_back(27);
      clusterSizes.push_back(39);
      clusterSizes.push_back(75);
      clusterSizes.push_back(195);
*/
      std::vector<ElementID> clusterElements;
      clusterElements.push_back(Cu_EL);
      clusterElements.push_back(Au_EL);
      std::vector<ElementID> ionElements;
      ionElements.push_back(Ar_EL);
      ionElements.push_back(Xe_EL);
      std::vector<Float> ionEnergies;
      ionEnergies.push_back(100*eV);
      ionEnergies.push_back(200*eV);
      ionEnergies.push_back(300*eV);
      ionEnergies.push_back(400*eV);
      ionEnergies.push_back(500*eV);
      mdbuilder::bomb_MetalCluster_on_Polyethylene_with_Ions(/*8,12,15*/6,9,15,
                                                             clusterSizes,
                                                             clusterElements,
                                                             ionElements,
                                                             ionEnergies,
                                                             numberOfImpacts);
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl = mdbuilder::build_Fullerite_C60(2,2,3);

      yaatk::text_ofstream fomde("Fullerite.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      mdtk::SimLoop sl = mdbuilder::build_Fulleride_C60(2,2,3,Cu_EL,3);

      yaatk::text_ofstream fomde("Fulleride.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();

      mdtk::SimLoop sl;
      mdbuilder::initialize_simloop(sl);

//    sl = mdbuilder::build_Fullerite_C60(2,2,3);
//    sl = mdbuilder::build_FCC_lattice(10,10,10,Cu_EL);
//    sl = mdbuilder::build_Graphite_lattice(12,14,3);
      sl = mdbuilder::build_Polyethylene_lattice_with_folds(4,6,10);

      mdbuilder::quench(sl,1.0*K,200*ps,0.01*ps,"_tmp-quench-to-1K");
      mdbuilder::heatUp(sl,300.0*K,false,200*ps,2.0*ps,"_tmp-heatup-to-300K");

      yaatk::text_ofstream fomde("300K.mde");
      sl.saveToMDE(fomde);
      fomde.close();
    }
    if (0)
    {
      glLoadIdentity();
      std::vector<ElementID> ionElements;
      ionElements.push_back(Ar_EL);
      ionElements.push_back(Xe_EL);
      std::vector<Float> ionEnergies;
      ionEnergies.push_back(100*eV);
//      ionEnergies.push_back(200*eV);
//      ionEnergies.push_back(300*eV);
      ionEnergies.push_back(400*eV);
      mdbuilder::build_FCC_metal_bombardment_with_ions(ionElements,
                                                       ionEnergies,
                                                       128);
    }
    if (0)
    {
      glLoadIdentity();
      {
        std::vector<Float> fullereneEnergies;
        fullereneEnergies.push_back(100*eV);
        fullereneEnergies.push_back(200*eV);
        fullereneEnergies.push_back(300*eV);
        fullereneEnergies.push_back(400*eV);
        mdbuilder::build_FCC_metal_bombardment_with_C60(fullereneEnergies,
                                                        128,
                                                        12,12,12);
      }
      glLoadIdentity();
      {
        std::vector<ElementID> ionElements;
        ionElements.push_back(Ar_EL);
        ionElements.push_back(Xe_EL);
        ionElements.push_back(Cu_EL);
        std::vector<Float> ionEnergies;
        ionEnergies.push_back(100*eV);
        ionEnergies.push_back(200*eV);
        ionEnergies.push_back(300*eV);
        ionEnergies.push_back(400*eV);
        mdbuilder::build_fullerite_bombardment_with_ions(ionElements,
                                                         ionEnergies,
                                                         128,
                                                         3,3,3);
      }
    }
    if (0)
    {
      glLoadIdentity();
      {
        std::vector<Float> impactEnergies;
        impactEnergies.push_back(25*eV);
        impactEnergies.push_back(50*eV);
        impactEnergies.push_back(100*eV);
        impactEnergies.push_back(200*eV);
        mdbuilder::build_metal_C60_mixing(impactEnergies);
      }
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
    printf("Depth=%d Stencil=%d Accum=%d\n", z, s, a);
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
    if (!strcmp(argv[argi],"--cluster-size") || !std::strcmp(argv[argi],"-s"))
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
    }

    if (!strcmp(argv[argi],"--number-of-impacts") || !std::strcmp(argv[argi],"-i"))
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
\n\
      --help     display this help and exit\n\
      --version  output version information and exit\n\
      --cluster-size  <cluster size>  generate experiment for a specified cluster size\n\
      --number-of-impacts  <number of impacts>  generate specific number of impacts for each experiment (default : 1024)\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
      return 0;
    }
  }

  if (clusterSize == 0)
  {
    std::cerr << "You should specify cluster size with --cluster-size option. Run with -h option for details.\n";
    return -1;
  }

  TRACE(clusterSize);
  TRACE(numberOfImpacts);

  srand(12345);

#ifdef MDBUILDER_OSMESA
  TRACE("Using OSMesa");
  buildWithOSMesa();
  return 0;
#else
  TRACE("Using FLTK OpenGL");
  Fl_Window win(500, 300, "MDBuilder");
  MDBuilderWindow mygl(10, 10, win.w()-20, win.h()-20);
  win.end();
  win.resizable(mygl);
  win.show();
  return(Fl::run());
#endif
}
