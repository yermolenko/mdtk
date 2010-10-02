/*
   The MainWindow class for the molecular dynamics trajectory viewer.

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

#include "MainWindow.hpp"
#include <sstream>
#include <fstream>

#define MAX_LOG_BUFFER_LEN 100000

#include <mdtk/tools.hpp>

namespace xmde
{

  using namespace mdtk;

  int MainWindow::handle(int e)
  {
    switch(e)
      {
      default:
	return Fl_Window::handle(e);
      }
  } 

  void
  MainWindow::current_atomindex_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr = (MainWindow*)(w->parent()->parent()->parent());
    int new_index = int(((Fl_Counter *)w)->value());
    MainWindow_Ptr->setAtomViewIndex(new_index/*-1*/);
  } 

  void
  MainWindow::current_stateindex_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr = (MainWindow*)(w->parent()->parent()->parent());
    int new_index = int(((Fl_Counter *)w)->value());
  
    MainWindow_Ptr->loadNewState(new_index);
  } 

  void	MainWindow::out(std::string s)
  {
    const char * LogMsg = s.c_str();
	
    int dlen = strlen(LogMsg);
    if (log_pos+dlen < MAX_LOG_BUFFER_LEN-1)
      {
	strcat(log_buffer,LogMsg);
	log_pos += strlen(LogMsg);
	atom_info->value(log_buffer);
	atom_info->position(log_pos,log_pos);
	atom_info->redraw();
      }   	
  }
  
  void	MainWindow::clear_out()
  {
    log_pos = 0;
    strcpy(log_buffer,"");
    atom_info->value(log_buffer);
    atom_info->position(log_pos,log_pos);
    atom_info->redraw();
  }
  
  void	MainWindow::setAtomViewIndex(int index)
  {
    atoms_view_box->selectAtom(index);
    current_atomindex->value(index/*+1*/);

    clear_out();

    std::ostringstream os;
    os << "Atom " << index << ":" << std::endl; 
  
    mdtk::Atom &atom = *(atoms_view_box->getAtoms()->operator[](index));
  
    os TRACESS(int(atom.ID));
    os TRACESS(atom.Z/mdtk::e);
    os TRACESS(atom.M/mdtk::amu);
    os TRACESS(atom.V);
    os TRACESS(atom.coords/mdtk::Ao);
//    os TRACESS(atom.an);
//    os TRACESS(atom.apply_barrier);
    os TRACESS(atom.apply_PBC);
    os TRACESS(atom.apply_ThermalBath);
//    os TRACESS(atom.ejected);
    os TRACESS(atom.tag);
    os TRACESS(atom.globalIndex);
  
    out(os.str());  
  
  }

  MainWindow::MainWindow(std::string &bsf,std::vector<std::string>& fileList,
    VisBox* avb, bool instantAnimate):
    Fl_Window(700,495,"MDTK Trajectory Viewer"),atoms_view_box(avb)
  {
    base_state_filename = bsf;
    stateList = fileList;
    stateIndex = 0;
  
    log_buffer = new char[MAX_LOG_BUFFER_LEN];
    log_buffer[0]='\0';
    log_pos = 0;


    tabs = new Fl_Tabs(5,5,690,485);
    atoms_view_group = new Fl_Group(5,25,690,485,"Controls and Info");
    atoms_view_group->begin();

    current_stateindex	= new Fl_Counter(415,450,150,30,"Viewing state #");
    current_stateindex->align(FL_ALIGN_TOP);
    current_stateindex->lstep(10);
    current_stateindex->step(1);
    current_stateindex->precision(0);
    current_stateindex->bounds(0,stateList.size()-1);
    current_stateindex->value(0);
    current_stateindex->callback((Fl_Callback*)current_stateindex_cb);
 
    btn_rescale = new Fl_Light_Button(580,380,105,30,"Auto rescale");
    btn_rescale->callback((Fl_Callback*)btn_rescale_cb);
    btn_rescale->value(/*atoms_view_box->GetAllowRescale()*/false);

    btn_animate = new Fl_Light_Button(580,450,105,30,"Animate");
    btn_animate->callback((Fl_Callback*)btn_animate_cb);
    btn_animate->value(instantAnimate);

    btn_scale_up = new Fl_Button(580,320,105,30,"Scale Up");
    btn_scale_up->tooltip(btn_scale_up_tooltip);
    btn_scale_up->callback((Fl_Callback*)btn_scale_up_cb);

    btn_scale_down = new Fl_Button(580,350,105,30,"Scale Down");
    btn_scale_down->tooltip(btn_scale_down_tooltip);
    btn_scale_down->callback((Fl_Callback*)btn_scale_down_cb);

//    atoms_view_box = new VisBox(15,35,500,500,base_state_filename,fileList[stateIndex]);

    btn_unified_atoms = new Fl_Light_Button(160,320,165,30,
					    "Unicolored atoms");
    btn_unified_atoms->tooltip(btn_unified_atoms_tooltip);
    btn_unified_atoms->callback((Fl_Callback*)btn_unified_atoms_cb);
    btn_unified_atoms->value(!atoms_view_box->GetNativeVertexColorsState());

    btn_atoms_color = new Fl_Button(160,350,165,30,"Atoms color");
    btn_atoms_color->tooltip(btn_atoms_color_tooltip);
    btn_atoms_color->callback((Fl_Callback*)btn_atoms_color_cb);

    btn_bg_color = new Fl_Button(160,380,165,30,"Background color");
    btn_bg_color->tooltip(btn_bg_color_tooltip);
    btn_bg_color->callback((Fl_Callback*)btn_bg_color_cb);

    atom_quality = new Fl_Slider(350, 320, 25, 125, "Render\nquality");
    atom_quality->labelsize(12);
    atom_quality->minimum(16);
    atom_quality->maximum(3);
    atom_quality->value(14);
    atom_quality->step(1);
    atom_quality->type(FL_VERT_NICE_SLIDER);
    atom_quality->callback((Fl_Callback*)atom_quality_cb);

    animate_delay = new Fl_Counter(620, 422, 65, 26, "Delay");
    animate_delay->labelsize(12);
    animate_delay->minimum(0.1);
    animate_delay->maximum(5.0);
    animate_delay->value(0.1);
    animate_delay->step(0.1);
    animate_delay->align(FL_ALIGN_LEFT);
    animate_delay->type(FL_SIMPLE_COUNTER);

    btn_show_axes = new Fl_Light_Button(160,420,165,30,"Show axes");
    btn_show_axes->tooltip(btn_show_axes_tooltip);
    btn_show_axes->callback((Fl_Callback*)btn_show_axes_cb);
    btn_show_axes->value(atoms_view_box->GetAxesState());

    btn_show_bath = new Fl_Light_Button(160,450,165,30,"Show thermal bath");
    btn_show_bath->tooltip(btn_show_bath_tooltip);
    btn_show_bath->callback((Fl_Callback*)btn_show_bath_cb);
    btn_show_bath->value(atoms_view_box->GetBathState());

/*
    btn_fixed_lights = new Fl_Light_Button(20,450,125,30,"Fixed Lights");
    btn_fixed_lights->tooltip(btn_fixed_lights_tooltip);
    btn_fixed_lights->callback((Fl_Callback*)btn_fixed_lights_cb);
    btn_fixed_lights->value(1);
*/

    Fl_Box* lightDirBox = new Fl_Box(FL_UP_FRAME,15,320,135,160,"Light Direction");
    lightDirBox->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

    light_x_dir	= new Fl_Slider(30, 345, 25, 115,"X");
    light_x_dir->labelsize(12);
    light_x_dir->minimum(-1.0);
    light_x_dir->maximum(+1.0);
    light_x_dir->value(1.0);
    light_x_dir->step(0.1);
    light_x_dir->type(FL_VERT_NICE_SLIDER);
    light_x_dir->callback((Fl_Callback*)light_x_dir_cb);

    light_y_dir	= new Fl_Slider(70, 345, 25, 115,"Y");
    light_y_dir->labelsize(12);
    light_y_dir->minimum(-1.0);
    light_y_dir->maximum(+1.0);
    light_y_dir->value(1.0);
    light_y_dir->step(0.1);
    light_y_dir->type(FL_VERT_NICE_SLIDER);
    light_y_dir->callback((Fl_Callback*)light_y_dir_cb);

    light_z_dir	= new Fl_Slider(110, 345, 25, 115,"Z");
    light_z_dir->labelsize(12);
    light_z_dir->minimum(-1.0);
    light_z_dir->maximum(+1.0);
    light_z_dir->value(1.0);
    light_z_dir->step(0.1);
    light_z_dir->type(FL_VERT_NICE_SLIDER);
    light_z_dir->callback((Fl_Callback*)light_z_dir_cb);

    btn_quick_save_image = new Fl_Button(415,320,125,30,"Save video");
    btn_quick_save_image->callback((Fl_Callback*)btn_quick_save_image_cb);

    btn_save_image = new Fl_Button(415,350,125,30,"Save image");
    btn_save_image->tooltip(btn_save_image_tooltip);
    btn_save_image->callback((Fl_Callback*)btn_save_image_cb);

    btn_save_mde = new Fl_Button(415,380,125,30,"Save in.mde.gz");
    btn_save_mde->callback((Fl_Callback*)btn_save_mde_cb);

    Fl_Box* rollBox = new Fl_Box(FL_EMBOSSED_FRAME,415,35,270,280,"Rotate and clip");
    rollBox->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);

    roll_x = new Fl_Roller(455, 60, 20, 235, "X");
    roll_x->tooltip(roll_x_tooltip);
    roll_x->type(0);
    roll_x->labelsize(12);
    roll_x->minimum(-360);
    roll_x->maximum(360);
    roll_x->value(0);
    roll_x->step(1);
    roll_x->callback((Fl_Callback*)roll_x_cb);

    val_xmin = new Fl_Slider(430, 60, 25, 235, "min");
    val_xmin->type(FL_VERT_NICE_SLIDER);
    val_xmin->labelsize(12);
    val_xmin->minimum(101);
    val_xmin->maximum(-1);
    val_xmin->value(0);
    val_xmin->step(0.1);
    val_xmin->when(FL_WHEN_RELEASE);
    val_xmin->callback((Fl_Callback*)val_xmin_cb);

    val_xmax = new Fl_Slider(475, 60, 25, 235, "max");
    val_xmax->type(FL_VERT_NICE_SLIDER);
    val_xmax->labelsize(12);
    val_xmax->minimum(101);
    val_xmax->maximum(-1);
    val_xmax->value(100);
    val_xmax->step(0.1);
    val_xmax->when(FL_WHEN_RELEASE);
    val_xmax->callback((Fl_Callback*)val_xmin_cb);

    roll_y = new Fl_Roller(540, 60, 20, 235, "Y");
    roll_y->tooltip(roll_y_tooltip);
    roll_y->type(0);
    roll_y->labelsize(12);
    roll_y->minimum(-360);
    roll_y->maximum(360);
    roll_y->value(0);
    roll_y->step(1);
    roll_y->callback((Fl_Callback*)roll_y_cb);

    val_ymin = new Fl_Slider(515, 60, 25, 235, "min");
    val_ymin->type(FL_VERT_NICE_SLIDER);
    val_ymin->labelsize(12);
    val_ymin->minimum(101);
    val_ymin->maximum(-1);
    val_ymin->value(0);
    val_ymin->step(0.1);
    val_ymin->when(FL_WHEN_RELEASE);
    val_ymin->callback((Fl_Callback*)val_xmin_cb);

    val_ymax = new Fl_Slider(560, 60, 25, 235, "max");
    val_ymax->type(FL_VERT_NICE_SLIDER);
    val_ymax->labelsize(12);
    val_ymax->minimum(101);
    val_ymax->maximum(-1);
    val_ymax->value(100);
    val_ymax->step(0.1);
    val_ymax->when(FL_WHEN_RELEASE);
    val_ymax->callback((Fl_Callback*)val_xmin_cb);

    roll_z = new Fl_Roller(625, 60, 20, 235, "Z");
    roll_z->tooltip(roll_z_tooltip);
    roll_z->type(0);
    roll_z->labelsize(12);
    roll_z->minimum(-360);
    roll_z->maximum(360);
    roll_z->value(0);
    roll_z->step(1);
    roll_z->callback((Fl_Callback*)roll_z_cb);

    val_zmin = new Fl_Slider(600, 60, 25, 235, "min");
    val_zmin->type(FL_VERT_NICE_SLIDER);
    val_zmin->labelsize(12);
    val_zmin->minimum(101);
    val_zmin->maximum(-1);
    val_zmin->value(0);
    val_zmin->step(0.1);
    val_zmin->when(FL_WHEN_RELEASE);
    val_zmin->callback((Fl_Callback*)val_xmin_cb);

    val_zmax = new Fl_Slider(645, 60, 25, 235, "max");
    val_zmax->type(FL_VERT_NICE_SLIDER);
    val_zmax->labelsize(12);
    val_zmax->minimum(101);
    val_zmax->maximum(-1);
    val_zmax->value(100);
    val_zmax->step(0.1);
    val_zmax->when(FL_WHEN_RELEASE);
    val_zmax->callback((Fl_Callback*)val_xmin_cb);

    new Fl_Box(FL_UP_FRAME,15,35,385,235+45,NULL);

    atom_info = new Fl_Multiline_Output(25,55,365,215,"Selected Atom info");
    atom_info->textcolor(FL_BLUE);
    atom_info->align(FL_ALIGN_TOP);
    atom_info->value(""); 

    current_atomindex	= new Fl_Counter(25,275,220,30,/*"Selected atom index : "*/NULL);
    current_atomindex->align(FL_ALIGN_LEFT);
    current_atomindex->lstep(100);
    current_atomindex->step(1);
    current_atomindex->precision(0);
    current_atomindex->bounds(0,atoms_view_box->getAtomsCount()-1);
    current_atomindex->value(0);
    current_atomindex->callback((Fl_Callback*)current_atomindex_cb);

    btn_show_selected = new Fl_Light_Button(255,275,130,30,"Show selected");
    btn_show_selected->callback((Fl_Callback*)btn_show_selected_cb);
    btn_show_selected->value(atoms_view_box->GetShowSelectedState());

    atoms_view_group->end();


    about_group = new Fl_Group(5,25,690,485,"About");
    about_group->begin();
    product_info = 
      std::string("Molecular dynamics trajectory viewer")+"\nusing "+
      "MDTK" + " " +
      mdtk::release_info.PRODUCT_VERSION + 
      "\n\n" + 
      "Copyright (C) 2003-2010 Oleksandr Yermolenko\n <oleksandr.yermolenko@@gmail.com>\n\n" +
      "Run the program with --version or --help options for details."
      ;
    about_box = new Fl_Box(5,25,690,485,product_info.c_str());
//    about_box->box(FL_EMBOSSED_BOX);
    about_box->labelfont(FL_HELVETICA_BOLD);
    about_group->end();

//    info_group = new Fl_Group(5,25,765,540,"Info");
//    info_group->begin();
//    info_group->end();

    tabs->end();
    end();

    atoms_view_box->show();
/*
    glw = new Fl_Gl_Window(15,35,500,500,"AAA");
    glw->show();
*/
    show();
    atoms_view_box->redraw();
    redraw();
    atoms_view_box->invalidate();
    atoms_view_box->hide();
    atoms_view_box->show();
	
//    initEnergyChart();
    loadNewState(0);

    atoms_view_box->SetAllowRescale(false);

    /*
      roll_x->value(-90.0);
      roll_x_cb(roll_x, NULL);
    */
    callback(window_cb);
    if (btn_animate->value()) Fl::add_timeout(1.0, timer_callback);
  }

  MainWindow::~MainWindow()
  {
  }

  void	MainWindow::redrawGL()
  {
    atoms_view_box->redraw();
  }

  void MainWindow::btn_save_image_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    VisBox* VisBox_Ptr;
    VisBox_Ptr = MainWindow_Ptr->atoms_view_box;

    char *tmp_filename = fl_file_chooser
      (
       "Choose a file to save image to ...",
       "Windows Bitmap Files (*.bmp)",
       0,0
       );
    if (tmp_filename && (!fl_filename_isdir(tmp_filename)))
      {
	VisBox_Ptr->QSaveImageToFile(tmp_filename);

	VisBox_Ptr->make_current();

	fl_message("It seems that file was successfully saved !");
      }

  }


  void MainWindow::quickSaveBitmap()
  {
    char tmp_filename[4096];
	
    sprintf(tmp_filename,"%010d.bmp",/*MainWindow_Ptr->*/stateIndex);
	
    if (!fl_filename_isdir(tmp_filename))
      {
	atoms_view_box->QSaveImageToFile(tmp_filename);
      }
  }


  void MainWindow::btn_quick_save_image_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    VisBox* VisBox_Ptr;
    VisBox_Ptr = MainWindow_Ptr->atoms_view_box;

    yaatk::mkdir("video");
    for(size_t i = 0; i < MainWindow_Ptr->stateList.size()/*-1*/; i++)
      {
//	MainWindow_Ptr->loadNewState(i);
        MainWindow_Ptr->current_stateindex->value(i);
        current_stateindex_cb(MainWindow_Ptr->current_stateindex,NULL);
        yaatk::chdir("video");
        while (!Fl::ready()) {};
	Fl::flush();
        while (!Fl::ready()) {};
        MainWindow_Ptr->quickSaveBitmap();
//        sleep(1);
        yaatk::chdir("..");
      }
    fl_message("Sequence of images saved in the \"video\" directory!");
  }

  void MainWindow::btn_save_mde_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    VisBox* VisBox_Ptr;
    VisBox_Ptr = MainWindow_Ptr->atoms_view_box;

    char *tmp_filename = fl_file_chooser
      (
       "Choose a file to save MDE in ...",
       "mde (*.mde)",
       0,0
       );
    if (tmp_filename && (!fl_filename_isdir(tmp_filename)))
      {
	VisBox_Ptr->SaveToMDE(tmp_filename);
      }  
  }



  void MainWindow::btn_scale_up_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());
    VisBox* VisBox_Ptr
      = MainWindow_Ptr->atoms_view_box;

    if (1.2*VisBox_Ptr->GetScale()<=VisBox_Ptr->GetMaxScale())
      VisBox_Ptr->SetScale(1.2*VisBox_Ptr->GetScale());
  }

  void MainWindow::btn_scale_down_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());
    VisBox* VisBox_Ptr = MainWindow_Ptr->atoms_view_box;

    VisBox_Ptr->SetScale(VisBox_Ptr->GetScale()/1.2);
  }

  void MainWindow::light_x_dir_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double x_dir;
    x_dir = ((Fl_Slider *)w)->value();

    MainWindow_Ptr->atoms_view_box->
      SetLightXDir(x_dir);
  }

  void MainWindow::light_y_dir_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double y_dir;
    y_dir = ((Fl_Slider *)w)->value();

    MainWindow_Ptr->atoms_view_box->
      SetLightYDir(y_dir);
  }

  void MainWindow::light_z_dir_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double z_dir;
    z_dir = ((Fl_Slider *)w)->value();

    MainWindow_Ptr->atoms_view_box->
      SetLightZDir(z_dir);
  }

  void MainWindow::roll_x_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double new_rot_x,old_rot_x;
    new_rot_x = ((Fl_Roller *)w)->value();
    old_rot_x = MainWindow_Ptr->atoms_view_box->Get_old_rot_x();

    MainWindow_Ptr->atoms_view_box->
      RollAround(new_rot_x - old_rot_x,1,0,0);
  }

  void MainWindow::roll_y_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double new_rot_y,old_rot_y;
    new_rot_y = ((Fl_Roller *)w)->value();
    old_rot_y = MainWindow_Ptr->atoms_view_box->Get_old_rot_y();

    MainWindow_Ptr->atoms_view_box->
      RollAround(new_rot_y - old_rot_y,0,1,0);
  }

  void MainWindow::roll_z_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    double new_rot_z,old_rot_z;
    new_rot_z = ((Fl_Roller *)w)->value();
    old_rot_z = MainWindow_Ptr->atoms_view_box->Get_old_rot_z();

    MainWindow_Ptr->atoms_view_box->
      RollAround(new_rot_z - old_rot_z,0,0,1);
  }

  void MainWindow::val_xmin_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    MainWindow_Ptr->atoms_view_box->
      ReArrange(
		MainWindow_Ptr->val_xmin->value()/100.0,
		MainWindow_Ptr->val_xmax->value()/100.0,
		MainWindow_Ptr->val_ymin->value()/100.0,
		MainWindow_Ptr->val_ymax->value()/100.0,
		MainWindow_Ptr->val_zmin->value()/100.0,
		MainWindow_Ptr->val_zmax->value()/100.0
		);
  }

  void MainWindow::atom_quality_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    MainWindow_Ptr->atoms_view_box->SetAtomsQuality(
                     int(MainWindow_Ptr->atom_quality->value()));
  }

  void MainWindow::btn_fixed_lights_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetFixedLightsState(v);
  }

  void MainWindow::btn_show_axes_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetAxesState(v);
    MainWindow_Ptr->atoms_view_box->SetBarrierState(v);
  }

  void MainWindow::btn_show_bath_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetBathState(v);
  }

  void MainWindow::btn_show_selected_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetShowSelectedState(v);
  }
  void MainWindow::btn_rescale_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetAllowRescale(v);
  }
  void MainWindow::btn_animate_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    
//    double delay = MainWindow_Ptr->animate_delay->value();
    
    if (v) Fl::repeat_timeout(0.5, timer_callback);

//    MainWindow_Ptr->atoms_view_box->SetAllowRescale(v);
  }

  void
MainWindow::timer_callback(void*)
{
  puts("TICK");

  if (MainWindow_GlobalPtr && MainWindow_GlobalPtr->btn_animate->value())
  {
    xmde::MainWindow* MainWindow_Ptr;
    MainWindow_Ptr = MainWindow_GlobalPtr;

    xmde::VisBox* VisBox_Ptr;
    VisBox_Ptr = MainWindow_Ptr->atoms_view_box;

    int oldval = MainWindow_Ptr->current_stateindex->value();
    int maxval = MainWindow_Ptr->current_stateindex->maximum();
    MainWindow_Ptr->current_stateindex->value((oldval+1)%(maxval+1));
    current_stateindex_cb(MainWindow_Ptr->current_stateindex,NULL);

    double delay = MainWindow_Ptr->animate_delay->value();

    if (MainWindow_Ptr->btn_animate->value()) Fl::repeat_timeout(delay, timer_callback);
  }
}

  void MainWindow::btn_unified_atoms_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());

    bool v = ((Fl_Light_Button *)w)->value();
    MainWindow_Ptr->atoms_view_box->SetNativeVertexColorsState(!v);
  }

  void MainWindow::btn_bg_color_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());
    Color c = MainWindow_Ptr->atoms_view_box->GetBGColor();

    unsigned char r,g,b;
    AnalyseRGB(c,r,g,b);

    if (fl_color_chooser("Choose Background Color",r,g,b))
      {
	c = CombineRGB(r,g,b);
	MainWindow_Ptr->atoms_view_box->SetBGColor(c);
      }
  }

  void MainWindow::btn_atoms_color_cb(Fl_Widget *w, void *)
  {
    MainWindow* MainWindow_Ptr;
    MainWindow_Ptr =
      (MainWindow*)(w->parent()->parent()->parent());
    Color c = MainWindow_Ptr->atoms_view_box->GetVertexColor();

    unsigned char r,g,b;
    AnalyseRGB(c,r,g,b);

    if (fl_color_chooser("Choose Vertex Color",r,g,b))
      {
	c = CombineRGB(r,g,b);
	MainWindow_Ptr->atoms_view_box->SetVertexColor(c);
      }
  }

  void MainWindow::window_cb(Fl_Widget* widget, void*)
  {
//    if (fl_choice("Do you really want to exit?","No","Yes",NULL)==1)
    {
      ((Fl_Window*)widget)->hide();
      exit(0);
    }
  }

  const char* MainWindow::
  btn_fixed_lights_tooltip = "Toggle Fixed/Free Lights";
  const char* MainWindow::
  btn_show_axes_tooltip = "Show/Hide Axes";
  const char* MainWindow::
  btn_show_bath_tooltip = "Show/Hide Thermal Bath";
  const char* MainWindow::
  btn_unified_atoms_tooltip = "Toggle Colored/Unified Atoms";
  const char* MainWindow::
  roll_x_tooltip = "Rolls around X-axis";
  const char* MainWindow::
  roll_y_tooltip = "Rolls around Y-axis";
  const char* MainWindow::
  roll_z_tooltip = "Rolls around Z-axis";

  const char* MainWindow::
  btn_bg_color_tooltip = "Choose Background Color";
  const char* MainWindow::
  btn_atoms_color_tooltip = "Choose Unified Atoms Color";


  const char* MainWindow::
  btn_save_image_tooltip = "Save image to file (*.bmp)";

  const char* MainWindow::
  btn_scale_up_tooltip = "Scale = Scale*1.2";
  const char* MainWindow::
  btn_scale_down_tooltip = "Scale = Scale/1.2";




  void  MainWindow::loadNewState(int index)
  {
    //  if (index > 0) atoms_view_box->SetAllowRescale(false);
    stateIndex = index;
    label((std::string("MDTK Trajectory Viewer [Control Window] - ")+stateList[stateIndex]).c_str());
    atoms_view_box -> label((std::string("MDTK 3D View - ")+stateList[stateIndex]).c_str());
    atoms_view_box -> loadNewML(base_state_filename,stateList[stateIndex]);
    setAtomViewIndex(int(current_atomindex->value())/*-1*/);
  }  
/*
  void MainWindow::initEnergyChart()
  {
    mdtk::SimLoop *ml;
    ml = new mdtk::SimLoop;
    ml->initNLafterLoading = false;
    for(size_t i = 0; i < stateList.size();i++)
      {
	using mdtk::Exception;

	YAATK_IFSTREAM_CREATE_ZIPPED(std::ifstream,fi,base_state_filename.c_str()); 

	ml->loadFromStream(fi);
	{
	  YAATK_IFSTREAM_CREATE_ZIPPED(std::ifstream,fixva,stateList[i].c_str()); 
	  ml->loadFromStreamXVA(fixva);
	  YAATK_IFSTREAM_CLOSE_ZIPPED(fixva,stateList[i].c_str()); 
	}  

	YAATK_IFSTREAM_CLOSE_ZIPPED(fi,base_state_filename.c_str()); 
 	  
	ev.push_back(ml->getEnergyCur()/mdtk::eV);  	  
 	  
      }  
    delete ml;
  }  
*/
}

xmde::MainWindow* MainWindow_GlobalPtr = NULL;

