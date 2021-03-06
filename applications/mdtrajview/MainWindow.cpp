/*
   The MainWindow class for the molecular dynamics trajectory viewer.

   Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
   2012, 2013, 2015 Oleksandr Yermolenko
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

#include "MainWindow.hpp"
#include <sstream>
#include <fstream>
#include <ctime>

#define MAX_LOG_BUFFER_LEN 100000

#include <mdtk/tools.hpp>

namespace xmde
{

using namespace mdtk;

int
MainWindow::handle(int e)
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
MainWindow::set_atom_properties_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  mdtk::Atom* a = MainWindow_Ptr->renderBox->getSelectedAtomPtr();

  a->coords.x = MainWindow_Ptr->atom_coords_x->value()*Ao;
  a->coords.y = MainWindow_Ptr->atom_coords_y->value()*Ao;
  a->coords.z = MainWindow_Ptr->atom_coords_z->value()*Ao;

  a->V.x = MainWindow_Ptr->atom_v_x->value();
  a->V.y = MainWindow_Ptr->atom_v_y->value();
  a->V.z = MainWindow_Ptr->atom_v_z->value();

  MainWindow_Ptr->renderBox->saveSelectedAtomProperies();

  MainWindow_Ptr->renderBox->reArrange(-1,101,-1,101,-1,101);
/*
  MainWindow_Ptr->renderBox->
    reArrange(
      MainWindow_Ptr->val_xmin->value()/100.0,
      MainWindow_Ptr->val_xmax->value()/100.0,
      MainWindow_Ptr->val_ymin->value()/100.0,
      MainWindow_Ptr->val_ymax->value()/100.0,
      MainWindow_Ptr->val_zmin->value()/100.0,
      MainWindow_Ptr->val_zmax->value()/100.0
      );
*/
  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::current_stateindex_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr = (MainWindow*)(w->parent()->parent()->parent());
  int new_index = int(((Fl_Counter *)w)->value());
  
  MainWindow_Ptr->loadNewSnapshot(new_index);
} 

void
MainWindow::out(std::string s)
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
  
void
MainWindow::clear_out()
{
  log_pos = 0;
  strcpy(log_buffer,"");
  atom_info->value(log_buffer);
  atom_info->position(log_pos,log_pos);
  atom_info->redraw();
}
  
void
MainWindow::setAtomViewIndex(int index)
{
  renderBox->selectedAtomIndex = index;
  current_atomindex->value(index/*+1*/);

  mdtk::Atom* a = renderBox->getSelectedAtomPtr();

  atom_coords_x->value(a->coords.x/Ao);
  atom_coords_y->value(a->coords.y/Ao);
  atom_coords_z->value(a->coords.z/Ao);

  atom_v_x->value(a->V.x);
  atom_v_y->value(a->V.y);
  atom_v_z->value(a->V.z);

  clear_out();

  std::ostringstream os;
  os << "Atom " << index << ":" << std::endl; 
  
  mdtk::Atom &atom = renderBox->getAtoms()->operator[](index);

  using namespace mdtk;
  
  os TRACESS_NO_ENDL(int(atom.ID),"; ");
  os TRACESS_NO_ENDL(atom.Z/e,"; ");
  os TRACESS(atom.M/amu);
  os TRACESS(atom.V);
  os TRACESS(atom.coords/Ao);
  os TRACESS(atom.PBC_count);
//    os TRACESS(atom.an);
//    os TRACESS(atom.apply_barrier);
  os TRACESS(atom.PBCEnabled());
  os TRACESS(atom.lateralPBCEnabled());
  os TRACESS(atom.apply_ThermalBath);
//    os TRACESS(atom.ejected);
  os TRACESS(atom.tagbits);
  os TRACESS(atom.globalIndex);
  
  out(os.str());  

  renderBox->redraw();
}

MainWindow::MainWindow(VisBox* avb, bool instantAnimate):
  Fl_Window(1000,550,"MDTK Trajectory Viewer"),
  renderBox(avb),
  stateList(),
  stateIndex(0)
{
  log_buffer = new char[MAX_LOG_BUFFER_LEN];
  log_buffer[0]='\0';
  log_pos = 0;

  Fl_Tabs*
    tabs = new Fl_Tabs(5,5,990,540);
  Fl_Group*
    atoms_view_group = new Fl_Group(5,25,990,540,
				    "Controls and Info");
  atoms_view_group->begin();

  current_stateindex = new Fl_Counter(415,450,150,30,
				      "Viewing state #");
  current_stateindex->align(FL_ALIGN_TOP);
  current_stateindex->lstep(10);
  current_stateindex->step(1);
  current_stateindex->precision(0);
//  current_stateindex->bounds(0,stateList.size()-1);
  current_stateindex->value(0);
  current_stateindex->callback((Fl_Callback*)current_stateindex_cb);
 
  {
    Fl_Light_Button* t
      = new Fl_Light_Button(570,360,115,20,
			    "Auto rescale");
    t->callback((Fl_Callback*)btn_rescale_cb);
    t->value(false);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,460,105,20,
			    "Unfold PBC");
    t->callback(btn_bool_toggle_cb,
		&renderBox->unfoldPBC);
    t->value(renderBox->unfoldPBC);
  }

  btn_animate = new Fl_Light_Button(580,450,105,30,
				    "Animate");
  btn_animate->callback((Fl_Callback*)btn_animate_cb);
  btn_animate->value(instantAnimate);

  btn_simulate = new Fl_Light_Button(580,500,105,30,
				    "Simulate");
  btn_simulate->callback((Fl_Callback*)btn_simulate_cb);
  btn_simulate->value(false);

  checkbtn_quench = new Fl_Check_Button(475,500,105,30,
                                        "Quench");
  checkbtn_quench->value(false);

  {
    Fl_Button* t
      = new Fl_Button(570,320,115,20,
		      "Scale Up");
    t->tooltip(btn_scale_up_tooltip);
    t->callback((Fl_Callback*)btn_scale_up_cb);
  }

  {
    Fl_Button* t
      = new Fl_Button(570,340,115,20,
		      "Scale Down");
    t->tooltip(btn_scale_down_tooltip);
    t->callback((Fl_Callback*)btn_scale_down_cb);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(15,480,135,20,
			    "Multicolor atoms");
    t->tooltip(btn_colored_atoms_tooltip);
    t->callback(btn_bool_toggle_cb,
		&renderBox->nativeVertexColors);
    t->value(renderBox->nativeVertexColors);
  }

  {
    Fl_Button* t
      = new Fl_Button(15,500,135,20,
		      "Atoms color");
    t->tooltip(btn_atoms_color_tooltip);
    t->callback((Fl_Callback*)btn_atoms_color_cb);
  }

  {
    Fl_Button* t
      = new Fl_Button(15,520,135,20,
		      "BG color");
    t->tooltip(btn_bg_color_tooltip);
    t->callback((Fl_Callback*)btn_bg_color_cb);
  }

  {
    Fl_Slider* t
      = new Fl_Slider(160, 320, 25, 185,
		      "Render\nquality");
    t->labelsize(12);
    t->minimum(16);
    t->maximum(3);
    t->value(14);
    t->step(1);
    t->type(FL_VERT_NICE_SLIDER);
    t->callback(set_int_and_invalidate_cb,
		&renderBox->atomsQuality);
  }

  animate_delay = new Fl_Counter(620, 422, 65, 26, 
				 "Delay");
  animate_delay->labelsize(12);
  animate_delay->minimum(0.01);
  animate_delay->maximum(5.0);
  animate_delay->value(0.1);
  animate_delay->step(0.01);
  animate_delay->align(FL_ALIGN_LEFT);
  animate_delay->type(FL_SIMPLE_COUNTER);

  {
    Fl_Box* t 
      = new Fl_Box(FL_EMBOSSED_FRAME,
		   695,35,290,140,
		   "Collision tree / Motion trace");
    t->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(300,320,105,20,
			    "Atoms");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showAtoms);
    t->value(renderBox->showAtoms);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(695+10, 40+15, 135, 20,
			    "Show trace");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCTree);
    t->value(renderBox->showCTree);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(695+10+135, 40+15, 135, 20,
			    "All times at once");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCTreeAllTimes);
    t->value(renderBox->showCTreeAllTimes);
  }

  {
    Fl_Counter* t
      = new Fl_Counter(695+10+135, 40+15+20+5, 135, 20,
		       "Energy threshold, eV: ");
    t->labelsize(12);
    t->minimum(0.1);
    t->maximum(900);
    t->value(renderBox->energyThresholdCTree);
    t->lstep(5);
    t->step(0.1);
    t->align(FL_ALIGN_LEFT);
//   t->type(FL_SIMPLE_COUNTER);
    t->callback(set_double_cb,
		&renderBox->energyThresholdCTree);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(695+10, 40+15+20+5+20+5, 270, 20,
			    "Interpolate using lines");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCTreeConnected);
    t->value(renderBox->showCTreeConnected);
  }
    
  {
    Fl_Light_Button* t
      = new Fl_Light_Button(695+10, 40+15+20+5+20+5+20+20, 270, 20,
			    "'Phantom' atoms in keyframes");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCTreeAtoms);
    t->value(renderBox->showCTreeAtoms);
  }

  {
    Fl_Counter* t;
    t = new Fl_Counter(695+10+135, 40+15+20+5+20+5+20, 135, 20,
		       "Narrowness of lines: ");
    t->labelsize(12);
    t->minimum(1);
    t->maximum(900);
    t->value(renderBox->downscaleCTree);
    t->lstep(10);
    t->step(1);
    t->align(FL_ALIGN_LEFT);
    t->callback(set_double_cb,
		&renderBox->downscaleCTree);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,340,105,20,
			    "Axes");
    t->tooltip(btn_show_axes_tooltip);
    t->callback(btn_bool_toggle_cb,
		&renderBox->showAxes);
    t->value(renderBox->showAxes);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,360,210,20,
			    "Thermal bath borders");
    t->tooltip(btn_show_bath_tooltip);
    t->callback(btn_bool_toggle_cb,
		&renderBox->showBath);
    t->value(renderBox->showBath);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,380,210,20,
			    "Thermal bath sketch");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showBathSketch);
    t->value(renderBox->showBathSketch);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,320,105,20,
			    "Bonds");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showBonds);
    t->value(renderBox->showBonds);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,400,105,20,
			    "Custom1");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCustom1);
    t->value(renderBox->showCustom1);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195,420,105,20,
			    "Custom2");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCustom2);
    t->value(renderBox->showCustom2);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(195+105,400,105,20,
			    "Custom3");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showCustom3);
    t->value(renderBox->showCustom3);
  }

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(300,340,105,20,
			    "Tiny Atoms");
    t->callback(btn_bool_toggle_cb,
		&renderBox->tinyAtoms);
    t->value(renderBox->tinyAtoms);
  }

  {
    Fl_Box* t 
      = new Fl_Box(FL_UP_FRAME,15,320,135,160,
		   "Light Direction");
    t->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
  }

  {
    Fl_Slider* t
      = new Fl_Slider(30, 345, 25, 115,"X");
    t->labelsize(12);
    t->minimum(-1.0);
    t->maximum(+1.0);
    t->value(renderBox->light0_dir[0]);
    t->step(0.1);
    t->type(FL_VERT_NICE_SLIDER);
    t->callback(set_float_cb,
		renderBox->light0_dir+0);
  }

  {
    Fl_Slider* t
      = new Fl_Slider(70, 345, 25, 115,"Y");
    t->labelsize(12);
    t->minimum(-1.0);
    t->maximum(+1.0);
    t->value(renderBox->light0_dir[1]);
    t->step(0.1);
    t->type(FL_VERT_NICE_SLIDER);
    t->callback(set_float_cb,
		renderBox->light0_dir+1);
  }

  {
    Fl_Slider* t
      = new Fl_Slider(110, 345, 25, 115,"Z");
    t->labelsize(12);
    t->minimum(-1.0);
    t->maximum(+1.0);
    t->value(renderBox->light0_dir[2]);
    t->step(0.1);
    t->type(FL_VERT_NICE_SLIDER);
    t->callback(set_float_cb,
		renderBox->light0_dir+2);
  }

  {
    Fl_Button* t
      = new Fl_Button(415,320,145,20,
		      "Save video");
    t->callback((Fl_Callback*)btn_quick_save_image_cb);
  }

  {
    Fl_Button* t
      = new Fl_Button(415,340,145,20,
		      "Save image");
    t->tooltip(btn_save_image_tooltip);
    t->callback((Fl_Callback*)btn_save_image_cb);
  }

  {
    Fl_Button* t
      = new Fl_Button(415,360,145,20,
		      "Save hi-res image");
    t->callback((Fl_Callback*)btn_save_tiled_image_cb);
  }

  {
    Fl_Button* t
      = new Fl_Button(415,380,145,20,
		      "Save state (rotated)");
    t->callback((Fl_Callback*)btn_save_state_cb, (void*)0);
  }

  {
    Fl_Button* t
      = new Fl_Button(415,400,145,20,
		      "Save state");
    t->callback((Fl_Callback*)btn_save_state_cb, (void*)1);
  }

  {
    Fl_Box* t
      = new Fl_Box(FL_EMBOSSED_FRAME,415,35,270,280,
		   "Rotate and clip");
    t->align(FL_ALIGN_TOP | FL_ALIGN_INSIDE);
  }

  roll_x = new Fl_Roller(455, 80, 20, 215, "X");
  roll_x->tooltip(roll_x_tooltip);
  roll_x->type(0);
  roll_x->labelsize(12);
  roll_x->minimum(0);
  roll_x->maximum(360);
  roll_x->value(0);
  roll_x->step(1);
  roll_x->callback((Fl_Callback*)roll_x_cb);

  val_xmin = new Fl_Slider(430, 80, 25, 215, "min");
  val_xmin->type(FL_VERT_NICE_SLIDER);
  val_xmin->labelsize(12);
  val_xmin->minimum(101);
  val_xmin->maximum(-1);
  val_xmin->value(0);
  val_xmin->step(0.1);
  val_xmin->when(FL_WHEN_RELEASE);
  val_xmin->callback((Fl_Callback*)val_xminmax_cb);

  val_xmax = new Fl_Slider(475, 80, 25, 215, "max");
  val_xmax->type(FL_VERT_NICE_SLIDER);
  val_xmax->labelsize(12);
  val_xmax->minimum(101);
  val_xmax->maximum(-1);
  val_xmax->value(100);
  val_xmax->step(0.1);
  val_xmax->when(FL_WHEN_RELEASE);
  val_xmax->callback((Fl_Callback*)val_xminmax_cb);

  roll_y = new Fl_Roller(540, 80, 20, 215, "Y");
  roll_y->tooltip(roll_y_tooltip);
  roll_y->type(0);
  roll_y->labelsize(12);
  roll_y->minimum(0);
  roll_y->maximum(360);
  roll_y->value(0);
  roll_y->step(1);
  roll_y->callback((Fl_Callback*)roll_y_cb);

  val_ymin = new Fl_Slider(515, 80, 25, 215, "min");
  val_ymin->type(FL_VERT_NICE_SLIDER);
  val_ymin->labelsize(12);
  val_ymin->minimum(101);
  val_ymin->maximum(-1);
  val_ymin->value(0);
  val_ymin->step(0.1);
  val_ymin->when(FL_WHEN_RELEASE);
  val_ymin->callback((Fl_Callback*)val_xminmax_cb);

  val_ymax = new Fl_Slider(560, 80, 25, 215, "max");
  val_ymax->type(FL_VERT_NICE_SLIDER);
  val_ymax->labelsize(12);
  val_ymax->minimum(101);
  val_ymax->maximum(-1);
  val_ymax->value(100);
  val_ymax->step(0.1);
  val_ymax->when(FL_WHEN_RELEASE);
  val_ymax->callback((Fl_Callback*)val_xminmax_cb);

  roll_z = new Fl_Roller(625, 80, 20, 215, "Z");
  roll_z->tooltip(roll_z_tooltip);
  roll_z->type(0);
  roll_z->labelsize(12);
  roll_z->minimum(0);
  roll_z->maximum(360);
  roll_z->value(0);
  roll_z->step(1);
  roll_z->callback((Fl_Callback*)roll_z_cb);

  val_zmin = new Fl_Slider(600, 80, 25, 215, "min");
  val_zmin->type(FL_VERT_NICE_SLIDER);
  val_zmin->labelsize(12);
  val_zmin->minimum(101);
  val_zmin->maximum(-1);
  val_zmin->value(0);
  val_zmin->step(0.1);
  val_zmin->when(FL_WHEN_RELEASE);
  val_zmin->callback((Fl_Callback*)val_xminmax_cb);

  val_zmax = new Fl_Slider(645, 80, 25, 215, "max");
  val_zmax->type(FL_VERT_NICE_SLIDER);
  val_zmax->labelsize(12);
  val_zmax->minimum(101);
  val_zmax->maximum(-1);
  val_zmax->value(100);
  val_zmax->step(0.1);
  val_zmax->when(FL_WHEN_RELEASE);
  val_zmax->callback((Fl_Callback*)val_xminmax_cb);

  {
    Fl_Button* t
      = new Fl_Button(430, 55, 25*2+20, 20,
		      "45\xB0 rot");
    t->callback(btn_view_cb,roll_x);
  }

  {
    Fl_Button* t
      = new Fl_Button(515, 55, 25*2+20, 20,
		      "45\xB0 rot");
    t->callback(btn_view_cb,roll_y);
  }

  {
    Fl_Button* t
      = new Fl_Button(600, 55, 25*2+20, 20,
		      "45\xB0 rot");
    t->callback(btn_view_cb,roll_z);
  }

  {
    new Fl_Box(FL_UP_FRAME,15,35,385+5,235+45,NULL);
  }

  atom_info = new Fl_Multiline_Output(25,55,370,165,
				      "Selected Atom info");
  atom_info->textcolor(FL_BLUE);
  atom_info->align(FL_ALIGN_TOP);
  atom_info->value(""); 

  atom_coords_x = new Fl_Value_Input(25+40,235,110,20,"Pos:");
  atom_coords_x->range(-10000,+10000);
  atom_coords_x->precision(4);
  atom_coords_x->align(FL_ALIGN_LEFT);
  atom_coords_x->callback(set_atom_properties_cb);

  atom_coords_y = new Fl_Value_Input(25+40+110*1,235,110,20,"");
  atom_coords_y->range(-10000,+10000);
  atom_coords_y->precision(4);
  atom_coords_y->align(FL_ALIGN_LEFT);
  atom_coords_y->callback(set_atom_properties_cb);

  atom_coords_z = new Fl_Value_Input(25+40+110*2,235,110,20,"");
  atom_coords_z->range(-10000,+10000);
  atom_coords_z->precision(4);
  atom_coords_z->align(FL_ALIGN_LEFT);
  atom_coords_z->callback(set_atom_properties_cb);

  atom_v_x = new Fl_Value_Input(25+40,235+20,110,20,"Vel:");
  atom_v_x->range(-1e10,+1e10);
  atom_v_x->precision(2);
  atom_v_x->align(FL_ALIGN_LEFT);
  atom_v_x->callback(set_atom_properties_cb);

  atom_v_y = new Fl_Value_Input(25+40+110*1,235+20,110,20,"");
  atom_v_y->range(-1e10,+1e10);
  atom_v_y->precision(2);
  atom_v_y->align(FL_ALIGN_LEFT);
  atom_v_y->callback(set_atom_properties_cb);

  atom_v_z = new Fl_Value_Input(25+40+110*2,235+20,110,20,"");
  atom_v_z->range(-1e10,+1e10);
  atom_v_z->precision(2);
  atom_v_z->align(FL_ALIGN_LEFT);
  atom_v_z->callback(set_atom_properties_cb);

  current_atomindex = new Fl_Counter(25,285,220,20,NULL);
  current_atomindex->align(FL_ALIGN_LEFT);
  current_atomindex->lstep(100);
  current_atomindex->step(1);
  current_atomindex->precision(0);
//  current_atomindex->bounds(0,renderBox->getAtomsCount()-1);
  current_atomindex->value(0);
  current_atomindex->callback((Fl_Callback*)current_atomindex_cb);

  {
    Fl_Light_Button* t
      = new Fl_Light_Button(255,285,140,20,
			    "Show selected");
    t->callback(btn_bool_toggle_cb,
		&renderBox->showSelected);
    t->value(renderBox->showSelected);
  }

  atoms_view_group->end();

  {
    Fl_Group* t
      = new Fl_Group(5,25,990,540,"About");
    t->begin();
    product_info = 
      std::string("Molecular dynamics trajectory viewer")+"\nusing "+
      "MDTK" + " " +
      mdtk::release_info.PRODUCT_VERSION + 
      "\n\n" + 
      "Copyright (C) 2003-2013 Oleksandr Yermolenko\n <oleksandr.yermolenko@@gmail.com>\n\n" +
      "Run the program with --version or --help options for details."
      ;
    
    {
      Fl_Box* t
	= new Fl_Box(5,25,990,540,product_info.c_str());
      t->labelfont(FL_HELVETICA_BOLD);
    }
    t->end();
  }
  
  tabs->end();
  end();

  renderBox->show();

  show();
  renderBox->redraw();
  redraw();
  renderBox->invalidate();
  renderBox->hide();
  renderBox->show();

  updateStateList();

  renderBox->allowRescale = false;
//  renderBox->reArrange(-1,101,-1,101,-1,101);
  renderBox->redraw();

  callback(window_cb);
  { // rotate around x axis by 6*45=270 degrees
    for(size_t ri = 0; ri < 6; ri++)
    {
      renderBox->make_current();
      btn_view_cb(roll_x,roll_x);
      renderBox->redraw();
    }
  }
  if (btn_animate->value()) Fl::add_timeout(1.0, timer_callback);
}

void
MainWindow::updateStateList()
{
  stateList.clear();

  MDTrajectory::const_iterator t = renderBox->mdt.begin();
  while (t != renderBox->mdt.end())
  {
    stateList.push_back(t->second.name);
    ++t;
  }

  current_stateindex->bounds(0,stateList.size()-1);

  current_atomindex->bounds(0,renderBox->getAtomsCount()-1);

  if (stateList.size() > 0)
    loadNewSnapshot(0);

  REQUIRE(renderBox->selectedAtomIndex < renderBox->getAtomsCount());
  setAtomViewIndex(renderBox->selectedAtomIndex);
}

MainWindow::~MainWindow()
{
}

void
MainWindow::redrawGL()
{
  renderBox->redraw();
}

void
MainWindow::btn_save_image_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  VisBox* VisBox_Ptr;
  VisBox_Ptr = MainWindow_Ptr->renderBox;

  char fname[1000];
#ifdef MDTRAJVIEW_PNG
  sprintf(fname,"%010d-small.png",
          int(MainWindow_Ptr->current_stateindex->value()));
  char *tmp_filename = fl_file_chooser
    (
      "Choose a file to save image to ...",
      "PNG files (*.png)",
      fname,0
      );
#else
  sprintf(fname,"%010d-small.bmp",
          int(MainWindow_Ptr->current_stateindex->value()));
  char *tmp_filename = fl_file_chooser
    (
      "Choose a file to save image to ...",
      "Windows Bitmap Files (*.bmp)",
      fname,0
      );
#endif
  if (tmp_filename && (!fl_filename_isdir(tmp_filename)))
  {
    VisBox_Ptr->saveImageToFile(tmp_filename);

    VisBox_Ptr->make_current();

//    fl_message("It seems that file was successfully saved !");
  }
}

void
MainWindow::btn_save_tiled_image_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  VisBox* VisBox_Ptr;
  VisBox_Ptr = MainWindow_Ptr->renderBox;

  char fname[1000];
#ifdef MDTRAJVIEW_PNG
  sprintf(fname,"%010d.png",
          int(MainWindow_Ptr->current_stateindex->value()));
  char *tmp_filename = fl_file_chooser
    (
      "Choose a file to save image to ...",
      "PNG files (*.png)",
      fname,0
      );
#else
  sprintf(fname,"%010d.bmp",
          int(MainWindow_Ptr->current_stateindex->value()));
  char *tmp_filename = fl_file_chooser
    (
      "Choose a file to save image to ...",
      "Windows Bitmap Files (*.bmp)",
      fname,0
      );
#endif
  if (tmp_filename && (!fl_filename_isdir(tmp_filename)))
  {
    VisBox_Ptr->saveTiledImageToFile(tmp_filename);

    VisBox_Ptr->make_current();

//    fl_message("It seems that file was successfully saved !");
  }
}

void
MainWindow::quickSaveBitmap()
{
  char tmp_filename[4096];

#ifdef MDTRAJVIEW_PNG
  sprintf(tmp_filename,"%010d__%.5f_ps__.png",stateIndex,renderBox->ml_->simTime/ps);
#else
  sprintf(tmp_filename,"%010d__%.5f_ps__.bmp",stateIndex,renderBox->ml_->simTime/ps);
#endif

  if (!fl_filename_isdir(tmp_filename))
  {
    renderBox->saveImageToFile(tmp_filename);
  }
}

void
MainWindow::btn_quick_save_image_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  VisBox* VisBox_Ptr;
  VisBox_Ptr = MainWindow_Ptr->renderBox;

//  char fname[1000];
//  sprintf(fname,"mdtk-video-%ld",std::time(NULL));
  char *tmp_filename = fl_dir_chooser
    (
      "Choose a directory to save video to ...",
//      fname,0
      "",0
      );
//  if (tmp_filename) yaatk::mkdir(tmp_filename);
  if (tmp_filename && (fl_filename_isdir(tmp_filename)))
  {
    std::string trajdir = yaatk::getcwd();
    yaatk::chdir(tmp_filename);

    char videoDirName[1024];
    bool videoDirExists = true;
    int videoDirIndex = 0;
    while (videoDirExists)
    {
      sprintf(videoDirName,"video%02d",videoDirIndex);
      if (fl_filename_isdir(videoDirName))
        videoDirIndex++;
      else
        videoDirExists = false;
      REQUIRE(videoDirIndex < 100);
    }

    yaatk::mkdir(videoDirName);
    yaatk::chdir(trajdir.c_str());
    for(size_t i = 0; i < MainWindow_Ptr->stateList.size()/*-1*/; i++)
    {
//	MainWindow_Ptr->loadNewState(i);
      MainWindow_Ptr->current_stateindex->value(i);
      current_stateindex_cb(MainWindow_Ptr->current_stateindex,NULL);
      yaatk::chdir(tmp_filename);
      yaatk::chdir(videoDirName);
      while (!Fl::ready()) {};
      Fl::flush();
      while (!Fl::ready()) {};
      MainWindow_Ptr->quickSaveBitmap();
//        sleep(1);
      yaatk::chdir(trajdir.c_str());
    }
  }
}

void
MainWindow::btn_save_mde_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  VisBox* VisBox_Ptr;
  VisBox_Ptr = MainWindow_Ptr->renderBox;

  char *tmp_filename = fl_file_chooser
    (
      "Choose a file to save MDE in ...",
      "mde (*.mde)",
      0,0
      );
  if (tmp_filename && (!fl_filename_isdir(tmp_filename)))
  {
    VisBox_Ptr->saveToMDE(tmp_filename);
  }  
}

void
MainWindow::btn_save_state_cb(Fl_Widget *w, void *p)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  VisBox* VisBox_Ptr;
  VisBox_Ptr = MainWindow_Ptr->renderBox;

  bool discardRotation = (bool) p;

  char *tmp_filename = fl_dir_chooser("Choose a directory to save files to ...","",0);
  if (tmp_filename && (fl_filename_isdir(tmp_filename)))
  {
    yaatk::ChDir cd(tmp_filename,false);
    const char *id = fl_input("Id of the new state:","new000");
    if (id)
      VisBox_Ptr->saveState(id, discardRotation);
    else
      VisBox_Ptr->saveState("new", discardRotation);
  }
}

void
MainWindow::btn_view_cb(Fl_Widget *w, void *p)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  Fl_Roller* roll = (Fl_Roller*) p;
  double v = roll->value();
  v -= 45; if (v <= 0) v = 360+v;
  roll->value(v);
  roll->callback()(roll,NULL);
}

void
MainWindow::btn_scale_up_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());
  VisBox* VisBox_Ptr
    = MainWindow_Ptr->renderBox;

  if (1.2*VisBox_Ptr->scale<=VisBox_Ptr->maxScale)
    VisBox_Ptr->scale = 1.2*VisBox_Ptr->scale;

  VisBox_Ptr->redraw();
}

void
MainWindow::btn_scale_down_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());
  VisBox* VisBox_Ptr = MainWindow_Ptr->renderBox;

  VisBox_Ptr->scale = VisBox_Ptr->scale/1.2;
  VisBox_Ptr->redraw();
}

void
MainWindow::btn_bool_toggle_cb(Fl_Widget *w, void *pbool)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  bool v = ((Fl_Light_Button *)w)->value();
  *((bool*)pbool) = v;
  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::set_double_cb(Fl_Widget *w, void *pdouble)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  double v = ((Fl_Slider *)w)->value();
  *((double*)pdouble) = v;

  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::set_float_cb(Fl_Widget *w, void *pfloat)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  float v = ((Fl_Slider *)w)->value();
  *((float*)pfloat) = v;

  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::set_int_cb(Fl_Widget *w, void *pint)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  int v = ((Fl_Slider *)w)->value();
  *((int*)pint) = v;

  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::set_int_and_invalidate_cb(Fl_Widget *w, void *pint)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  int v = ((Fl_Slider *)w)->value();
  *((int*)pint) = v;

  MainWindow_Ptr->renderBox->invalidate();
  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::roll_x_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  double new_rot_x,old_rot_x;
  new_rot_x = ((Fl_Roller *)w)->value();
  old_rot_x = MainWindow_Ptr->renderBox->old_rot_x;

  MainWindow_Ptr->renderBox->
    rollAround(new_rot_x - old_rot_x,1,0,0);
}

void
MainWindow::roll_y_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  double new_rot_y,old_rot_y;
  new_rot_y = ((Fl_Roller *)w)->value();
  old_rot_y = MainWindow_Ptr->renderBox->old_rot_y;

  MainWindow_Ptr->renderBox->
    rollAround(new_rot_y - old_rot_y,0,1,0);
}

void
MainWindow::roll_z_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  double new_rot_z,old_rot_z;
  new_rot_z = ((Fl_Roller *)w)->value();
  old_rot_z = MainWindow_Ptr->renderBox->old_rot_z;

  MainWindow_Ptr->renderBox->
    rollAround(new_rot_z - old_rot_z,0,0,1);
}

void
MainWindow::val_xminmax_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  MainWindow_Ptr->renderBox->
    reArrange(
      MainWindow_Ptr->val_xmin->value()/100.0,
      MainWindow_Ptr->val_xmax->value()/100.0,
      MainWindow_Ptr->val_ymin->value()/100.0,
      MainWindow_Ptr->val_ymax->value()/100.0,
      MainWindow_Ptr->val_zmin->value()/100.0,
      MainWindow_Ptr->val_zmax->value()/100.0
      );
}

void
MainWindow::btn_rescale_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  bool v = ((Fl_Light_Button *)w)->value();
  MainWindow_Ptr->renderBox->allowRescale = v;
  MainWindow_Ptr->renderBox->reArrange(-1,101,-1,101,-1,101);
  MainWindow_Ptr->renderBox->redraw();
}

void
MainWindow::btn_animate_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  bool v = ((Fl_Light_Button *)w)->value();
    
//    double delay = MainWindow_Ptr->animate_delay->value();
    
  if (v) Fl::repeat_timeout(0.5, timer_callback);

//    MainWindow_Ptr->renderBox->SetAllowRescale(v);
}

void
MainWindow::btn_simulate_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());

  bool v = ((Fl_Light_Button *)w)->value();

  if (v)
  {
    bool quench = MainWindow_Ptr->checkbtn_quench->value();
    MainWindow_Ptr->renderBox->loadDataFromSimulation(quench);
  }
}

void
MainWindow::addMDSnapshot(const MDSnapshot& s)
{
  renderBox->mdt[s.time] = s;

  int lastStateIndex = current_stateindex->value();
  updateStateList();
  loadNewSnapshot(stateList.size()-1);
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
    VisBox_Ptr = MainWindow_Ptr->renderBox;

    int oldval = MainWindow_Ptr->current_stateindex->value();
    int maxval = MainWindow_Ptr->current_stateindex->maximum();
    MainWindow_Ptr->current_stateindex->value((oldval+1)%(maxval+1));
    current_stateindex_cb(MainWindow_Ptr->current_stateindex,NULL);

    double delay = MainWindow_Ptr->animate_delay->value();

    if (MainWindow_Ptr->btn_animate->value()) Fl::repeat_timeout(delay, timer_callback);
  }
}

void
MainWindow::btn_bg_color_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());
  Color c = MainWindow_Ptr->renderBox->bgColor;

  unsigned char r,g,b;
  analyseRGB(c,r,g,b);

  if (fl_color_chooser("Choose Background Color",r,g,b))
  {
    c = combineRGB(r,g,b);
    MainWindow_Ptr->renderBox->bgColor = c;
    MainWindow_Ptr->renderBox->redraw();
  }
}

void
MainWindow::btn_atoms_color_cb(Fl_Widget *w, void *)
{
  MainWindow* MainWindow_Ptr;
  MainWindow_Ptr =
    (MainWindow*)(w->parent()->parent()->parent());
  Color c = MainWindow_Ptr->renderBox->vertexColor;

  unsigned char r,g,b;
  analyseRGB(c,r,g,b);

  if (fl_color_chooser("Choose Vertex Color",r,g,b))
  {
    c = combineRGB(r,g,b);
    MainWindow_Ptr->renderBox->vertexColor = c;
    MainWindow_Ptr->renderBox->redraw();
  }
}

void
MainWindow::window_cb(Fl_Widget* widget, void*)
{
//    if (fl_choice("Do you really want to exit?","No","Yes",NULL)==1)
  {
    ((Fl_Window*)widget)->hide();
    exit(0);
  }
}

const char* MainWindow::
btn_show_axes_tooltip = "Show/Hide Axes";
const char* MainWindow::
btn_show_bath_tooltip = "Show/Hide Thermal Bath";
const char* MainWindow::
btn_colored_atoms_tooltip = "Toggle Colored/Unified Atoms";
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

void
MainWindow::loadNewSnapshot(int index)
{
  //  if (index > 0) renderBox->SetAllowRescale(false);
  stateIndex = index;
  label((std::string("MDTK Trajectory Viewer [Control Window] - ")+stateList[stateIndex]).c_str());
  renderBox -> label((std::string("MDTK 3D View - ")+stateList[stateIndex]).c_str());
  renderBox -> loadNewSnapshot(stateIndex);

  if (stateIndex != current_stateindex->value())
    current_stateindex->value(stateIndex);

  setAtomViewIndex(int(current_atomindex->value())/*-1*/);
}  

}

xmde::MainWindow* MainWindow_GlobalPtr = NULL;

