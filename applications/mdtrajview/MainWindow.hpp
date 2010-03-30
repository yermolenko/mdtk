/*
   The MainWindow class for the molecular dynamics trajectory viewer
   (header file)

   Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009 Oleksandr
   Yermolenko <oleksandr.yermolenko@gmail.com>

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

#ifndef	mde_MainWindow_h
#define	mde_MainWindow_h

#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Tabs.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Roller.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Chart.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Counter.H>
#include <FL/Fl_Value_Output.H>
#include <FL/Fl_Multiline_Output.H> 

#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Radio_Light_Button.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Chart.H>


#include <FL/Fl_Gl_Window.H>
#include <FL/Enumerations.H>
#include <FL/gl.h>
#include <FL/fl_ask.H>

#include <FL/Fl_Color_Chooser.H>
#include <FL/Fl_File_Chooser.H>


#include "bmpImage.hpp"
#include "VisBox.hpp"

#include <vector>

namespace xmde
{

  class MainWindow : public Fl_Window
  {
    void quickSaveBitmap();

    std::string product_info;

    Fl_Light_Button
      *btn_show_axes,
      *btn_show_bath,
      *btn_show_selected,
      *btn_rescale,
      *btn_animate,
      *btn_unified_atoms,
      *btn_fixed_lights;

    Fl_Button  *btn_quick_save_image;
    Fl_Button  *btn_save_image;
    Fl_Button  *btn_save_mde;

    Fl_Button  *btn_bg_color;
    Fl_Button  *btn_atoms_color;

    Fl_Button  *btn_scale_up;
    Fl_Button  *btn_scale_down;

    Fl_Slider* light_x_dir;
    Fl_Slider* light_y_dir;
    Fl_Slider* light_z_dir;


    Fl_Counter* current_atomindex;

    Fl_Counter* current_stateindex;

    std::vector<double> ev;
  
    static	const	char
      *btn_show_axes_tooltip,
      *btn_show_bath_tooltip,
      *btn_unified_atoms_tooltip,
      *btn_fixed_lights_tooltip,
      *roll_x_tooltip,
      *roll_y_tooltip,
      *roll_z_tooltip,
      *btn_atoms_color_tooltip,
      *btn_bg_color_tooltip,
      *btn_save_image_tooltip,
      *btn_scale_up_tooltip,
      *btn_scale_down_tooltip;
					
    VisBox  *atoms_view_box;//,*tmp_gb;

//    Fl_Gl_Window* glw;


    Fl_Tabs     *tabs;
    Fl_Group    *atoms_view_group;
    Fl_Group    *about_group;
//    Fl_Group    *info_group;
    Fl_Box      *about_box;

    Fl_Multiline_Output* atom_info; 
  
    Fl_Roller			*roll_x;
    Fl_Roller			*roll_y;
    Fl_Roller			*roll_z;

    Fl_Slider       *val_xmin,
      *val_xmax,
      *val_ymin,
      *val_ymax,
      *val_zmin,
      *val_zmax;

    Fl_Counter       *animate_delay;
                  
    Fl_Slider       *atom_quality;

    static void current_atomindex_cb(Fl_Widget *, void *);
    static void current_stateindex_cb(Fl_Widget *, void *);


    static void light_x_dir_cb(Fl_Widget *, void *);
    static void light_y_dir_cb(Fl_Widget *, void *);
    static void light_z_dir_cb(Fl_Widget *, void *);


    static void roll_x_cb(Fl_Widget *, void *);
    static void roll_y_cb(Fl_Widget *, void *);
    static void roll_z_cb(Fl_Widget *, void *);

    static void val_xmin_cb(Fl_Widget *, void *);

    static void atom_quality_cb(Fl_Widget *, void *);

    static void btn_show_axes_cb(Fl_Widget *, void *);
    static void btn_show_bath_cb(Fl_Widget *, void *);
    static void btn_show_selected_cb(Fl_Widget *, void *);
    static void btn_rescale_cb(Fl_Widget *, void *);
    static void timer_callback(void *);
    static void btn_animate_cb(Fl_Widget *, void *);
    static void btn_unified_atoms_cb(Fl_Widget *, void *);
    static void btn_fixed_lights_cb(Fl_Widget *, void *);

    static void btn_bg_color_cb(Fl_Widget *, void *);
    static void btn_atoms_color_cb(Fl_Widget *, void *);

    static void btn_save_image_cb(Fl_Widget *, void *);
    static void btn_quick_save_image_cb(Fl_Widget *, void *);
    static void btn_save_mel_cb(Fl_Widget *, void *);
    static void btn_save_mde_cb(Fl_Widget *, void *);
  
    static void btn_scale_up_cb(Fl_Widget *, void *);
    static void btn_scale_down_cb(Fl_Widget *, void *);

    static void window_cb(Fl_Widget *, void *);

  public:
    int  handle(int);
    std::string base_state_filename;
    MainWindow(std::string&,std::vector<std::string>&,VisBox* ,bool);
    ~MainWindow();
    void redrawGL();
    void	setAtomViewIndex(int index);
  
  
    std::vector<std::string> stateList;
    int stateIndex;
    void loadNewState(int index);
  
  
    char *  log_buffer;
    int log_pos;
  
    void out(std::string );
    void clear_out();
  };

}

extern xmde::MainWindow* MainWindow_GlobalPtr;

#endif

