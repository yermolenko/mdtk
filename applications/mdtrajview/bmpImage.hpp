/*
   The BMP image class header file.

   Copyright (C) 2004, 2009, 2010, 2012 Oleksandr Yermolenko
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

#ifndef grctk_bmp_image_h
#define grctk_bmp_image_h

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace grctk 
{

  class bmpImage
  {
    unsigned char *data_;
    unsigned long data_size;
    unsigned long w_,h_;
    void write_word(const void *ptr, std::FILE *stream);
    void write_dword(const void *ptr, std::FILE *stream);
    void write_byte(const void *ptr, std::FILE *stream);
  public:
    bmpImage(unsigned long w, unsigned long h, unsigned char *data);
    ~bmpImage();
    int SaveToFile(const char* filename);
    unsigned long getPixel(unsigned long i, unsigned long j)
      {
	unsigned long c;
	std::memcpy(&c,&(data_[(i*w_+j)*3+0]),3);
/*
	unsigned long c = 0;
	c |= data_[(i*w_+j)*3+0] << 8;
	c |= data_[(i*w_+j)*3+1] << 16;
	c |= data_[(i*w_+j)*3+2] << 32;
*/
	return c;
      }
    void setPixel(unsigned long i, unsigned long j, unsigned long c)
      {
	std::memcpy(&(data_[(i*w_+j)*3+0]),&c,3);
/*
	data_[(i*w_+j)*3+0] = (c & (1<<8-1)) >> 0;
	data_[(i*w_+j)*3+1] = (c & (1<<16-1)) >> 8;
	data_[(i*w_+j)*3+2] = (c & (1<<32-1)) >> 16;
*/
      }
    unsigned char* getRawDataPtr() {return data_;}
  };

}

#endif

