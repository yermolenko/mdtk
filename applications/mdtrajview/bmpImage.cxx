/*
   The BMP image class.

   Copyright (C) 2004, 2009 Oleksandr Yermolenko
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

#include "bmpImage.hpp"

namespace grctk
{

  void bmpImage::write_word(const void *ptr, FILE *stream)
  {
    fwrite(ptr,2,1,stream);	
  }
  void bmpImage::write_dword(const void *ptr, FILE *stream)
  {
    fwrite(ptr,4,1,stream);	
  }
  void bmpImage::write_byte(const void *ptr, FILE *stream)
  {
    fwrite(ptr,1,1,stream);	
  }

  bmpImage::bmpImage(unsigned long w, unsigned long h,
		     unsigned char *data)
  {
    unsigned long i;
    w_ = w;
    h_ = h;
    data_size = w_*h_*3;
//std::printf("SIZE=%lu %lu %lu",w_,h_,data_size);
    data_ = new unsigned char[data_size]; 
    for(i = 0; i < data_size;i++)
      data_[i] = data[i];
  }

  bmpImage::~bmpImage()
  {
    delete [] data_;
  }

  int bmpImage::SaveToFile(const char* filename)
  {
    FILE	*fo;
    if ((fo = fopen(filename, "wb")) == NULL)
      {
	fprintf(stderr, "Cannot open output file.\n");
	return -1;
      }
	
    // signature
    write_byte("B",fo);
    write_byte("M",fo);

    ///////////////////////////////////
    // file size
    ///////////////////////////////////
    unsigned long bmpdata_size = (((w_*3)%4==0)?(w_*3):(((w_*3)/4+1)*4))*h_;
    unsigned long file_size = bmpdata_size+54;
    write_dword(&file_size,fo);

    // reserved
    //	unsigned long reserved = 0;
    //	write_dword(&reserved,fo);
	
    ///////////////////////////////////
    // bmp data offset
    ///////////////////////////////////
    //	unsigned long bmpdata_offset = 0;
    //	write_dword(&bmpdata_offset,fo);

    unsigned char zero = 0;
    unsigned char z36 = 0x36;
    write_byte(&zero,fo);
    write_byte(&zero,fo);
    write_byte(&zero,fo);
    write_byte(&zero,fo);
    write_byte(&z36,fo);
    write_byte(&zero,fo);
    write_byte(&zero,fo);
    write_byte(&zero,fo);

	
    // length of bmp header
    unsigned long bmp_header = 0x28;
    write_dword(&bmp_header,fo);
	
    // width and height in pixels
    write_dword(&w_,fo);
    write_dword(&h_,fo);
	
    // number of planes
    unsigned short number_of_planes = 1;
    write_word(&number_of_planes,fo);
	
    // bpp
    unsigned short bpp = 24;
    write_word(&bpp,fo);
	
    // compression
    unsigned long compression = 0;
    write_dword(&compression,fo);
	
    //////// BMP data size ///////////
    write_dword(&bmpdata_size,fo);
    //////////////////////////////////

    ///// Horizontal resolution expressed in pixel per meter.	
    unsigned long hres = 2834;
    write_dword(&hres,fo);

    ///// Vertical resolution expressed in pixel per meter.	
    unsigned long vres = 2834;
    write_dword(&vres,fo);

    unsigned long colors = 0;
    write_dword(&colors,fo);
	
    unsigned long imp_colors = 0;//16777216;
    write_dword(&imp_colors,fo);

    // pallete
	
    // bitmap data
    unsigned long actual_bmpdata_size = 0;;
    unsigned long i;	
    for(i = 0; i < h_;i++)
    {
      for(unsigned long j = 0; j < w_;j++)
        for(unsigned long b = 0; b < 3;b++)
          {write_byte(&data_[(i*w_+j)*3+b],fo);actual_bmpdata_size++;}
      unsigned char zero = 0;
      if ((w_*3)%4!=0)
        for(unsigned long zi = 1;zi <= ((w_*3)/4+1)*4-w_*3; zi++)
          {write_byte(&zero,fo);actual_bmpdata_size++;}
    }
    if (actual_bmpdata_size != bmpdata_size) {fprintf(stderr,"Internal error when writing BMP\n");};
/*
    for(i = 0; i < data_size;i++)
      write_byte(&data_[i],fo);
*/	
    fclose(fo);	

    return 0;
  }

}

