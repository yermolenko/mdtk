/*
   Yet another auxiliary toolkit.

   Copyright (C) 2003, 2005, 2006, 2009, 2010 Oleksandr Yermolenko
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

#include "mdtk/yaatk.hpp"
#include "mdtk/Exception.hpp"

#include <zlib.h>
#include <cstring>

#include <sstream>

namespace yaatk
{

/*extern*/ int yaatk_extraID = 0;

#define MDTK_GZ_BUFFER_SIZE 10000

int
zip_file(const char *zipName, const char* unzipName)
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  gzFile   unzipped = gzopen(unzipName,"rb");
//  if (unzipped == 0) TRACE(unzipName);
  REQUIRE(unzipped != 0);
  gzFile   zipped   = gzopen(zipName,"wb");
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while((unzippedFileSize = gzread(unzipped,buf,MDTK_GZ_BUFFER_SIZE)) > 0)
  {
    REQUIRE(unzippedFileSize != -1);
    int bytesWritten    = gzwrite(zipped,buf,unzippedFileSize);
    REQUIRE(unzippedFileSize == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(unzipped);  
  gzclose(zipped);
  return 0;
}

int
zip_stringstream(const char *zipName, std::stringstream &uzs)
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  gzFile   zipped   = gzopen(zipName,"wb");
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while((uzs.read(buf,MDTK_GZ_BUFFER_SIZE),unzippedFileSize = uzs.gcount()) > 0)
  {
    REQUIRE(unzippedFileSize != -1);
    int bytesWritten    = gzwrite(zipped,buf,unzippedFileSize);
    REQUIRE(unzippedFileSize == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(zipped);
  return 0;
}

std::string extractDir(std::string trajNameFinal)
{
    std::string mde_dirname = trajNameFinal;
    {
      int i;
      for(i = mde_dirname.size()-1; i >= 0; i--)
        if (mde_dirname[i] == DIR_DELIMIT_CHAR) break;
      i++;
      mde_dirname.resize(i);
    }  
    return mde_dirname;
}  

std::string extractLastItem(std::string trajNameFinal)
{
    std::string mde_dirname = trajNameFinal;
    std::string res;
//    std::string x;
    if (mde_dirname[mde_dirname.size()-1]==DIR_DELIMIT_CHAR)
      mde_dirname.resize(mde_dirname.size()-1);
    {
      int i;
      for(i = mde_dirname.size()-1; i >= 0; i--)
        if (mde_dirname[i] == DIR_DELIMIT_CHAR) break;
      i++;
      res = mde_dirname.substr(i,mde_dirname.size()-i);
    }  
    return res;
}  

  // ZippedStreams zippedStreams; //deprecated


#define  MDTK_GZ_BUFFER_SIZE 10000

void unzip_file(const char *zipName, const char* unzipName)
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
//  TRACE(zipName);
  gzFile zipped   = gzopen(zipName,"rb");
  REQUIRE(zipped != 0);
  FILE*  unzipped = fopen(unzipName,"wb");
//  if (unzipped == 0) TRACE(unzipName);
  REQUIRE(unzipped != 0);
  int unzippedFileSize;
  while ((unzippedFileSize = gzread(zipped,buf,MDTK_GZ_BUFFER_SIZE)) > 0)
  {
    int bytesWritten = fwrite(buf,unzippedFileSize,1,unzipped);
    REQUIRE(1 == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(zipped);
  fclose(unzipped);  
}

void unzip_stringstream(const char *zipName, std::stringstream& os)
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  gzFile zipped   = gzopen(zipName,"rb");
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while ((unzippedFileSize = gzread(zipped,buf,MDTK_GZ_BUFFER_SIZE)) > 0)
  {
    os.write(buf,unzippedFileSize);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(zipped);
}


#define  MDTK_FILECMP_BUFFER_SIZE 10000

bool
isIdentical(const std::string& file1,const std::string& file2)
{
  using mdtk::Exception;
  char buf1[MDTK_FILECMP_BUFFER_SIZE];
  char buf2[MDTK_FILECMP_BUFFER_SIZE];
  gzFile  f1   = gzopen(file1.c_str(),"rb");
  if (f1 == 0) return false;
  gzFile  f2   = gzopen(file2.c_str(),"rb");
  if (f2 == 0) {gzclose(f1);return false;}
  int f1FileSize, f2FileSize;
  do
  {
    f1FileSize = gzread(f1,buf1,MDTK_FILECMP_BUFFER_SIZE);
    f2FileSize = gzread(f2,buf2,MDTK_FILECMP_BUFFER_SIZE);
    REQUIRE(f1FileSize != -1 && f2FileSize != -1);
    if (f1FileSize != f2FileSize) {gzclose(f1);gzclose(f2);return false;}
    if (f1FileSize == 0 && f2FileSize == 0) break;
    if (std::memcmp(buf1,buf2,f1FileSize)) {gzclose(f1);gzclose(f2);return false;}
  }while(1);
  REQUIRE(f1FileSize != -1 && f2FileSize != -1);
  gzclose(f1);
  gzclose(f2);  
  return true;
}  


}


