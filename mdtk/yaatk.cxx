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
  char cmd[2000];sprintf(cmd,"gzip -c >%s",zipName);
#ifndef __WIN32__
  FILE* zipped   = popen(cmd,"w");
#else
  FILE* zipped   = _popen(cmd,"wb");
#endif
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while((uzs.read(buf,MDTK_GZ_BUFFER_SIZE),unzippedFileSize = uzs.gcount()) > 0)
  {
    REQUIRE(unzippedFileSize != -1);
    int bytesWritten    = fwrite(buf,1,unzippedFileSize,zipped);
    REQUIRE(unzippedFileSize == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
#ifndef __WIN32__
  pclose(zipped);
#else
  _pclose(zipped);
#endif
  return 0;
}

Stream::ZipInvokeInfo Stream::zipInvokeInfoGlobal=ZipInvokeInfo("gzip_internal",".gz");

std::vector<Stream::ZipInvokeInfo> Stream::zipInvokeInfoList = Stream::initZipInvokeInfoList();

std::vector<Stream::ZipInvokeInfo> 
Stream::initZipInvokeInfoList()
{
  std::vector<ZipInvokeInfo> v;
  v.push_back(ZipInvokeInfo("gzip_internal",".gz"));
  v.push_back(ZipInvokeInfo("gzip",".gz"));
  v.push_back(ZipInvokeInfo("bzip2",".bz2"));
  v.push_back(ZipInvokeInfo("xz",".xz"));
  return v;
}

std::string
Stream::getZippedExt()
{
  return zipInvokeInfo.extension;
}

void
Stream::guessZipType()
{
  for(size_t i = 0; i < zipInvokeInfoList.size(); i++)
  {
    const ZipInvokeInfo& z = zipInvokeInfoList[i];
    if (z.command != "nozip" && 
	filename.find(z.extension) == filename.size()-z.extension.size())
      {
	zipInvokeInfo = z;
	filename = filename.substr(0,filename.size()-z.extension.size());
	return;
      }
  }
  zipInvokeInfo = ZipInvokeInfo("nozip","");
}

Stream::Stream(std::string fname,bool isOutput,bool isBinary)
      :std::stringstream(
isBinary?
(std::stringstream::in | std::stringstream::out | std::stringstream::binary)
:
(std::stringstream::in | std::stringstream::out)
)
      ,filename(fname),output(isOutput),opened(false),zipInvokeInfo(zipInvokeInfoGlobal)
{
  if (!output) guessZipType();
  TRACE(zipInvokeInfo.command);
  open();
}

void Stream::open()
{
  if (!opened)
  {
    if (!output) 
      opened = !unZipMe();
    else
      opened = true;
  }
}

void
Stream::close()
{
  if (opened)
  {
    if (output) 
      opened = zipMe(); 
    else
      opened = false;
  }
}

int 
Stream::zipMe()
{
  if (zipInvokeInfo.command=="gzip_internal") return zipMe_internal();
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  FILE* zipped;
  if (zipInvokeInfo.command!="nozip")
  {
    char cmd[2000];sprintf(cmd,"%s -c >%s",zipInvokeInfo.command.c_str(),getZippedFileName().c_str());
#ifndef __WIN32__
    zipped   = popen(cmd,"w");
#else
    zipped   = _popen(cmd,"wb");
#endif
  }
  else
    zipped   = fopen(getZippedFileName().c_str(),"wb");
    
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while((read(buf,MDTK_GZ_BUFFER_SIZE),unzippedFileSize = gcount()) > 0)
  {
    REQUIRE(unzippedFileSize != -1);
    int bytesWritten    = fwrite(buf,1,unzippedFileSize,zipped);
    REQUIRE(unzippedFileSize == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
  if (zipInvokeInfo.command!="nozip")
  {
#ifndef __WIN32__
    int pclose_status = pclose(zipped);
#else
    int pclose_status = _pclose(zipped);
#endif
    REQUIRE(!pclose_status);
  }
  else
  {
    int fclose_status = fclose(zipped);
    REQUIRE(!fclose_status);
  }
  return 0;
}

int
Stream::unZipMe()
{
  if (zipInvokeInfo.command=="gzip_internal") return unZipMe_internal();
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  FILE* zipped;
  if (zipInvokeInfo.command!="nozip")
  {
    char cmd[2000];sprintf(cmd,"%s -dc %s",zipInvokeInfo.command.c_str(),getZippedFileName().c_str());
#ifndef __WIN32__
    zipped   = popen(cmd,"r");
#else
    zipped   = _popen(cmd,"rb");
#endif
  }
  else
    zipped   = fopen(getZippedFileName().c_str(),"rb");

  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while ((unzippedFileSize = fread(buf,1,MDTK_GZ_BUFFER_SIZE,zipped)) > 0)
  {
    write(buf,unzippedFileSize);
  }
  
  REQUIRE(unzippedFileSize == 0);
  if (zipInvokeInfo.command!="nozip")
  {
#ifndef __WIN32__
    int pclose_status = pclose(zipped);
#else
    int pclose_status = _pclose(zipped);
#endif
  REQUIRE(!pclose_status);
  }
  else
  {
    int fclose_status = fclose(zipped);
    REQUIRE(!fclose_status);
  }

  return 0;
}

int
Stream::zipMe_internal()
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  gzFile   zipped   = gzopen(getZippedFileName().c_str(),"wb");
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while((read(buf,MDTK_GZ_BUFFER_SIZE),unzippedFileSize = gcount()) > 0)
  {
    REQUIRE(unzippedFileSize != -1);
    int bytesWritten    = gzwrite(zipped,buf,unzippedFileSize);
    REQUIRE(unzippedFileSize == bytesWritten);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(zipped);
  return 0;
}

int
Stream::unZipMe_internal()
{
  using mdtk::Exception;
  char buf[MDTK_GZ_BUFFER_SIZE];
  gzFile zipped   = gzopen(getZippedFileName().c_str(),"rb");
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while ((unzippedFileSize = gzread(zipped,buf,MDTK_GZ_BUFFER_SIZE)) > 0)
  {
    write(buf,unzippedFileSize);
  }
  REQUIRE(unzippedFileSize == 0);
  gzclose(zipped);
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
  char cmd[2000];sprintf(cmd,"gzip -dc %s",zipName);
#ifndef __WIN32__
  FILE* zipped   = popen(cmd,"r");
#else
  FILE* zipped   = _popen(cmd,"rb");
#endif
  REQUIRE(zipped != 0);
  int unzippedFileSize;
  while ((unzippedFileSize = fread(buf,1,MDTK_GZ_BUFFER_SIZE,zipped)) > 0)
  {
    os.write(buf,unzippedFileSize);
  }
  
  REQUIRE(unzippedFileSize == 0);
#ifndef __WIN32__
  pclose(zipped);
#else
  _pclose(zipped);
#endif
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


