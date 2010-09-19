/*
   Yet another auxiliary toolkit (header file).

   Copyright (C) 2003, 2005, 2006, 2009 Oleksandr Yermolenko
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

#ifndef yaatk_hpp
#define yaatk_hpp

#include <iostream>
#include <iomanip>
#include <vector>
#include <mdtk/Exception.hpp>
#include <mdtk/config.hpp>

#include <stdio.h>


#ifndef __WIN32__
#include <unistd.h>
#include <errno.h>
#define _open open
#define _close close
#define _O_RDONLY O_RDONLY
#define _O_BINARY 0
#define _filelength(HANDLE) lseek(HANDLE,0,SEEK_END)
#endif

#ifdef __WIN32__
  #include <direct.h>
#else
  #include <sys/stat.h>
#endif

#include <fcntl.h>
#ifdef __WIN32__
#include <io.h>
#else
#endif  

enum YAATK_FSTREAM_MODE {YAATK_FSTREAM_TEXT = 0, YAATK_FSTREAM_BIN = 1};

namespace yaatk
{

#define MPI_TEST_SUCCESS(x) \
{ \
  int mpi_err_code = (x);\
  if (mpi_err_code!=MPI_SUCCESS) \
  { \
    int rank; \
    std::cerr << "*** SIMPLE_MPI_ERROR; At" << #x << std::endl << std::flush; \
    int local_rc = MPI_Comm_rank(MPI_COMM_WORLD,&rank); \
    if (local_rc == MPI_SUCCESS) \
      std::cerr << "*** MPI_ERROR; rank = " << rank << ". At" << #x << std::endl << std::flush; \
    else \
      std::cerr << "*** MPI_ERROR; At" << #x << std::endl << std::flush; \
    throw MPI_Exception(mpi_err_code); \
  } \
};

#define REQUIREM(cond,msg) \
{ \
  if (!(cond)) \
  { \
    std::cerr << (msg) << std::endl << std::flush; \
    throw mdtk::Exception(msg); \
  }  \
}

#define REQUIRE(cond) \
{ \
  if (!(cond)) \
  { \
    std::cerr<< "Assertion "  << (#cond) << " FAILED !!!" << std::endl << std::flush; \
    throw mdtk::Exception((#cond)); \
  }  \
}



#define TRACE_EMPH(x) std::cout << #x << " : " << "<<<##  " << (x) << "" << std::endl;


#define TRACE(x) std::cout << #x << " : " << (x) << std::endl;
#define ERRTRACE(x) std::cerr << #x << " : " << (x) << std::endl;
#ifdef MDE_PARALLEL
#define PTRACE(x) if (comm_rank==0) std::cout << #x << " : " << (x) << std::endl; 
#define PTRACE_SIMPLE(x) if (comm_rank==0) std::cout << x; 
#else
#define PTRACE(x) std::cout << #x << " : " << (x) << std::endl; 
#define PTRACE_SIMPLE(x) std::cout << x; 
#endif

#ifdef MDE_PARALLEL
#define PLOG(x) if (comm_rank==0) std::cout << (x); 
#else
#define PLOG(x) std::cout << (x); 
#endif

#ifdef MDE_PARALLEL
#define PVLOG(x) if (comm_rank==0 && verboseTrace) std::cout << (x) << std::flush; 
#else
#define PVLOG(x) if (verboseTrace) std::cout << (x) << std::flush; 
#endif



#define TRACESS(x) << #x << " : " << (x) << std::endl;


#define YAATK_BIN_WRITE(FSTREAM_INST,VAR_INST) \
  FSTREAM_INST.write((char*)&(VAR_INST),sizeof(VAR_INST))

#define YAATK_BIN_READ(FSTREAM_INST,VAR_INST) \
  FSTREAM_INST.read((char*)&(VAR_INST),sizeof(VAR_INST))

#define YAATK_FSTREAM_WRITE(FSTREAM_INST,VAR_INST,SMODE) \
  if (SMODE == YAATK_FSTREAM_TEXT) {FSTREAM_INST << VAR_INST << "\n";} else YAATK_BIN_WRITE(FSTREAM_INST, VAR_INST);

#define YAATK_FSTREAM_WRITE_NONL(FSTREAM_INST,VAR_INST,SMODE) \
  if (SMODE == YAATK_FSTREAM_TEXT) {FSTREAM_INST << VAR_INST;} else YAATK_BIN_WRITE(FSTREAM_INST, VAR_INST);

#define YAATK_FSTREAM_READ(FSTREAM_INST,VAR_INST,SMODE) \
  if (SMODE == YAATK_FSTREAM_TEXT) {FSTREAM_INST >> VAR_INST;} else YAATK_BIN_READ (FSTREAM_INST, VAR_INST);


#define YAATK_FSTREAM_CREATE(FSTREAM_CLASS,FSTREAM_INST,FSTREAM_FILENAME) \
  FSTREAM_CLASS FSTREAM_INST(FSTREAM_FILENAME); \
  REQUIRE(FSTREAM_INST!=0);

#define YAATK_FSTREAM_CREATE_OPT(FSTREAM_CLASS,FSTREAM_INST,FSTREAM_FILENAME,FSTREAM_OPT) \
  FSTREAM_CLASS FSTREAM_INST(FSTREAM_FILENAME,FSTREAM_OPT); \
  REQUIRE(FSTREAM_INST!=0);

#define YAATK_FSTREAM_CLOSE(FSTREAM_INST) \
  REQUIRE(FSTREAM_INST!=0); \
  FSTREAM_INST.close();




inline
void mkdir(const char *name)
{
#ifdef __WIN32__
  ::_mkdir(name);
#else
  ::mkdir(name,S_IRWXU);
#endif  
}

inline
void chdir(const char *name)
{
#ifdef __WIN32__
  ::_chdir(name);
#else
  ::chdir(name);
#endif  
}

#define DIR_DELIMIT_CHAR '/'
#define DIR_DELIMIT_STR "/"


int
zip_stringstream(const char *zipName, std::stringstream& uzs);
void 
unzip_stringstream(const char *zipName, std::stringstream& os);



void unzip_file(const char *zipName, const char* unzipName);

std::string extractDir(std::string trajNameFinal);
std::string extractLastItem(std::string trajNameFinal);


#define YAATK_UNZIP_FILE(FILENAME_TO_UNZIP) \
{ \
  static char param1[1024]; \
  sprintf(param1,"%s."YAATK_ZIP_EXT,FILENAME_TO_UNZIP); \
  yaatk::unzip_file(param1,yaatk::zippedStreams.getZippedNameTMP(FILENAME_TO_UNZIP)); \
  { \
    bool unzipped_ok = false; \
    { \
      char unzipped_name[1024]; \
      sprintf(unzipped_name,"%s",yaatk::zippedStreams.getZippedNameTMP(FILENAME_TO_UNZIP)); \
      int handle; \
      handle = _open(unzipped_name, _O_RDONLY | _O_BINARY); \
      if (handle == -1) TRACE(errno); \
      if ((handle != -1) && (_filelength(handle) > 0)) unzipped_ok = true; \
      _close(handle); \
    }  \
  }  \
}  

#define YAATK_IFSTREAM_CREATE_ZIPPED(FSTREAM_CLASS,FSTREAM_VAR,ZIPPED_FILENAME_TO_OPEN) \
  yaatk::zippedStreams.setZippedNameTMP(ZIPPED_FILENAME_TO_OPEN); \
  YAATK_UNZIP_FILE(ZIPPED_FILENAME_TO_OPEN); \
  YAATK_FSTREAM_CREATE(FSTREAM_CLASS,FSTREAM_VAR,yaatk::zippedStreams.getZippedNameTMP(ZIPPED_FILENAME_TO_OPEN));

#define YAATK_IFSTREAM_CREATE_ZIPPED_OPT(FSTREAM_CLASS,FSTREAM_VAR,ZIPPED_FILENAME_TO_OPEN,FSTREAM_OPT) \
  yaatk::zippedStreams.setZippedNameTMP(ZIPPED_FILENAME_TO_OPEN); \
  YAATK_UNZIP_FILE(ZIPPED_FILENAME_TO_OPEN); \
  YAATK_FSTREAM_CREATE_OPT(FSTREAM_CLASS,FSTREAM_VAR,yaatk::zippedStreams.getZippedNameTMP(ZIPPED_FILENAME_TO_OPEN),FSTREAM_OPT);



#define YAATK_IFSTREAM_CLOSE_ZIPPED(FSTREAM_VAR,ZIPPED_FILENAME_TO_OPEN) \
{ \
  YAATK_FSTREAM_CLOSE(FSTREAM_VAR); \
    { \
      static char file_to_delete[1024]; \
      sprintf(file_to_delete,"%s",yaatk::zippedStreams.getZippedNameTMP(ZIPPED_FILENAME_TO_OPEN)); \
      int removeSuccess = std::remove(file_to_delete); \
      if (removeSuccess != 0) TRACE(file_to_delete); \
      REQUIRE(removeSuccess == 0); \
    }  \
  yaatk::zippedStreams.afterClose(ZIPPED_FILENAME_TO_OPEN); \
}



struct ZippedStream
{
  std::string filename;
  std::string tmpname;
  ZippedStream(std::string filename_, std::string tmpname_)
  :filename(filename_), tmpname(tmpname_)
  {
  }  
};  

#define ZIPPED_TMP_FOLDER "_zipped_tmp"

extern int yaatk_extraID;

class ZippedStreams
{
  std::vector<ZippedStream> streams;
public:  
  void  setZippedNameTMP(std::string filename)
  {
    char tmp_[1024];
    int i;
    for(i = filename.size()-1; i >= 0; i--)
      if (filename.c_str()[i] == '/' || filename.c_str()[i] == '\\')
      { i++;
        break;
      }  
    if (i < 0) i = 0;
//    sprintf(tmp_,ZIPPED_TMP_FOLDER"\\""%s.%d",filename.substr(i,filename.size()-i).c_str(),streams.size());
    int extraID = 0;
#ifdef MDE_PARALLEL
    extraID = mdtk::comm_rank;
    TRACE(extraID);
#endif
    sprintf(tmp_,ZIPPED_TMP_FOLDER"/""%s.%lu.%d",filename.substr(i,filename.size()-i).c_str(),(unsigned long)streams.size(),extraID+yaatk_extraID);
    streams.push_back(ZippedStream(filename,std::string(tmp_)));  
  }  
  const char* getZippedNameTMP(std::string filename)
  {
    using mdtk::Exception;
    size_t i;
    for(i = 0; i < streams.size(); i++)
      if (streams[i].filename == filename) break;
    REQUIRE(i < streams.size());
    return streams[i].tmpname.c_str();
  }
  void afterClose(std::string filename)
  {
    using mdtk::Exception;
    size_t i;
    for(i = 0; i < streams.size(); i++)
      if (streams[i].filename == filename) break;
    REQUIRE(i < streams.size());
    streams.erase(streams.begin()+i);
  }  
  ZippedStreams():streams()
  {
    mkdir(ZIPPED_TMP_FOLDER);
  }  
};

extern ZippedStreams zippedStreams;

bool
isIdentical(const std::string& file1,const std::string& file2);


int zip_file(const char *zipName, const char* unzipName);



#define YAATK_ZIP_EXT "gz"



#define YAATK_ZIP_FILE(FILENAME_TO_ZIP) \
{ \
  static char file_zipped[1024]; \
  sprintf(file_zipped,"%s."YAATK_ZIP_EXT,FILENAME_TO_ZIP); \
  if (!yaatk::zip_file(file_zipped,FILENAME_TO_ZIP)) \
  { \
  int removeSuccess = std::remove(FILENAME_TO_ZIP); \
  if (removeSuccess != 0) TRACE(FILENAME_TO_ZIP); \
  REQUIRE(removeSuccess == 0); \
  } \
}


}

#endif



