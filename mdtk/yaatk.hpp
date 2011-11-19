/*
   Yet another auxiliary toolkit (header file).

   Copyright (C) 2003, 2005, 2006, 2009, 2010, 2011 Oleksandr Yermolenko
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

#include <cstdio>

#include <sstream>


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
  {                                                                     \
  if (SMODE == YAATK_FSTREAM_TEXT) {FSTREAM_INST >> VAR_INST;} else YAATK_BIN_READ (FSTREAM_INST, VAR_INST); \
  if (FSTREAM_INST.fail()) throw mdtk::Exception("Error in reading variable "#VAR_INST); \
  }

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

inline
std::string getcwd()
{
  const int maxlen = 1000;
  char dir[maxlen];
  char* getcwd_retval;
#ifdef __WIN32__
  getcwd_retval = ::_getcwd(dir,maxlen);
#else
  getcwd_retval = ::getcwd(dir,maxlen);
#endif
  REQUIRE(getcwd_retval);
  return std::string(dir);
}

inline
void remove(const char *name)
{
  std::remove(name);
}

inline
void rename(const char *oldname, const char *newname)
{
  std::rename(oldname,newname);
}

#define DIR_DELIMIT_CHAR '/'
#define DIR_DELIMIT_STR "/"


  class Stream : public std::stringstream
  {
    int zipMe();
    int unZipMe();
    int zipMe_internal();
    int unZipMe_internal();
    std::string filename;
    bool output;
    bool opened;
  public:
    bool isOpened() {return opened;}
    struct ZipInvokeInfo
    {
      std::string command;
      std::string extension;
      ZipInvokeInfo(std::string c, std::string e)
	:command(c),extension(e){}
    };
    static std::vector<ZipInvokeInfo> zipInvokeInfoList;
    static std::vector<ZipInvokeInfo> initZipInvokeInfoList();
    std::string getZippedExt();
    std::string getFileName() {return filename;}
    void guessZipTypeByExtension();
    void guessZipTypeByPresence();
    void setFileName(std::string fname) {filename = fname;}
    void setOutputMode(bool om) {output = om;}
    std::string getZippedFileName() {return filename+getZippedExt();}
    static ZipInvokeInfo zipInvokeInfoGlobal;
    ZipInvokeInfo zipInvokeInfo;
    Stream(std::string fname,bool isOutput,bool isBinary);
    virtual ~Stream() { close(); }
    void open();
    void close();
  };



  class binary_fstream : public Stream
  {
  public:
    binary_fstream(std::string fname,bool isOutput)
      :Stream(fname,isOutput,true) {}
    virtual ~binary_fstream() {}
  };

  class binary_ifstream : public binary_fstream
  {
  public:
    binary_ifstream(std::string fname)
      :binary_fstream(fname,false) {}
    virtual ~binary_ifstream() {}
  };

  class binary_ofstream : public binary_fstream
  {
  public:
    binary_ofstream(std::string fname)
      :binary_fstream(fname,true) {}
    virtual ~binary_ofstream() {}
  };


  class text_fstream : public Stream
  {
  public:
    text_fstream(std::string fname,bool isOutput)
      :Stream(fname,isOutput,false) {}
    virtual ~text_fstream() {}
  };

  class text_ifstream : public text_fstream
  {
  public:
    text_ifstream(std::string fname)
      :text_fstream(fname,false) {}
    virtual ~text_ifstream() {}
  };

  class text_ofstream : public text_fstream
  {
  public:
    text_ofstream(std::string fname)
      :text_fstream(fname,true) {}
    virtual ~text_ofstream() {}
  };

inline
bool exists(std::string filename)
{
  bool retval = false;
  try
  {
    yaatk::binary_ifstream fi(filename);
    if (fi.isOpened()) {fi.close(); retval = true;}
  }catch (...) {;}
  return retval;
}

std::string extractDir(std::string trajNameFinal);
std::string extractLastItem(std::string trajNameFinal);

bool
isIdentical(const std::string& file1,const std::string& file2);

}

#endif



