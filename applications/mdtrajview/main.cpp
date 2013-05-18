/*
   mdtrajview (the molecular dynamics trajectory viewer)

   Copyright (C) 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
   2012, 2013 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include "VisBox.hpp"
#include "MainWindow.hpp"

#include "mdtk/Atom.hpp"
#include "mdtk/config.hpp"
#include "mdtk/consts.hpp"
#include "mdtk/SimLoop.hpp"
#include "mdtk/SimLoopSaver.hpp"

#include <fstream>
#include <algorithm>
#include <cstring>

using namespace std;
using namespace mdtk;

void
findIntermediateStates(std::string trajDir,std::vector<std::string>& states)
{
  std::vector<std::string> fileList = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < fileList.size(); ++i)
  {
    const std::string& filename = fileList[i];
    if (filename.find("mde0") == 0 && filename.find(".xva") != std::string::npos)
      states.push_back(filename);
  }

  sort(states.begin(),states.end());
}

void
findBaseFiles(std::string trajDir,std::vector<std::string>& states)
{
  std::vector<std::string> fileList = yaatk::listFiles(yaatk::getcwd());

  for(size_t i = 0; i < fileList.size(); ++i)
  {
    const std::string& filename = fileList[i];
    if (filename.find(".mde") != std::string::npos)
      states.push_back(filename);
  }
}

int main(int argc, char *argv[])
{
  bool instantAnimate = false;
  std::vector<std::string> fileList;
  std::string baseFile;
  bool loadPartialSnapshots = false;
  bool xvasSpecified = false;
  bool basefilesOnly = false;

  bool newFileFormats = true;
  if (yaatk::exists("in.mde") ||
      yaatk::exists("mde_init") ||
      yaatk::exists("simloop.conf") ||
      yaatk::exists("mde_final"))
  {
    newFileFormats = false;
  }

  for(int argi = 1; argi < argc; ++argi)
  {
    if ((argv[argi][0] != '-'))
    {
      baseFile = argv[argi];
      argi++;
      while (argi < argc && argv[argi][0] != '-')
      {
        fileList.push_back(argv[argi]);
        xvasSpecified = true;
        argi++;
      }
      if (argi == argc)
        break;
    }

    if (yaatk::isOption(argv[argi],"partial-snapshots",'s'))
    {
      loadPartialSnapshots = true;
    }

    if (yaatk::isOption(argv[argi],"instant-animation",'a'))
    {
      instantAnimate = true;
    }

    if (yaatk::isOption(argv[argi],"legacy-file-formats",'l'))
    {
      newFileFormats = false;
    }

    if (yaatk::isOption(argv[argi],"new-file-formats",'n'))
    {
      newFileFormats = true;
    }

    if (yaatk::isOption(argv[argi],"version"))
    {
      std::cout << "mdtrajview (Molecular dynamics trajectory viewer) ";
      mdtk::release_info.print();
      return 0;
    }

    if (yaatk::isOption(argv[argi],"help",'h'))
    {
      std::cout << "\
Usage: mdtrajview [OPTION] [base file with complete info] [XVA files with incomplete info]... \n\
Visualizes molecular dynamics trajectory previously simulated by mdtrajsim.\n\
Operates on results of mdtrajsim's run in the current directory.\n\
\n\
      -a, --instant-animation    start animation immediately\n\
      -s, --partial-snapshots    try to load partial snapshots file (snapshots.conf)\n\
      -l, --legacy-file-formats  use legacy input file formats (heuristics is by default)\n\
      -n, --new-file-formats     use new input file formats (heuristics is by default)\n\
      -h, --help                 display this help and exit\n\
      --version                  output version information and exit\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
      return 0;
    }
  }

  if (!newFileFormats)
  {

  if (!xvasSpecified && baseFile == "")
  {
    std::vector<std::string> basefiles;
    basefiles.push_back("in.mde");
    findBaseFiles("./",basefiles);
    basefiles.push_back("mde_init");
    basefiles.push_back("simloop.conf");
    basefiles.push_back("mde_final");
    for(size_t i = 0; i < basefiles.size(); i++)
    {
      if (yaatk::exists(basefiles[i]))
      {baseFile = basefiles[i];break;}
      else
      {baseFile = "";}
    }
    findIntermediateStates("./",fileList);
  }

  if (baseFile == "")
  {
    std::cerr << "Can't find output of successful mdtrajsim run in the current directory. But you can specify files directly. Run with -h option for details.\n";
    return 1;
  }

  if (fileList.size() > 0 && fileList[0].find(".xva") == std::string::npos)
  {
    basefilesOnly = true;
    std::vector<std::string> fileList_tmp = fileList;
    fileList.clear();
    fileList.push_back(baseFile);
    for(size_t i = 0; i < fileList_tmp.size(); ++i)
      fileList.push_back(fileList_tmp[i]);
  }

  TRACE(loadPartialSnapshots);
  TRACE(xvasSpecified);
  TRACE(baseFile);
  for(size_t i = 0; i < fileList.size(); i++)
  {
    TRACE(fileList[i]);
  }

  std::vector<std::string> fileList_tmp = fileList;
  if (fileList.size() > 5000)
  {
    fileList.clear();
    for(size_t i = 0; i < fileList_tmp.size(); ++i)
    {
      if (i % 10 == 0) fileList.push_back(fileList_tmp[i]);
    }
  }

//  const size_t maxsize = 2000;
//  if (fileList.size() > maxsize) fileList.resize(maxsize);
  }

  xmde::VisBox avb(15,35,500,500);
  avb.set_non_modal();

  if (!newFileFormats)
  {
  if (!basefilesOnly)
    avb.loadDataFromFiles(baseFile,fileList,loadPartialSnapshots);
  else
    avb.loadDataFromMDLoopStates(fileList);
  }
  else
  {
    std::vector<std::string> states_ng;
    if (baseFile == "")
      states_ng = mdtk::SimLoopSaver::listIds();
    else
      states_ng.push_back(baseFile);

    std::vector<std::string> xvas;
    if (fileList.size() == 0)
      findIntermediateStates("./",xvas);
    if (fileList.size() > 0 && fileList[0].find(".xva") == std::string::npos)
    {
      basefilesOnly = true;
      for(size_t i = 0; i < fileList.size(); ++i)
        states_ng.push_back(fileList[i]);
    }
    else
    {
      for(size_t i = 0; i < fileList.size(); ++i)
        xvas.push_back(fileList[i]);
    }

    avb.loadDataFromFilesOfNewFileFormat(states_ng,xvas,loadPartialSnapshots);
  }

  xmde::MainWindow w(&avb, instantAnimate);
  MainWindow_GlobalPtr = &w;

  Fl::run();

  return 0;
}

