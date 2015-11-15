/* 
   mdtrajsim (molecular dynamics trajectory simulator)

   Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012,
   2013, 2015 Oleksandr Yermolenko <oleksandr.yermolenko@gmail.com>

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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <algorithm>

#ifdef MPIBATCH
#include "mpi.h"
#endif

#include <mdtk/SimLoop.hpp>
#include <mdtk/SimLoopSaver.hpp>
#include <mdtk/SnapshotList.hpp>

#include "../common.h"

#include "tools.hpp"

using namespace mdtk;

class CustomSimLoop : public SimLoop
{
public:
  SnapshotList snapshotList;
  CustomSimLoop();
  bool isItTimeToSave(Float interval);
  void doBeforeIteration();
  void doAfterIteration();
};

CustomSimLoop::CustomSimLoop()
  :SimLoop(),snapshotList()
{
  verboseTrace = true;
}

bool
CustomSimLoop::isItTimeToSave(Float interval)
{
  return (simTime == 0.0 || 
          int(simTime/interval) != int((simTime - dt_prev)/interval));
}

void
CustomSimLoop::doBeforeIteration()
{
  if (simTime < 2.0*ps) simTimeSaveTrajInterval = 1.0*ps;
  else simTimeSaveTrajInterval = 5.0*ps;

  if (simTime < 2.0*ps)
  {
    if (isItTimeToSave(5e-16*5*4))
      snapshotList.getSnapshot(*this);
  }
  else
  {
    if (simTime < 4.0*ps)
    {
      if (isItTimeToSave(5e-16*25*4))
        snapshotList.getSnapshot(*this);
    }
    else
    {
      if (isItTimeToSave(5e-16*50*4))
        snapshotList.getSnapshot(*this);
    }
  }

  if (iteration%iterationFlushStateInterval == 0 && iteration != 0)
  {
    snapshotList.writestate();
  }
}

void
CustomSimLoop::doAfterIteration()
{
  return;
}

bool
isAlreadyFinished()
{
  {
    if (yaatk::exists("completed.ok")) return true;
  }
  {
    if (yaatk::exists("completed.error")) return true;
  }
  {
    std::ifstream ifinal("in.mde.after_crash");
    if (ifinal)
      return true;
    else
      ifinal.close();
  }
  return false;
}

int runTraj(std::string inputFilesId = "")
{
  if (isAlreadyFinished() || isLockedByOthers())
    return 0;

  if (!yaatk::DataState::isClean())
    return 1;

  TRACE(mdtk::buildID);

  try
  {
    CustomSimLoop mdloop;
    mdtk::SimLoopSaver mds(mdloop);

    setupPotentials(mdloop);

    if (inputFilesId != "")
    {
      if (mds.load(inputFilesId) & mdtk::SimLoopSaver::LOADED_R)
      {
        PRINT("Input files loaded\n");
      }
      else
      {
        REQUIRE(yaatk::exists(inputFilesId.c_str()));
        yaatk::text_ifstream fi(inputFilesId.c_str());
        mdloop.loadFromStream(fi);
        fi.close();
      }

      mds.write();
    }
    else
    {
      if (mds.listIterations().size() > 0)
      {
        mds.loadIterationLatest();
      }
      else
      {
        std::vector<std::string> ids = mds.listIds();
        if (ids.size() > 0)
          mds.load(ids[0]);
      }
    }

    if (yaatk::exists("snapshots.conf"))
      mdloop.snapshotList.loadstate();

    mdloop.iterationFlushStateInterval = 1000;
    mdloop.execute();
    mds.write();
    mdloop.snapshotList.writestate();

    if (mdloop.simTime >= mdloop.simTimeFinal) // is simulation really finished ?
    {
      yaatk::text_ofstream fo2("completed.ok");
      fo2.close();

      mds.removeIterations();
    }
  }
  catch(mdtk::Exception& e)
  {
    std::cerr << "MDTK exception in the main thread: "
              << e.what() << std::endl;
    std::cerr << "Probably wrong input file format or no input files" << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "Unknown exception in the main thread." << std::endl;
    return 2;
  }

  return 0;
}

int runImpactSequence(
  std::string inputFilesId,
  int startFromImpact,
  int finalImpact = std::numeric_limits<int>::max())
{
  int retcode = 0;

  CustomSimLoop mdloop;
  mdtk::SimLoopSaver mds(mdloop);

  int loadRetCode = mds.load(inputFilesId);
  REQUIRE(loadRetCode & mdtk::SimLoopSaver::LOADED_R);

  double bombX, bombY;

  yaatk::binary_ifstream ionpos_bin("ionpos.bin");
  REQUIRE(ionpos_bin);
  int recordLength = ionpos_bin.getDataLength();
  REQUIRE(recordLength % (sizeof(bombX) + sizeof(bombY)) == 0);
  int impactCount = recordLength/(sizeof(bombX) + sizeof(bombY));
  TRACE(impactCount);

  if (impactCount > 0 && impactCount-1 < finalImpact)
    finalImpact = impactCount-1;
  TRACE(finalImpact);

  for(int impact = 0; impact <= finalImpact; ++impact)
  {
    YAATK_BIN_READ(ionpos_bin,bombX);
    REQUIRE(ionpos_bin);
    YAATK_BIN_READ(ionpos_bin,bombY);
    REQUIRE(ionpos_bin);

    Atom& projectile = mdloop.atoms.back();
    projectile.coords.x = bombX;
    projectile.coords.y = bombY;

    if (impact < startFromImpact)
      continue;

    {
      char trajDirName[1024];
      sprintf(trajDirName,"%08d",impact);

      {
        yaatk::ChDir cd(trajDirName,false);

        if (isAlreadyFinished())
          continue;
        else
          placeMyLock();

#ifdef MPIBATCH
        cout << "Process #" << comm_rank() << " of " << comm_size() << " tries ";
#else
        cout << "I am trying ";
#endif
        cout << "to take trajectory " << yaatk::getcwd() << " : ";

        if (isLockedByOthers())
        {
          cout << "LOCKED" << std::endl;
          removeMyLock();
          continue;
        }

        if (!yaatk::DataState::isClean())
        {
          cout << "MAY BE CORRUPTED" << std::endl;
          removeMyLock();
          continue;
        }

        cout << "OK" << std::endl;
      }

      {
        yaatk::ChDir cd(trajDirName);

        {
          yaatk::DataState ds;

          yaatk::binary_ofstream ionpos_used_bin("ionpos_used.bin");
          yaatk::text_ofstream ionpos_used_txt("ionpos_used.txt");

          YAATK_BIN_WRITE(ionpos_used_bin,bombX);
          YAATK_BIN_WRITE(ionpos_used_bin,bombY);

          ionpos_used_txt << bombX << " "
                          << bombY << "\n";


          ionpos_used_bin.close();
          ionpos_used_txt.close();
        }

/*
  mdloop.iteration = 0;
  mdloop.simTime = 0.0*ps;
  mdloop.simTimeFinal = 10.0*ps;
  mdloop.simTimeSaveTrajInterval = 0.1*ps;
*/

        mdloop.simTimeSaveTrajInterval = 1e6*ps;

        mds.write();

        retcode |= runTraj();

        mds.removeIterations(false,true);

        removeMyLock();
      }
      {
        yaatk::ChDir cd(trajDirName);
        if (isAlreadyFinished())
        {
          placeMyLock();
          if (!isLockedByOthers())
          {
            shrinkLogFiles();
            removeXVAfiles();
          }
          removeMyLock();
        }
      }
    }
  }

  ionpos_bin.close();

  return retcode;
}

int
main(int argc, char *argv[])
{
#ifdef MPIBATCH
  MPI_TEST_SUCCESS(MPI_Init(&argc,&argv));

  std::ostringstream oss_stdout_redir_name;
  oss_stdout_redir_name << yaatk::getDateTime_YYYYMMDD_HHMMSS() << "-";
  oss_stdout_redir_name << "stdout.txt-mpibatch";
  oss_stdout_redir_name << getProcessID();
  yaatk::StreamToFileRedirect* cout_redir
    = new yaatk::StreamToFileRedirect(std::cout,oss_stdout_redir_name.str());

  std::ostringstream oss_stderr_redir_name;
  oss_stderr_redir_name << yaatk::getDateTime_YYYYMMDD_HHMMSS() << "-";
  oss_stderr_redir_name << "stderr.txt-mpibatch";
  oss_stderr_redir_name << getProcessID();
  yaatk::StreamToFileRedirect* cerr_redir
    = new yaatk::StreamToFileRedirect(std::cerr,oss_stderr_redir_name.str());

  if (comm_rank() == 0)
  {
    std::cout << "mdtrajsim (Molecular dynamics trajectory simulator), mpibatch version ";
    mdtk::release_info.print();
    std::cout << "-------------------------" << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  std::cout << "MPI process #" << comm_rank() << " of " << comm_size() << " started. "
            << "It's name is \"" << comm_name() << "\"" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  SleepForSeconds(10*comm_rank());
#endif

  myLockFilename += getProcessID();

  int startFromImpact = 0;
  int maxFinalImpact = std::numeric_limits<int>::max();
#ifdef MPIBATCH
  bool batch = true;
#else
  bool batch = false;
#endif
  bool commonUsage = false; // by default perform experiment-specific simulation

  std::string inputFilesId = "base";

  if (yaatk::exists("mdtrajsim.option.batch"))
  {
    commonUsage = false;
    batch = true;
  }

  for(int argi = 1; argi < argc; ++argi)
  {
    if ((argv[argi][0] != '-'))
    {
      inputFilesId = argv[argi];
      argi++;
      while (argi < argc && argv[argi][0] != '-')
      {
        argi++;
      }
      if (argi == argc)
        break;
    }

    if (yaatk::isOption(argv[argi],"common-usage",'c'))
    {
      commonUsage = true;
    }

    if (yaatk::isOption(argv[argi],"start-from-impact",'i'))
    {
      commonUsage = false;

      argi++;

      if (!(argi < argc))
      {
        std::cerr << "You should specify index of the first impact to simulate, e.g. --start-from-impact 0\n";
        return -1;
      }
      std::istringstream iss(argv[argi]);
      iss >> startFromImpact;
      if (!(startFromImpact >= 0))
      {
        std::cerr << "Wrong starting impact number\n";
        return -1;
      }
    }

    if (yaatk::isOption(argv[argi],"max-final-impact",'m'))
    {
      commonUsage = false;

      argi++;

      if (!(argi < argc))
      {
        std::cerr << "You should specify maximum index of the last impact to simulate, e.g. --max-final-impact 100\n";
        return -1;
      }
      std::istringstream iss(argv[argi]);
      iss >> maxFinalImpact;
      if (!(maxFinalImpact >= 0))
      {
        std::cerr << "Wrong maximum index of the last impact\n";
        return -1;
      }
    }

    if (yaatk::isOption(argv[argi],"batch",'b'))
    {
      commonUsage = false;
      batch = true;
    }

    if (yaatk::isOption(argv[argi],"version"))
    {
      std::cout << "mdtrajsim (Molecular dynamics trajectory simulator) ";
      mdtk::release_info.print();
      return 0;
    }

    if (yaatk::isOption(argv[argi],"help",'h'))
    {
      std::cout << "\
Usage: mdtrajsim [OPTION]... [FILE]\n\
Simulates molecular dynamics trajectory described by the files with FILE id\n\
\n\
Common options:\n\
      -c, --common-usage           force common usage\n\
      -h, --help                   display this help and exit\n\
      --version                    output version information and exit\n\
Experiment-specific options:\n\
      -i, --start-from-impact <i>  start simulation from the i-th impact\n\
      -m, --max-final-impact <m>   stop simulation after finishing m-th impact\n\
      -b, --batch                  perform batch execution of experiments in the current directory\n\
\n\
Report bugs to <oleksandr.yermolenko@gmail.com>\n\
";
      return 0;
    }
  }

  int retcode = 0;

  if (commonUsage)
  {
    PRINT("Performing simple simulation.\n");
    TRACE(inputFilesId);
    retcode = runTraj(inputFilesId);
  }
  else
  {
    PRINT("Peforming simulation of impact sequence(s).\n");
    TRACE(startFromImpact);
    TRACE(maxFinalImpact);
    TRACE(batch);
    TRACE(inputFilesId);

    try
    {
      if (batch)
      {
        for(int finalImpact = 0; finalImpact <= maxFinalImpact; finalImpact+=1)
        {
          std::vector<std::string> expDirs
            = yaatk::listDirectories(yaatk::getcwd());
          for(size_t i = 0; i < expDirs.size(); ++i)
          {
            const std::string& expname = expDirs[i];
            if (expname.find("_on_PE_") != std::string::npos)
            {
              TRACE(expname);
              yaatk::ChDir cd(expname, false);
              std::vector<std::string> haloDirs
                = yaatk::listDirectories(yaatk::getcwd());
              for(size_t j = 0; j < haloDirs.size(); ++j)
              {
                const std::string& haloname = haloDirs[j];
                if (haloname.find("halo-") != std::string::npos)
                {
                  TRACE(haloname);
                  yaatk::ChDir cd(haloname, false);
                  retcode |= runImpactSequence("../base", startFromImpact, finalImpact);
                }
              }
            }
          }
        }
      }
      else
        retcode |= runImpactSequence(inputFilesId, startFromImpact, maxFinalImpact);
    }
    catch(mdtk::Exception& e)
    {
      std::cerr << "MDTK exception in the main thread: "
                << e.what() << std::endl;
      std::cerr << "Probably wrong input file format or no input files" << std::endl;
      retcode = 1;
    }
    catch(...)
    {
      std::cerr << "Unknown exception in the main thread." << std::endl;
      retcode = 2;
    }
  }

#ifdef MPIBATCH
  if (cout_redir != NULL)
    delete cout_redir;
  if (cerr_redir != NULL)
    delete cerr_redir;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_TEST_SUCCESS(MPI_Finalize());
#endif

  return retcode;
}

