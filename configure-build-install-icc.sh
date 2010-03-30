#!/bin/sh
#
# Script for performing out-of-source build and installation of MDTK.
# (Example of using Intel C++ compiler).
#
# Copyright (C) 2009 Oleksandr Yermolenko
# <oleksandr.yermolenko@gmail.com>
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.
#

mkdir _build-icc
cd _build-icc
  source /opt/intel/Compiler/11.1/064/bin/intel64/iccvars_intel64.sh
  LANG=C CC=/opt/intel/Compiler/11.1/064/bin/intel64/icc CXX=/opt/intel/Compiler/11.1/064/bin/intel64/icpc cmake -D BUILD_INTELCC=1 ..
  LANG=C make
  LANG=C make install
cd ..
