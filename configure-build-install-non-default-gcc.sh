#!/bin/sh
#
# Script for performing out-of-source build and installation of MDTK.
# (Example of using non-default compiler).
#
# Copyright (C) 2009 Oleksandr Yermolenko
# <oleksandr.yermolenko@gmail.com>
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.
#

mkdir _build-non-default-gcc
cd _build-non-default-gcc
  CC=gcc-4.4 CXX=g++-4.4 cmake ..
  make
  make install
cd ..
