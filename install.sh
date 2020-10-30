#!/usr/bin/env bash

# This script will help run compile IDP-denovo

echo Working in the idpdenovo directory:
cd "$(dirname "$0")"
echo $PWD
cd plugins/gmap-2014-12-24/
echo Configuring and compiling gmap installation in:
echo $PWD
./configure --prefix=$PWD/../BIN
make
make check
make install
echo Finished compiling GMAP
cd ../../
echo Compiling idpdenovo in directory:
echo $PWD
make clean
make
