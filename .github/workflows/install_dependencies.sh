#!/bin/bash

set -euxo pipefail

# Don't install dependencies if they already exist
test -d /usr/local/etc/HAVE-DEPENDENCIES-1 && exit 0

# Install cmake
# We need a modern cmake to build AMReX
mkdir cmake
pushd cmake
wget https://github.com/Kitware/CMake/releases/download/v3.25.2/cmake-3.25.2-linux-x86_64.tar.gz
tar xzf cmake-3.25.2-linux-x86_64.tar.gz
rsync -r cmake-3.25.2-linux-x86_64/ /usr/local
popd

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
mkdir adios2
pushd adios2
wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.8.3.tar.gz
tar xzf v2.8.3.tar.gz
cd ADIOS2-2.8.3
mkdir build
cd build
cmake ..
make -j$(nproc)
sudo make -j$(nproc) install
popd

# Install NSIMD
# NSIMD allows writing explicitly SIMD-vectorized code
# Note: These instructions assume that the system has x86_64 CPUs
mkdir nsimd
pushd nsimd
wget https://github.com/agenium-scale/nsimd/archive/refs/tags/v3.0.1.tar.gz
tar xzf v3.0.1.tar.gz
cd nsimd-3.0.1
mkdir build
cd build
cmake                                           \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo           \
    -DCMAKE_C_COMPILER=gcc                      \
    -DCMAKE_CXX_COMPILER=g++                    \
    -Dsimd=SSE2                                 \
    -DCMAKE_INSTALL_PREFIX=/usr/local           \
    ..
make -j$(nproc)
sudo make -j$(nproc) install
popd

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - Depends on ADIOS2
mkdir openPMD-api
pushd openPMD-api
wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.14.5.tar.gz
tar xzf 0.14.5.tar.gz
cd openPMD-api-0.14.5
mkdir build
cd build
# We need to disable testing because there is a bug in the enclosed `catch2` library
cmake .. -DBUILD_TESTING=OFF
make -j$(nproc)
sudo make -j$(nproc) install
sudo chmod a+rx /usr/local/bin/openpmd-pipe
popd

# Install Silo
# Silo defines a standard for laying out AMR data in a file
mkdir silo
pushd silo
wget https://github.com/LLNL/Silo/releases/download/v4.11/silo-4.11.tar.gz
tar xzf silo-4.11.tar.gz
cd silo-4.11
mkdir build
cd build
../configure                                                                                            \
    --disable-fortran                                                                                   \
    --enable-optimization                                                                               \
    --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib \
    --prefix=/usr/local
make -j$(nproc)
sudo make -j$(nproc) install
popd

# Install yaml-cpp
# yaml-cpp allows reading and writing YAML files
mkdir yaml-cpp
pushd yaml-cpp
wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz
tar xzf yaml-cpp-0.6.3.tar.gz
cd yaml-cpp-yaml-cpp-0.6.3
mkdir build
cd build
cmake ..
make -j$(nproc)
sudo make -j$(nproc) install
popd

# Install AMReX
# AMReX provides adaptive mesh refinement
mkdir amrex
pushd amrex
wget https://github.com/AMReX-Codes/amrex/archive/23.02.tar.gz
tar xzf 23.02.tar.gz
cd amrex-23.02
mkdir build
cd build
cmake                                           \
    -DCMAKE_BUILD_TYPE=Debug                    \
    -DAMReX_OMP=ON                              \
    -DAMReX_PARTICLES=ON                        \
    -DCMAKE_INSTALL_PREFIX=/usr/local           \
    ..
make -j$(nproc)
sudo make -j$(nproc) install
popd

# Remember that we installed the dependencies
mkdir -p /usr/local/etc
: >/usr/local/etc/HAVE-DEPENDENCIES-1
