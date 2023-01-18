# How to build this Docker image:
#     docker build --tag einsteintoolkit/cactus-amrex .
#     docker push einsteintoolkit/cactus-amrex

FROM ubuntu:20.04

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
# TODO: Add CVS to check out Lorene
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        ca-certificates \
        g++ \
        gfortran \
        git \
        hdf5-tools \
        libboost-all-dev \
        libfftw3-dev \
        libgit2-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libudev-dev \
        perl \
        pkg-config \
        python2 \
        python3 \
        python3-pip \
        python3-requests \
        rsync \
        subversion \
        vim \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# jinja2 >= 3.1 removes jinja2.Markup, causing failures:
# https://github.com/bokeh/bokeh/pull/11174
RUN pip install jinja2==3.0.3 && \
    pip install bokeh==2.0.1 && \
    pip install matplotlib && \
    pip install requests && \
    pip install pygit2==1.0.3

# Install cmake
# We need a modern cmake to build AMReX
RUN mkdir dist && \
    (cd dist && \
    wget https://github.com/Kitware/CMake/releases/download/v3.21.2/cmake-3.21.2-linux-x86_64.tar.gz && \
    tar xzf cmake-3.21.2-linux-x86_64.tar.gz && \
    rsync -r cmake-3.21.2-linux-x86_64/ /usr/local && \
    true) && \
    rm -rf dist

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.8.3.tar.gz && \
    tar xzf v2.8.3.tar.gz && \
    cd ADIOS2-2.8.3 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install NSIMD
# NSIMD allows writing explicitly SIMD-vectorized code
# Note: This assumes that the system has x86_64 CPUs
RUN mkdir src && \
    (cd src && \
    wget https://github.com/agenium-scale/nsimd/archive/refs/tags/v3.0.1.tar.gz && \
    tar xzf v3.0.1.tar.gz && \
    cd nsimd-3.0.1 && \
    mkdir build && \
    cd build && \
    cmake \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_C_COMPILER=gcc \
        -DCMAKE_CXX_COMPILER=g++ \
        -Dsimd=SSE2 \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - Depends on ADIOS2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.14.5.tar.gz && \
    tar xzf 0.14.5.tar.gz && \
    cd openPMD-api-0.14.5 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install Silo
# openPMD-api defines a standard for laying out AMR data in a file
RUN mkdir src && \
    (cd src && \
    wget https://github.com/LLNL/Silo/releases/download/v4.11/silo-4.11.tar.gz && \
    tar xzf silo-4.11.tar.gz && \
    cd silo-4.11 && \
    mkdir build && \
    cd build && \
    ../configure \
        --disable-fortran \
        --enable-optimization \
        --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib \
        --prefix=/usr/local && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install ssht
# ssht provides spin-weighted spherical harmonicsa
# TODO: Use "cmake install" instead of "cp"
# RUN mkdir src && \
#     (git clone https://github.com/astro-informatics/ssht.git && \
#     mkdir ssht/build && \
#     cd ssht/build && \
#     cmake .. && \
#     make install) && \
#     rm -fr src
RUN mkdir src && \
    (cd src && \
    wget https://github.com/astro-informatics/ssht/archive/v1.5.1.tar.gz && \
    tar xzf v1.5.1.tar.gz && \
    cd ssht-1.5.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install yaml-cpp
# yaml-cpp allows reading and writing YAML files
RUN mkdir src && \
    (cd src && \
    wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz && \
    tar xzf yaml-cpp-0.6.3.tar.gz && \
    cd yaml-cpp-yaml-cpp-0.6.3 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install AMReX
# AMReX provides adaptive mesh refinement
# - Enable Fortran for `docker/Dockerfile`
# - Keep AMReX source tree around for debugging
# - Install this last because it changes most often
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/22.11.tar.gz && \
    tar xzf 22.11.tar.gz && \
    cd amrex-22.11 && \
    mkdir build && \
    cd build && \
    cmake \
        -DCMAKE_BUILD_TYPE=Debug \
        -DAMReX_OMP=ON \
        -DAMReX_PARTICLES=ON \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true)
