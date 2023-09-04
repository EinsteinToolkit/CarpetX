# How to build this Docker image:

#     docker build --file carpetx-cuda.dockerfile --tag einsteintoolkit/carpetx:cuda-real64 .
#     docker push einsteintoolkit/carpetx:cuda-real64

#     docker build --build-arg real_precision=real32 --file carpetx-cuda.dockerfile --tag einsteintoolkit/carpetx:cuda-real32 .
#     docker push einsteintoolkit/carpetx:cuda-real32

FROM nvidia/cuda:12.2.0-devel-ubuntu22.04

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        ca-certificates \
        cmake \
        cvs \
        g++ \
        gfortran \
        git \
        hdf5-tools \
        less \
        libblosc-dev \
        libboost-all-dev \
        libbz2-dev \
        libfftw3-dev \
        libgit2-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libudev-dev \
        meson \
        ninja-build \
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

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.9.1.tar.gz && \
    tar xzf v2.9.1.tar.gz && \
    cd ADIOS2-2.9.1 && \
    cmake -B build -G Ninja -S . \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_TESTING=OFF \
        -DADIOS2_BUILD_EXAMPLES=OFF \
        -DADIOS2_USE_Fortran=OFF \
        && \
    cmake --build build && \
    cmake --install build && \
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
    cmake -B build -G Ninja -S . \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -Dsimd=AVX2 \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - depends on ADIOS2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.15.2.tar.gz && \
    tar xzf 0.15.2.tar.gz && \
    cd openPMD-api-0.15.2 && \
    cmake -B build -G Ninja -S . \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_TESTING=OFF \
        -DBUILD_EXAMPLES=OFF \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install RePrimAnd
# RePrimAnd is a physics package for nuclear equations of state
RUN mkdir src && \
    (cd src && \
    wget https://github.com/wokast/RePrimAnd/archive/refs/tags/v1.6.tar.gz && \
    tar xzf v1.6.tar.gz && \
    cd RePrimAnd-1.6 && \
    meson build --buildtype=release --prefix=/usr/local && \
    ninja -C build && \
    ninja -C build install && \
    true) && \
    rm -rf src

# Install Silo
# Silo defines a standard for laying out AMR data in a file
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
        --prefix=/usr/local \
        && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install ssht
# ssht provides spin-weighted spherical harmonics
RUN mkdir src && \
    (cd src && \
    wget https://github.com/astro-informatics/ssht/archive/v1.5.2.tar.gz && \
    tar xzf v1.5.2.tar.gz && \
    cd ssht-1.5.2 && \
    cmake -B build -G Ninja -S . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install yaml-cpp
# yaml-cpp allows reading and writing YAML files
RUN mkdir src && \
    (cd src && \
    wget https://github.com/jbeder/yaml-cpp/archive/yaml-cpp-0.6.3.tar.gz && \
    tar xzf yaml-cpp-0.6.3.tar.gz && \
    cd yaml-cpp-yaml-cpp-0.6.3 && \
    cmake -B build -G Ninja -S .  \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DYAML_CPP_BUILD_TESTS=OFF \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ASDF
# ASDF is an I/O library like HDF5
# - depends on yaml-cpp
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/asdf-cxx/archive/refs/tags/version/7.2.1.tar.gz && \
    tar xzf 7.2.1.tar.gz && \
    cd asdf-cxx-version-7.2.1 && \
    cmake -B build -G Ninja -S . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install SimulationIO
# SimulationIO is an I/O library like HDF5
# - depends on asdf-cxx
# - depends on yaml-cpp
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/SimulationIO/archive/refs/tags/version/9.0.3.tar.gz && \
    tar xzf 9.0.3.tar.gz && \
    cd SimulationIO-version-9.0.3 && \
    cmake -B build -G Ninja -S . -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX=/usr/local && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

ARG real_precision=real64

# Install AMReX
# AMReX provides adaptive mesh refinement
# - Enable Fortran for `docker/Dockerfile`
# - Install this last because it changes most often
# Should we keep the AMReX source tree around for debugging?
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/23.08.tar.gz && \
    tar xzf 23.08.tar.gz && \
    cd amrex-23.08 && \
    case $real_precision in \
        real32) precision=SINGLE;; \
        real64) precision=DOUBLE;; \
        *) exit 1;; \
    esac && \
    cmake -B build -G Ninja -S . \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DAMReX_CUDA_ARCH=7.0 \
        -DAMReX_FORTRAN=OFF \
        -DAMReX_FORTRAN_INTERFACES=OFF \
        -DAMReX_GPU_BACKEND=CUDA \
        -DAMReX_OMP=ON \
        -DAMReX_PARTICLES=ON \
        -DAMReX_PRECISION="$precision" \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src
