# How to build this Docker image:

#     docker build --file carpetx-rocm.dockerfile --tag einsteintoolkit/carpetx:rocm-real64 .
#     docker push einsteintoolkit/carpetx:rocm-real64

#     docker build --build-arg real_precision=real32 --file carpetx-rocm.dockerfile --tag einsteintoolkit/carpetx:rocm-real32 .
#     docker push einsteintoolkit/carpetx:rocm-real32

FROM rocm/dev-ubuntu-22.04:5.6

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        build-essential \
        ca-certificates \
        cvs \
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
        rocprim-dev \
        rocrand-dev \
        rsync \
        subversion \
        vim \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# # jinja2 >= 3.1 removes jinja2.Markup, causing failures:
# # https://github.com/bokeh/bokeh/pull/11174
# RUN pip install jinja2==3.0.3 && \
#     pip install bokeh==2.0.1 && \
#     pip install matplotlib && \
#     pip install requests && \
#     pip install pygit2==1.0.3

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
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.9.0.tar.gz && \
    tar xzf v2.9.0.tar.gz && \
    cd ADIOS2-2.9.0 && \
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
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.15.1.tar.gz && \
    tar xzf 0.15.1.tar.gz && \
    cd openPMD-api-0.15.1 && \
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
# ssht provides spin-weighted spherical harmonics
RUN mkdir src && \
    (cd src && \
    wget https://github.com/astro-informatics/ssht/archive/v1.5.2.tar.gz && \
    tar xzf v1.5.2.tar.gz && \
    cd ssht-1.5.2 && \
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

ARG real_precision=real64

# Install AMReX
# AMReX provides adaptive mesh refinement
# - Enable Fortran for `docker/Dockerfile`
# - Install this last because it changes most often
# Should we  keep AMReX source tree around for debugging?
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/23.08.tar.gz && \
    tar xzf 23.08.tar.gz && \
    cd amrex-23.08 && \
    mkdir build && \
    cd build && \
    case $real_precision in \
        real32) precision=SINGLE;; \
        real64) precision=DOUBLE;; \
        *) exit 1;; \
    esac && \
    cmake \
        -DAMReX_GPU_BACKEND=HIP \
        -DAMReX_AMD_ARCH=gfx90a \
        -DAMReX_OMP=OFF \
        -DAMReX_PARTICLES=ON \
        -DAMReX_PRECISION="$precision" \
        -DBUILD_SHARED_LIBS=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_C_COMPILER=/opt/rocm/llvm/bin/clang \
        -DCMAKE_CXX_COMPILER=/opt/rocm/llvm/bin/clang++ \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_PREFIX_PATH='/opt/rocm/lib/cmake/AMDDeviceLibs;/opt/rocm/lib/cmake/amd_comgr;/opt/rocm/lib/cmake/hip;/opt/rocm/lib/cmake/hiprand;/opt/rocm/lib/cmake/hsa-runtime64;/opt/rocm/lib/cmake/rocprim;/opt/rocm/lib/cmake/rocrand' \
        -DMPI_C_ADDITIONAL_INCLUDE_DIRS='/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi;/usr/lib/x86_64-linux-gnu/openmpi/include' \
        -DMPI_C_LIB_NAMES='mpi;open-rte;open-pal;hwloc' \
        -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS='/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi;/usr/lib/x86_64-linux-gnu/openmpi/include' \
        -DMPI_CXX_LIB_NAMES='mpi_cxx;mpi;open-rte;open-pal;hwloc' \
        -DMPI_hwloc_LIBRARY=/lib/x86_64-linux-gnu/libhwloc.so \
        -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so \
        -DMPI_mpi_cxx_LIBRARY=/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so \
        -DMPI_open-pal_LIBRARY=/lib/x86_64-linux-gnu/libopen-pal.so \
        -DMPI_open-rte_LIBRARY=/lib/x86_64-linux-gnu/libopen-rte.so \
        .. && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src
