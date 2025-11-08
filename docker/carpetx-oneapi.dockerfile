# How to build this Docker image:

#     docker build --file carpetx-oneapi.dockerfile --tag einsteintoolkit/carpetx:oneapi-real64 .
#     docker push einsteintoolkit/carpetx:oneapi-real64

#     docker build --build-arg real_precision=real32 --file carpetx-oneapi.dockerfile --tag einsteintoolkit/carpetx:oneapi-real32 .
#     docker push einsteintoolkit/carpetx:oneapi-real32

# FROM intel/oneapi-basekit:2025.2.0-0-devel-ubuntu24.04
FROM intel/oneapi-basekit:2025.2.2-0-devel-ubuntu24.04
# [AMReX linker step is broken]
# FROM intel/oneapi-basekit:2025.3.0-0-devel-ubuntu24.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANGUAGE=en_US.en \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 28DA432DAAC8BAEA && \
    apt-get update --allow-insecure-repositories && \
    apt-get --yes --no-install-recommends install \
        ca-certificates \
        clang-format \
        cmake \
        cvs \
        diffutils \
        g++ \
        gcc \
        gdb \
        gfortran \
        git \
        hdf5-filter-plugin \
        hdf5-filter-plugin-blosc-serial \
        hdf5-filter-plugin-zfp-serial \
        hdf5-plugin-lzf \
        hdf5-tools \
        hwloc-nox \
        language-pack-en \
        less \
        libblosc-dev \
        libblosc2-dev \
        libboost-all-dev \
        libbz2-dev \
        libfftw3-dev \
        libgit2-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        liblz4-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libprotobuf-dev \
        libreadline-dev \
        libtool \
        libudev-dev \
        libyaml-cpp-dev \
        libzfp-dev \
        libzstd-dev \
        locales \
        m4 \
        meson \
        ninja-build \
        nlohmann-json3-dev \
        numactl \
        perl \
        pkgconf \
        protobuf-compiler \
        python3 \
        python3-numpy \
        python3-pip \
        python3-requests \
        rsync \
        subversion \
        vim \
        wget \
        zlib1g-dev \
        && \
    rm -rf /var/lib/apt/lists/*

# Remove troublesome libraries
RUN find /opt/intel -name 'impi.pc' -delete && \
    find /opt/intel -name 'libhwloc.a' -delete && \
    find /opt/intel -name 'libhwloc.so' -delete && \
    find /opt/intel -name 'libmpi*a' -delete && \
    find /opt/intel -name 'libmpi*so' -delete && \
    find /opt/intel -name 'mpi.h' -delete && \
    find /opt/intel -name 'mpicc' -delete && \
    find /opt/intel -name 'mpicxx' -delete && \
    find /opt/intel -name 'mpiexec' -delete && \
    find /opt/intel -name 'mpirun' -delete

# # Install blosc2
# # blosc2 is a compression library, comparable to zlib
# RUN mkdir src && \
#     (cd src && \
#     wget https://github.com/Blosc/c-blosc2/archive/refs/tags/v2.18.0.tar.gz && \
#     tar xzf v2.18.0.tar.gz && \
#     cd c-blosc2-2.18.0 && \
#     cmake -B build -G Ninja \
#         -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#         -DCMAKE_INSTALL_PREFIX=/usr/local \
#         -DBUILD_BENCHMARKS=OFF \
#         -DBUILD_EXAMPLES=OFF \
#         -DBUILD_FUZZERS=OFF \
#         -DBUILD_STATIC=OFF \
#         -DBUILD_TESTS=OFF \
#         -DPREFER_EXTERNAL_LZ4=ON  \
#         -DPREFER_EXTERNAL_ZLIB=ON \
#         -DPREFER_EXTERNAL_ZSTD=ON \
#         && \
#     cmake --build build && \
#     cmake --install build && \
#     true) && \
#     rm -rf src

# Install MGARD
# MGARD is a lossy compression library
# Note: -DMGARD_ENABLE_SYCL=ON does not work
RUN mkdir src && \
    (cd src && \
    wget https://github.com/CODARcode/MGARD/archive/refs/tags/1.5.2.tar.gz && \
    tar xzf 1.5.2.tar.gz && \
    cd MGARD-1.5.2 && \
    cmake -B build -G Ninja \
        -DBUILD_TESTING=OFF \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_PREFIX_PATH=/usr/local \
        -DMGARD_ENABLE_OPENMP=ON \
        -DMGARD_ENABLE_SERIAL=ON \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
# - depends on blosc2
# - depends on MGARD
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.10.2.tar.gz && \
    tar xzf v2.10.2.tar.gz && \
    cd ADIOS2-2.10.2 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_TESTING=OFF \
        -DADIOS2_BUILD_EXAMPLES=OFF \
        -DADIOS2_Blosc2_PREFER_SHARED=ON \
        -DADIOS2_USE_BZip2=ON \
        -DADIOS2_USE_Blosc2=ON \
        -DADIOS2_USE_Fortran=OFF \
        -DADIOS2_USE_HDF5=ON \
        -DADIOS2_USE_MGARD=ON \
        -DADIOS2_USE_ZFP=ON \
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
    wget https://github.com/eschnett/asdf-cxx/archive/refs/tags/version/7.3.2.tar.gz && \
    tar xzf 7.3.2.tar.gz && \
    cd asdf-cxx-version-7.3.2 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
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
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -Dsimd=AVX2 \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# COPY patches/openPMD-api.patch /cactus/patches/

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - depends on ADIOS2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.16.1.tar.gz && \
    tar xzf 0.16.1.tar.gz && \
    cd openPMD-api-0.16.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_EXAMPLES=OFF \
        -DBUILD_TESTING=OFF \
        -DopenPMD_BUILD_SHARED_LIBS=ON \
        -DopenPMD_USE_MPI=ON \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install RePrimAnd
# RePrimAnd is a physics package for nuclear equations of state
RUN mkdir src && \
    (cd src && \
    wget https://github.com/wokast/RePrimAnd/archive/refs/tags/v1.7.tar.gz && \
    tar xzf v1.7.tar.gz && \
    cd RePrimAnd-1.7 && \
    meson build --buildtype=release --prefix=/usr/local && \
    ninja -C build && \
    ninja -C build install && \
    true) && \
    rm -rf src

# Install Silo
# Silo defines a standard for laying out AMR data in a file
RUN mkdir src && \
    (cd src && \
    wget https://github.com/LLNL/Silo/releases/download/4.11.1/silo-4.11.1.tar.xz && \
    tar xJf silo-4.11.1.tar.xz && \
    cd silo-4.11.1 && \
    mkdir build && \
    cd build && \
    ../configure \
        --disable-fortran \
        --disable-static \
        --enable-optimization \
        --enable-shared \
        --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/include,/usr/lib/x86_64-linux-gnu/hdf5/serial/lib \
        --prefix=/usr/local \
        && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    true) && \
    rm -rf src

# Install Conduit
# Conduit defines a standard for laying out AMR data in a file.
# - depends on Silo
# TODO:
# - enable CUDA? HIP?
# -DADIOS_DIR=/usr/local   # conduit doesn't find ADIOS because there is no FindADIOS.cmake
# -DZFP_DIR=/usr           # conduit doesn't find zfp because the library directory is wrong
RUN mkdir src && \
    (cd src && \
    wget https://github.com/LLNL/conduit/releases/download/v0.9.5/conduit-v0.9.5-src-with-blt.tar.gz && \
    tar xzf conduit-v0.9.5-src-with-blt.tar.gz && \
    cd conduit-v0.9.5 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DCONDUIT_ENABLE_TESTS=OFF \
        -DENABLE_COVERAGE=OFF \
        -DENABLE_DOCS=OFF \
        -DENABLE_EXAMPLES=OFF \
        -DENABLE_FORTRAN=OFF \
        -DENABLE_MPI=ON \
        -DENABLE_OPENMP=ON \
        -DENABLE_PYTHON=ON \
        -DENABLE_RELAY_WEBSERVER=OFF \
        -DENABLE_TESTS=OFF \
        -DENABLE_UTILS=ON \
        -DHDF5_DIR=/usr \
        -DSILO_DIR=/usr/local \
        -DZLIB_DIR=/usr \
        src \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install SimulationIO
# SimulationIO is an I/O library like HDF5
# - depends on asdf-cxx
# - depends on yaml-cpp
# Currently disabling ASDF because there is a confusion with C++11/17 standards
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/SimulationIO/archive/refs/tags/version/9.0.3.tar.gz && \
    tar xzf 9.0.3.tar.gz && \
    cd SimulationIO-version-9.0.3 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DENABLE_ASDF_CXX=OFF \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ssht
# ssht provides spin-weighted spherical harmonics
RUN mkdir src && \
    (cd src && \
    wget https://github.com/astro-informatics/ssht/archive/v1.5.2.tar.gz && \
    tar xzf v1.5.2.tar.gz && \
    cd ssht-1.5.2 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_TESTING=OFF \
        && \
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
    wget https://github.com/AMReX-Codes/amrex/archive/25.11.tar.gz && \
    tar xzf 25.11.tar.gz && \
    rm -rf /opt/intel/oneapi/mpi && \
    cd amrex-25.11 && \
    case $real_precision in \
        real32) precision=SINGLE;; \
        real64) precision=DOUBLE;; \
        *) exit 1;; \
    esac && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_C_COMPILER=icx \
        -DCMAKE_CXX_COMPILER=icpx \
        -DCMAKE_PREFIX_PATH='/opt/intel/oneapi' \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DAMReX_FORTRAN=OFF \
        -DAMReX_FORTRAN_INTERFACES=OFF \
        -DAMReX_GPU_BACKEND=SYCL \
        -DAMReX_INTEL_ARCH=intel_gpu_pvc \
        -DAMReX_OMP=OFF \
        -DAMReX_PARTICLES=ON \
        -DAMReX_PRECISION="$precision" \
        -DMPI_CXX_ADDITIONAL_INCLUDE_DIRS='/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi;/usr/lib/x86_64-linux-gnu/openmpi/include' \
        -DMPI_CXX_LIB_NAMES='mpi_cxx;mpi;open-rte;open-pal;hwloc' \
        -DMPI_C_ADDITIONAL_INCLUDE_DIRS='/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi;/usr/lib/x86_64-linux-gnu/openmpi/include' \
        -DMPI_C_LIB_NAMES='mpi;open-rte;open-pal;hwloc' \
        -DMPI_hwloc_LIBRARY=/lib/x86_64-linux-gnu/libhwloc.so \
        -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so \
        -DMPI_mpi_cxx_LIBRARY=/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so \
        -DMPI_open-pal_LIBRARY=/lib/x86_64-linux-gnu/libopen-pal.so \
        -DMPI_open-rte_LIBRARY=/lib/x86_64-linux-gnu/libopen-rte.so \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Find libraries in /usr/local/lib64
RUN echo /usr/local/lib64 >/etc/ld.so.conf.d/usr-local-lib64.conf && \
    ldconfig
