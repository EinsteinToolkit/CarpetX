# How to build this Docker image:

#     docker build --file carpetx-cpu.dockerfile --tag einsteintoolkit/carpetx:cpu-real64 .
#     docker push einsteintoolkit/carpetx:cpu-real64

#     docker build --build-arg real_precision=real32 --file carpetx-cpu.dockerfile --tag einsteintoolkit/carpetx:cpu-real32 .
#     docker push einsteintoolkit/carpetx:cpu-real32

# noble is ubuntu:24.04
# FROM amd64/ubuntu:noble-20240801
# FROM amd64/ubuntu:noble-20240904.1
FROM amd64/ubuntu:noble-20241011

ENV DEBIAN_FRONTEND=noninteractive \
    LANGUAGE=en_US.en \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

RUN mkdir /cactus
WORKDIR /cactus

# Install system packages
# - Boost on Ubuntu requires OpenMPI
#        elfutils
#        python2
RUN apt-get update && \
    apt-get --yes --no-install-recommends install \
        bzip2 \
        ca-certificates \
        clang-format \
        cmake \
        cvs \
        diffutils \
        elfutils \
        g++ \
        gcc \
        gdb \
        gfortran \
        git \
        hdf5-tools \
        hwloc-nox \
        language-pack-en \
        less \
        libblosc-dev \
        libboost-all-dev \
        libbz2-dev \
        libfftw3-dev \
        libgit2-dev \
        libgsl-dev \
        libhdf5-dev \
        libhwloc-dev \
        libiberty-dev \
        liblzma-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpapi-dev \
        libpetsc-real-dev \
        libprotobuf-dev \
        libtool \
        libudev-dev \
        libyaml-cpp-dev \
        libzstd-dev \
        locales \
        m4 \
        make \
        meson \
        ninja-build \
        numactl \
        papi-tools \
        patch \
        perl \
        pkgconf \
        protobuf-compiler \
        python3 \
        python3-pip \
        python3-requests \
        rsync \
        subversion \
        tar \
        vim \
        wget \
        xz-utils \
        zlib1g-dev \
        zstd \
        && \
    rm -rf /var/lib/apt/lists/*

# Install HPCToolkit
# Install this first because it is expensive to build
# Try to reuse build tools from Ubuntu, but do not use any libraries because HPC Toolkit is a bit iffy to install.
RUN mkdir src && \
    (cd src && \
    wget https://github.com/spack/spack/archive/refs/tags/v0.22.2.tar.gz && \
    tar xzf v0.22.2.tar.gz && \
    export SPACK_ROOT="$(pwd)/spack-0.22.2" && \
    mkdir -p "${HOME}/.spack" && \
    echo 'config: {install_tree: {root: /spack}}' >"${HOME}/.spack/config.yaml" && \
    . ${SPACK_ROOT}/share/spack/setup-env.sh && \
    spack external find \
        autoconf \
        automake \
        cmake \
        diffutils \
        elfutils \
        gmake \
        libtool \
        m4 \
        meson \
        ninja \
        numactl \
        perl \
        pkgconf \
        python \
    && \
    spack install --fail-fast hpctoolkit ~viewer && \
    spack view --dependencies no hardlink /hpctoolkit hpctoolkit && \
    true) && \
    rm -rf src "${HOME}/.spack"

# Install blosc2
# blosc2 is a compression library, comparable to zlib
RUN mkdir src && \
    (cd src && \
    wget https://github.com/Blosc/c-blosc2/archive/refs/tags/v2.15.1.tar.gz && \
    tar xzf v2.15.1.tar.gz && \
    cd c-blosc2-2.15.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_BENCHMARKS=OFF \
        -DBUILD_EXAMPLES=OFF \
        -DBUILD_FUZZERS=OFF \
        -DBUILD_STATIC=OFF \
        -DBUILD_TESTS=OFF \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install MGARD
# MGARD is a lossy compression library
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
        -DADIOS2_USE_Blosc2=ON \
        -DADIOS2_USE_Fortran=OFF \
        -DADIOS2_USE_MGARD=ON \
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

# Install openPMD-api
# openPMD-api defines a standard for laying out AMR data in a file
# - depends on ADIOS2
RUN mkdir src && \
    (cd src && \
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.16.0.tar.gz && \
    tar xzf 0.16.0.tar.gz && \
    cd openPMD-api-0.16.0 && \
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
    wget https://github.com/AMReX-Codes/amrex/archive/24.11.tar.gz && \
    tar xzf 24.11.tar.gz && \
    cd amrex-24.11 && \
    case $real_precision in \
        real32) precision=SINGLE;; \
        real64) precision=DOUBLE;; \
        *) exit 1;; \
    esac && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DAMReX_FORTRAN=OFF \
        -DAMReX_FORTRAN_INTERFACES=OFF \
        -DAMReX_OMP=ON \
        -DAMReX_PARTICLES=ON \
        -DAMReX_PRECISION="$precision" \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Find libraries in /usr/local/lib64
RUN echo /usr/local/lib64 >/etc/ld.so.conf.d/usr-local-lib64.conf && \
    ldconfig
