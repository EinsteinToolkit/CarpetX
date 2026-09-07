# How to build this Docker image:

#     docker build --file carpetx-native-arm64v8-cpu.dockerfile --tag einsteintoolkit/carpetx:native-arm64v8-cpu-real64 .
#     docker push einsteintoolkit/carpetx:native-arm64v8-cpu-real64

#     docker build --build-arg real_precision=real32 --file carpetx-native-arm64v8-cpu.dockerfile --tag einsteintoolkit/carpetx:native-arm64v8-cpu-real32 .
#     docker push einsteintoolkit/carpetx:native-arm64v8-cpu-real32

# noble is ubuntu:24.04
# FROM arm64v8/ubuntu:noble-20250714
# FROM arm64v8/ubuntu:noble-20250805
# FROM arm64v8/ubuntu:noble-20251001
# FROM arm64v8/ubuntu:noble-20251013
# FROM arm64v8/ubuntu:noble-20260113
# FROM arm64v8/ubuntu:noble-20260210.1
# FROM arm64v8/ubuntu:noble-20260324
# FROM arm64v8/ubuntu:noble-20260410
# FROM arm64v8/ubuntu:noble-20260509.1
# FROM arm64v8/ubuntu:resolute-20260610
FROM arm64v8/ubuntu:resolute-20260811.1

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
        curl \
        cvs \
        diffutils \
        elfutils \
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
        libiberty-dev \
        liblz4-dev \
        liblzma-dev \
        libopenblas-dev \
        libopenmpi-dev \
        libpapi-dev \
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
        make \
        meson \
        ninja-build \
        nlohmann-json3-dev \
        numactl \
        papi-tools \
        patch \
        perl \
        pkgconf \
        protobuf-compiler \
        python3 \
        python3-numpy \
        python3-pip \
        python3-requests \
        python3-venv \
        rsync \
        subversion \
        tar \
        unzip \
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
    wget https://github.com/spack/spack/archive/refs/tags/v1.2.2.tar.gz && \
    tar xzf v1.2.2.tar.gz && \
    export SPACK_ROOT="$(pwd)/spack-1.2.2" && \
    mkdir -p "${HOME}/.spack" && \
    echo 'config: {install_tree: {root: /spack}}' >"${HOME}/.spack/config.yaml" && \
    . ${SPACK_ROOT}/share/spack/setup-env.sh && \
    wget https://github.com/spack/spack-packages/archive/refs/tags/v2026.06.0.tar.gz && \
    tar xzf v2026.06.0.tar.gz && \
    spack repo add --scope site spack-packages-2026.06.0/repos/spack_repo/builtin && \
    spack external find \
        autoconf \
        automake \
        cmake \
        curl \
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

# Install MGARD
# MGARD is a lossy compression library
RUN mkdir src && \
    (cd src && \
    wget https://github.com/CODARcode/MGARD/archive/refs/tags/1.6.0.tar.gz && \
    tar xzf 1.6.0.tar.gz && \
    cd MGARD-1.6.0 && \
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

# Install SZ3
# SZ3 is a lossy compression library
RUN mkdir src && \
    (cd src && \
    wget https://github.com/szcompressor/SZ3/releases/download/v3.3.2/SZ3-v3.3.2.zip && \
    unzip SZ3-v3.3.2.zip && \
    cd SZ3-master && \
    cmake -B build -G Ninja \
        -DBUILD_H5Z_FILTER=ON \
        -DBUILD_MDZ=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_PREFIX_PATH=/usr/local \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
# - depends on blosc2
# - depends on MGARD
#
#     wget https://github.com/GTkorvo/dill/commit/94b4c437182318a0d446fb6f82511fac0e4f2516.patch && \
#     (cd thirdparty/dill/dill && patch -p1 <../94b4c437182318a0d446fb6f82511fac0e4f2516.patch) && \
#
#     wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.12.0.tar.gz && \
#     tar xzf v2.12.0.tar.gz && \
#     cd ADIOS2-2.12.0 && \
#
#    wget https://github.com/ornladios/ADIOS2/archive/45035d24b4505b241037e2a1b35a4bbf1af782a8.tar.gz && \
#    tar xzf 45035d24b4505b241037e2a1b35a4bbf1af782a8.tar.gz && \
#    cd ADIOS2-45035d24b4505b241037e2a1b35a4bbf1af782a8 && \
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.12.1.tar.gz && \
    tar xzf v2.12.1.tar.gz && \
    cd ADIOS2-2.12.1 && \
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
        -DADIOS2_USE_SZ3=ON \
        -DADIOS2_USE_ZFP=ON && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ASDF
# ASDF is an I/O library like HDF5
# - depends on yaml-cpp
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/asdf-cxx/archive/refs/tags/version/8.0.0.tar.gz && \
    tar xzf 8.0.0.tar.gz && \
    cd asdf-cxx-version-8.0.0 && \
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
# Note: This assumes that the system has aarch64 (ARM) CPUs
RUN mkdir src && \
    (cd src && \
    wget https://github.com/agenium-scale/nsimd/archive/refs/tags/v3.0.1.tar.gz && \
    tar xzf v3.0.1.tar.gz && \
    cd nsimd-3.0.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -Dsimd=aarch64 \
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
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.17.1.tar.gz && \
    tar xzf 0.17.1.tar.gz && \
    cd openPMD-api-0.17.1 && \
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
    meson setup build --buildtype=release --prefix=/usr/local -Dcpp_std=c++14 && \
    ninja -C build && \
    ninja -C build install && \
    true) && \
    rm -rf src

# Install Silo
# Silo defines a standard for laying out AMR data in a file
RUN mkdir src && \
    (cd src && \
    wget https://github.com/LLNL/Silo/releases/download/4.12.1/Silo-4.12.1.tar.xz && \
    tar xJf Silo-4.12.1.tar.xz && \
    cd Silo-4.12.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DSILO_BUILD_FOR_BSD_LICENSE=ON \
        -DSILO_ENABLE_BROWSER=OFF \
        -DSILO_ENABLE_FORTRAN=OFF \
        -DSILO_ENABLE_HDF5=ON \
        -DSILO_ENABLE_JSON=OFF \
        -DSILO_ENABLE_PYTHON_MODULE=OFF \
        -DSILO_ENABLE_SHARED=ON \
        -DSILO_ENABLE_SILEX=OFF \
        -DSILO_ENABLE_SILOCK=ON \
        -DSILO_ENABLE_TESTS=OFF \
        -DSILO_HDF5_SZIP_DIR=/usr/local \
        && \
    cmake --build build && \
    cmake --install build && \
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
    wget https://github.com/LLNL/conduit/releases/download/v0.9.8/conduit-v0.9.8-src-with-blt.tar.gz && \
    tar xzf conduit-v0.9.8-src-with-blt.tar.gz && \
    cd conduit-v0.9.8 && \
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
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
        -DENABLE_ASDF_CXX=ON \
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

# # Install Fuka
# RUN mkdir src
# WORKDIR src
# RUN git clone --branch ET_2025_05 https://bitbucket.org/fukaws/fuka
# WORKDIR fuka
# # Ensure this macro is defined at build time to enable multi-threaded importers
# RUN sed -i -E 's%// #define DEFAULT_KAD_MEM%#define DEFAULT_KAD_MEM%' include/memory.hpp
# WORKDIR build_release
# RUN env HOME_KADATH=/cactus/src/fuka \
#     cmake -B build -G Ninja \
#         -DCMAKE_BUILD_TYPE=RelWithDebInfo \
#         -DCMAKE_INSTALL_PREFIX=/usr/local \
#         -DGRHAYL_EOS=OFF \
#         -DPAR_VERSION=ON
# RUN cmake --build build
# WORKDIR ..
# RUN cp lib/libkadath.a /usr/local/lib
# RUN ls -l /cactus/src/fuka/include
# RUN false
# WORKDIR ../..
# RUN rm -rf src
# RUN echo
# RUN ls -l /usr/local/include
# RUN ls -l /usr/local/lib
# RUN false

ARG real_precision=real64

# Install AMReX
# AMReX provides adaptive mesh refinement
# - Enable Fortran for `docker/Dockerfile`
# - Install this last because it changes most often
# Should we keep the AMReX source tree around for debugging?
RUN mkdir src && \
    (cd src && \
    wget https://github.com/AMReX-Codes/amrex/archive/26.01.tar.gz && \
    tar xzf 26.01.tar.gz && \
    cd amrex-26.01 && \
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
