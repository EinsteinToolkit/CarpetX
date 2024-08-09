# How to build this Docker image:

#     docker build --file carpetx-rocm.dockerfile --tag einsteintoolkit/carpetx:rocm-real64 .
#     docker push einsteintoolkit/carpetx:rocm-real64

#     docker build --build-arg real_precision=real32 --file carpetx-rocm.dockerfile --tag einsteintoolkit/carpetx:rocm-real32 .
#     docker push einsteintoolkit/carpetx:rocm-real32

# FROM rocm/dev-ubuntu-22.04:6.1.2
FROM rocm/dev-ubuntu-24.04:6.2

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
        ca-certificates \
        cmake \
        cvs \
        diffutils \
        g++ \
        gcc \
        gdb \
        gfortran \
        git \
        hdf5-tools \
        hiprand-dev \
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
        libopenblas-dev \
        libopenmpi-dev \
        libpetsc-real-dev \
        libtool \
        libudev-dev \
        libyaml-cpp-dev \
        libzstd-dev \
        locales \
        m4 \
        meson \
        ninja-build \
        numactl \
        perl \
        pkgconf \
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

# # Install HPCToolkit
# # Install this first because it is expensive to build
# # Installing HPCToolkit without rocm since the install fails when this option is enabled.
# RUN mkdir src && \
#     (cd src && \
#     wget https://github.com/spack/spack/archive/refs/tags/v0.21.0.tar.gz && \
#     tar xzf v0.21.0.tar.gz && \
#     export SPACK_ROOT="$(pwd)/spack-0.21.0" && \
#     mkdir -p "${HOME}/.spack" && \
#     echo 'config: {install_tree: {root: /spack}}' >"${HOME}/.spack/config.yaml" && \
#     . ${SPACK_ROOT}/share/spack/setup-env.sh && \
#     spack external find \
#         autoconf \
#         automake \
#         cmake \
#         diffutils \
#         elfutils \
#         gmake \
#         libtool \
#         m4 \
#         meson \
#         ninja \
#         numactl \
#         perl \
#         pkgconf \
#         python \
#     && \
#     spack install --fail-fast hpctoolkit ~viewer && \
#     spack view --dependencies no hardlink /hpctoolkit hpctoolkit && \
#     true) && \
#     rm -rf src "${HOME}/.spack"

# Install blosc2
# blosc2 is a compression library, comparable to zlib
RUN mkdir src && \
    (cd src && \
    wget https://github.com/Blosc/c-blosc2/archive/refs/tags/v2.15.0.tar.gz && \
    tar xzf v2.15.0.tar.gz && \
    cd c-blosc2-2.15.0 && \
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

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
RUN mkdir src && \
    (cd src && \
    wget https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.10.1.tar.gz && \
    tar xzf v2.10.1.tar.gz && \
    cd ADIOS2-2.10.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_TESTING=OFF \
        -DADIOS2_BUILD_EXAMPLES=OFF \
        -DADIOS2_Blosc2_PREFER_SHARED=ON \
        -DADIOS2_USE_Blosc2=ON \
        -DADIOS2_USE_Fortran=OFF \
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
    wget https://github.com/openPMD/openPMD-api/archive/refs/tags/0.15.2.tar.gz && \
    tar xzf 0.15.2.tar.gz && \
    cd openPMD-api-0.15.2 && \
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
    wget https://github.com/AMReX-Codes/amrex/archive/24.08.tar.gz && \
    tar xzf 24.08.tar.gz && \
    cd amrex-24.08 && \
    case $real_precision in \
        real32) precision=SINGLE;; \
        real64) precision=DOUBLE;; \
        *) exit 1;; \
    esac && \
    env LDFLAGS=-Wl,-rpath,/opt/rocm/lib \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_CXX_COMPILER=/opt/rocm/llvm/bin/clang++ \
        -DCMAKE_C_COMPILER=/opt/rocm/llvm/bin/clang \
        -DCMAKE_PREFIX_PATH='/opt/rocm/lib/cmake/AMDDeviceLibs;/opt/rocm/lib/cmake/amd_comgr;/opt/rocm/lib/cmake/hip;/opt/rocm/lib/cmake/hiprand;/opt/rocm/lib/cmake/hsa-runtime64;/opt/rocm/lib/cmake/rocprim;/opt/rocm/lib/cmake/rocrand' \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DAMReX_AMD_ARCH=gfx90a \
        -DAMReX_FORTRAN=OFF \
        -DAMReX_FORTRAN_INTERFACES=OFF \
        -DAMReX_GPU_BACKEND=HIP \
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
