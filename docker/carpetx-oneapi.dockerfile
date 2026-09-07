# How to build this Docker image:

#     docker build --file carpetx-oneapi.dockerfile --tag einsteintoolkit/carpetx:oneapi-real64 .
#     docker push einsteintoolkit/carpetx:oneapi-real64

#     docker build --build-arg real_precision=real32 --file carpetx-oneapi.dockerfile --tag einsteintoolkit/carpetx:oneapi-real32 .
#     docker push einsteintoolkit/carpetx:oneapi-real32

# This is `carpetx-native-oneapi.dockerfile` with every MPI-using library moved
# onto the **MPI ABI** (MPI-5.0 §20.2) instead of a directly-linked MPI. Two
# packages provide that ABI over the distribution's MPI:
#
#   mpi_abi_wrapper  the C bindings: `mpi.h`, `libmpi_abi`, `mpicc`, `mpicxx`
#   mpif             the Fortran bindings: `mpif.h`, `mpi.mod`, `mpi_f08.mod`,
#                    `libmpif`, `mpifort` -- needed because §20.4 leaves
#                    `MPI_Fint` and the Fortran modules out of the ABI entirely
#
# Both install into `/usr/local`, so `/usr/local/bin/{mpicc,mpicxx,mpifort}` is
# the toolchain every stanza below is pointed at, and neither `libmpi.so.40`
# (Open MPI) nor `libmpi.so.12` (Intel MPI) may appear in the `ldd` output of
# anything this file installs. Which MPI actually
# runs a binary is then a *run-time* choice: `libmpi_abi` dlopens the wrapper,
# and another ABI-providing implementation can be put in front of it without
# recompiling. `/usr/local/bin/mpiexec` is a **forwarder**, not a launcher:
# starting a job stays the wrapped MPI's business, so it resolves that MPI's own
# `mpiexec` -- `$MPI_ABI_WRAPPER_MPIEXEC`, else the path found beside
# `MPI_C_COMPILER` at configure time -- and passes every argument through. That
# is the same run-time indirection `libmpi_abi` uses to find `libmpiwrapper`, so
# the launcher follows the library when a binary is re-pointed.
#
# Consequences for the apt list below:
# - `libhdf5-dev` and `libpetsc-real-dev` are **gone**. Both are built here
#   instead, over the ABI. PETSc drags in a whole parallel stack when installed
#   from apt -- parallel HDF5, hypre, MUMPS, ScaLAPACK, PtScotch, SuperLU-DIST,
#   CombBLAS, FFTW3-MPI -- every one of which links `libmpi.so.40`; dropping it
#   removes all of them at once, and the PETSc built here has none of them.
# - `libaec-dev` is **added**: it provides the szip encoder that Debian's HDF5
#   was built with, which the HDF5 built here would otherwise lose.
# - `libcurl4-openssl-dev` and `libjpeg-dev` are **added** because they were
#   only ever present as transitive dependencies of `libpetsc-real-dev`, and
#   dropping it took them with it. The Cactus option list points LIBCURL_DIR
#   and LIBJPEG_DIR at /usr, so their absence stops the CST outright:
#   "LIBJPEG not found at /usr". Naming them explicitly is right regardless --
#   depending on PETSc to supply an image's libjpeg was always an accident.
# - `python3-dev` is **added** for the same kind of reason. Conduit is built
#   with `-DENABLE_PYTHON=ON`, so its `find_package(Python3 ... Development
#   NumPy)` needs `Python.h` and the numpy headers; those arrived only as a
#   transitive dependency of `libboost-all-dev`, by way of
#   `libboost-python1.83-dev`. Naming it explicitly is right regardless, and is
#   required outright in the images below that narrow Boost.
# - `libopenmpi-dev` is **gone**. The MPI being wrapped here is Intel MPI, out
#   of the base image, and there is to be exactly one MPI in the image.
# - `libboost-all-dev` becomes **`libboost-filesystem-dev`**, and this is not
#   cosmetic: `libboost-all-dev` depends on `libboost-mpi-dev`, which depends on
#   `libopenmpi-dev`, so leaving it in would reinstall Ubuntu's Open MPI through
#   the back door and defeat the line above. `libboost-filesystem-dev` reaches
#   no MPI package at all (checked with `apt-cache depends --recurse`).
#   **Not plain `libboost-dev`**, which is the obvious narrowing and does not
#   work: RePrimAnd's `meson.build` asks for `dependency('boost')`, and Meson's
#   Boost detection wants to see boost *libraries* even when the dependency is
#   header-only -- with headers alone it reports "Run-time dependency Boost
#   found: NO (tried system)" and RePrimAnd stops. Measured on Ubuntu 26.04:
#   `libboost-dev` installs 0 boost libraries and fails; `libboost-filesystem-dev`
#   installs 6 and gives "Boost found: YES 1.90.0 (/usr)".
#   `filesystem` is the right component to keep rather than an arbitrary one:
#   `desert-arm64v8.cfg` sets `BOOST_LIBS="boost_filesystem"`, so Cactus does
#   link it. Nothing else does -- CarpetX's own Boost use is header-only
#   (`boost::math` in `repos/CarpetX/Algo/src/roots.hxx`) and RePrimAnd's
#   `dependency('boost')` names no modules -- so for those it only has to be
#   present to be found.
#
# The matching Cactus option list needs `MPI_DIR = /usr/local`,
# `HDF5_DIR = /usr/local` (with `hdf5_hl_fortran`/`hdf5_fortran`, the CMake
# spelling, rather than Debian's `hdf5hl_fortran`), and `PETSC_DIR = /usr/local`.
#
# oneAPI-specific notes:
# - **the wrapped MPI is Intel MPI**, at `$I_MPI_ROOT` =
#   `/opt/intel/oneapi/mpi/2021.18`, which the base image ships. This reverses
#   `carpetx-native-oneapi.dockerfile`'s policy: that file deletes Intel MPI outright
#   (a `find /opt/intel -name 'libmpi*so' -delete` stanza, and a second
#   `rm -rf /opt/intel/oneapi/mpi` inside the AMReX stanza) in order to use
#   Ubuntu's Open MPI instead. Under the ABI there is no reason to throw the
#   vendor MPI away: it is wrapped like any other, and a run-time swap can still
#   put a different implementation in front of it. Both deletion stanzas are
#   therefore gone -- except for the `libhwloc.a`/`libhwloc.so` half of the
#   first, which is about Intel's hwloc colliding with the system one and has
#   nothing to do with MPI.
# - mpi_abi_wrapper is configured with `mpiicx`, Intel MPI's *Intel-compiler*
#   wrapper, rather than its GNU-compiler `mpicc`, so that the wrapper is built
#   by the same `icx` that compiles everything else in this image.
# - **the base image advertises Intel MPI's headers to every compiler
#   invocation**, through `C_INCLUDE_PATH` and `CPLUS_INCLUDE_PATH`. That
#   shadows the ABI's own `mpi.h` everywhere, and no compiler flag can undo it,
#   because /usr/local/include is a built-in system directory and `-I` on such a
#   directory is ignored. The wrapped MPI's include directory is therefore moved
#   aside once mpi_abi_wrapper has been built; see the stanza that does it.
# - the base image moves from `intel/oneapi-basekit` to
#   `intel/oneapi-toolkit`. basekit is dead-ended at 2025.3.2 with no Ubuntu 26
#   tag; oneapi-toolkit is its successor and bundles Intel MPI, MKL and the
#   compilers together. **Watch the AMReX link step**: `carpetx-native-oneapi.dockerfile`
#   is deliberately pinned back to basekit 2025.2.2 because every later basekit
#   broke it. Whether that is fixed in the 2026.1.0 toolkit is unknown; if it
#   recurs, fall back to basekit 2025.2.2 and keep everything else here.
# - the apt list is `carpetx-cpu.dockerfile`'s (minus `libopenmpi-dev`,
#   with Boost narrowed), not `carpetx-native-oneapi.dockerfile`'s. That file drops a
#   dozen packages -- `make`, `patch`, `tar`, `xz-utils`, `bzip2`, ... -- because
#   it builds neither HPCToolkit nor PETSc from source. The PETSc stanza below
#   does, and `--download-*` needs those tools.
# - AMReX's nine hard-wired `-DMPI_*` variables are **gone**. They named Open
#   MPI's `libmpi`, `libmpi_cxx`, `libopen-rte` and `libopen-pal` by absolute
#   path, which is both the wrong MPI now and exactly the kind of direct link
#   this file forbids. If `icpx` turns out not to accept what
#   `mpi_abi_wrapper`'s `mpicc -showme:link` reports, the ABI spelling is far
#   simpler, the ABI having no `open-rte`/`open-pal`/`mpi_cxx` split:
#   `-DMPI_C_LIB_NAMES=mpi_abi`,
#   `-DMPI_mpi_abi_LIBRARY=/usr/local/lib/libmpi_abi.so` and
#   `-DMPI_C_ADDITIONAL_INCLUDE_DIRS=/usr/local/include`.

# The basekit images are `carpetx-native-oneapi.dockerfile`'s history; that file is
# pinned to 2025.2.2 because the AMReX linker step is broken on every later
# basekit. basekit itself is now dead-ended -- 2025.3.2 is its last release and
# it has no Ubuntu 26 tag -- so oneapi-toolkit is the only way forward.
# FROM intel/oneapi-basekit:2025.2.0-0-devel-ubuntu24.04
# FROM intel/oneapi-basekit:2025.2.2-0-devel-ubuntu24.04
# FROM intel/oneapi-basekit:2025.3.0-0-devel-ubuntu24.04
# FROM intel/oneapi-basekit:2025.3.1-0-devel-ubuntu24.04
# FROM intel/oneapi-basekit:2025.3.2-0-devel-ubuntu24.04
# FROM intel/oneapi-toolkit:2026.0.0-devel-ubuntu24.04
# FROM intel/oneapi-toolkit:2026.0.1-devel-ubuntu24.04
FROM intel/oneapi-toolkit:2026.1.0-devel-ubuntu26.04

# The base image puts Intel MPI's `bin` second on PATH, six entries ahead of
# /usr/local/bin, so a bare `mpiexec`, `mpirun`, `mpicc` or `mpicxx` means the
# *wrapped* MPI's rather than the ABI prefix's. Today the two launchers are the
# same binary -- `/usr/local/bin/mpiexec` is a forwarder, and the launcher it
# forwards to is the one found beside `MPI_C_COMPILER` at configure time, i.e.
# Intel's -- so this is not a correctness bug now. It matters for two other
# reasons:
# - the run-time swap this whole file exists to enable. Point
#   `MPI_ABI_WRAPPER_LIB` at another implementation and the forwarder follows it
#   through `MPI_ABI_WRAPPER_MPIEXEC`, while a bare Intel `mpiexec` does not.
#   Launcher and library must match, and a mismatch does not fail loudly -- it
#   silently starts N one-rank jobs. Cactus's own runscript (`desert.run`) calls
#   `mpiexec` unqualified, so this is the ordering it will get.
# - `mpicc`/`mpicxx`. Intel MPI's wrappers stop working further down anyway,
#   once its `include` directory is moved aside, so a bare `mpicc` resolving to
#   them would fail confusingly rather than compile against the ABI.
# The CUDA variant needs the same line for the same reason; the ROCm and CPU
# variants do not, having /usr/local/bin ahead of their MPI already.
ENV PATH=/usr/local/bin:${PATH}

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
        bison \
        bzip2 \
        ca-certificates \
        clang-format \
        cmake \
        curl \
        cvs \
        diffutils \
        elfutils \
        flex \
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
        libaec-dev \
        libblosc-dev \
        libblosc2-dev \
        libboost-filesystem-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libfftw3-dev \
        libgit2-dev \
        libgsl-dev \
        libhwloc-dev \
        libiberty-dev \
        libjpeg-dev \
        liblz4-dev \
        liblzma-dev \
        libopenblas-dev \
        libpapi-dev \
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
        python3-dev \
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

# # Install HPCToolkit
# # Install this first because it is expensive to build
# # HPCToolkit is not built here, matching `carpetx-native-oneapi.dockerfile`.
# RUN mkdir src && \
#     (cd src && \
#     wget https://github.com/spack/spack/archive/refs/tags/v1.2.2.tar.gz && \
#     tar xzf v1.2.2.tar.gz && \
#     export SPACK_ROOT="$(pwd)/spack-1.2.2" && \
#     mkdir -p "${HOME}/.spack" && \
#     echo 'config: {install_tree: {root: /spack}}' >"${HOME}/.spack/config.yaml" && \
#     . ${SPACK_ROOT}/share/spack/setup-env.sh && \
#     wget https://github.com/spack/spack-packages/archive/refs/tags/v2026.06.0.tar.gz && \
#     tar xzf v2026.06.0.tar.gz && \
#     spack repo add --scope site spack-packages-2026.06.0/repos/spack_repo/builtin && \
#     spack external find \
#         autoconf \
#         automake \
#         cmake \
#         curl \
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

# Remove troublesome libraries
# Only the hwloc half of `carpetx-native-oneapi.dockerfile`'s stanza survives: Intel's
# static and shared hwloc collide with the system one. Every `libmpi*`, `mpi.h`,
# `mpicc`, `mpicxx`, `mpiexec` and `mpirun` deletion is gone -- that MPI is the
# one being wrapped below.
RUN find /opt/intel -name 'libhwloc.a' -delete && \
    find /opt/intel -name 'libhwloc.so' -delete

# Put the ABI's libraries ahead of Intel MPI's on the linker's search path
# `/usr/local/lib` is not in ld's default list on Debian, while the base image
# puts Intel MPI's `lib` on LIBRARY_PATH, so a bare `-lmpi` would find Intel
# MPI's `libmpi.so`. Nothing here asks for `-lmpi` -- the ABI is `-lmpi_abi` --
# but the ordering costs nothing and removes the possibility.
# The include side of the same problem is *not* solvable this way; see the
# stanza after mpi_abi_wrapper below.
ENV LIBRARY_PATH=/usr/local/lib:${LIBRARY_PATH}

# Install mpi_abi_wrapper
# mpi_abi_wrapper wraps an MPI library and provides the MPI ABI
# The MPI being wrapped is Intel MPI from the base image, named explicitly:
# left to itself find_package(MPI) would search PATH, and `mpiicx` is Intel
# MPI's Intel-compiler wrapper, which is what the rest of this image compiles
# with. MPI_HOME points at the same prefix so that any nested find_package(MPI)
# inside the wrapper's own build agrees.
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/mpi_abi_wrapper/archive/refs/tags/v1.2.0.tar.gz && \
    tar xzf v1.2.0.tar.gz && \
    cd mpi_abi_wrapper-1.2.0 && \
    cmake -B build -G Ninja \
        -DBUILD_TESTING=OFF \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_PREFIX_PATH=/usr/local \
        -DMPI_C_COMPILER="${I_MPI_ROOT}/bin/mpiicx" \
        -DMPI_HOME="${I_MPI_ROOT}" \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Hide the wrapped MPI's headers from everything below
# On `carpetx-cpu.dockerfile` nothing has to be done here: no other
# `mpi.h` on that image sits in a directory the compiler searches unasked. This
# image is different -- the base exports Intel MPI's include directory in
# `C_INCLUDE_PATH` and `CPLUS_INCLUDE_PATH` -- and that **cannot be out-ranked
# from the command line**. GCC searches `-I` and CPATH first, then `-isystem`
# and C_INCLUDE_PATH/CPLUS_INCLUDE_PATH, then its built-in directories; but
# /usr/local/include *is* a built-in directory, and "if a standard system
# include directory is also specified with -I, the -I option is ignored" -- it
# stays at the end of the chain. Measured, all three ways: with Intel MPI's
# directory on C_INCLUDE_PATH, plain compilation, `-I/usr/local/include` and
# `CPATH=/usr/local/include` all pick Intel MPI's header. Only `-isystem` or
# removing the directory changes the answer, and `-isystem` would have to be
# threaded through every stanza's flags.
#
# The symptom if this is skipped: mpif compiles against Intel MPI's `mpi.h` and
# fails on `MPI_ABI_VERSION undeclared` and on
# `sizeof(MPI_Status) == 8 * sizeof(MPI_Fint)` -- Intel MPI is MPICH-derived,
# with MPI_STATUS_SIZE 5, where the ABI mandates 8. HDF5, ADIOS2, openPMD, Silo
# and Conduit would all have hit it next.
#
# So the headers are moved aside, which makes the invariant this file rests on
# -- nothing below compiles against the MPI being wrapped -- a property of the
# filesystem rather than of a search path. It is version-independent, unlike
# editing `2021.18` out of the two variables. **Only the headers move**:
# `libmpi.so` stays where it is, because `libmpiwrapper` linked it just above
# and dlopens it at run time. `mpif.h` and the `mpi*.mod` files go along with
# them, which is intended -- mpif supplies its own, and a Fortran thorn finding
# Intel MPI's would be the same bug one language over.
#
# Placed after mpi_abi_wrapper, the one stanza that legitimately needs these
# headers, and before everything that must not see them.
RUN mv "${I_MPI_ROOT}/include" "${I_MPI_ROOT}/include.wrapped-by-mpi-abi"

# Install mpif
# mpif provides the Fortran bindings for the MPI ABI
# - depends on mpi_abi_wrapper
# The ABI is specified for C only: MPI-5.0 §20.4 leaves `MPI_Fint` and the
# Fortran modules out of it deliberately, so `mpi.h` and `libmpi_abi` give a
# Fortran caller nothing. mpif supplies `mpif.h`, `use mpi` and `use mpi_f08`
# on top of them, and is what makes HDF5's Fortran interface, PETSc's Fortran
# stubs and Cactus's own Fortran thorns possible over the ABI.
# MPI_C_COMPILER is named explicitly rather than left to find_package(MPI):
# CMake's FindMPI searches PATH, finds the wrapped MPI's own `mpicc` -- Intel
# MPI's, in this image -- and mpif then stops with "MPI_C_LIBRARIES ... names no
# libmpi_abi". MPI_HOME is what mpifort records as its default MPI prefix
# (`mpifort -showme:mpiprefix`).
# CMAKE_Fortran_COMPILER is pinned, which matters only in this image: the oneAPI
# base puts `ifx` on PATH alongside apt's gfortran, and CMake is otherwise free
# to pick either. `ifx` is named so that this image has one Fortran toolchain
# throughout, matching the `icx`/`icpx` the rest of it uses. HDF5 below carries
# the same pin and PETSc inherits it through mpifort, so `libmpif`, HDF5's
# Fortran interface and PETSc's Fortran stubs all agree.
# **The Cactus option list must follow**: Fortran `.mod` files are
# compiler-specific, so `desert-oneapi.cfg`'s `F90 = gfortran` has to become
# `F90 = ifx`, or `use mpi` in a thorn gets a module file it cannot read.
RUN mkdir src && \
    (cd src && \
    wget https://github.com/eschnett/mpif/archive/refs/tags/v1.0.1.tar.gz && \
    tar xzf v1.0.1.tar.gz && \
    cd mpif-1.0.1 && \
    cmake -B build -G Ninja \
        -DBUILD_SHARED_LIBS=ON \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_Fortran_COMPILER=ifx \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_HOME=/usr/local \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# From here on, every find_package(MPI) must resolve to the ABI rather than to
# the system Open MPI on PATH. The stanzas below say so explicitly with
# -DMPI_C_COMPILER, but that does not reach a find_package(MPI) issued from
# *inside* another package's config file -- and HDF5's `hdf5-config.cmake` does
# exactly that, to rebuild the MPI::MPI_C target its exported targets depend
# on. A consumer that merely links HDF5 therefore acquires an -lmpi against the
# system MPI, silently: measured on SZ3's H5Z filter and on Silo, both of which
# came out with a direct NEEDED on libmpi.so.40. MPI_HOME is a documented
# FindMPI hint (`MPI_HINT_DIRS`, FindMPI.cmake) read from the environment, so
# it reaches those nested calls where -D cannot.
ENV MPI_HOME=/usr/local

# Install HDF5
# HDF5 is the I/O library nearly everything below reads and writes through
# - depends on mpi_abi_wrapper (parallel I/O) and mpif (the Fortran interface)
# Built here rather than taken from apt because Debian's parallel flavour
# (`libhdf5-openmpi-dev`) links libmpi.so.40, and its serial flavour has no
# MPI-IO at all. Cactus wants one HDF5 with all three language interfaces:
# - HDF5_ALLOW_UNSUPPORTED is required because HDF5 refuses to combine the C++
#   interface with parallel I/O. Debian's `libhdf5-openmpi-cpp` is built the
#   same way; the C++ interface simply has no parallel entry points of its own.
# - HDF5 2.x renamed two of these options from their 1.14 spellings, and both
#   renames fail quietly rather than loudly: `ALLOW_UNSUPPORTED` became
#   `HDF5_ALLOW_UNSUPPORTED` (that one does at least stop the configure), and
#   `HDF5_ENABLE_Z_LIB_SUPPORT` became `HDF5_ENABLE_ZLIB_SUPPORT` -- the old
#   name is merely reported as an unused variable, and the build then finishes
#   with no zlib and hence no `deflate` filter at all.
# - the CMake build spells the Fortran libraries `hdf5_fortran` /
#   `hdf5_hl_fortran`, where Debian spells the latter `hdf5hl_fortran`. The
#   Cactus option list has to follow this spelling.
# - the `hdf5-filter-plugin*` and `hdf5-plugin-lzf` apt packages install their
#   dynamically-loaded filters under /usr/lib/*/hdf5/{serial,openmpi}/plugins,
#   built against Debian's `libhdf5_serial`. This build looks in
#   /usr/local/hdf5/lib/plugin instead, so it does not load them -- which is
#   the right outcome, since loading them would pull a second HDF5 with a
#   different soname into the process. Those packages remain useful to apt's
#   own /usr/bin/h5dump; a file needing one of those filters has to be read
#   with that, not with /usr/local/bin/h5dump.
RUN mkdir src && \
    (cd src && \
    wget https://github.com/HDFGroup/hdf5/releases/download/2.2.0/hdf5-2.2.0.tar.gz && \
    tar xzf hdf5-2.2.0.tar.gz && \
    cd hdf5-2.2.0 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_Fortran_COMPILER=ifx \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DBUILD_STATIC_LIBS=OFF \
        -DBUILD_TESTING=OFF \
        -DHDF5_ALLOW_UNSUPPORTED=ON \
        -DHDF5_BUILD_CPP_LIB=ON \
        -DHDF5_BUILD_EXAMPLES=OFF \
        -DHDF5_BUILD_FORTRAN=ON \
        -DHDF5_BUILD_HL_LIB=ON \
        -DHDF5_BUILD_TOOLS=ON \
        -DHDF5_ENABLE_PARALLEL=ON \
        -DHDF5_ENABLE_SZIP_ENCODING=ON \
        -DHDF5_ENABLE_SZIP_SUPPORT=ON \
        -DHDF5_ENABLE_ZLIB_SUPPORT=ON \
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
        -DMPI_Fortran_COMPILER=/usr/local/bin/mpifort \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install PETSc
# PETSc provides the linear and nonlinear solvers CarpetX/PDESolvers uses
# - depends on mpi_abi_wrapper, mpif, HDF5, and a BLAS/LAPACK
#
# The external packages below are the set an HPC site would normally enable,
# and the reason they are here is CarpetX/PDESolvers: it assembles the whole
# refinement hierarchy into one flat `MatCreateAIJ` and then chooses no
# preconditioner in the source at all, so whatever PETSc defaults to (ILU
# serially, block-Jacobi/ILU in parallel) is what ends up solving an elliptic
# problem that wants multigrid. Every package here is reachable from
# `PDESolvers::petsc_options` alone, with no code change:
#
#   hypre         BoomerAMG, the standard algebraic multigrid for exactly this
#                 (`-pc_type hypre -pc_hypre_type boomeramg`). This is the one
#                 that matters; the rest are supporting cast.
#   MUMPS         parallel sparse direct (`-pc_type lu
#                 -pc_factor_mat_solver_type mumps`) -- the fallback when AMG
#                 stalls, and the usual coarse-grid solver
#   SuperLU_DIST  the other parallel direct solver, same role
#   SuperLU       serial direct
#   SuiteSparse   UMFPACK / CHOLMOD / KLU, serial direct
#   ScaLAPACK     required by MUMPS
#   METIS         graph partitioning
#   PT-Scotch     parallel partitioning and ordering, for MUMPS and SuperLU_DIST
#
# - **PT-Scotch rather than ParMETIS**, the other obvious choice for that last
#   slot: ParMETIS is licensed for non-commercial research use only, so it must
#   not go into an image that gets pushed to a registry. PT-Scotch is CeCILL-C,
#   and its `libptscotchparmetisv3` answers the ParMETIS v3 API that MUMPS and
#   SuperLU_DIST actually ask for. Debian's PETSc makes the same substitution
#   for the same reason.
# - **--download-* rather than a stanza each.** PETSc builds these with the
#   compilers it was given, so they come out over the ABI automatically, and
#   their versions stay matched to the PETSc release, which is the combination
#   PETSc tests. Wiring eight packages up by hand and then persuading PETSc to
#   accept them is the part that goes wrong.
# - PT-Scotch generates its parsers at build time, which is why `flex` and
#   `bison` are in the apt list above. Without them configure dies deep in the
#   download stage with "PTScotch needs flex installed".
# - --with-hdf5-dir draws a warning -- "Using version 2.2.0 of package HDF5,
#   PETSc is tested with 1.14". Configure accepts it; nothing in Cactus uses
#   PetscViewerHDF5, so this is availability rather than something relied on.
#   Drop --with-hdf5-dir if it ever turns into a hard error.
# - PETSc's configure compiles and runs small MPI programs. Intel MPI does not
#   refuse to run as root the way Open MPI does, so the OMPI_ALLOW_RUN_AS_ROOT
#   pair `carpetx-cpu.dockerfile` needs here is not required; it is left
#   in place anyway, harmless, in case the wrapped MPI is ever changed back.
# - --with-mpiexec names the ABI prefix's own `mpiexec`, which forwards to the
#   wrapped MPI's launcher rather than being one. Naming the forwarder rather
#   than /usr/bin/mpiexec keeps the recorded launcher correct if this image is
#   later re-pointed at another MPI; PETSc would otherwise bake in whatever it
#   found on PATH.
# - the stock build installs `libpetsc.so` and `include/petsc.h` under one
#   prefix, which is exactly the layout `ExternalLibraries-PETSc/src/detect.sh`
#   looks for -- unlike Debian's `libpetsc_real.so` under /usr/lib/petsc.
RUN mkdir src && \
    (cd src && \
    export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 && \
    wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.25.5.tar.gz && \
    tar xzf petsc-3.25.5.tar.gz && \
    cd petsc-3.25.5 && \
    ./configure \
        --prefix=/usr/local \
        --download-hypre \
        --download-metis \
        --download-mumps \
        --download-ptscotch \
        --download-scalapack \
        --download-suitesparse \
        --download-superlu \
        --download-superlu_dist \
        --with-blaslapack-lib=-lopenblas \
        --with-cc=/usr/local/bin/mpicc \
        --with-cxx=/usr/local/bin/mpicxx \
        --with-debugging=0 \
        --with-fc=/usr/local/bin/mpifort \
        --with-hdf5-dir=/usr/local \
        --with-mpiexec=/usr/local/bin/mpiexec \
        --with-shared-libraries=1 \
        --with-x=0 \
        COPTFLAGS=-O2 CXXOPTFLAGS=-O2 FOPTFLAGS=-O2 \
        && \
    make && \
    make install && \
    true) && \
    rm -rf src

# Install MGARD
# MGARD is a lossy compression library
# Note: -DMGARD_ENABLE_SYCL=ON does not work
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
# - BUILD_H5Z_FILTER links HDF5, which drags MPI in behind it; hence the
#   compilers, even though SZ3 itself knows nothing about MPI
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
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install ADIOS2
# ADIOS2 is a parallel I/O library, comparable to HDF5
# - depends on blosc2
# - depends on MGARD
# - depends on mpi_abi_wrapper: ADIOS2 is built against the MPI ABI, not against
#   the system MPI directly. ADIOS2 inserts its own `cmake/` at the front of
#   CMAKE_MODULE_PATH and its `cmake/FindMPI.cmake` defers to CMake's bundled
#   module, so pointing CMAKE_MODULE_PATH at mpi_abi_wrapper's FindMPI shim
#   would be shadowed. Instead we name the wrappers, which CMake's own FindMPI
#   interrogates with `-showme:compile` / `-showme:link` -- flags that
#   mpi_abi_wrapper's `bin/mpicc` answers the way a real mpicc would.
# - ADIOS2_USE_MPI=ON rather than the AUTO default, so that a failure to find
#   the MPI ABI is an error here instead of a silently serial ADIOS2.
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
        -DADIOS2_USE_MPI=ON \
        -DADIOS2_USE_SZ3=ON \
        -DADIOS2_USE_ZFP=ON \
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
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
# Note: This assumes that the system has x86_64 CPUs
RUN mkdir src && \
    (cd src && \
    wget https://github.com/agenium-scale/nsimd/archive/refs/tags/v3.0.1.tar.gz && \
    tar xzf v3.0.1.tar.gz && \
    cd nsimd-3.0.1 && \
    cmake -B build -G Ninja \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
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
# - depends on mpi_abi_wrapper, and this is not optional: openPMD calls
#   find_package(MPI) itself, and left alone it finds the system Open MPI and
#   then fails to link ADIOS2 -- "undefined reference to
#   adios2::ADIOS::ADIOS(ompi_communicator_t*)", because ADIOS2's MPI_Comm is
#   the ABI's. The subtler half is HDF5: it is *also* mandatory that
#   find_package(HDF5) resolve to the ABI-built HDF5 in /usr/local, since
#   openPMD hands an MPI_Comm to H5Pset_fapl_mpio. Mixing there does not fail
#   to link -- C linkage hides the type mismatch -- it silently passes an ABI
#   handle to a library expecting an ompi_communicator_t*.
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
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
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
# - Silo has no MPI of its own, but it links HDF5 and so inherits HDF5's MPI;
#   the compilers are named for the same reason as in the SZ3 stanza
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
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Install Conduit
# Conduit defines a standard for laying out AMR data in a file.
# - depends on Silo
# - depends on mpi_abi_wrapper: `libconduit_relay_mpi*` and
#   `libconduit_blueprint_mpi` take an MPI_Comm across their public interface.
#   HDF5_DIR moves from /usr to /usr/local for the same reason -- Conduit's
#   relay writes HDF5, and the two must agree on what an MPI_Comm is.
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
        -DHDF5_DIR=/usr/local \
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
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
# - depends on mpi_abi_wrapper: AMReX is where CarpetX's own communication
#   happens, so its MPI_Comm has to be the same one Cactus calls MPI_Init on.
#   AMReX_MPI defaults to ON and is left that way; only the compilers change.
#   The nine `-DMPI_*` variables `carpetx-native-oneapi.dockerfile` sets here, and its
#   `rm -rf /opt/intel/oneapi/mpi`, are both gone: see the oneAPI notes in the
#   header.
# - `/usr/local` is appended to CMAKE_PREFIX_PATH. That variable is set here as
#   an explicit `/opt/intel/oneapi`-only list, which would otherwise hide the
#   ABI's own CMake package from any nested find_package.
# - **AMReX_SYCL_SPLIT_KERNEL=OFF is what makes this link finish.** The option
#   defaults to ON and adds `-fsycl-device-code-split=per_kernel`
#   (`Tools/CMake/AMReXSYCL.cmake`), which asks `sycl-post-link` for one device
#   image per kernel -- and AMReX generates a very large number of them through
#   templates. Observed on the 2026.1.0 toolkit: after `[121/122]`, i.e. at the
#   final link of `libamrex_3d.so`, `sycl-post-link` sat at 100% CPU with
#   12.7 GB resident and climbing, with no end in sight.
#
#   This is the "AMReX linker step is broken" note that pins
#   `carpetx-native-oneapi.dockerfile` to basekit 2025.2.2. That file's `FROM`
#   history marks 2025.3.0, 2025.3.1, 2025.3.2 and toolkit 2026.0.0 as broken;
#   2026.1.0 makes five releases, so the regression is not one to wait out.
#   Turning the split off is the one lever that acts directly on the tool that
#   is blowing up, rather than around it.
#
#   The cost is at run time, not build time: without the per-kernel split there
#   is one larger device image, so the first run JITs more than it strictly
#   needs. If a future toolkit fixes `sycl-post-link`, this line is the thing to
#   remove; if instead this proves insufficient, the fallback is
#   basekit 2025.2.2 -- it ships Intel MPI 2021.16 and exports the same
#   `C_INCLUDE_PATH`, so every other stanza in this file works there unchanged.
#
#   Unrelated but worth knowing while reading the link line: AMReX enables
#   `-flink-huge-device-code` when `CMAKE_BUILD_TYPE MATCHES "Debug"`, and
#   CMake's `MATCHES` is a regex, so `RelWithDebInfo` matches it by substring.
#   That large-device-code workaround is therefore already in effect here, by
#   accident rather than by choice.
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
        -DCMAKE_C_COMPILER=icx \
        -DCMAKE_CXX_COMPILER=icpx \
        -DCMAKE_PREFIX_PATH='/opt/intel/oneapi;/usr/local' \
        -DCMAKE_INSTALL_PREFIX=/usr/local \
        -DBUILD_SHARED_LIBS=ON \
        -DAMReX_FORTRAN=OFF \
        -DAMReX_FORTRAN_INTERFACES=OFF \
        -DAMReX_GPU_BACKEND=SYCL \
        -DAMReX_INTEL_ARCH=intel_gpu_pvc \
        -DAMReX_OMP=OFF \
        -DAMReX_PARTICLES=ON \
        -DAMReX_PRECISION="$precision" \
        -DAMReX_SYCL_SPLIT_KERNEL=OFF \
        -DMPI_C_COMPILER=/usr/local/bin/mpicc \
        -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
        && \
    cmake --build build && \
    cmake --install build && \
    true) && \
    rm -rf src

# Find libraries in /usr/local/lib64
RUN echo /usr/local/lib64 >/etc/ld.so.conf.d/usr-local-lib64.conf && \
    ldconfig
