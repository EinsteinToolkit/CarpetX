#include "driver.hxx"
#include "mpi_types.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <cctk_Parameters.h>

#include <AMReX_Functional.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Orientation.H>
#include <AMReX_ParReduce.H>

#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace CarpetX {

template <typename T, int D> MPI_Datatype reduction_mpi_datatype() {
  static MPI_Datatype datatype = MPI_DATATYPE_NULL;
  if (datatype == MPI_DATATYPE_NULL) {
    MPI_Type_contiguous(sizeof(reduction<T, D>) / sizeof(T),
                        mpi_datatype<T>::value, &datatype);
    char name[MPI_MAX_OBJECT_NAME];
    int namelen;
    MPI_Type_get_name(mpi_datatype<T>::value, name, &namelen);
    std::ostringstream buf;
    buf << "reduction<" << name << "," << D << ">";
    std::string newname = buf.str();
    MPI_Type_set_name(datatype, newname.c_str());
    MPI_Type_commit(&datatype);
  }
  return datatype;
}

namespace {
template <typename T>
void mpi_reduce_typed(const void *restrict x0, void *restrict y0,
                      const int *restrict length) {
  const T *restrict x = static_cast<const T *>(x0);
  T *restrict y = static_cast<T *>(y0);
#pragma omp simd
  for (int i = 0; i < *length; ++i)
    y[i] += x[i];
}

void mpi_reduce(void *restrict x, void *restrict y, int *restrict length,
                MPI_Datatype *restrict datatype) {
  // Analyze MPI datatype
  int num_integers, num_addresses, num_datatypes, combiner;
  MPI_Type_get_envelope(*datatype, &num_integers, &num_addresses,
                        &num_datatypes, &combiner);
  assert(combiner == MPI_COMBINER_CONTIGUOUS);
  assert(num_integers == 1);
  assert(num_datatypes == 1);
  std::vector<int> integers(num_integers);
  std::vector<MPI_Aint> addresses(num_addresses);
  std::vector<MPI_Datatype> datatypes(num_datatypes);
  MPI_Type_get_contents(*datatype, num_integers, num_addresses, num_datatypes,
                        integers.data(), addresses.data(), datatypes.data());
  MPI_Datatype inner_datatype = datatypes.at(0);
  if (inner_datatype == MPI_FLOAT)
    return mpi_reduce_typed<reduction<float, dim> >(x, y, length);
  if (inner_datatype == MPI_DOUBLE)
    return mpi_reduce_typed<reduction<double, dim> >(x, y, length);
  if (inner_datatype == MPI_LONG_DOUBLE)
    return mpi_reduce_typed<reduction<long double, dim> >(x, y, length);
  CCTK_ERROR("Unsupported MPI datatype");
}
} // namespace

MPI_Op reduction_mpi_op() {
  static MPI_Op op = MPI_OP_NULL;
  if (op == MPI_OP_NULL)
    MPI_Op_create(mpi_reduce, 1 /*commutes*/, &op);
  return op;
}

////////////////////////////////////////////////////////////////////////////////

namespace {

// Number of active octants (out of 8) surrounding a grid point, i.e. eight
// times the volume weight with which the point enters a reduction.
//
// For vertex-centred directions, points on the boundary of the valid region
// and points at refinement boundaries are counted with a reduced weight.
// `masked` reports whether a point is covered by a finer level.
//
// This is shared between the CPU and the GPU code path so that both compute
// the same weights.
template <typename Mask>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
active_octants(const vect<int, dim> &ipos, const vect<int, dim> &indextype,
               const vect<int, dim> &imin, const vect<int, dim> &imax,
               const Mask &masked) {
  constexpr vect<vect<int, dim>, dim> di = {vect<int, dim>::unit(0),
                                            vect<int, dim>::unit(1),
                                            vect<int, dim>::unit(2)};
  constexpr vect<vect<vect<int, dim>, dim>, 2> dirs = {-di, +di};

  // The octants adjacent to each of the six faces
  constexpr unsigned char faces[2][dim] = {
      {0b01010101, 0b00110011, 0b00001111},
      {0b10101010, 0b11001100, 0b11110000}};

  const vect<vect<int, dim>, 2> ibnd = {imin, imax - 1};

  // For vertex-centred grids, ensure that points at the outer boundary are
  // counted with a weight of 1/2
  unsigned char outer_active = 0b11111111;
  for (int f = 0; f < 2; ++f)
    for (int d = 0; d < dim; ++d)
      if (indextype[d] == 0 && ipos[d] == ibnd[f][d])
        outer_active &= ~faces[f][d];

  // For vertex-centred grids, ensure that points at refinement boundaries are
  // counted with a weight of 1/2
  unsigned char inner_active;
  if (!masked(ipos)) {
    inner_active = 0b11111111;
  } else {
    inner_active = 0b00000000;
    for (int f = 0; f < 2; ++f)
      for (int d = 0; d < dim; ++d)
        if (indextype[d] == 0 &&
            !(ipos[d] != ibnd[f][d] && masked(ipos + dirs[f][d])))
          inner_active |= faces[f][d];
  }

  assert((~outer_active & ~inner_active & 0b11111111) == 0);

  const unsigned char active = outer_active & inner_active;
  int count = 0;
  for (int b = 0; b < 8; ++b)
    count += (active >> b) & 1;
  return count;
}

// Like `active_octants`, but for a grid function array on the GPU, where the
// bounds of the valid region are not available directly: the array covers the
// fab box, which is the valid box grown by `ng`.
//
// We access the bounds via `lbound`/`ubound` instead of `Array4::begin`/`end`
// because AMReX 26.02 removes those members. `ubound` is inclusive.
template <typename T, typename Mask>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
active_octants_in_fab(const amrex::Array4<const T> &restrict vars,
                      const vect<int, dim> &ipos,
                      const vect<int, dim> &indextype, const vect<int, dim> &ng,
                      const Mask &masked) {
  const amrex::Dim3 lo = amrex::lbound(vars);
  const amrex::Dim3 hi = amrex::ubound(vars);
  const vect<int, dim> imin{lo.x + ng[0], lo.y + ng[1], lo.z + ng[2]};
  const vect<int, dim> imax{hi.x + 1 - ng[0], hi.y + 1 - ng[1],
                            hi.z + 1 - ng[2]};
  return active_octants(ipos, indextype, imin, imax, masked);
}

template <typename T>
reduction<T, dim>
reduce_array(const amrex::Array4<const T> &restrict vars, const int n,
             const vect<int, dim> &tmin, const vect<int, dim> &tmax,
             const vect<int, dim> &indextype, const vect<int, dim> &imin,
             const vect<int, dim> &imax,
             const amrex::Array4<const int> *restrict const finemask,
             const vect<T, dim> &x0, const vect<T, dim> &dx) {
  const auto masked = [&](const vect<int, dim> &ipos) {
    return finemask && (*finemask)(ipos[0], ipos[1], ipos[2]);
  };

  const CCTK_REAL dV = prod(dx);

  // Use per-loop reduction objects to reduce round-off error
  reduction<T, dim> redk;
  // TODO: use loop.hxx code to loop over grid
  for (int k = tmin[2]; k < tmax[2]; ++k) {
    reduction<T, dim> redj;
    for (int j = tmin[1]; j < tmax[1]; ++j) {
      reduction<T, dim> redi;
      for (int i = tmin[0]; i < tmax[0]; ++i) {
        const vect<int, dim> ipos = {i, j, k};

        const int nactive = active_octants(ipos, indextype, imin, imax, masked);
        if (nactive > 0) {
          const T W = nactive / T(8);

          const vect<T, dim> x = x0 + ipos * dx;
          redi += reduction<T, dim>(x, W * dV, vars(i, j, k, n));
        }
      }
      redj += redi;
    }
    redk += redj;
  }
  return redk;
}
} // namespace

reduction<CCTK_REAL, dim> reduce(int gi, int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;

  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  assert(group.grouptype == CCTK_GF);

  reduction<CCTK_REAL, dim> red;
  // TODO: Parallelize over patches and levels
  for (auto &restrict patchdata : ghext->patchdata) {
    for (auto &restrict leveldata : patchdata.leveldata) {
      const auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const amrex::MultiFab &mfab = *groupdata.mfab.at(tl);
      std::unique_ptr<amrex::iMultiFab> finemask_imfab;

      warn_if_invalid(groupdata, vi, tl, make_valid_int(),
                      []() { return "Before reduction"; });

      const vect<int, dim> indextype = groupdata.indextype;

      const auto &restrict geom = patchdata.amrcore->Geom(leveldata.level);
      const CCTK_REAL *restrict const x01 = geom.ProbLo();
      const CCTK_REAL *restrict const dx1 = geom.CellSize();
      const vect<CCTK_REAL, dim> dx = {dx1[0], dx1[1], dx1[2]};
      const vect<CCTK_REAL, dim> x0v = {x01[0], x01[1], x01[2]};
      const auto x0 = x0v + indextype * dx / 2;

      const int fine_level = leveldata.level + 1;
      if (fine_level < int(patchdata.leveldata.size())) {
        const auto &restrict fine_leveldata =
            patchdata.leveldata.at(fine_level);
        const auto &restrict fine_groupdata = *fine_leveldata.groupdata.at(gi);
        const amrex::MultiFab &fine_mfab = *fine_groupdata.mfab.at(tl);

        const amrex::IntVect reffact{2, 2, 2};

        finemask_imfab = std::make_unique<amrex::iMultiFab>(makeFineMask(
            mfab, fine_mfab.boxArray(), reffact, geom.periodicity(),
            /*coarse value*/ 0, /* fine value */ 1));
      }

#ifndef AMREX_USE_GPU
      // CPU

      auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
      // TODO: check that multi-threading actually helps (and we are
      // not dominated by memory latency)
      // TODO: document required version of OpenMP to use custom reductions
#ifdef __NVCOMPILER
#pragma omp parallel
      {
        auto &outer = red;
        reduction<CCTK_REAL, dim> red;
#else
#pragma omp parallel reduction(reduction : red)
#endif
        for (amrex::MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
          const amrex::Box &bx = mfi.tilebox(); // current tile (without ghosts)
          const vect<int, dim> tmin{bx.smallEnd(0), bx.smallEnd(1),
                                    bx.smallEnd(2)};
          const vect<int, dim> tmax{bx.bigEnd(0) + 1, bx.bigEnd(1) + 1,
                                    bx.bigEnd(2) + 1};
          const amrex::Box &vbx =
              mfi.validbox(); // interior region (without ghosts)
          const vect<int, dim> imin{vbx.smallEnd(0), vbx.smallEnd(1),
                                    vbx.smallEnd(2)};
          const vect<int, dim> imax{vbx.bigEnd(0) + 1, vbx.bigEnd(1) + 1,
                                    vbx.bigEnd(2) + 1};

          const amrex::Array4<const CCTK_REAL> &vars = mfab.array(mfi);

          std::unique_ptr<amrex::Array4<const int> > finemask;
          if (finemask_imfab) {
            finemask = std::make_unique<amrex::Array4<const int> >(
                finemask_imfab->array(mfi));
            // Ensure the mask has the correct size
            assert(finemask->begin.x == vars.begin.x);
            assert(finemask->begin.y == vars.begin.y);
            assert(finemask->begin.z == vars.begin.z);
            assert(finemask->end.x == vars.end.x);
            assert(finemask->end.y == vars.end.y);
            assert(finemask->end.z == vars.end.z);
          }

          red += reduce_array(vars, vi, tmin, tmax, indextype, imin, imax,
                              finemask.get(), x0, dx);
        }
#ifdef __NVCOMPILER
#pragma omp critical(CarpetX_reduce)
        outer += red;
      }
#endif

#else
      // GPU
      //
      // Reduce over all boxes of the level in a single kernel launch.
      // `amrex::ParReduce` iterates over all local boxes of the `MultiFab`
      // and passes the local box index as first argument.

      using tuple_type = reduction<CCTK_REAL, dim>::tuple_type;

      const amrex::IntVect nghosts = mfab.nGrowVect();
#ifdef CCTK_DEBUG
      // We derive the valid box from the fab box below, which requires that
      // every fab is grown uniformly
      for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi)
        assert(mfi.fabbox() == amrex::grow(mfi.validbox(), nghosts));
#endif
      const vect<int, dim> ng{nghosts[0], nghosts[1], nghosts[2]};

      const amrex::MultiArray4<const CCTK_REAL> data_ma = mfab.const_arrays();
      const bool have_mask = bool(finemask_imfab);
      const amrex::MultiArray4<const int> mask_ma =
          have_mask ? finemask_imfab->const_arrays()
                    : amrex::MultiArray4<const int>();

      const CCTK_REAL dV = prod(dx);
      const int vi1 = vi;
      const vect<int, dim> indextype1 = indextype;
      const vect<CCTK_REAL, dim> x01 = x0;
      const vect<CCTK_REAL, dim> dx1 = dx;

      // First pass: the sums and the extrema
      reduction<CCTK_REAL, dim> level_red = amrex::ParReduce(
          amrex::TypeList<
              amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpSum,
              amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpMax,
              amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
              amrex::ReduceOpSum, amrex::ReduceOpSum>{},
          amrex::TypeList<CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL,
                          CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL, CCTK_REAL,
                          CCTK_REAL>{},
          mfab,
          amrex::IntVect(0), // interior region, without ghosts
          [=] CCTK_DEVICE(const int b, const int i, const int j,
                          const int k) noexcept -> tuple_type {
            const amrex::Array4<const CCTK_REAL> &restrict vars = data_ma[b];
            const auto masked = [=](const vect<int, dim> &ipos) {
              return have_mask && mask_ma[b](ipos[0], ipos[1], ipos[2]);
            };

            const vect<int, dim> ipos{i, j, k};
            const int nactive =
                active_octants_in_fab(vars, ipos, indextype1, ng, masked);
            if (nactive == 0)
              return tuple_type(reduction<CCTK_REAL, dim>());
            const CCTK_REAL W = nactive / CCTK_REAL(8);

            const vect<CCTK_REAL, dim> x = x01 + ipos * dx1;
            return tuple_type(
                reduction<CCTK_REAL, dim>(x, W * dV, vars(i, j, k, vi1)));
          });

      // Second pass: the locations of the extrema. AMReX's reduction
      // operators cannot express argmin/argmax, so we reduce the smallest
      // linear index among the points attaining the extremum, and convert
      // that index back to a position below.
      //
      // The extrema are those of this patch and level; `reduction::operator+`
      // keeps the locations belonging to the smaller minimum resp. the larger
      // maximum when the levels, patches, and MPI processes are combined.
      //
      // Points attaining the extremum are not unique. Which one is found is
      // unspecified, as it is for the CPU code path, where the order in which
      // the OpenMP threads combine their results is unspecified as well.
      if (level_red.vol > 0) {
        // The index space holding all valid points of this level
        const amrex::Box lvlbox = mfab.boxArray().minimalBox();
        const amrex::Dim3 blo = amrex::lbound(lvlbox);
        const amrex::Long nx = lvlbox.length(0);
        const amrex::Long ny = lvlbox.length(1);
        constexpr amrex::Long noindex = std::numeric_limits<amrex::Long>::max();

        // If the data contain nans then no point compares equal to the
        // extremum, and the locations remain unknown (nan)
        const CCTK_REAL lmin = level_red.min;
        const CCTK_REAL lmax = level_red.max;

        const amrex::GpuTuple<amrex::Long, amrex::Long> locs = amrex::ParReduce(
            amrex::TypeList<amrex::ReduceOpMin, amrex::ReduceOpMin>{},
            amrex::TypeList<amrex::Long, amrex::Long>{}, mfab,
            amrex::IntVect(0), // interior region, without ghosts
            [=] CCTK_DEVICE(const int b, const int i, const int j,
                            const int k) noexcept
                -> amrex::GpuTuple<amrex::Long, amrex::Long> {
              const amrex::Array4<const CCTK_REAL> &restrict vars = data_ma[b];
              const auto masked = [=](const vect<int, dim> &ipos) {
                return have_mask && mask_ma[b](ipos[0], ipos[1], ipos[2]);
              };

              const vect<int, dim> ipos{i, j, k};
              const int nactive =
                  active_octants_in_fab(vars, ipos, indextype1, ng, masked);
              if (nactive == 0)
                return {noindex, noindex};

              const amrex::Long idx =
                  (amrex::Long(k - blo.z) * ny + (j - blo.y)) * nx +
                  (i - blo.x);
              const CCTK_REAL val = vars(i, j, k, vi1);
              return {val == lmin ? idx : noindex, val == lmax ? idx : noindex};
            });

        const auto decode = [&](const amrex::Long idx) {
          if (idx == noindex)
            return vect<CCTK_REAL, dim>::pure(0.0 / 0.0);
          const vect<int, dim> ipos{int(idx % nx) + blo.x,
                                    int(idx / nx % ny) + blo.y,
                                    int(idx / (nx * ny)) + blo.z};
          return x01 + ipos * dx1;
        };
        level_red.minloc = decode(amrex::get<0>(locs));
        level_red.maxloc = decode(amrex::get<1>(locs));
      }

      red += level_red;

#endif
    }
  }

  MPI_Datatype datatype = reduction_mpi_datatype<CCTK_REAL, dim>();
  MPI_Op op = reduction_mpi_op();
  MPI_Allreduce(MPI_IN_PLACE, &red, 1, datatype, op, MPI_COMM_WORLD);

  return red;
}

} // namespace CarpetX
