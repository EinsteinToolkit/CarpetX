#ifndef CARPETX_CARPETX_BOUNDARIES_HXX
#define CARPETX_CARPETX_BOUNDARIES_HXX

#include "driver.hxx"
#include "loop_device.hxx"

#ifdef CCTK_DEBUG
#include <atomic>
#include <string>
#endif

namespace CarpetX {

namespace boundaries_detail {

// TODO: Move these functions to loop.hxx

template <typename F>
CCTK_ATTRIBUTE_NOINLINE __attribute__((__flatten__, __hot__)) void
loop_region(const F &f, const Arith::vect<int, dim> &imin,
            const Arith::vect<int, dim> &imax) {
  assert(all(imin < imax));

  const amrex::Box box(amrex::IntVect(imin[0], imin[1], imin[2]),
                       amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1));
  amrex::ParallelFor(box, [=] CCTK_DEVICE(const int i, const int j, const int k)
                              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                const Arith::vect<int, dim> p{i, j, k};
                                f(p);
                              });
}

} // namespace boundaries_detail
using namespace boundaries_detail;

////////////////////////////////////////////////////////////////////////////////

#ifdef CCTK_DEBUG

// INSTRUMENT (BUGFIX_TODO.md R2 / B10), debug builds only and OFF unless
// `CARPETX_LOG_BC_PASS` is set in the environment.
//
// B2's claim is that the two boundary passes are DISJOINT: pass 1
// (`skip_interpatch_corners`) skips every interpatch x outer corner, and
// pass 2 (`interpatch_corners_only`) writes exactly those and nothing else.
// The two witnesses B2 inherited cannot say that. CapyrX's `n_outer_skipped`
// (`interpolate.cxx:720`) is the INTERPOLATOR's count of `|O| + |C|` and B2
// does not change it; the post-pass-2 `contains_nan()` sweep
// (`schedule.cxx`) is an absence-of-NaN, i.e. a negative result with no
// positive control. This counts cells.
//
//   CARPETX_LOG_BC_PASS=1   per-`apply_boundary_conditions` totals
//   CARPETX_LOG_BC_PASS=2   the above, plus one line per face region
//
// Read it at OMP_NUM_THREADS=1: the counters are atomic and correct at any
// thread count, but the per-region lines are emitted from inside an
// `omp parallel` region and interleave.
//
// UNITS: cells, summed over the boxes of one MultiFab. A patch split into
// several AMReX boxes has its own ghost cells per box, so this is a count of
// (box, cell) pairs and only equals the unique cell count at one box per
// patch -- which every rig on BUGFIX_TODO.md's table pins.
struct bc_pass_census_t {
  std::atomic<long long> corner_cells_skipped;
  std::atomic<long long> corner_cells_kept;
  std::atomic<long long> other_cells_skipped;
  std::atomic<long long> other_cells_kept;
  std::atomic<long long> corner_regions_skipped;
  std::atomic<long long> corner_regions_kept;
  std::atomic<long long> other_regions_skipped;
  std::atomic<long long> other_regions_kept;
};
extern bc_pass_census_t bc_pass_census;

// 0 = off, 1 = totals, 2 = totals + per-region lines. Read from the
// environment once, on first use.
int bc_pass_census_level();
void bc_pass_census_reset();
void bc_pass_census_report(const std::string &groupname, int patch, int level,
                           bc_pass_t bc_pass);

#endif // CCTK_DEBUG

struct BoundaryCondition {
  const GHExt::PatchData::LevelData::GroupData &groupdata;
  const GHExt::PatchData &patchdata;
  const amrex::Geometry &geom;
  amrex::FArrayBox &dest;

  // Which subset of the faces/edges/corners this pass is responsible for. Set
  // once by the constructor; `apply_on_face` is a member template that takes no
  // function arguments (see the explicit instantiations in
  // `boundaries_impl_*.cxx`), so a runtime member is the only way to reach it
  // without changing 26 instantiation signatures. Nothing selects anything but
  // `all` yet.
  const bc_pass_t bc_pass;

  // Interior of the domain: Do not set any points in this region
  Arith::vect<int, dim> imin, imax;
  Arith::vect<CCTK_REAL, dim> xmin, xmax, dx;

  // Destination region
  Arith::vect<int, dim> dmin, dmax;

  Loop::GF3D2layout layout;
  CCTK_REAL *restrict destptr;

  BoundaryCondition(const GHExt::PatchData::LevelData::GroupData &groupdata,
                    amrex::FArrayBox &dest,
                    bc_pass_t bc_pass = bc_pass_t::all);

  BoundaryCondition(const BoundaryCondition &) = delete;
  BoundaryCondition(BoundaryCondition &&) = delete;
  BoundaryCondition &operator=(const BoundaryCondition &) = delete;
  BoundaryCondition &operator=(BoundaryCondition &&) = delete;

  void apply() const;

  template <int NI, int NJ, int NK> void apply_on_face() const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI>
  void apply_on_face_symbcx(const Arith::vect<int, dim> &bmin,
                            const Arith::vect<int, dim> &bmax) const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
            symmetry_t SCJ, boundary_t BCJ>
  void apply_on_face_symbcxy(const Arith::vect<int, dim> &bmin,
                             const Arith::vect<int, dim> &bmax) const;

  template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
            symmetry_t SCJ, boundary_t BCJ, symmetry_t SCK, boundary_t BCK>
  void apply_on_face_symbcxyz(const Arith::vect<int, dim> &bmin,
                              const Arith::vect<int, dim> &bmax) const;
};

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_BOUNDARIES_HXX
