#include "driver.hxx"
#include "interp.hxx"
#include "mpi_types.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <defs.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Functions.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Table.h>

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

namespace CarpetX {
using Arith::pown;

namespace {

// Interpolate a grid function at one point, dimensionally recursive
template <typename T, int order, int centering> struct interpolator {
  static constexpr vect<bool, dim> indextype{(centering & 0b100) != 0,
                                             (centering & 0b010) != 0,
                                             (centering & 0b001) != 0};

  const GridDescBase &grid;
#ifdef CCTK_DEBUG
  const int gi;
#endif
  const int vi;
#ifdef CCTK_DEBUG
  const int patch;
  const int level;
#endif
  const amrex::Array4<const T> &vars;
  const vect<int, dim> &derivs;
  // Allow outer boundaries as interpolation sources
  const vect<vect<bool, dim>, 2> allowed_boundaries;

  interpolator(const GridDescBase &grid, CCTK_ATTRIBUTE_UNUSED const int gi,
               const int vi, const int patch, const int level,
               const amrex::Array4<const T> &vars, const vect<int, dim> &derivs,
               const vect<vect<bool, dim>, 2> &allowed_boundaries)
      : grid(grid),
#ifdef CCTK_DEBUG
        gi(gi),
#endif
        vi(vi),
#ifdef CCTK_DEBUG
        patch(patch), level(level),
#endif
        vars(vars), derivs(derivs), allowed_boundaries(allowed_boundaries) {
  }

  static constexpr T eps() {
    using std::pow;
    return pow(std::numeric_limits<T>::epsilon(), T(3) / 4);
  }

  // TODO: Check whether interpolated variables are valid

  // Base case: only access a grid point
  template <int dir>
  std::enable_if_t<(dir == -1), T>
  interpolate(const vect<int, dim> &i, const vect<CCTK_REAL, dim> &di) const {
    const amrex::IntVect j(i[0] + vars.begin.x, i[1] + vars.begin.y,
                           i[2] + vars.begin.z);
#ifdef CCTK_DEBUG
    assert(vars.contains(j[0], j[1], j[2]));
#endif
    const T val = vars(j, vi);
#ifdef CCTK_DEBUG
    using std::isfinite;
    if (!(isfinite(val))) {
      std::cerr << "!isfinite gi=" << gi
                << " groupname=" << CCTK_FullGroupName(gi) << " vi=" << vi
                << " i=" << i << " di=" << di << " val=" << val << "\n";
      for (int c = -1; c <= +1; ++c)
        for (int b = -1; b <= +1; ++b)
          for (int a = -1; a <= +1; ++a)
            if (vars.contains(j[0] + a, j[1] + b, j[2] + c))
              std::cerr << "  val[" << a << "," << b << "," << c
                        << "]=" << vars(j[0] + a, j[1] + b, j[2] + c, vi)
                        << "\n";
    }
    assert(isfinite(val));
#endif
    return val;
  }

  // General case: interpolate in one direction, then recurse
  template <int dir>
  std::enable_if_t<(dir >= 0), T>
  interpolate(const vect<int, dim> &i, const vect<CCTK_REAL, dim> &di) const {
    static_assert(dir < dim);
    const auto DI = vect<int, dim>::unit(dir);

    // Ignore the centering for interpolation
    // switch ((centering >> (2 - dir)) & 1)

    const T x = di[dir] - order / T(2);
    // #ifdef CCTK_DEBUG
    //     using std::fabs;
    //     assert(fabs(x) <= T(0.5) + eps());
    // #endif

    switch (order) {
    case 0: {
      const T y0 = interpolate<dir - 1>(i, di);
      switch (derivs[dir]) {
      case 0:
        return y0;
      case 1:
        return 0;
      case 2:
        return 0;
      }
    }
    case 1: {
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      switch (derivs[dir]) {
      case 0:
        return (1 / T(2) - x) * y0 + (1 / T(2) + x) * y1;
      case 1:
        return (-y0 + y1) / grid.dx[dir];
      case 2:
        return 0;
      }
    }
    case 2: {
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (-1 / T(2) * x + 1 / T(2) * pown(x, 2)) * y0 +
               (1 - pown(x, 2)) * y1 +
               (1 / T(2) * x + 1 / T(2) * pown(x, 2)) * y2;
      case 1:
        return ((-1 / T(2) + x) * y0 - 2 * x * y1 + (1 / T(2) + x) * y2) /
               grid.dx[dir];
      case 2:
        return (y0 - 2 * y1 + y2) / pown(grid.dx[dir], 2);
      }
    }
    case 3: {
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      const T y3 = interpolate<dir - 1>(i + 3 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (-1 / T(16) + 1 / T(24) * x + 1 / T(4) * pown(x, 2) -
                1 / T(6) * pown(x, 3)) *
                   y0 +
               (9 / T(16) - 9 / T(8) * x - 1 / T(4) * pown(x, 2) +
                1 / T(2) * pown(x, 3)) *
                   y1 +
               (9 / T(16) + 9 / T(8) * x - 1 / T(4) * pown(x, 2) -
                1 / T(2) * pown(x, 3)) *
                   y2 +
               (-1 / T(16) - 1 / T(24) * x + 1 / T(4) * pown(x, 2) +
                1 / T(6) * pown(x, 3)) *
                   y3;
      case 1:
        return ((1 / T(24) + 1 / T(2) * x - 1 / T(2) * pown(x, 2)) * y0 +
                (-9 / T(8) - 1 / T(2) * x + 3 / T(2) * pown(x, 2)) * y1 +
                (9 / T(8) - 1 / T(2) * x - 3 / T(2) * pown(x, 2)) * y2 +
                (-1 / T(24) + 1 / T(2) * x + 1 / T(2) * pown(x, 2)) * y3) /
               grid.dx[dir];
      case 2:
        return ((1 / T(2) - x) * y0 + (-1 / T(2) + 3 * x) * y1 +
                (-1 / T(2) - 3 * x) * y2 + (1 / T(2) + x) * y3) /
               pown(grid.dx[dir], 2);
      }
    }
    case 4: {
      const T y0 = interpolate<dir - 1>(i, di);
      const T y1 = interpolate<dir - 1>(i + DI, di);
      const T y2 = interpolate<dir - 1>(i + 2 * DI, di);
      const T y3 = interpolate<dir - 1>(i + 3 * DI, di);
      const T y4 = interpolate<dir - 1>(i + 4 * DI, di);
      switch (derivs[dir]) {
      case 0:
        return (1 / T(12) * x - 1 / T(24) * pown(x, 2) -
                1 / T(12) * pown(x, 3) + 1 / T(24) * pown(x, 4)) *
                   y0 +
               (-2 / T(3) * x + 2 / T(3) * pown(x, 2) + 1 / T(6) * pown(x, 3) -
                1 / T(6) * pown(x, 4)) *
                   y1 +
               (1 - 5 / T(4) * pown(x, 2) + 1 / T(4) * pown(x, 4)) * y2 +
               (2 / T(3) * x + 2 / T(3) * pown(x, 2) - 1 / T(6) * pown(x, 3) -
                1 / T(6) * pown(x, 4)) *
                   y3 +
               (-1 / T(12) * x - 1 / T(24) * pown(x, 2) +
                1 / T(12) * pown(x, 3) + 1 / T(24) * pown(x, 4)) *
                   y4;
      case 1:
        return ((1 / T(12) - 1 / T(12) * x - 1 / T(4) * pown(x, 2) +
                 1 / T(6) * pown(x, 3)) *
                    y0 +
                (-2 / T(3) + 4 / T(3) * x + 1 / T(2) * pown(x, 2) -
                 2 / T(3) * pown(x, 3)) *
                    y1 +
                (-5 / T(2) * x + pown(x, 3)) * y2 +
                (2 / T(3) + 4 / T(3) * x - 1 / T(2) * pown(x, 2) -
                 2 / T(3) * pown(x, 3)) *
                    y3 +
                (-1 / T(12) - 1 / T(12) * x + 1 / T(4) * pown(x, 2) +
                 1 / T(6) * pown(x, 3)) *
                    y4) /
               grid.dx[dir];
      case 2:
        return ((-1 / T(12) - 1 / T(2) * x + 1 / T(2) * pown(x, 2)) * y0 +
                (4 / T(3) + x - 2 * pown(x, 2)) * y1 +
                (-5 / T(2) + 3 * pown(x, 2)) * y2 +
                (4 / T(3) - x - 2 * pown(x, 2)) * y3 +
                (-1 / T(12) + 1 / T(2) * x + 1 / T(2) * pown(x, 2)) * y4) /
               pown(grid.dx[dir], 2);
      }
    }
    default:
      assert(0);
    } // switch order
  }

  template <typename Particles>
  void interpolate3d(const Particles &particles,
                     std::vector<T> &varresult) const {
    const auto x0 = grid.x0 + (2 * grid.lbnd - !indextype) * grid.dx / 2;
    // const auto x1 = x0 + (grid.lsh - 1 - indextype) * grid.dx;
    const auto dx = grid.dx;

    assert(vars.end.x - vars.begin.x == grid.lsh[0] - indextype[0]);
    assert(vars.end.y - vars.begin.y == grid.lsh[1] - indextype[1]);
    assert(vars.end.z - vars.begin.z == grid.lsh[2] - indextype[2]);

    // We assume that the input is synchronized, i.e. that all ghost
    // zones are valid, but all outer boundaries are invalid.

    // The point must lie inside the domain. At outer boundaries the
    // point may be in the boundary region, but at ghost boundaries
    // the point cannot be in the ghost region. We define as "ghost"
    // region here the interpolation stencil size which is `order /
    // 2`.
    // const auto x0_allowed =
    //     x0 +
    //     (!allowed_boundaries[0] * grid.nghostzones + (order - 1) / T(2)) * dx
    //     - eps();
    // const auto x1_allowed =
    //     x1 -
    //     (!allowed_boundaries[1] * grid.nghostzones + (order - 1) / T(2)) * dx
    //     + eps();

    // The allowed index range is [i0, i1)
    const auto i0_allowed = !allowed_boundaries[0] * grid.nghostzones;
    const auto i1_allowed =
        grid.lsh - (!allowed_boundaries[1] * grid.nghostzones + order);

    const int np = int(varresult.size());

#pragma omp parallel for simd
    for (int n = 0; n < np; ++n) {
      const vect<T, dim> x{particles[n].rdata(0), particles[n].rdata(1),
                           particles[n].rdata(2)};
      // // Ensure the point is inside the domain
      // assert(all(x >= x0_allowed && x <= x1_allowed));

      // Find stencil anchor (i.e. the leftmost stencil point)
      const auto qi = (x - x0) / dx;
      const auto lrint1 = [](auto a) {
        using std::lrint;
        return int(lrint(a));
      };
      auto i = fmap(lrint1, qi - order / T(2));
      auto di = qi - i;
      // Consistency check
      assert(all(i >= 0 && i + order < grid.lsh));

      // Push point away from boundaries if they are just a little outside
      for (int d = 0; d < dim; ++d) {
        if (i[d] + order / 2 < i0_allowed[d] &&
            di[d] - order / T(2) >= +T(0.5) - eps()) {
          i[d] += 1;
          di[d] -= 1;
        }
        if (i[d] >= i1_allowed[d] && di[d] - order / T(2) <= -T(0.5) + eps()) {
          i[d] -= 1;
          di[d] += 1;
        }
      }

      // Avoid points on boundaries
      const bool is_allowed = all(i >= i0_allowed && i < i1_allowed);

#ifdef CCTK_DEBUG
      // mp_slave_2.md §4 instrumentation: log the donor (patch, anchor index)
      // chosen for every interpolation query point, so a failing query's
      // source index (idata(1), stable across both SyncGroupsByDirI passes)
      // can be cross-referenced against CapyrX_MultiPatch's victim dump to
      // pin whether the donor cell is a genuine ghost zone (i < nghostzones)
      // or an interior overlap-band cell that slave_overlap can mutate
      // (mp_slave_2.md §5). Opt-in via env var to avoid flooding every debug
      // run.
      //
      // BUGFIX_TODO.md [O8] / B10, added for step A9: four more fields, so that
      // the question "does this anchor's tensor-product support {i..i+order}
      // reach the DONOR patch's outer-boundary ghosts?" is answered from the
      // dump instead of from a hand derivation of the donor's shape. The
      // unapplied form of exactly such a derivation is what forced the
      // retraction of the whole coordinate-staleness result (mp_slave_5.md),
      // so the bound is now printed by the same code that enforces it.
      //
      //   index  -- GridDescBase::component IS mfp.index() (schedule.cxx:154).
      //             It therefore joins to the `index=` field of CapyrX's
      //             VICTIM/GHOSTSKIP lines, NOT to their `component=`, which is
      //             a per-level tile counter in active_levels_t::loop_parallel
      //             / loop_serially (schedule.cxx:868, :893). Do not confuse
      //             the two: with tiling enabled one `index` can carry several
      //             `component`s.
      //   lsh    -- this box's allocated shape, the upper bound the support has
      //             to be tested against (i1_allowed above subtracts `order`
      //             AND nghostzones from it, so i1_allowed is NOT that bound).
      //   bbox   -- true where the box face touches the patch's AMReX domain,
      //             i.e. where a ghost plane is O, I or C rather than G.
      //   allowed-- patch_allowed_boundaries as this interpolator received it:
      //             the per-face anchor policy, which is `is_outer_boundary`
      //             on a domain face (CapyrX multipatch.cxx:75-91 negated at
      //             CapyrX interpolate.cxx:516) and is forced true on a
      //             non-domain face. A face bounds an O region iff
      //             bbox && allowed; bbox && !allowed is an interpatch face
      //             (I/C); !bbox is an intra-patch box ghost (G).
      //
      // Fields are APPENDED, so a parser anchored on the old form sees the new
      // lines as malformed rather than silently mis-parsing them; A8's
      // a8_report.py was widened to accept both.
      //
      // READING THIS STREAM: `OMP_NUM_THREADS=1` AND `MPIEXEC=none`, and they
      // are two different preconditions, not one restated (BUGFIX_TODO.md B10).
      // This is a bare `std::cerr` chain, one `<<` per field, emitted from
      // inside an `omp parallel for`: at 16 threads the lines INTERLEAVE, and
      // the damage is invisible to `wc -l` because every thread still writes
      // its own newline -- 19.8M of 21.2M lines malformed at 16 threads and 0
      // at 1 (`[P26]`, `[P32]`).  Separately, at ONE thread, mpiexec's stderr
      // forwarder DROPS BYTES from the head of a line under an instrument
      // flood (`[P31]`), which no amount of atomicity here would fix.  B10
      // considered building each line in an `ostringstream` and emitting it
      // with a single `<<`; rejected, because it removes only the first hazard,
      // leaves both operational gates in place unchanged, and would make every
      // stream A8 and A9 have already measured incomparable with the next one.
      // `tools/run_split.sh` pins both.
      //
      // COST WHEN ON, measured (evidence/fix/b10/before/h_report.txt): 149664
      // lines on `color.par` and 378624 on `color_ghost.par`, at
      // `cctk_itlast = 0`.  It shares one environment variable with seven other
      // blocks, so it cannot be switched on alone (`[P27]`).
      {
        static const bool log_donors = std::getenv("CAPYRX_LOG_DONORS") != nullptr;
        if (log_donors) {
          std::cerr << "DONOR patch=" << patch << " level=" << level
                    << " gi=" << gi << " groupname=" << CCTK_FullGroupName(gi)
                    << " vi=" << vi << " n=" << particles[n].idata(1)
                    << " i=(" << i[0] << "," << i[1] << "," << i[2] << ")"
                    << " nghostzones=(" << grid.nghostzones[0] << ","
                    << grid.nghostzones[1] << "," << grid.nghostzones[2] << ")"
                    << " is_allowed=" << is_allowed
                    << " index=" << grid.component
                    << " lsh=(" << grid.lsh[0] << "," << grid.lsh[1] << ","
                    << grid.lsh[2] << ")"
                    << " bbox_lo=(" << grid.bbox[0][0] << "," << grid.bbox[0][1]
                    << "," << grid.bbox[0][2] << ")"
                    << " bbox_hi=(" << grid.bbox[1][0] << "," << grid.bbox[1][1]
                    << "," << grid.bbox[1][2] << ")"
                    << " allowed_lo=(" << allowed_boundaries[0][0] << ","
                    << allowed_boundaries[0][1] << ","
                    << allowed_boundaries[0][2] << ")"
                    << " allowed_hi=(" << allowed_boundaries[1][0] << ","
                    << allowed_boundaries[1][1] << ","
                    << allowed_boundaries[1][2] << ")"
                    << "\n";
        }
      }
#endif

      if (!is_allowed) {
        CCTK_VERROR("Interpolation anchor is not allowed, as it lies outside "
                    "of the interior region: "
                    "patch = %d "
                    "n = %d "
                    "i = (%d, %d, %d) "
                    "i0_allowed = (%d, %d, %d) "
                    "i1_allowed = (%d, %d, %d) "
                    "x = (%f, %f, %f).",
                    grid.patch, n, i[0], i[1], i[2], i0_allowed[0],
                    i0_allowed[1], i0_allowed[2], i1_allowed[0], i1_allowed[1],
                    i1_allowed[2], x[0], x[1], x[2]);
      }

      assert(is_allowed);

      const T res = !is_allowed ? -2 : interpolate<dim - 1>(i, di);

      varresult[n] = res;
    }
  }
};

} // namespace

int InterpLocalUniform(int /*N_dims*/, int /*param_table_handle*/,
                       /***** coordinate system *****/
                       const CCTK_REAL /*coord_origin*/[],
                       const CCTK_REAL /*coord_delta*/[],
                       /***** interpolation points *****/
                       int /*N_interp_points*/, int /*interp_coords_type_code*/,
                       const void *const /*interp_coords*/[],
                       /***** input arrays *****/
                       int /*N_input_arrays*/,
                       const CCTK_INT /*input_array_dims*/[],
                       const CCTK_INT /*input_array_type_codes*/[],
                       const void *const /*input_arrays*/[],
                       /***** output arrays *****/
                       int /*N_output_arrays*/,
                       const CCTK_INT /*output_array_type_codes*/[],
                       void *const /*output_arrays*/[]) {
  CCTK_ERROR("Dummy InterpLocalUniform function called");
}

extern "C" CCTK_INT CarpetX_InterpGridArrays(
    cGH const *const cctkGH, int const N_dims, int const local_interp_handle,
    int const param_table_handle, int const coord_system_handle,
    int const N_interp_points, int const interp_coords_type_code,
    void const *const coords[], int const N_input_arrays,
    CCTK_INT const input_array_variable_indices[], int const N_output_arrays,
    CCTK_INT const output_array_type_codes[], void *const output_arrays[]) {
  /* TODO: verify that the interface with SymmetryInterpolate can be simply
     copied from Carpet like below */
  //  if (CCTK_IsFunctionAliased("SymmetryInterpolate")) {
  //    return SymmetryInterpolate(
  //        cctkGH, N_dims, local_interp_handle, param_table_handle,
  //        coord_system_handle, N_interp_points, interp_coords_type_code,
  //        coords, N_input_arrays, input_array_variable_indices,
  //        N_output_arrays, output_array_type_codes, output_arrays);
  //  } else {
  return CarpetX_DriverInterpolate(
      cctkGH, N_dims, local_interp_handle, param_table_handle,
      coord_system_handle, N_interp_points, interp_coords_type_code, coords,
      N_input_arrays, input_array_variable_indices, N_output_arrays,
      output_array_type_codes, output_arrays);
  //  }
}

extern "C" CCTK_INT CarpetX_DriverInterpolate(
    CCTK_POINTER_TO_CONST const cctkGH, CCTK_INT const N_dims,
    CCTK_INT const local_interp_handle, CCTK_INT const param_table_handle,
    CCTK_INT const coord_system_handle, CCTK_INT const N_interp_points,
    CCTK_INT const interp_coords_type_code,
    CCTK_POINTER_TO_CONST const coords[], CCTK_INT const N_input_arrays,
    CCTK_INT const input_array_variable_indices[],
    CCTK_INT const N_output_arrays, CCTK_INT const output_array_type_codes[],
    CCTK_POINTER const output_arrays[]) {
  DECLARE_CCTK_PARAMETERS;

  // We do not support local interpolators yet
  const int carpetx_interp_handle = CCTK_InterpHandle("CarpetX");
  assert(carpetx_interp_handle >= 0);
  if (carpetx_interp_handle != local_interp_handle) {
    CCTK_VERROR("Incorrect local interpolator handle provided, only 'CarpetX' "
                "is allowed: %d != %d",
                local_interp_handle, carpetx_interp_handle);
  }

  // This verifies that the order in param_table_handle matches the order of the
  // runtime parameter from CarpetX.
  //
  // amr_interlude.md AMR-D9, fixed in step AMR-B7.  These two checks used to be
  // bare `assert`s.  `NDEBUG` is defined in NEITHER of the two configurations
  // this branch is built in (`[P199]`), so they abort the OPTIMIZED build as
  // well, with SIGABRT and no message -- and what they catch is a plain
  // parameter mistake, which a caller cannot diagnose from a bare SIGABRT.
  // `PunctureTracker::interp_order` defaults to 1 and `interpolation_order`
  // defaults to 1, so the two agree until a parameter file sets one of them,
  // and then every tracked run dies without saying why.
  //
  // The table's "order" is read here and nowhere else: the interpolation
  // dispatches on `interpolation_order` (the two `switch (interpolation_order)`
  // below), so a mismatch cannot be honoured, only reported.
  CCTK_INT order;
  int n_elems = Util_TableGetInt(param_table_handle, &order, "order");
  {
    const cFunctionData *const caller =
        CCTK_ScheduleQueryCurrentFunction(static_cast<const cGH *>(cctkGH));
    const char *const thorn = caller ? caller->thorn : "<not in a schedule bin>";
    const char *const routine =
        caller ? caller->routine : "<not in a schedule bin>";
    if (n_elems != 1)
      CCTK_VERROR("DriverInterpolate: the caller (%s::%s) did not put exactly "
                  "one integer \"order\" into its parameter table "
                  "(Util_TableGetInt returned %d).  CarpetX requires that key, "
                  "and requires it to equal CarpetX::interpolation_order, "
                  "which is %d.",
                  thorn, routine, n_elems, int(interpolation_order));
    if (order != interpolation_order)
      CCTK_VERROR("DriverInterpolate: interpolation order mismatch.  The "
                  "caller (%s::%s) asked for order=%d in its parameter table, "
                  "but CarpetX::interpolation_order is %d.  CarpetX "
                  "interpolates at CarpetX::interpolation_order and ignores "
                  "the parameter table's value, so it refuses here rather "
                  "than silently interpolating at an order nobody asked for.  "
                  "Set the calling thorn's own order parameter and "
                  "CarpetX::interpolation_order equal.  Note that "
                  "CarpetX::interpolation_order is STEERABLE=always while a "
                  "caller's order parameter usually is not, so steering it "
                  "mid-run breaks this equality too.",
                  thorn, routine, int(order), int(interpolation_order));
  }

  std::vector<CCTK_INT> varinds;
  varinds.resize(N_output_arrays);
  n_elems = Util_TableGetIntArray(param_table_handle, N_output_arrays,
                                  varinds.data(), "operand_indices");
  if (n_elems == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert(N_input_arrays == N_output_arrays);
    for (int i = 0; i < N_input_arrays; i++) {
      varinds.at(i) = input_array_variable_indices[i];
    }
  } else if (n_elems == N_output_arrays) {
    for (int i = 0; i < n_elems; i++) {
      varinds.at(i) = input_array_variable_indices[varinds.at(i)];
    }
  } else {
    CCTK_VERROR("TableGetIntArray failed with error code %d", n_elems);
  }

  std::vector<CCTK_INT> operations;
  operations.resize(N_output_arrays, 0);
  n_elems = Util_TableGetIntArray(param_table_handle, N_output_arrays,
                                  operations.data(), "operation_codes");
  if (n_elems == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    assert(N_input_arrays == N_output_arrays);
  } else if (n_elems != N_output_arrays) {
    CCTK_ERROR("TableGetIntArray failed.");
  }

  const CCTK_POINTER resultptrs = (CCTK_POINTER)output_arrays;
  const bool allow_boundaries = true;
  CarpetX_Interpolate(
      cctkGH, N_interp_points, static_cast<const CCTK_REAL *>(coords[0]),
      static_cast<const CCTK_REAL *>(coords[1]),
      static_cast<const CCTK_REAL *>(coords[2]), N_output_arrays,
      varinds.data(), operations.data(), allow_boundaries, resultptrs);

  return 0;
}

CarpetX::InterpolationSetup::InterpolationSetup(
    CCTK_ATTRIBUTE_UNUSED const cGH *restrict const cctkGH,
    const CCTK_INT npoints, const CCTK_REAL *restrict const globalsx,
    const CCTK_REAL *restrict const globalsy,
    const CCTK_REAL *restrict const globalsz,
    const bool require_level0_donors)
    : npoints(npoints), require_level0_donors(require_level0_donors) {
  DECLARE_CCTK_PARAMETERS;
  assert(in_global_mode(cctkGH));

  static const bool have_MultiPatch_GlobalToLocal2 =
      CCTK_IsFunctionAliased("MultiPatch_GlobalToLocal2");

  // Convert global to patch-local coordinates
  // TODO: Call this only if there is a non-trivial patch system
  // TODO: Don't call this for multipatch interpolation, precalculate this
  // instead
  std::vector<CCTK_INT> patches(npoints);
  std::vector<CCTK_REAL> localsx(npoints);
  std::vector<CCTK_REAL> localsy(npoints);
  std::vector<CCTK_REAL> localsz(npoints);

  if (have_MultiPatch_GlobalToLocal2) {
    MultiPatch_GlobalToLocal2(npoints, globalsx, globalsy, globalsz,
                              patches.data(), localsx.data(), localsy.data(),
                              localsz.data());
  } else {
    for (int n = 0; n < npoints; ++n) {
      patches.at(n) = 0;
    }

    for (int n = 0; n < npoints; ++n) {
      localsx.at(n) = globalsx[n];
    }

    for (int n = 0; n < npoints; ++n) {
      localsy.at(n) = globalsy[n];
    }

    for (int n = 0; n < npoints; ++n) {
      localsz.at(n) = globalsz[n];
    }
  }

  // Apply symmetries to coordinates
  assert(!reflection_x);
  assert(!reflection_y);
  assert(!reflection_upper_x);
  assert(!reflection_upper_y);
  assert(!reflection_upper_z);

  if (reflection_z) {
    symmetry_reflected_z.resize(npoints);
    assert(ghext->num_patches() == 1);
    constexpr int patch = 0;
    const amrex::Geometry &geom = ghext->patchdata.at(patch).amrcore->Geom(0);
    const CCTK_REAL *restrict const xmin = geom.ProbLo();
#pragma omp simd
    for (int n = 0; n < npoints; ++n) {
      const auto refl = localsz[n] < xmin[2];
      symmetry_reflected_z[n] = refl;
      if (refl) {
        localsz[n] = 2 * xmin[2] - localsz[n];
      }
    }
  }

  // Project particles into the domain for AMReX's distribution
  // AMReX silently drops particles that are outside the domain. We
  // can't have this. We thus push them back into the domain. Of
  // course, these modified coordinates are not useful for
  // interpolating, so we have both the `local` (true) and the `pos`
  // (AMReX) coordinates.
  std::vector<CCTK_REAL> posx(npoints);
  std::vector<CCTK_REAL> posy(npoints);
  std::vector<CCTK_REAL> posz(npoints);

#pragma omp simd
  for (int n = 0; n < npoints; ++n) {
    const int patch = patches.at(n);
    const amrex::Geometry &geom = ghext->patchdata.at(patch).amrcore->Geom(0);
    const CCTK_REAL *restrict const xmin = geom.ProbLo();
    const CCTK_REAL *restrict const xmax = geom.ProbHi();
    const CCTK_REAL *restrict const dx = geom.CellSize();
    using std::clamp;
    // Push the particle at least 1/2 grid spacing into the domain
    // TODO: push by less because 1/2 is too much if there are many AMR levels
    posx[n] = clamp(localsx[n], xmin[0] + dx[0] / 2, xmax[0] - dx[0] / 2);
    posy[n] = clamp(localsy[n], xmin[1] + dx[1] / 2, xmax[1] - dx[1] / 2);
    posz[n] = clamp(localsz[n], xmin[2] + dx[2] / 2, xmax[2] - dx[2] / 2);
  }

  // Create particle containers
  std::vector<PinnedParticleTile> pinned_particle_tiles(ghext->num_patches());
  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    PinnedParticleTile &pinned_particle_tile = pinned_particle_tiles.at(patch);
    // here the two slots represents components in the structure-of-arrays (SoA)
    // layout
    pinned_particle_tile.define(0, 0);
  }

  // Set particle positions
  // TODO: parallelize this loop
  const int proc = amrex::ParallelDescriptor::MyProc();
  for (int n = 0; n < npoints; ++n) {
    const int patch = patches.at(n);
    amrex::Particle<3, 2> p;
    p.id() = Particle::NextID();
    p.cpu() = proc;
    p.pos(0) = posx[n]; // AMReX distribution position
    p.pos(1) = posy[n];
    p.pos(2) = posz[n];
    p.rdata(0) = localsx[n]; // actual particle coordinate
    p.rdata(1) = localsy[n];
    p.rdata(2) = localsz[n];
    p.idata(0) = proc; // source process
    p.idata(1) = n;    // source index
    pinned_particle_tiles.at(patch).push_back(p);
  }

  containers.resize(ghext->num_patches());
  int owned_patches = 0;
  for (int patch = 0; patch < ghext->num_patches(); ++patch) {
    const PinnedParticleTile &pinned_particle_tile =
        pinned_particle_tiles.at(patch);

    const auto &restrict patchdata = ghext->patchdata.at(patch);
    containers.at(patch) = patchdata.amrcore.get();
    const int level = 0;
    const auto &restrict leveldata = patchdata.leveldata.at(level);
    const amrex::MFIter mfi(*leveldata.fab);
    // The mfi can be invalid if the number of processes does not evenly divide
    // the number of blocks
    if (!mfi.isValid()) {
      // This process owns no level-0 box of this patch, so it has nowhere to
      // put the query points it collected for it.  Dropping them is harmless
      // only if there are none: a dropped point is never inserted into any
      // container, so nobody interpolates it, nobody answers it, and this
      // process's receive buffer comes back short.  That is caught -- but it
      // is caught 500 lines below, by a size comparison whose whole diagnosis
      // is the words "Internal error", on a rank that cannot say which patch
      // it was or how many points it lost.
      //
      // So refuse here, where the three numbers still exist.  The predicate is
      // exact rather than conservative: it cannot fire on any configuration
      // that works today, because a process with no points for this patch
      // loses nothing by skipping it.
      //
      // There is deliberately NO parameter to downgrade this to a warning.  A
      // hatch would be a lie: the run does not survive the drop either way,
      // and a warning would only move the abort back to the site that cannot
      // name anything.
      //
      // WHEN THIS FIRES.  `amr.refine_grid_layout` (CarpetX::refine_grid_layout,
      // default yes) makes AMReX chop each patch's level-0 box array until
      // every process has a box, so a freshly decomposed grid satisfies this
      // by construction.  A RECOVERED grid does not go through that path at
      // all: `RecoverGridStructure` restores the box array out of the
      // checkpoint, so a checkpoint written at N processes carries a box array
      // built for N and nothing rebuilds it for M.
      const int np = int(pinned_particle_tile.numParticles());
      if (np > 0) {
        // CCTK_VERROR discards buffered stdout, and the box census a reader
        // wants next to this message is on stdout.
        std::fflush(nullptr);
        CCTK_VERROR(
            "This process (%d of %d) owns no level-0 box of patch %d, but it "
            "holds %d of its own %lld interpolation query point(s) for that "
            "patch. Those points would be dropped here and answered by "
            "nobody. Patch %d's level-0 box array has %d box(es) for %d "
            "process(es). CarpetX::refine_grid_layout = yes makes AMReX chop "
            "every patch's level-0 box array until each process holds a box, "
            "and it is what a freshly decomposed grid relies on; a RECOVERED "
            "grid does not, because the box array is restored from the "
            "checkpoint and is therefore the array that was built for the "
            "process count the checkpoint was written at. Recover at that "
            "process count, or start from initial data.",
            proc, amrex::ParallelDescriptor::NProcs(), patch, np,
            (long long)npoints, patch, leveldata.fab->size(),
            amrex::ParallelDescriptor::NProcs());
      }
      continue;
    }
    ++owned_patches;

    ParticleTile &particle_tile = containers.at(patch).GetParticles(
        level)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

    const auto old_np = particle_tile.numParticles();
    const auto new_np = old_np + pinned_particle_tile.numParticles();
    particle_tile.resize(new_np);
    amrex::copyParticles(particle_tile, pinned_particle_tile, 0, old_np,
                         pinned_particle_tile.numParticles());
  }

  // Say so when the refusal above did NOT fire, because a guard that speaks
  // only when it refuses turns every quiet run into an unmeasured one: there
  // is then no way to read "every process owns a box of every patch" apart
  // from "this binary has no such check".
  //
  // TWO INDEPENDENT ONE-SHOTS, and that is not tidiness.  With a single flag
  // the announcement is spent on whichever call comes first, and on a patch
  // system the first call is routinely the EMPTY one -- so the run carries
  // "nothing to check" from a process that holds two hundred thousand query
  // points a few lines later, and never carries the informative line at all.
  // That is exactly what happened to the C-AMR2 clearance line, where it went
  // unnoticed long enough to be written into three reports before anyone
  // measured it.  Two flags cost nothing and cannot be spent on each other.
  //
  // WHAT THIS CHANNEL IS AND IS NOT.  It reaches a log only from the ROOT
  // process: the flesh `freopen`s stdout to the null device on every non-root
  // process unless the run is given `-r`
  // (`flesh/src/main/CommandLine.c:783`), so a `CCTK_VINFO` from process 1 is
  // discarded before `CCTK_VInfo` is even reached.  The per-process channel is
  // the refusal above, which is a `CCTK_VERROR` and therefore goes to stderr,
  // which the flesh does not redirect.  So: this line witnesses that the check
  // RUNS; what witnesses that it HELD on every process is that the run
  // survived, because the refusal has no warn-only mode.
  {
    static bool announced_holds = false;
    static bool announced_empty = false;
    if (npoints == 0) {
      if (!announced_empty) {
        announced_empty = true;
        CCTK_VINFO("Level-0 patch ownership has nothing to check on this "
                   "process (%d of %d): it holds no interpolation query point "
                   "at all, over %d patch(es)",
                   proc, amrex::ParallelDescriptor::NProcs(),
                   ghext->num_patches());
      }
    } else if (!announced_holds) {
      announced_holds = true;
      CCTK_VINFO("Level-0 patch ownership holds: this process (%d of %d) owns "
                 "a level-0 box of %d of the %d patch(es) and none of its "
                 "%lld query point(s) is dropped",
                 proc, amrex::ParallelDescriptor::NProcs(), owned_patches,
                 ghext->num_patches(), (long long)npoints);
    }
  }

  // Send particles to interpolation points
  for (auto &container : containers) {
#ifdef CCTK_DEBUG
    const int patch = int(&container - containers.data());

    std::size_t old_nparticles = 0;
    std::set<int> oldids;
    {
      const auto &levels = container.GetParticles();
      for (const auto &level : levels) {
        const int lev = int(&level - levels.data());
        for (amrex::ParConstIter<3, 2> pti(container, lev); pti.isValid();
             ++pti) {
          const auto &particles = pti.GetArrayOfStructs();
          const int component = MFPointer(pti).index();
          for (const auto &particle : particles)
            oldids.insert(particle.id());
          old_nparticles += particles.size();
        }
      }
      const MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
      MPI_Allreduce(MPI_IN_PLACE, &old_nparticles, 1,
                    mpi_datatype<std::size_t>::value, MPI_SUM, comm);
    }
#endif

    container.Redistribute();

#ifdef CCTK_DEBUG
    std::size_t new_nparticles = 0;
    std::set<int> newids;
    {
      const auto &levels = container.GetParticles();
      for (const auto &level : levels) {
        const int lev = int(&level - levels.data());
        for (amrex::ParConstIter<3, 2> pti(container, lev); pti.isValid();
             ++pti) {
          const int component = MFPointer(pti).index();
          const auto &particles = pti.GetArrayOfStructs();
          for (const auto &particle : particles)
            newids.insert(particle.id());
          new_nparticles += particles.size();
        }
      }
      const MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
      MPI_Allreduce(MPI_IN_PLACE, &new_nparticles, 1,
                    mpi_datatype<std::size_t>::value, MPI_SUM, comm);
    }
    if (new_nparticles != old_nparticles) {
      for (const auto oldid : oldids)
        if (!newids.count(oldid))
          CCTK_VWARN(CCTK_WARN_ALERT, "old id %d not present in new ids",
                     oldid);
      for (const auto newid : newids)
        if (!oldids.count(newid))
          CCTK_VWARN(CCTK_WARN_ALERT, "new id %d not present in old ids",
                     newid);
      CCTK_VERROR(
          "We lost interpolation points on patch %d. Before redistributing: "
          "%zu particles, after redistributing: %zu particles",
          patch, old_nparticles, new_nparticles);
    }
#endif
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// C-AMR2 -- NO INTERPATCH GHOST POINT MAY BE ANSWERED FROM A `level > 0` BOX.
//
//   C-AMR2.  No patch's interpatch ghost point may be answered from a
//   `level > 0` box.
//
// This is the second half of the contract that replaces the blanket refusal of
// mesh refinement on a multi-patch grid.  The first half (C-AMR, in
// `fillpatch.cxx`) is about the PROLONGATION: it forbids reading coarse data
// from outside the owning patch.  This half is about the DONOR of an
// interpatch seam value, and the two are independent -- a configuration can
// satisfy C-AMR with five coarse cells to spare and violate this one.
//
// THE MECHANISM.  `InterpolationSetup`'s constructor inserts every query
// particle at level 0 and then calls `Redistribute()`, whose `Where()` loop
// runs FINEST LEVEL FIRST.  A particle that lands inside a refined box is
// therefore assigned to the refined level, and `Interpolate` answers it from
// prolongated fine data instead of evolved coarse data.  Nothing is NaN and
// nothing aborts: the values are finite, they are simply drawn from a
// different level than the seam fill was designed around, and the switch
// happens mid-run, at a moment set by wherever the refinement boxes have
// travelled to.  That is what makes it worth a refusal rather than a warning:
// a seam whose donor changes character during a run is not a discretisation
// the rest of the scheme can be reasoned about.
//
// WHY THE CALLER DECLARES THIS AND THE DRIVER DOES NOT INFER IT.  This class
// has two kinds of caller and only one of them is C-AMR2's subject:
//
//   - an interpatch seam fill, whose query points ARE the ghost points of a
//     patch boundary.  For those, a `level > 0` answer is the violation.
//   - `CarpetX_Interpolate` and everything that funnels into it -- puncture
//     tracking, wave extraction, any thorn calling `CCTK_InterpGridArrays`.
//     For those a `level > 0` answer is CORRECT and usually the point: a
//     tracked puncture sits inside its own refinement box by construction.
//
// The driver cannot tell the two apart from the coordinates, so it does not
// try.  `require_level0_donors` is passed by the caller, defaults to `false`,
// and only a caller that knows its points are interpatch ghosts sets it.  A
// default of `true` would have been fail-safe in the abstract and wrong here:
// it would silently impose a contract on callers that never asked for one, and
// refusing a configuration that works is the one thing this guard may not do.
//
// WHY IT RUNS AS A PRE-PASS AND NOT INSIDE THE INTERPOLATION LOOP.  Two
// reasons, and the second is the load-bearing one.  (a) A refusal should
// happen before the work, not after it.  (b) The interpolation kernel below
// carries a bare `assert(all(i >= 0 && i + order < grid.lsh))` inside an
// `omp parallel for`, live in the optimized build, and on a violating
// multi-patch configuration that assert has been observed to fire
// NON-DETERMINISTICALLY -- 13 of 24 attempts on one geometry, from inside the
// cache rebuild that the regrid repair triggers during initialisation.  A
// guard that ran after the loop would be shadowed by that abort in half the
// runs and would report a different thing in the other half.  This pass reads
// only particle counts, so it is deterministic and it always comes first.
//
// COST.  One `ParConstIter` sweep over levels `>= 1` only, per `Interpolate`
// call, and only when a caller has asked for the check.  At
// `max_num_levels = 1` there is no such level and the loop body never runs.
//
// THE COUNTS ARE RANK-LOCAL.  Which rank answers a given query point is
// decided by `Redistribute`, so a violation may be visible to one process and
// not to another.  Every process evaluates the predicate and any one of them
// stops the job; the message says whose count it is.
//
////////////////////////////////////////////////////////////////////////////////

void CarpetX::InterpolationSetup::RefuseAboveLevel0Donors(
    const CCTK_INT nvars, const CCTK_INT *restrict const varinds) const {
  DECLARE_CCTK_PARAMETERS;

  if (!require_level0_donors)
    return;

  long long answered_above0 = 0;
  int maxlevel = 0, nrefined_patches = 0;
  int first_patch = -1, first_level = -1;
  CCTK_REAL first_x = 0, first_y = 0, first_z = 0;

  for (const auto &patchdata : ghext->patchdata) {
    const int patch = patchdata.patch;
    const int nlevels = int(patchdata.leveldata.size());
    if (nlevels > 1)
      ++nrefined_patches;
    for (int level = 1; level < nlevels; ++level) {
      for (amrex::ParConstIter<3, 2> pti(containers.at(patch), level);
           pti.isValid(); ++pti) {
        const int np = pti.numParticles();
        if (np <= 0)
          continue;
        answered_above0 += np;
        if (level > maxlevel)
          maxlevel = level;
        if (first_patch < 0) {
          // `rdata(0..2)` is the TRUE patch-local query coordinate;
          // `pos(0..2)` is the clamped position `Redistribute` located the
          // particle with, which is not the coordinate the answer is for.
          const auto &particles = pti.GetArrayOfStructs();
          first_patch = patch;
          first_level = level;
          first_x = particles[0].rdata(0);
          first_y = particles[0].rdata(1);
          first_z = particles[0].rdata(2);
        }
      }
    }
  }

  // Nothing above level 0.  Report it once per process, but only once a
  // refined level exists at all -- otherwise the line would say nothing,
  // every run would carry it, and it would stop being evidence.
  //
  // A ZERO POINT COUNT GETS ITS OWN SENTENCE, and that is not fussiness: a
  // `patch_system = "Cartesian"` run reaches here with `npoints == 0`, and a
  // line reading "C-AMR2 holds" over an empty set is a check that passed on
  // nothing.  Say which of the two happened.
  static bool announced_holds = false;
  if (answered_above0 == 0) {
    if (nrefined_patches > 0 && !announced_holds) {
      announced_holds = true;
      if (npoints == 0)
        CCTK_VINFO("C-AMR2 has nothing to check on this process: it holds no "
                   "interpatch query point at all, while %d patch(es) carry a "
                   "refined level. This is what a single-patch grid looks "
                   "like here",
                   nrefined_patches);
      else
        CCTK_VINFO("C-AMR2 holds: all %lld of this process's interpatch query "
                   "points are answered from level 0, with %d patch(es) "
                   "carrying a refined level",
                   (long long)npoints, nrefined_patches);
    }
    return;
  }

  const char *const var0 = nvars > 0 ? CCTK_FullVarName(int(varinds[0])) : NULL;
  const bool warn_only = CCTK_EQUALS(multipatch_amr_contract, "warn");

  // `CCTK_VERROR` discards buffered stdout, and the level and box census a
  // reader wants next to this message is on stdout.
  std::fflush(nullptr);

  static bool announced_violated = false;
  if (warn_only && announced_violated)
    return;
  announced_violated = true;

  char msg[2400];
  std::snprintf(
      msg, sizeof msg,
      "C-AMR2 is VIOLATED. %lld of this process's %lld interpatch query "
      "points (filling %s%s) are answered from a level > 0 box; the deepest "
      "is level %d, and the first of them is ANSWERED ON patch %d at level "
      "%d, at that patch's local coordinate (%.17g, %.17g, %.17g) -- which is "
      "the patch that OWNS the point and carries the refined level, not the "
      "patch whose ghost zone asked for it. Those values are finite and "
      "nothing "
      "will abort: they are prolongated fine data where the seam fill expects "
      "evolved coarse data, and which of the two a given seam point gets "
      "changes during the run as the refinement boxes move. A refined region "
      "has reached the zone from which another patch's interpatch ghost "
      "points are drawn -- that zone extends inward of the patch interface by "
      "the patch overlap plus the ghost width, so it can be met without any "
      "refinement box being near the interface itself. Either keep the "
      "refined region out of it (BoxInBox::radius_*, BoxInBox::position_*), "
      "or move the zone by raising the patch system's resolution or its inner "
      "boundary -- a resolution change must keep every patch's cell count a "
      "multiple of CarpetX::blocking_factor_{x,y,z} in every direction, so "
      "the reachable values are quantised and the first arithmetically "
      "sufficient one is usually not reachable. "
      "CarpetX::multipatch_amr_contract = \"warn\" downgrades this to a "
      "warning, for diagnosis only.",
      answered_above0, (long long)npoints, var0 ? var0 : "unknown variable",
      nvars > 1 ? " and others" : "", maxlevel, first_patch, first_level,
      double(first_x), double(first_y), double(first_z));

  if (warn_only)
    CCTK_VWARN(CCTK_WARN_ALERT, "%s", msg);
  else
    CCTK_VERROR("%s", msg);
}

void CarpetX::InterpolationSetup::Interpolate(
    CCTK_ATTRIBUTE_UNUSED const cGH *restrict const cctkGH,
    const CCTK_INT nvars, const CCTK_INT *restrict const varinds,
    const CCTK_INT *restrict const operations,
    const std::vector<Arith::vect<Arith::vect<bool, 3>, 2> >
        &allowed_boundaries, //  [patch][face][direction]
    const CCTK_POINTER resultptrs_) const {
  DECLARE_CCTK_PARAMETERS;

  // C-AMR2, before any work: a caller whose query points are interpatch ghosts
  // gets them answered from level 0 or not at all.
  RefuseAboveLevel0Donors(nvars, varinds);

  // Define result variables
  const int nprocs = amrex::ParallelDescriptor::NProcs();
  std::vector<std::vector<CCTK_REAL> > results(nprocs); // [nprocs]

  // Interpolate
  constexpr int tl = 0;
  struct givi_t {
    int gi, vi;
  };
  std::vector<givi_t> givis(nvars);
  for (int v = 0; v < nvars; ++v) {
    int gi = CCTK_GroupIndexFromVarI(varinds[v]);
    assert(gi >= 0);
    assert(gi < CCTK_NumGroups());
    int vi = varinds[v] - CCTK_FirstVarIndexI(gi);
    assert(vi >= 0);
    assert(vi < CCTK_NumVarsInGroupI(gi));
    givis.at(v) = {gi, vi};
  }

  for (const auto &patchdata : ghext->patchdata) {
    const int patch = patchdata.patch;
    for (const auto &leveldata : patchdata.leveldata) {
      const int level = leveldata.level;

      // TODO: use OpenMP
      for (amrex::ParConstIter<3, 2> pti(containers.at(patch), level);
           pti.isValid(); ++pti) {
        const MFPointer mfp(pti);
        const GridDesc grid(leveldata, mfp);
        // const int component = mfp.index();

        // Derive per-box stencil-anchor permissions. bbox[f][d] is true when
        // the box face touches the patch's AMReX domain boundary (both physical
        // outer boundaries and inter-patch boundaries). Interior intra-patch
        // faces (bbox=false) allow anchoring unconditionally, because AMReX
        // guarantees their ghost zones are filled by the time this runs. The
        // force_conservative_intrapatch escape hatch that used to relax that
        // for a caller running before AMReX's fill-patch pass went with the
        // caller itself in BUGFIX_TODO.md step B3; see interp.hxx.
        vect<vect<bool, dim>, 2> patch_allowed_boundaries;
        for (int f = 0; f < 2; ++f)
          for (int d = 0; d < dim; ++d)
            patch_allowed_boundaries[f][d] =
                grid.bbox[f][d] ? allowed_boundaries.at(patch)[f][d] : true;

        const int np = pti.numParticles();
        const auto &particles = pti.GetArrayOfStructs();

        std::vector<std::vector<CCTK_REAL> > varresults(nvars);

        // TODO: Don't re-calculate interpolation coefficients for each
        // variable
        for (int v = 0; v < nvars; ++v) {
          const int gi = givis.at(v).gi;
          const int vi = givis.at(v).vi;
          const auto &restrict groupdata = *leveldata.groupdata.at(gi);
          const int centering = groupdata.indextype[0] * 0b100 +
                                groupdata.indextype[1] * 0b010 +
                                groupdata.indextype[2] * 0b001;
          assert(all(groupdata.nghostzones == grid.nghostzones));
          const amrex::Array4<const CCTK_REAL> &vars =
              groupdata.mfab.at(tl)->array(pti);
          vect<int, dim> derivs;
          int op = operations[v];
          while (op > 0) {
            const int dir = op % 10 - 1;
            if (dir >= 0) {
              assert(dir >= 0 && dir < dim);
              ++derivs[dir];
            }
            op /= 10;
          }
          auto &varresult = varresults.at(v);
          varresult.resize(np);

          switch (centering) {
          case 0b000: {
            // Vertex centering

            switch (interpolation_order) {
            case 0: {
              const interpolator<CCTK_REAL, 0, 0b000> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 1: {
              const interpolator<CCTK_REAL, 1, 0b000> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 2: {
              const interpolator<CCTK_REAL, 2, 0b000> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 3: {
              const interpolator<CCTK_REAL, 3, 0b000> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 4: {
              const interpolator<CCTK_REAL, 4, 0b000> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            default:
              CCTK_VERROR("Interpolation order %d for centering [%d,%d,%d] not "
                          "yet supported",
                          int(interpolation_order), groupdata.indextype[0],
                          groupdata.indextype[1], groupdata.indextype[2]);
            } // switch interpolation_order
            break;
          } // case 0b000

          case 0b111: {
            // Cell centering

            switch (interpolation_order) {
            case 0: {
              const interpolator<CCTK_REAL, 0, 0b111> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 1: {
              const interpolator<CCTK_REAL, 1, 0b111> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 2: {
              const interpolator<CCTK_REAL, 2, 0b111> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 3: {
              const interpolator<CCTK_REAL, 3, 0b111> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            case 4: {
              const interpolator<CCTK_REAL, 4, 0b111> interp{
                  grid,  gi,   vi,     patch,
                  level, vars, derivs, patch_allowed_boundaries};
              interp.interpolate3d(particles, varresult);
              break;
            }
            default:
              CCTK_VERROR("Interpolation order %d for centering [%d,%d,%d] not "
                          "yet supported",
                          int(interpolation_order), groupdata.indextype[0],
                          groupdata.indextype[1], groupdata.indextype[2]);
            } // switch interpolation_order
            break;
          } // case 0b111

          default:
            CCTK_VERROR("Centering [%d,%d,%d] not yet supported",
                        groupdata.indextype[0], groupdata.indextype[1],
                        groupdata.indextype[2]);
          } // switch centering

        } // for var

        for (int n = 0; n < np; ++n) {
          const int proc = particles[n].idata(0);
          const int id = particles[n].idata(1);
          auto &result = results.at(proc);
          result.push_back(id);
          for (int v = 0; v < nvars; ++v)
            result.push_back(varresults.at(v).at(n));
        }
      }
    }
  }

  // Collect particles back
  const MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
  const MPI_Datatype datatype = mpi_datatype<CCTK_REAL>::value;

  std::vector<int> sendcounts(nprocs);
  std::vector<int> senddispls(nprocs);
  int total_sendcount = 0;
  for (int p = 0; p < nprocs; ++p) {
    const auto &result = results.at(p);
    sendcounts.at(p) = result.size();
    senddispls.at(p) = total_sendcount;
    total_sendcount += sendcounts.at(p);
  }
  std::vector<int> recvcounts(nprocs);
  MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT,
               comm);
  std::vector<int> recvdispls(nprocs);
  int total_recvcount = 0;
  for (int p = 0; p < nprocs; ++p) {
    recvdispls.at(p) = total_recvcount;
    total_recvcount += recvcounts.at(p);
  }

  std::vector<CCTK_REAL> sendbuf(total_sendcount);
  for (int p = 0; p < nprocs; ++p) {
    // TODO: Don't copy, store data here right away
    assert(p >= 0);
    assert(p < int(results.size()));
    const auto &result = results.at(p);
    assert(sendcounts.at(p) == int(result.size()));
    assert(p >= 0);
    assert(p < int(senddispls.size()));
    assert(senddispls.at(p) >= 0);
    assert(senddispls.at(p) + sendcounts.at(p) <= int(sendbuf.size()));
    std::copy(result.begin(), result.end(), sendbuf.data() + senddispls.at(p));
  }
  std::vector<CCTK_REAL> recvbuf(total_recvcount);
  MPI_Alltoallv(sendbuf.data(), sendcounts.data(), senddispls.data(), datatype,
                recvbuf.data(), recvcounts.data(), recvdispls.data(), datatype,
                comm);
#ifdef CCTK_DEBUG
  // Check consistency of received ids
  std::vector<bool> idxs(npoints, false);
  for (int n = 0; n < npoints; ++n) {
    const int offset = (nvars + 1) * n;
    const int idx = int(recvbuf.at(offset));
    assert(!idxs.at(idx));
    idxs.at(idx) = true;
  }
  for (int n = 0; n < npoints; ++n)
    assert(idxs.at(n));
#endif

  // Set result
  CCTK_REAL *const restrict *const restrict resultptrs =
      static_cast<CCTK_REAL *const *>(resultptrs_);
  if (int(recvbuf.size()) != (nvars + 1) * npoints)
    CCTK_ERROR("Internal error");
  for (int n = 0; n < npoints; ++n) {
    const int offset = (nvars + 1) * n;
    const int idx = int(recvbuf.at(offset));
    for (int v = 0; v < nvars; ++v)
      resultptrs[v][idx] = recvbuf.at(offset + 1 + v);
  }

  // Apply symmetries to interpolated values
  assert(!reflection_x);
  assert(!reflection_y);
  assert(!reflection_upper_x);
  assert(!reflection_upper_y);
  assert(!reflection_upper_z);

  if (reflection_z) {
    // The code below is only valid for Psi4
    assert(nvars == 2);
    assert(varinds[0] == CCTK_VarIndex("Weyl::Psi4re"));
    assert(varinds[1] == CCTK_VarIndex("Weyl::Psi4im"));
    // l^a = et^a + er^a
    // n^a = et^a - er^a
    // m^a = etheta^a + i ephi^a
    // Psi4 = C_abcd m-bar^b n^b m-bar^c n^d
    for (int n = 0; n < npoints; ++n) {
      if (symmetry_reflected_z[n]) {
        resultptrs[0][n] = -resultptrs[0][n];
        resultptrs[1][n] = +resultptrs[1][n];
      }
    }
  }
}

extern "C" void CarpetX_Interpolate(const CCTK_POINTER_TO_CONST cctkGH_,
                                    const CCTK_INT npoints,
                                    const CCTK_REAL *restrict const globalsx,
                                    const CCTK_REAL *restrict const globalsy,
                                    const CCTK_REAL *restrict const globalsz,
                                    const CCTK_INT nvars,
                                    const CCTK_INT *restrict const varinds,
                                    const CCTK_INT *restrict const operations,
                                    const CCTK_INT allow_boundaries,
                                    const CCTK_POINTER resultptrs_) {
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);

  const InterpolationSetup setup(cctkGH, npoints, globalsx, globalsy, globalsz);

  // Replicate original behaviour:
  // allow_boundaries=true  -> allow stencil on all faces
  // allow_boundaries=false -> forbid stencil on all bbox faces
  const int npatches = ghext->num_patches();
  const vect<vect<bool, dim>, 2> uniform{
      {{bool(allow_boundaries), bool(allow_boundaries), bool(allow_boundaries)},
       {bool(allow_boundaries), bool(allow_boundaries),
        bool(allow_boundaries)}}};
  const std::vector<vect<vect<bool, dim>, 2> > allowed_boundaries(npatches,
                                                                  uniform);
  setup.Interpolate(cctkGH, nvars, varinds, operations, allowed_boundaries,
                    resultptrs_);
}

} // namespace CarpetX
