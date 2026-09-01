/*
 * BUGFIX_TODO.md step D1 (C12) -- the failing-before test for CarpetX's outer
 * boundary kernel.
 *
 * WHAT IT MEASURES.  Fill an affine field
 *
 *     f(x,y,z) = c0 + cx*x + cy*y + cz*z
 *
 * in the grid INTERIOR, sync, and compare every outer-boundary point against
 * the analytic value.  `linear extrapolation` reproduces an affine field
 * exactly by construction (`f(0) = 2 f(h) - f(2h)`), and `neumann` and
 * `reflection` reproduce it exactly in any direction in which it is constant.
 * `TestBoundaries_ExtrapolationParamCheck` refuses every configuration in
 * which that is not true of ALL SIX faces, so on any rig that starts, the
 * correct answer is "the analytic field, everywhere", and the tolerance can be
 * -- and by default is -- exactly ZERO.
 *
 * WHY THE OFFSET TABLE.  The pre-D1 kernel took one finite difference along
 * the diagonal and rescaled it by `sqrt(sum(pow2(dst-src)) / sum(pow2(delta)))`
 * in INTEGER arithmetic.  That is exact when every boundary offset is equal
 * (`sqrt(2/2) = 1`, `sqrt(8/2) = 2`) and wrong otherwise, so the error is a
 * function of the per-direction offset triple and of nothing else.  Reporting
 * max|error| alone would hide that.  The table keyed on
 * `(|off_x|, |off_y|, |off_z|)` is what makes the signature visible, and it is
 * `multipatch_case.md` M5(b)'s falsifier: if the equal-offset classes are NOT
 * exact before the fix, the reading of the kernel that D1 is built on is wrong.
 *
 * SINGLE RANK.  The table is accumulated in host memory and is not reduced
 * across ranks, so the routine refuses `nProcs > 1` rather than reporting a
 * rank-local table as if it were global.  The defect is not rank dependent;
 * there is nothing to gain from the extra machinery.
 */

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loop.hxx"

#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>

namespace TestBoundaries {

namespace {

// One accumulator bucket per distinct offset triple.
struct class_stats {
  long long npoints = 0;
  CCTK_REAL max_abs_error = 0;
  // The point at which `max_abs_error` was seen, for the report.
  std::array<int, 3> worst_I = {0, 0, 0};
  std::array<CCTK_REAL, 3> worst_X = {0, 0, 0};
  CCTK_REAL worst_got = 0, worst_want = 0;
};

// key = (|off_x|, |off_y|, |off_z|), zero in a direction that is not a
// boundary direction for that point.  `std::map` so the report comes out
// sorted, and because the number of keys is bounded by
// `(2*ghost_size+1)^3` and checked against `extrapolation_max_classes`.
std::map<std::array<int, 3>, class_stats> classes;

// Total points classified.  A pass on zero of them is a green light for
// nothing (BUGFIX_TODO.md Traps), so the report checks it.
long long npoints_seen = 0;

/*
 * CarpetX's boundary parameters are PRIVATE, so SHARES cannot reach them.
 * `CCTK_ParameterGet` can, and `CarpetX/TestProlongate` already takes exactly
 * this route for `prolongation_order` / `prolongation_type`, with the same
 * comment about the parameter being private on purpose.  Reading the driver's
 * real values is the point: a copy of the configuration kept in this thorn's
 * own parameters could disagree with what the driver is doing, and the guard
 * below would then be validating a rig that is not running.
 */
const char *carpetx_keyword(const char *const name) {
  int type = -1;
  const void *const p = CCTK_ParameterGet(name, "CarpetX", &type);
  if (!p)
    CCTK_VERROR("TestBoundaries: CarpetX has no parameter \"%s\"", name);
  if (type != PARAMETER_KEYWORD)
    CCTK_VERROR("TestBoundaries: CarpetX::%s is not a KEYWORD (type %d)", name,
                type);
  return *static_cast<const char *const *>(p);
}

CCTK_INT carpetx_int(const char *const name) {
  int type = -1;
  const void *const p = CCTK_ParameterGet(name, "CarpetX", &type);
  if (!p)
    CCTK_VERROR("TestBoundaries: CarpetX has no parameter \"%s\"", name);
  if (type != PARAMETER_INT && type != PARAMETER_BOOLEAN)
    CCTK_VERROR("TestBoundaries: CarpetX::%s is neither INT nor BOOLEAN "
                "(type %d)",
                name, type);
  return *static_cast<const CCTK_INT *>(p);
}

struct bc_config {
  // [face][dir], face 0 = lower
  const char *boundary[2][3];
  int reflection[2][3];
  int periodic[3];
  CCTK_INT ghost_size, ghost_size_d[3];
};

bc_config get_bc_config() {
  bc_config c;
  c.boundary[0][0] = carpetx_keyword("boundary_x");
  c.boundary[0][1] = carpetx_keyword("boundary_y");
  c.boundary[0][2] = carpetx_keyword("boundary_z");
  c.boundary[1][0] = carpetx_keyword("boundary_upper_x");
  c.boundary[1][1] = carpetx_keyword("boundary_upper_y");
  c.boundary[1][2] = carpetx_keyword("boundary_upper_z");
  c.reflection[0][0] = carpetx_int("reflection_x");
  c.reflection[0][1] = carpetx_int("reflection_y");
  c.reflection[0][2] = carpetx_int("reflection_z");
  c.reflection[1][0] = carpetx_int("reflection_upper_x");
  c.reflection[1][1] = carpetx_int("reflection_upper_y");
  c.reflection[1][2] = carpetx_int("reflection_upper_z");
  const CCTK_INT periodic = carpetx_int("periodic");
  c.periodic[0] = periodic && carpetx_int("periodic_x");
  c.periodic[1] = periodic && carpetx_int("periodic_y");
  c.periodic[2] = periodic && carpetx_int("periodic_z");
  c.ghost_size = carpetx_int("ghost_size");
  c.ghost_size_d[0] = carpetx_int("ghost_size_x");
  c.ghost_size_d[1] = carpetx_int("ghost_size_y");
  c.ghost_size_d[2] = carpetx_int("ghost_size_z");
  return c;
}

/*
 * Is the affine field with coefficient `coeff` in direction `d` reproduced
 * EXACTLY by the condition `keyword` on that face?  Returns nullptr when it
 * is, and otherwise the reason, for the refusal message.
 */
const char *why_not_exact(const char *keyword, const bool is_reflection,
                          const bool is_periodic, const CCTK_REAL coeff) {
  if (is_periodic)
    return "periodic: an affine field is not periodic, so the wrapped value is "
           "not the analytic one (and this test would be measuring the domain "
           "size, not the boundary kernel)";
  if (is_reflection)
    return coeff == 0
               ? nullptr
               : "reflection: exact only for a field that is EVEN about the "
                 "mirror plane. Set this direction's affine coefficient to 0";
  if (CCTK_Equals(keyword, "linear extrapolation"))
    return nullptr;
  if (CCTK_Equals(keyword, "neumann"))
    return coeff == 0 ? nullptr
                      : "neumann: copies the nearest interior value, which is "
                        "the analytic one only for a field that is CONSTANT in "
                        "this direction. Set this direction's affine "
                        "coefficient to 0";
  if (CCTK_Equals(keyword, "none"))
    return "none: nothing writes these boundary points, so with an "
           "interior-only fill they stay poison and any neighbouring edge or "
           "corner reads them";
  if (CCTK_Equals(keyword, "dirichlet"))
    return "dirichlet: writes a constant, which is not the affine field";
  if (CCTK_Equals(keyword, "robin"))
    return "robin: rescales by a ratio of radii, which does not reproduce an "
           "affine field";
  return "unrecognised boundary keyword";
}

} // namespace

extern "C" void TestBoundaries_ExtrapolationParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int nprocs = CCTK_nProcs(cctkGH);
  if (nprocs != 1)
    CCTK_VERROR("TestBoundaries: test_extrapolation = yes requires a single "
                "MPI rank (running on %d). The per-offset-class table is "
                "accumulated in host memory and is not reduced across ranks, "
                "and the defect it measures is not rank dependent. Rerun with "
                "one rank.",
                nprocs);

  const char *const dirnames[3] = {"x", "y", "z"};
  const CCTK_REAL coeff[3] = {affine_cx, affine_cy, affine_cz};
  const bc_config c = get_bc_config();

  for (int d = 0; d < 3; ++d) {
    for (int f = 0; f < 2; ++f) {
      const char *const keyword = c.boundary[f][d];
      const bool refl = bool(c.reflection[f][d]);
      const char *const why =
          why_not_exact(keyword, refl, bool(c.periodic[d]), coeff[d]);
      if (why)
        CCTK_VERROR(
            "TestBoundaries: the %s %s face is configured as \"%s\"%s with "
            "affine_c%s = %.17g, and that combination does not reproduce the "
            "affine test field exactly -- %s. This test compares every "
            "boundary point against the analytic field, so a face that cannot "
            "reproduce it would report an error that says nothing about the "
            "extrapolation kernel. Refusing rather than measuring the wrong "
            "thing.",
            f == 0 ? "lower" : "upper", dirnames[d], keyword,
            refl ? " (with reflection)" : "", dirnames[d], double(coeff[d]),
            why);
    }
  }

  CCTK_VINFO("TestBoundaries: affine outer-boundary test armed. "
             "f = %.17g + %.17g x + %.17g y + %.17g z; tolerance %.17g",
             double(affine_c0), double(affine_cx), double(affine_cy),
             double(affine_cz), double(extrapolation_tolerance));
}

extern "C" void TestBoundaries_WriteAffine(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_WriteAffine;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL c0 = affine_c0, cx = affine_cx, cy = affine_cy,
                  cz = affine_cz;

  // Everywhere in the ground truth. `p.x/y/z` are computed from the point's
  // GLOBAL index (`point_desc`: `x0 + (lbnd + I - 1/2) dx` for a
  // vertex-centred group), so they are correct in the boundary points too.
  Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    affine_exact(p.I) = c0 + cx * p.x + cy * p.y + cz * p.z;
  });

  // Interior only in the field under test: the sync that follows this routine
  // is then the ONLY writer of its boundary points.
  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    affine(p.I) = c0 + cx * p.x + cy * p.y + cz * p.z;
  });
}

extern "C" void TestBoundaries_ExtrapolationReset(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  classes.clear();
  npoints_seen = 0;
}

extern "C" void TestBoundaries_ExtrapolationReduce(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_ExtrapolationReduce;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL got = affine(p.I);
    const CCTK_REAL want = affine_exact(p.I);
    const CCTK_REAL err = got - want;
    affine_error(p.I) = err;

    using std::fabs;
    const CCTK_REAL aerr = fabs(err);

    /*
     * `p.NI[d]` is the outward boundary normal in direction `d` and is 0 where
     * the point is not in the outer boundary; `p.I0[d]` is the nearest
     * INTERIOR index in that direction (`Loop::GridDescBase::loop_box`).  They
     * are the same pair the boundary kernel calls `inormal` and `src`, so
     * `|I - I0|` is exactly the `|dst[d] - src[d]|` whose mis-scaling is the
     * defect.  Taking them from the loop rather than recomputing them from
     * `ghost_size` is deliberate: `boundary_box` derives them from `bbox`, so
     * this classification is right for a multi-box grid too, where a box that
     * does not touch the domain boundary contributes only interior points.
     */
    std::array<int, 3> key = {0, 0, 0};
    for (int d = 0; d < 3; ++d)
      if (p.NI[d] != 0)
        key[d] = p.I[d] > p.I0[d] ? p.I[d] - p.I0[d] : p.I0[d] - p.I[d];

#pragma omp critical(TestBoundaries_extrapolation_accum)
    {
      ++npoints_seen;
      class_stats &cs = classes[key];
      ++cs.npoints;
      if (aerr > cs.max_abs_error || cs.npoints == 1) {
        cs.max_abs_error = aerr;
        cs.worst_I = {p.I[0], p.I[1], p.I[2]};
        cs.worst_X = {p.x, p.y, p.z};
        cs.worst_got = got;
        cs.worst_want = want;
      }
    }
  });
}

extern "C" void TestBoundaries_ExtrapolationReport(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (int(classes.size()) > extrapolation_max_classes)
    CCTK_VERROR("TestBoundaries: %d distinct offset classes, more than "
                "extrapolation_max_classes = %d. Either the rig is far bigger "
                "than this instrument was meant for or the classification is "
                "wrong; refusing to print either way.",
                int(classes.size()), int(extrapolation_max_classes));

  CCTK_VINFO("TestBoundaries: outer-boundary extrapolation test, %lld points "
             "in %d offset classes",
             npoints_seen, int(classes.size()));
  const bc_config c = get_bc_config();
  CCTK_VINFO("  boundary_x/y/z       = \"%s\" / \"%s\" / \"%s\"",
             c.boundary[0][0], c.boundary[0][1], c.boundary[0][2]);
  CCTK_VINFO("  boundary_upper_x/y/z = \"%s\" / \"%s\" / \"%s\"",
             c.boundary[1][0], c.boundary[1][1], c.boundary[1][2]);
  CCTK_VINFO("  reflection lower/upper = %d%d%d / %d%d%d",
             c.reflection[0][0], c.reflection[0][1], c.reflection[0][2],
             c.reflection[1][0], c.reflection[1][1], c.reflection[1][2]);
  CCTK_VINFO("  periodic x/y/z = %d%d%d", c.periodic[0], c.periodic[1],
             c.periodic[2]);
  CCTK_VINFO("  ghost_size = %d, ghost_size_x/y/z = %d / %d / %d",
             int(c.ghost_size), int(c.ghost_size_d[0]),
             int(c.ghost_size_d[1]), int(c.ghost_size_d[2]));

  // Per-rank aggregates, keyed on how many directions are boundary directions:
  // 0 = interior, 1 = face, 2 = edge, 3 = corner.
  long long rank_npoints[4] = {0, 0, 0, 0};
  CCTK_REAL rank_maxerr[4] = {0, 0, 0, 0};
  // Split every non-interior class into "all its nonzero offsets are equal"
  // and "they are not" -- the pre-D1 kernel is exact on the former by
  // arithmetic accident and wrong on the latter.
  long long eq_npoints[2] = {0, 0};
  CCTK_REAL eq_maxerr[2] = {0, 0};

  CCTK_VINFO("  %-14s %10s %24s %26s %24s %24s", "offsets", "npoints",
             "max|error|", "worst point (i,j,k)", "got", "want");
  for (const auto &kv : classes) {
    const std::array<int, 3> &key = kv.first;
    const class_stats &cs = kv.second;

    int rank = 0, first = 0;
    bool all_equal = true;
    for (int d = 0; d < 3; ++d)
      if (key[d] != 0) {
        ++rank;
        if (first == 0)
          first = key[d];
        else if (key[d] != first)
          all_equal = false;
      }

    rank_npoints[rank] += cs.npoints;
    if (cs.max_abs_error > rank_maxerr[rank])
      rank_maxerr[rank] = cs.max_abs_error;
    if (rank > 0) {
      const int bucket = all_equal ? 0 : 1;
      eq_npoints[bucket] += cs.npoints;
      if (cs.max_abs_error > eq_maxerr[bucket])
        eq_maxerr[bucket] = cs.max_abs_error;
    }

    char offs[32];
    std::snprintf(offs, sizeof offs, "(%d,%d,%d)", key[0], key[1], key[2]);
    CCTK_VINFO("  %-14s %10lld %24.17g (%6d,%6d,%6d) %24.17g %24.17g", offs,
               cs.npoints, double(cs.max_abs_error), cs.worst_I[0],
               cs.worst_I[1], cs.worst_I[2], double(cs.worst_got),
               double(cs.worst_want));
  }

  static const char *const rank_label[4] = {"interior", "face", "edge",
                                            "corner"};
  for (int r = 0; r < 4; ++r)
    CCTK_VINFO("  %-8s npoints=%lld max|error|=%.17g", rank_label[r],
               rank_npoints[r], double(rank_maxerr[r]));
  CCTK_VINFO("  boundary, all offsets equal     npoints=%lld max|error|=%.17g",
             eq_npoints[0], double(eq_maxerr[0]));
  CCTK_VINFO("  boundary, offsets NOT all equal npoints=%lld max|error|=%.17g",
             eq_npoints[1], double(eq_maxerr[1]));

  CCTK_REAL maxerr = 0;
  for (const auto &kv : classes)
    if (kv.second.max_abs_error > maxerr)
      maxerr = kv.second.max_abs_error;

  /*
   * Flush before any of the refusals below.  `CCTK_VERROR` aborts, `stdout` is
   * a fully buffered FILE when it is redirected to a file, and an abort throws
   * the buffer away: the FAILING run -- the one whose table is the evidence --
   * lost all of it on two of the three D1 legs and the last 27 lines on the
   * third.  The table is the product of this routine; it has to reach the log
   * before the verdict does.
   */
  std::fflush(nullptr);

  // A pass on zero checked points is a green light for nothing
  // (BUGFIX_TODO.md Traps).
  if (npoints_seen == 0)
    CCTK_VERROR("TestBoundaries: the extrapolation test checked ZERO points. "
                "It cannot pass. Something is wrong with the schedule or with "
                "the grid.");
  if (rank_npoints[1] == 0)
    CCTK_VERROR("TestBoundaries: the extrapolation test found no FACE boundary "
                "points at all (%lld points seen). Nothing about the boundary "
                "kernel has been measured.",
                npoints_seen);

  if (maxerr > extrapolation_tolerance) {
    if (extrapolation_abort)
      CCTK_VERROR("TestBoundaries: FAILED. max|affine - affine_exact| = %.17g "
                  "> extrapolation_tolerance = %.17g. Every configured "
                  "boundary condition reproduces this affine field exactly, so "
                  "any nonzero error here is the boundary kernel's.",
                  double(maxerr), double(extrapolation_tolerance));
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestBoundaries: FAILED (extrapolation_abort = no). "
               "max|error| = %.17g > tolerance %.17g",
               double(maxerr), double(extrapolation_tolerance));
  } else {
    CCTK_VINFO("TestBoundaries: PASSED. max|affine - affine_exact| = %.17g <= "
               "extrapolation_tolerance = %.17g over %lld points",
               double(maxerr), double(extrapolation_tolerance), npoints_seen);
  }
}

} // namespace TestBoundaries
