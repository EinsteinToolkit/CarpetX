#ifndef CARPETX_CARPETX_BOUNDARIES_IMPL_HXX
#define CARPETX_CARPETX_BOUNDARIES_IMPL_HXX

#include "boundaries.hxx"

#include <array>
#include <functional>
#include <sstream>
#include <type_traits>

namespace CarpetX {

// TODO: The algorithm below applies symmetry and boundary conditions
// in all directions "at once", and thus has to be able to apply
// multiple boundary conditions simultaneously. This adds complexity
// and makes the code slow to build.
//
// An alternative would be to apply symmetry and boundary conditions
// sequentially, first on faces, then edges, then corners. On edges
// and corners we would choose one (preferred) symmetry or boundary
// condition to be applied there. This would simplify the code, reduce
// build times, and also improve run time. It would not be necessary
// any more to run on edges or corners; instead, one could extend some
// of the faces to cover edges and corners.

////////////////////////////////////////////////////////////////////////////////

constexpr int NEG = -1, INT = 0, POS = +1;

constexpr int maxncomps = 16;

template <int NI, int NJ, int NK>
void BoundaryCondition::apply_on_face() const {
  constexpr Arith::vect<int, dim> inormal{NI, NJ, NK};
  static_assert(!all(inormal == 0));

  using std::max, std::min;
  Arith::vect<int, dim> bmin, bmax;
  for (int d = 0; d < dim; ++d) {
    if (inormal[d] < 0) {
      // We have to fill the lower boundary
      bmin[d] = dmin[d];
      bmax[d] = min(dmax[d], imin[d]);
    } else if (inormal[d] > 0) {
      // We have to fill the upper boundary
      bmin[d] = max(dmin[d], imax[d]);
      bmax[d] = dmax[d];
    } else {
      // This direction is not a boundary
      bmin[d] = max(dmin[d], imin[d]);
      bmax[d] = min(dmax[d], imax[d]);
    }
  }

  // Do we actually own part of this boundary?
  bool npoints_is_zero = false;
  for (int d = 0; d < dim; ++d)
    npoints_is_zero |= bmax[d] <= bmin[d];
  if (npoints_is_zero)
    return;

  /*
   * BUGFIX_TODO.md step B2(b): the two partitioned passes are disjoint here.
   *
   * An "interpatch corner" is a region that is outside the patch domain in an
   * INTERPATCH direction and, simultaneously, in a direction carrying a real
   * outer boundary or symmetry condition. `MultiPatch_Interpolate` never fills
   * such a cell (it skips any ghost with `p.NI[d] != 0` on an outer face), and
   * the only honest source for it is the pure interpatch ghost one step inward
   * -- which does not exist until the interpolator has run. So pass 1 skips
   * these, the interpolator runs, and pass 2 writes them and nothing else.
   *
   * Two traps, both of which this had to be written around:
   *
   *  1. The dispatch below ERASES interpatch-ness: it maps every interpatch
   *     face to `symmetry_t::none` / `boundary_t::none`, and there is a
   *     `static_assert` downstream that it did. The classification therefore
   *     cannot be made anywhere below this point, and needs its own loop over
   *     the three directions.
   *
   *  2. `groupdata.boundaries` still carries the globally configured outer BC
   *     on interpatch faces (`get_group_boundaries` uses `get_symmetries(-1)`),
   *     and will until step B7. So an interpatch direction must be excluded
   *     from `has_outer_bc` EXPLICITLY -- the `!ip &&`. Without it an
   *     interpatch x interpatch edge reads as a "corner", pass 2 writes cells
   *     the interpolator owns, and this commit re-creates the double write it
   *     exists to remove.
   *
   * One classification is made here silently and is worth naming. A reflection
   * or periodic face carries `boundary_t::symmetry_boundary`, which is
   * `!= none`, so an interpatch x REFLECTION edge is counted as a corner and
   * deferred to pass 2. That is harmless -- it is still exactly one write, and
   * pass 2 runs where the interpatch leg is filled -- but the class is
   * `interpatch x (outer BC or symmetry)`, not `interpatch x outer`. (The
   * always-on `CCTK_VERROR` below forbids interpatch and `symmetry_boundary`
   * on the SAME face, which is a different thing.)
   *
   * A pure interpatch face, and an interpatch x interpatch edge or corner,
   * need no flag: the dispatch maps them to all-`none`, and the
   * `if constexpr (all(symmetries == none && boundaries == none))` early-out
   * in `apply_on_face_symbcxyz` already skips them.
   *
   * The whole block is inside `if (bc_pass != all)` so that a run which never
   * selects a partitioned pass -- every single-patch run, since the call site
   * gates the selection on `have_multipatch_boundaries` -- executes literally
   * none of it. That is what makes the inertness claim structural rather than
   * merely measured.
   */
  if (bc_pass != bc_pass_t::all) {
    bool has_interpatch = false, has_outer_bc = false;
    for (int d = 0; d < dim; ++d) {
      if (inormal[d] == 0)
        continue;
      const int f = inormal[d] > 0;
      const bool ip = patchdata.symmetries[f][d] == symmetry_t::interpatch;
      has_interpatch |= ip;
      has_outer_bc |= !ip && groupdata.boundaries[f][d] != boundary_t::none;
    }
    const bool is_interpatch_corner = has_interpatch && has_outer_bc;
    const bool skip =
        (bc_pass == bc_pass_t::skip_interpatch_corners && is_interpatch_corner) ||
        (bc_pass == bc_pass_t::interpatch_corners_only && !is_interpatch_corner);

#ifdef CCTK_DEBUG
    if (bc_pass_census_level() > 0) {
      long long ncells = 1;
      for (int d = 0; d < dim; ++d)
        ncells *= bmax[d] - bmin[d];
      auto &cells = is_interpatch_corner
                        ? (skip ? bc_pass_census.corner_cells_skipped
                                : bc_pass_census.corner_cells_kept)
                        : (skip ? bc_pass_census.other_cells_skipped
                                : bc_pass_census.other_cells_kept);
      auto &regions = is_interpatch_corner
                          ? (skip ? bc_pass_census.corner_regions_skipped
                                  : bc_pass_census.corner_regions_kept)
                          : (skip ? bc_pass_census.other_regions_skipped
                                  : bc_pass_census.other_regions_kept);
      cells += ncells;
      regions += 1;
      if (bc_pass_census_level() > 1) {
        std::ostringstream buf;
        buf << bc_pass;
#pragma omp critical
        CCTK_VINFO("BCPASSREGION group=%s patch=%d level=%d pass=%s "
                   "n=(%d,%d,%d) corner=%d skip=%d ncells=%lld "
                   "bmin=(%d,%d,%d) bmax=(%d,%d,%d)",
                   groupdata.groupname.c_str(), patchdata.patch,
                   groupdata.level, buf.str().c_str(), NI, NJ, NK,
                   int(is_interpatch_corner), int(skip), ncells, bmin[0],
                   bmin[1], bmin[2], bmax[0], bmax[1], bmax[2]);
      }
    }
#endif // CCTK_DEBUG

    if (skip)
      return;
  }

  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners multiple conditions will apply.
  if constexpr (NI == 0) {
    // interior
    apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::none>(bmin,
                                                                         bmax);
  } else {
    // face
    const int f = NI > 0;
    const symmetry_t symmetry_x = patchdata.symmetries[f][0];
    const boundary_t boundary_x = groupdata.boundaries[f][0];

    /*
     * Assert that that a face that has interpatch symmetry does not have any
     * kind of boundary applied to it. symmetry_boundary is reserved for
     * periodic/reflection faces. interpatch is patch-specific. Both appearing
     * at once is an internal inconsistency in the boundary setup.
     */
    if (symmetry_x == symmetry_t::interpatch &&
        boundary_x == boundary_t::symmetry_boundary) {
#pragma omp critical
      CCTK_VERROR(
          "Group '%s' patch %d: x-face f=%d is both symmetry_t::interpatch "
          "(patchdata.symmetries) and boundary_t::symmetry_boundary "
          "(groupdata.boundaries). This is an internal inconsistency: "
          "interpatch and symmetry_boundary are mutually exclusive.",
          groupdata.groupname.c_str(), patchdata.patch, f);
    }

    /*
     * get_group_boundaries() uses get_symmetries(-1). The -1 means "ignore
     * patches", so interpatch faces receive whatever outer BC is globally
     * configured in groupdata.boundaries.  We detect and suppress that here,
     * letting MultiPatch_Interpolate fill the interpatch ghost zone instead.
     */
    // B2(a) made this message pass-dependent: it used to say "suppressing"
    // unconditionally, which is now true for the two partitioned passes and
    // FALSE for `bc_pass_t::all`, where the configured outer BC really is
    // applied to the interpatch face (`main`'s behaviour, kept deliberately --
    // see the dispatch comment below).
#ifdef CCTK_DEBUG
    {
      if (symmetry_x == symmetry_t::interpatch &&
          boundary_x != boundary_t::none &&
          boundary_x != boundary_t::symmetry_boundary) {
        std::ostringstream buf;
        buf << bc_pass;
#pragma omp critical
        {
          CCTK_VINFO("apply_on_face: group '%s' patch %d x-face f=%d is "
                     "symmetry_t::interpatch with configured outer BC "
                     "(boundary=%d); bc_pass=%s => %s",
                     groupdata.groupname.c_str(), patchdata.patch, f,
                     static_cast<int>(boundary_x), buf.str().c_str(),
                     bc_pass == bc_pass_t::all
                         ? "APPLYING it (pre-B7 `all` behaviour)"
                         : "suppressing it; the ghost zone is filled by "
                           "MultiPatch_Interpolate");
        }
      }
    }
#endif // CCTK_DEBUG

    /*
     * The original condition required `boundary_x == none` for the interpatch
     * branch, but get_group_boundaries() (which is patch-agnostic) sets
     * boundary_x to the globally configured outer BC even for interpatch faces.
     * The original condition therefore NEVER matched for interpatch faces with
     * a configured outer BC, causing the outer BC to fire on interpatch ghost
     * cells instead of being suppressed.
     *
     * BUGFIX_TODO.md step B2(a): suppressing it UNCONDITIONALLY is wrong too.
     * A coarse temporary's interpatch faces are written by nobody -- the
     * interpolator only ever touches the real `mfab` -- and `mf_set_domain_bndry`
     * has just NaN-filled them, so `FillPatchInterp` would prolongate NaN. The
     * suppression is therefore tied to the pass: the two partitioned passes
     * (which run only where `MultiPatch_Interpolate` follows/precedes them) map
     * interpatch to none/none, and `bc_pass_t::all` keeps `main`'s behaviour of
     * writing the configured outer BC there -- finite garbage that something
     * later overwrites, rather than nothing at all.
     *
     * NOTE: after step B7 removes the outer BC from `groupdata.boundaries` on
     * interpatch faces, `boundary_x == none` holds there and the first
     * disjunct's second half makes this clause fire in EVERY pass again. That
     * is deliberate, and B8 is what refuses the configuration in which it
     * matters (more than one time level).
     */
    if ((symmetry_x == symmetry_t::none && boundary_x == boundary_t::none) ||
        (symmetry_x == symmetry_t::interpatch &&
         (boundary_x == boundary_t::none || bc_pass != bc_pass_t::all)) ||
        symmetry_x == symmetry_t::periodic) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::none>(
          bmin, bmax);
    } else if (symmetry_x == symmetry_t::reflection) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::reflection,
                           boundary_t::none>(bmin, bmax);
    } else if (boundary_x == boundary_t::dirichlet) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::dirichlet>(
          bmin, bmax);
    } else if (boundary_x == boundary_t::linear_extrapolation) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none,
                           boundary_t::linear_extrapolation>(bmin, bmax);
    } else if (boundary_x == boundary_t::neumann) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::neumann>(
          bmin, bmax);
    } else if (boundary_x == boundary_t::robin) {
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::robin>(
          bmin, bmax);
    } else {
#pragma omp critical
      CCTK_ERROR("internal error");
    }
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI>
void BoundaryCondition::apply_on_face_symbcx(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners multiple conditions will apply.
  if constexpr (NJ == 0) {
    // interior
    apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                          boundary_t::none>(bmin, bmax);
  } else {
    // face
    const int f = NJ > 0;
    const symmetry_t symmetry_y = patchdata.symmetries[f][1];
    const boundary_t boundary_y = groupdata.boundaries[f][1];

    // Hard invariant: interpatch and symmetry_boundary are mutually exclusive.
    if (symmetry_y == symmetry_t::interpatch &&
        boundary_y == boundary_t::symmetry_boundary) {
#pragma omp critical
      CCTK_VERROR(
          "Group '%s' patch %d: y-face f=%d is both symmetry_t::interpatch "
          "(patchdata.symmetries) and boundary_t::symmetry_boundary "
          "(groupdata.boundaries). This is an internal inconsistency: "
          "interpatch and symmetry_boundary are mutually exclusive.",
          groupdata.groupname.c_str(), patchdata.patch, f);
    }

    // B2(a): pass-dependent, see the x-face copy.
#ifdef CCTK_DEBUG
    if (symmetry_y == symmetry_t::interpatch &&
        boundary_y != boundary_t::none &&
        boundary_y != boundary_t::symmetry_boundary) {
      std::ostringstream buf;
      buf << bc_pass;
#pragma omp critical
      CCTK_VINFO("apply_on_face: group '%s' patch %d y-face f=%d is "
                 "symmetry_t::interpatch with configured outer BC "
                 "(boundary=%d); bc_pass=%s => %s",
                 groupdata.groupname.c_str(), patchdata.patch, f,
                 static_cast<int>(boundary_y), buf.str().c_str(),
                 bc_pass == bc_pass_t::all
                     ? "APPLYING it (pre-B7 `all` behaviour)"
                     : "suppressing it; the ghost zone is filled by "
                       "MultiPatch_Interpolate");
    }
#endif // CCTK_DEBUG

    // B2(a): the interpatch suppression is pass-dependent. See the long
    // comment on the x-face above.
    if ((symmetry_y == symmetry_t::none && boundary_y == boundary_t::none) ||
        (symmetry_y == symmetry_t::interpatch &&
         (boundary_y == boundary_t::none || bc_pass != bc_pass_t::all)) ||
        symmetry_y == symmetry_t::periodic) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::none>(bmin, bmax);
    } else if (symmetry_y == symmetry_t::reflection) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::reflection,
                            boundary_t::none>(bmin, bmax);
    } else if (boundary_y == boundary_t::dirichlet) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::dirichlet>(bmin, bmax);
    } else if (boundary_y == boundary_t::linear_extrapolation) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::linear_extrapolation>(bmin, bmax);
    } else if (boundary_y == boundary_t::neumann) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::neumann>(bmin, bmax);
    } else if (boundary_y == boundary_t::robin) {
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::robin>(bmin, bmax);
    } else {
#pragma omp critical
      CCTK_ERROR("internal error");
    }
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
          symmetry_t SCJ, boundary_t BCJ>
void BoundaryCondition::apply_on_face_symbcxy(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners multiple conditions will apply.
  if constexpr (NK == 0) {
    // interior
    apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                           boundary_t::none>(bmin, bmax);
  } else {
    // face
    const int f = NK > 0;
    const symmetry_t symmetry_z = patchdata.symmetries[f][2];
    const boundary_t boundary_z = groupdata.boundaries[f][2];

    // Hard invariant: interpatch and symmetry_boundary are mutually exclusive.
    if (symmetry_z == symmetry_t::interpatch &&
        boundary_z == boundary_t::symmetry_boundary) {
#pragma omp critical
      CCTK_VERROR(
          "Group '%s' patch %d: z-face f=%d is both symmetry_t::interpatch "
          "(patchdata.symmetries) and boundary_t::symmetry_boundary "
          "(groupdata.boundaries). This is an internal inconsistency: "
          "interpatch and symmetry_boundary are mutually exclusive.",
          groupdata.groupname.c_str(), patchdata.patch, f);
    }

    // B2(a): pass-dependent, see the x-face copy.
#ifdef CCTK_DEBUG
    if (symmetry_z == symmetry_t::interpatch &&
        boundary_z != boundary_t::none &&
        boundary_z != boundary_t::symmetry_boundary) {
      std::ostringstream buf;
      buf << bc_pass;
#pragma omp critical
      CCTK_VINFO("apply_on_face: group '%s' patch %d z-face f=%d is "
                 "symmetry_t::interpatch with configured outer BC "
                 "(boundary=%d); bc_pass=%s => %s",
                 groupdata.groupname.c_str(), patchdata.patch, f,
                 static_cast<int>(boundary_z), buf.str().c_str(),
                 bc_pass == bc_pass_t::all
                     ? "APPLYING it (pre-B7 `all` behaviour)"
                     : "suppressing it; the ghost zone is filled by "
                       "MultiPatch_Interpolate");
    }
#endif // CCTK_DEBUG

    // B2(a): the interpatch suppression is pass-dependent. See the long
    // comment on the x-face above.
    if ((symmetry_z == symmetry_t::none && boundary_z == boundary_t::none) ||
        (symmetry_z == symmetry_t::interpatch &&
         (boundary_z == boundary_t::none || bc_pass != bc_pass_t::all)) ||
        symmetry_z == symmetry_t::periodic) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::none>(bmin, bmax);
    } else if (symmetry_z == symmetry_t::reflection) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ,
                             symmetry_t::reflection, boundary_t::none>(bmin,
                                                                       bmax);
    } else if (boundary_z == boundary_t::dirichlet) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::dirichlet>(bmin, bmax);
    } else if (boundary_z == boundary_t::linear_extrapolation) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::linear_extrapolation>(bmin, bmax);
    } else if (boundary_z == boundary_t::neumann) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::neumann>(bmin, bmax);
    } else if (boundary_z == boundary_t::robin) {
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::robin>(bmin, bmax);
    } else {
#pragma omp critical
      CCTK_ERROR("internal error");
    }
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
          symmetry_t SCJ, boundary_t BCJ, symmetry_t SCK, boundary_t BCK>
void BoundaryCondition::apply_on_face_symbcxyz(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  constexpr Arith::vect<int, dim> inormal{NI, NJ, NK};
  static_assert(!all(inormal == 0));
  constexpr Arith::vect<boundary_t, dim> boundaries{BCI, BCJ, BCK};
  constexpr Arith::vect<symmetry_t, dim> symmetries{SCI, SCJ, SCK};

  // TODO: Move loop over components to the far outside

  const int ncomps = dest.nComp();
  if (CCTK_BUILTIN_EXPECT(ncomps > maxncomps, false))
    CCTK_VERROR("apply_on_face_symbcxyz Internal error: Found ncomps=%d, "
                "maxncomps=%d when applying "
                "boundary conditions",
                ncomps, maxncomps);
  const int cmin = 0;
  const int cmax = ncomps;

  // Periodic symmetries should have been translated to `none`
  static_assert(!any(symmetries == symmetry_t::periodic));

  /*
   * symmetry_t::interpatch must never reach this function. The fix in
   * apply_on_face{,_symbcx,_symbcxy} translates every interpatch face to
   * symmetry_t::none before the template chain reaches here.
   */
  static_assert(
      !any(symmetries == symmetry_t::interpatch),
      "symmetry_t::interpatch reached apply_on_face_symbcxyz. apply_on_face "
      "must always map interpatch to symmetry_t::none.");

  /*
   * boundary_t::symmetry_boundary is a sentinel type only. It must never be
   * forwarded for actual BC application.
   */
  static_assert(
      !any(boundaries == boundary_t::symmetry_boundary),
      "boundary_t::symmetry_boundary reached apply_on_face_symbcxyz. This "
      "sentinel type must not be used for actual BC application.");

  // mp_slave_7.md §5 / mp_slave_8.md instrumentation: bucket (b) (312
  // "corner" cells at an angular-interpatch face x the wedges' cartesian-
  // facing inner-radial face) empirically cleared when Fix #1 deferred
  // slave_overlap's write-back, even though mp_slave_3.md §5 showed those
  // cells never enter MultiPatch1_Interpolate's ghost-fill/slave-write loop
  // at all. This logs, for the CAPYRX_TESTMULTIPATCH::COLOR group only and
  // only for genuine 2-or-3-axis corners (bucket (b)'s cell class), which
  // branch below actually fires for this (NI,NJ,NK) face combination -- the
  // "NOOP" branch (all 3 axes none/none: nothing is ever written here by
  // apply_boundary_conditions, in EITHER BC pass) or the "WRITE" branch
  // (some axis has a real symmetry/boundary condition configured) -- plus
  // the per-axis symmetry_t/boundary_t values that decided it. This settles
  // whether bucket (b)'s cells are, per this function's own logic, ever
  // written by the "2nd BC pass" at all (mp_slave_3.md §5 assumed they are,
  // without checking).
#ifdef CCTK_DEBUG
  {
    constexpr int n_active_axes =
        int(NI != 0) + int(NJ != 0) + int(NK != 0);
    static const bool log_donors = std::getenv("CAPYRX_LOG_DONORS") != nullptr;
    if (log_donors && n_active_axes >= 2 &&
        groupdata.groupname == "CAPYRX_TESTMULTIPATCH::COLOR") {
      constexpr bool is_noop = all(symmetries == symmetry_t::none &&
                                   boundaries == boundary_t::none);
#pragma omp critical
      CCTK_VINFO(
          "CORNERBC branch=%s patch=%d level=%d NI=%d NJ=%d NK=%d "
          "sym=(%d,%d,%d) bnd=(%d,%d,%d) bmin=(%d,%d,%d) bmax=(%d,%d,%d)",
          is_noop ? "NOOP" : "WRITE", patchdata.patch, groupdata.level, NI, NJ,
          NK, static_cast<int>(SCI), static_cast<int>(SCJ),
          static_cast<int>(SCK), static_cast<int>(BCI), static_cast<int>(BCJ),
          static_cast<int>(BCK), bmin[0], bmin[1], bmin[2], bmax[0], bmax[1],
          bmax[2]);
    }
  }
#endif // CCTK_DEBUG

  if constexpr (all(symmetries == symmetry_t::none &&
                    boundaries == boundary_t::none)) {
    // If there are no boundary conditions to apply, then do
    // nothing.

    // do nothing

  } else {
    // This is the generic case for applying boundary conditions.

    /*
     * Detect the "corner cell catastrophe" scenario: a non-trivial BC is being
     * applied but at least one pass-through direction (symmetry=none,
     * boundary=none) has its bmin/bmax range extend into the ghost zone. This
     * happens for cells at the intersection of an interpatch face (one
     * direction) and an outer-BC face (another direction).
     *
     * In that situation the BC source for the pass-through direction is src[d]
     * = dst[d], which may lie in the interpatch ghost zone (not yet populated
     * by MultiPatch_Interpolate at the time apply_boundary_conditions runs).
     * The resulting BC value written to these corner cells is therefore
     * computed from uninitialized ghost data.
     *
     * IMPORTANT: CapyrX's MultiPatch_Interpolate skips corner cells that are on
     * any outer-boundary face (see loop_bnd skip logic in
     * CapyrX_MultiPatch/src/interpolate.cxx), so it never populates them.
     *
     * UNTIL BUGFIX_TODO.md step B2 this was reached by BOTH BC passes: the
     * first wrote these corner cells from a not-yet-populated interpatch ghost
     * zone, and SyncGroupsByDirI called apply_boundary_conditions a SECOND time
     * after MultiPatch_Interpolate to overwrite them from valid sources. Since
     * B2 the passes are disjoint -- pass 1 runs
     * `bc_pass_t::skip_interpatch_corners` and returns before reaching this
     * function for a corner, pass 2 runs `bc_pass_t::interpatch_corners_only`
     * -- so on the partitioned path this scenario is reached ONLY by pass 2,
     * and only after the interpatch sources exist. It is still reachable with
     * `bc_pass_t::all`, which is what every non-sync call site and every
     * `tl >= 1` uses.
     *
     * NOTE ON THE MESSAGE BELOW: the check that emits it is purely geometric --
     * it only tests whether a pass-through direction's bmin/bmax reaches into
     * the ghost zone -- so it says nothing about whether the sources are valid.
     * The pass is printed with it precisely so that the reader can tell the two
     * cases apart: `interpatch_corners_only` is the healthy case (sources
     * populated), `all` is the one where the old hazard survives.
     */
#ifdef CCTK_DEBUG
    {
      bool has_passthrough_in_ghost = false;

      for (int d = 0; d < dim; ++d) {
        if (symmetries[d] == symmetry_t::none &&
            boundaries[d] == boundary_t::none) {
          /*
           * Pass-through: src[d] = dst[d].  If the destination region extends
           * outside [imin,imax) in this direction the source is in the ghost
           * zone.
           */
          if (bmin[d] < imin[d] || bmax[d] > imax[d])
            has_passthrough_in_ghost = true;
        }
      }
      if (has_passthrough_in_ghost) {
        std::ostringstream passbuf;
        passbuf << bc_pass;
#pragma omp critical
        CCTK_VINFO("apply_on_face_symbcxyz: [group '%s' patch %d] Corner-cell "
                   "scenario: applying BC [NI=%d NJ=%d NK=%d] "
                   "on bmin=[%d,%d,%d] bmax=[%d,%d,%d] with imin=[%d,%d,%d] "
                   "imax=[%d,%d,%d] in bc_pass=%s. One or more pass-through "
                   "directions extend into the interpatch ghost zone. These "
                   "corner cells are SKIPPED by MultiPatch_Interpolate. In "
                   "bc_pass=interpatch_corners_only this is the healthy case: "
                   "the pass runs after MultiPatch_Interpolate, so the sources "
                   "are populated. In bc_pass=all the sources may not be, which "
                   "is the hazard BUGFIX_TODO.md step B2 partitions away on the "
                   "sync path. Geometric only; it does not by itself indicate a "
                   "bug.",
                   groupdata.groupname.c_str(), patchdata.patch, NI, NJ, NK,
                   bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2],
                   imin[0], imin[1], imin[2], imax[0], imax[1], imax[2],
                   passbuf.str().c_str());
      }
    }
#endif // CCTK_DEBUG

    Arith::vect<CCTK_REAL, maxncomps> dirichlet_values;
    for (int comp = 0; comp < ncomps; ++comp)
      dirichlet_values[comp] = groupdata.dirichlet_values.at(comp);

    Arith::vect<int, dim> neumann_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::neumann) {
        if (inormal[d] != 0) {
          neumann_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
          if (inormal[d] < 0)
            assert(neumann_source[d] < dmax[d]);
          else
            assert(neumann_source[d] >= dmin[d]);
        }
      } else {
        neumann_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<int, dim> linear_extrapolation_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::linear_extrapolation) {
        assert(inormal[d] != 0);
        linear_extrapolation_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
        if (inormal[d] < 0) {
          assert(linear_extrapolation_source[d] < dmax[d]);
          assert(linear_extrapolation_source[d] - inormal[d] < dmax[d]);
        } else {
          assert(linear_extrapolation_source[d] >= dmin[d]);
          assert(linear_extrapolation_source[d] - inormal[d] >= dmin[d]);
        }
      } else {
        linear_extrapolation_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<int, dim> robin_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::robin) {
        assert(inormal[d] != 0);
        robin_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
        if (inormal[d] < 0)
          assert(robin_source[d] < dmax[d]);
        else
          assert(robin_source[d] >= dmin[d]);
      } else {
        robin_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<CCTK_REAL, maxncomps> robin_values;
    for (int comp = 0; comp < ncomps; ++comp)
      robin_values[comp] = groupdata.robin_values.at(comp);

    Arith::vect<int, dim> reflection_offset;
    for (int d = 0; d < dim; ++d) {
      if (symmetries[d] == symmetry_t::reflection) {
        assert(inormal[d] != 0);
        reflection_offset[d] =
            inormal[d] < 0 ? 2 * imin[d] - groupdata.indextype.at(d)
                           : 2 * (imax[d] - 1) + groupdata.indextype.at(d);
        if (inormal[d] < 0)
          assert(reflection_offset[d] - bmin[d] < dmax[d]);
        else
          assert(reflection_offset[d] - (bmax[d] - 1) >= dmin[d]);
      } else {
        reflection_offset[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<CCTK_REAL, maxncomps> reflection_parities;
    for (int comp = 0; comp < ncomps; ++comp) {
      CCTK_REAL reflection_parity = +1;
      for (int d = 0; d < dim; ++d)
        if (symmetries[d] == symmetry_t::reflection)
          reflection_parity *= groupdata.parities.at(comp).at(d);
      using std::fabs;
      assert(fabs(reflection_parity) == 1);
      reflection_parities[comp] = reflection_parity;
    }

    // We cannot capture `destptr` directly (on Summit, with CUDA 11.5.2)
    // We cannot use a `restrict` declaration either.
    CCTK_REAL *const destptr1 = destptr;

    // mp_slave_7.md §5 / mp_slave_8.md instrumentation: per-cell src/val log
    // for the "WRITE" branch (the CORNERBC log above already reports whether
    // this face combination is NOOP or WRITE at the call level; this adds
    // the actual value for genuine 2+-axis corners in the color group, when
    // it does write).
#ifdef CCTK_DEBUG
    constexpr int n_active_axes = int(NI != 0) + int(NJ != 0) + int(NK != 0);
    static const bool log_donors_cells =
        std::getenv("CAPYRX_LOG_DONORS") != nullptr;
    const bool log_this_corner =
        log_donors_cells && n_active_axes >= 2 &&
        groupdata.groupname == "CAPYRX_TESTMULTIPATCH::COLOR";
    const int log_patch = patchdata.patch;
    const int log_level = groupdata.level;
#endif

    const auto kernel =
        [
#ifdef CCTK_DEBUG
            dmin = dmin, dmax = dmax, imin = imin, imax = imax,
            log_this_corner = log_this_corner, log_patch = log_patch,
            log_level = log_level,
#endif
            xmin = xmin, dx = dx, layout = layout, destptr = destptr1,
            //
            cmin, cmax, dirichlet_values, neumann_source,
            linear_extrapolation_source, robin_source, robin_values,
            reflection_offset,
            reflection_parities] CCTK_DEVICE(const Arith::vect<int, dim> &dst)
            __attribute__((__always_inline__, __flatten__)) {
              constexpr Arith::vect<int, dim> inormal{NI, NJ, NK};
              constexpr Arith::vect<boundary_t, dim> boundaries{BCI, BCJ, BCK};
              constexpr Arith::vect<symmetry_t, dim> symmetries{SCI, SCJ, SCK};

              // `src` is the point at which we are looking to determine
              // the boundary value
              Arith::vect<int, dim> src = dst;
              // `delta` (if nonzero) describes a second point at which
              // we are looking, so that we can calculate a gradient for
              // the boundary value
              Arith::vect<int, dim> delta{0, 0, 0};
              for (int d = 0; d < dim; ++d) {
                if (boundaries[d] == boundary_t::dirichlet) {
                  // do nothing
                } else if (boundaries[d] == boundary_t::linear_extrapolation) {
                  // Same slope:
                  //   f'(0)       = f'(h)
                  //   f(h) - f(0) = f(2h) - f(h)
                  //          f(0) = 2 f(h) - f(2h)
                  // f(0) is the boundary point
                  src[d] = linear_extrapolation_source[d];
                  delta[d] = -inormal[d];
                } else if (boundaries[d] == boundary_t::neumann) {
                  // Same value:
                  //   f(0) = f(h)
                  // f(0) is the boundary point
                  src[d] = neumann_source[d];
                } else if (boundaries[d] == boundary_t::robin) {
                  // Robin condition, specialized to 1/r fall-off:
                  //   f(r) = finf + C/r
                  // Determine fall-off constant `C`:
                  //   C = r * (f(r) - finf)
                  // Solve for value at boundary:
                  //   f(r+h) = finf + C / (r + h)
                  //          = finf + r / (r + h) * (f(r) - finf)
                  // Rewrite using Cartesian coordinates:
                  //   C = |x| * (f(x) - finf)
                  //   f(x') = finf + C / |x'|
                  //         = finf + |x| / |x'| * (f(x) - finf)
                  // f(x') is the boundary point
                  src[d] = robin_source[d];
                } else if (symmetries[d] == symmetry_t::reflection) {
                  src[d] = reflection_offset[d] - dst[d];
                } else if (symmetries[d] == symmetry_t::none &&
                           boundaries[d] == boundary_t::none) {
                  // this direction is not a boundary; do nothing
                } else {
                  // std::cerr << " dst=" << dst << " d=" << d
                  //           << " boundaries=" << boundaries
                  //           << " symmetries=" << symmetries <<
                  //           "\n";
                  assert(0);
                }
              }

      /*
       * For each direction where a non-trivial BC is actively applied
       * (non-dirichlet, non-none), the source coordinate must lie within the
       * domain interior [imin, imax). Sources for
       * Neumann/Robin/linear-extrapolation are always imin[d] or imax[d]-1.
       * Reflection maps dst into the interior. A failure here means the BC
       * stencil setup is internally broken (the source should never be in the
       * ghost zone for an actively-BC'd direction).
       */
#ifdef CCTK_DEBUG
              {
                assert(all(src >= dmin && src < dmax));

                for (int d_src = 0; d_src < dim; ++d_src) {
                  if (boundaries[d_src] != boundary_t::none &&
                      boundaries[d_src] != boundary_t::dirichlet) {
                    assert(src[d_src] >= imin[d_src] &&
                           src[d_src] < imax[d_src]);
                  }
                }
              }
#endif // CCTK_DEBUG

              for (int comp = cmin; comp < cmax; ++comp) {
                const CCTK_REAL dirichlet_value = dirichlet_values[comp];
                const CCTK_REAL robin_value = robin_values[comp];
                const CCTK_REAL reflection_parity = reflection_parities[comp];
                const Loop::GF3D2<CCTK_REAL> var(layout,
                                                 destptr + comp * layout.np);

#ifdef CCTK_DEBUG
                using std::isnan;
#endif

                CCTK_REAL val;
                if constexpr (any(boundaries == boundary_t::dirichlet)) {
                  val = dirichlet_value;
                } else {
                  val = var(src);

#ifdef CCTK_DEBUG
                  assert(!isnan(val));
#endif
                  if constexpr (any(boundaries == boundary_t::robin)) {
                    for (int d = 0; d < dim; ++d) {
                      if (boundaries[d] == boundary_t::robin) {
                        using std::sqrt;
                        // boundary point
                        const auto xb = xmin + dst * dx;
                        const auto rb = sqrt(sum(pow2(xb)));
                        // interior point
                        const auto xi = xmin + src * dx;
                        const auto ri = sqrt(sum(pow2(xi)));
                        const auto q = ri / rb;
                        val = robin_value + q * (val - robin_value);
                      }
                    }
                  }
                  if constexpr (any(boundaries ==
                                    boundary_t::linear_extrapolation)) {
                    // Calculate gradient
                    const CCTK_REAL grad = val - var(src + delta);
                    using std::sqrt;
                    val += sqrt(sum(pow2(dst - src)) / sum(pow2(delta))) * grad;
                  }
#ifdef CCTK_DEBUG
                  for (int d = 0; d < dim; ++d)
                    assert(dst[d] >= dmin[d] && dst[d] < dmax[d]);
#endif
                  if constexpr (any(symmetries == symmetry_t::reflection))
                    val *= reflection_parity;
                }
#ifdef CCTK_DEBUG
                assert(!isnan(val));
                // CCTK_VINFO is host-only; this kernel is also compiled for
                // the device (CCTK_DEVICE lambda), so the logging call must
                // not be compiled in the device pass.
#if !AMREX_DEVICE_COMPILE
                if (log_this_corner && comp == 0) {
#pragma omp critical
                  CCTK_VINFO("CORNERBC_WRITE patch=%d level=%d dst=(%d,%d,%d) "
                             "src=(%d,%d,%d) val=%.17g",
                             log_patch, log_level, dst[0], dst[1], dst[2],
                             src[0], src[1], src[2], double(val));
                }
#endif
#endif
                var.store(dst, val);
              }
            };

    // Note: Calling `loop_region` is much slower than calling `ParallelFor`
    // directly.
    // Maybe the `attribute(noinline)` is to blame?
    // loop_region(kernel, bmin, bmax);
    const amrex::Box box(amrex::IntVect(bmin[0], bmin[1], bmin[2]),
                         amrex::IntVect(bmax[0] - 1, bmax[1] - 1, bmax[2] - 1));
    amrex::ParallelFor(box,
                       [=] CCTK_DEVICE(const int i, const int j, const int k)
                           __attribute__((__always_inline__, __flatten__)) {
                             const Arith::vect<int, dim> p{i, j, k};
                             kernel(p);
                           });
  }
}

extern template void BoundaryCondition::apply_on_face<NEG, NEG, NEG>() const;
extern template void BoundaryCondition::apply_on_face<INT, NEG, NEG>() const;
extern template void BoundaryCondition::apply_on_face<POS, NEG, NEG>() const;
extern template void BoundaryCondition::apply_on_face<NEG, INT, NEG>() const;
extern template void BoundaryCondition::apply_on_face<INT, INT, NEG>() const;
extern template void BoundaryCondition::apply_on_face<POS, INT, NEG>() const;
extern template void BoundaryCondition::apply_on_face<NEG, POS, NEG>() const;
extern template void BoundaryCondition::apply_on_face<INT, POS, NEG>() const;
extern template void BoundaryCondition::apply_on_face<POS, POS, NEG>() const;

extern template void BoundaryCondition::apply_on_face<NEG, NEG, INT>() const;
extern template void BoundaryCondition::apply_on_face<INT, NEG, INT>() const;
extern template void BoundaryCondition::apply_on_face<POS, NEG, INT>() const;
extern template void BoundaryCondition::apply_on_face<NEG, INT, INT>() const;
// extern template void BoundaryCondition::apply_on_face<INT, INT, INT>() const;
extern template void BoundaryCondition::apply_on_face<POS, INT, INT>() const;
extern template void BoundaryCondition::apply_on_face<NEG, POS, INT>() const;
extern template void BoundaryCondition::apply_on_face<INT, POS, INT>() const;
extern template void BoundaryCondition::apply_on_face<POS, POS, INT>() const;

extern template void BoundaryCondition::apply_on_face<NEG, NEG, POS>() const;
extern template void BoundaryCondition::apply_on_face<INT, NEG, POS>() const;
extern template void BoundaryCondition::apply_on_face<POS, NEG, POS>() const;
extern template void BoundaryCondition::apply_on_face<NEG, INT, POS>() const;
extern template void BoundaryCondition::apply_on_face<INT, INT, POS>() const;
extern template void BoundaryCondition::apply_on_face<POS, INT, POS>() const;
extern template void BoundaryCondition::apply_on_face<NEG, POS, POS>() const;
extern template void BoundaryCondition::apply_on_face<INT, POS, POS>() const;
extern template void BoundaryCondition::apply_on_face<POS, POS, POS>() const;

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_BOUNDARIES_IMPL_HXX
