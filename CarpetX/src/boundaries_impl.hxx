#ifndef CARPETX_CARPETX_BOUNDARIES_IMPL_HXX
#define CARPETX_CARPETX_BOUNDARIES_IMPL_HXX

#include "boundaries.hxx"

#include <array>
#include <functional>
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
#ifdef CCTK_DEBUG
    {
      if (symmetry_x == symmetry_t::interpatch &&
          boundary_x != boundary_t::none &&
          boundary_x != boundary_t::symmetry_boundary) {
#pragma omp critical
        {
          CCTK_VINFO("apply_on_face: group '%s' patch %d x-face f=%d is "
                     "symmetry_t::interpatch; suppressing configured outer BC "
                     "(boundary=%d). Ghost zone will be filled by "
                     "MultiPatch_Interpolate.",
                     groupdata.groupname.c_str(), patchdata.patch, f,
                     static_cast<int>(boundary_x));
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
     * cells instead of being suppressed. Thus we treat ALL interpatch faces as
     * none/none regardless of groupdata.boundaries.
     */
    if ((symmetry_x == symmetry_t::none && boundary_x == boundary_t::none) ||
        symmetry_x == symmetry_t::interpatch ||
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

#ifdef CCTK_DEBUG
    if (symmetry_y == symmetry_t::interpatch &&
        boundary_y != boundary_t::none &&
        boundary_y != boundary_t::symmetry_boundary) {
#pragma omp critical
      CCTK_VINFO("apply_on_face: group '%s' patch %d y-face f=%d is "
                 "symmetry_t::interpatch; suppressing configured outer BC "
                 "(boundary=%d). Ghost zone will be filled by "
                 "MultiPatch_Interpolate.",
                 groupdata.groupname.c_str(), patchdata.patch, f,
                 static_cast<int>(boundary_y));
    }
#endif // CCTK_DEBUG

    if ((symmetry_y == symmetry_t::none && boundary_y == boundary_t::none) ||
        symmetry_y == symmetry_t::interpatch ||
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

#ifdef CCTK_DEBUG
    if (symmetry_z == symmetry_t::interpatch &&
        boundary_z != boundary_t::none &&
        boundary_z != boundary_t::symmetry_boundary) {
#pragma omp critical
      CCTK_VINFO("apply_on_face: group '%s' patch %d z-face f=%d is "
                 "symmetry_t::interpatch; suppressing configured outer BC "
                 "(boundary=%d). Ghost zone will be filled by "
                 "MultiPatch_Interpolate.",
                 groupdata.groupname.c_str(), patchdata.patch, f,
                 static_cast<int>(boundary_z));
    }
#endif // CCTK_DEBUG

    if ((symmetry_z == symmetry_t::none && boundary_z == boundary_t::none) ||
        symmetry_z == symmetry_t::interpatch ||
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
     * CapyrX_MultiPatch/src/interpolate.cxx). This first BC pass therefore
     * writes incorrect values to those corner cells because their interpatch
     * ghost sources are not yet populated. SyncGroupsByDirI corrects this with
     * a second apply_boundary_conditions call immediately after
     * MultiPatch_Interpolate (see schedule.cxx, the block beginning "Second BC
     * pass"). The warning below fires only during this first pass and is
     * expected; it does NOT indicate a permanent error.
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
#pragma omp critical
        CCTK_VINFO("apply_on_face_symbcxyz: [group '%s' patch %d] Corner-cell "
                   "scenario (first BC pass): applying BC [NI=%d NJ=%d NK=%d] "
                   "on bmin=[%d,%d,%d] bmax=[%d,%d,%d] with imin=[%d,%d,%d] "
                   "imax=[%d,%d,%d]. One or more pass-through directions "
                   "extend into the interpatch ghost zone, which is not yet "
                   "populated by MultiPatch_Interpolate. These corner cells "
                   "are SKIPPED by MultiPatch_Interpolate and will be "
                   "re-corrected by the second BC pass in SyncGroupsByDirI "
                   "(schedule.cxx). This message is expected during the first "
                   "pass and does not indicate a bug.",
                   groupdata.groupname.c_str(), patchdata.patch, NI, NJ, NK,
                   bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2],
                   imin[0], imin[1], imin[2], imax[0], imax[1], imax[2]);
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

    const auto kernel =
        [
#ifdef CCTK_DEBUG
            dmin = dmin, dmax = dmax, imin = imin, imax = imax,
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
