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

    if ((symmetry_x == symmetry_t::none && boundary_x == boundary_t::none) ||
        (symmetry_x == symmetry_t::interpatch &&
         boundary_x == boundary_t::none) ||
        symmetry_x == symmetry_t::periodic)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::none>(
          bmin, bmax);
    else if (symmetry_x == symmetry_t::reflection)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::reflection,
                           boundary_t::none>(bmin, bmax);
    else if (boundary_x == boundary_t::dirichlet)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::dirichlet>(
          bmin, bmax);
    else if (boundary_x == boundary_t::linear_extrapolation)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none,
                           boundary_t::linear_extrapolation>(bmin, bmax);
    else if (boundary_x == boundary_t::neumann)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::neumann>(
          bmin, bmax);
    else if (boundary_x == boundary_t::robin)
      apply_on_face_symbcx<NI, NJ, NK, symmetry_t::none, boundary_t::robin>(
          bmin, bmax);
    else
#pragma omp critical
      CCTK_ERROR("internal error");
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

    if ((symmetry_y == symmetry_t::none && boundary_y == boundary_t::none) ||
        (symmetry_y == symmetry_t::interpatch &&
         boundary_y == boundary_t::none) ||
        symmetry_y == symmetry_t::periodic)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::none>(bmin, bmax);
    else if (symmetry_y == symmetry_t::reflection)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::reflection,
                            boundary_t::none>(bmin, bmax);
    else if (boundary_y == boundary_t::dirichlet)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::dirichlet>(bmin, bmax);
    else if (boundary_y == boundary_t::linear_extrapolation)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::linear_extrapolation>(bmin, bmax);
    else if (boundary_y == boundary_t::neumann)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::neumann>(bmin, bmax);
    else if (boundary_y == boundary_t::robin)
      apply_on_face_symbcxy<NI, NJ, NK, SCI, BCI, symmetry_t::none,
                            boundary_t::robin>(bmin, bmax);
    else
#pragma omp critical
      CCTK_ERROR("internal error");
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

    if ((symmetry_z == symmetry_t::none && boundary_z == boundary_t::none) ||
        (symmetry_z == symmetry_t::interpatch &&
         boundary_z == boundary_t::none) ||
        symmetry_z == symmetry_t::periodic)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::none>(bmin, bmax);
    else if (symmetry_z == symmetry_t::reflection)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ,
                             symmetry_t::reflection, boundary_t::none>(bmin,
                                                                       bmax);
    else if (boundary_z == boundary_t::dirichlet)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::dirichlet>(bmin, bmax);
    else if (boundary_z == boundary_t::linear_extrapolation)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::linear_extrapolation>(bmin, bmax);
    else if (boundary_z == boundary_t::neumann)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::neumann>(bmin, bmax);
    else if (boundary_z == boundary_t::robin)
      apply_on_face_symbcxyz<NI, NJ, NK, SCI, BCI, SCJ, BCJ, symmetry_t::none,
                             boundary_t::robin>(bmin, bmax);
    else
#pragma omp critical
      CCTK_ERROR("internal error");
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
    CCTK_VERROR("Internal error: Found ncomps=%d, maxncomps=%d when applying "
                "boundary conditions",
                ncomps, maxncomps);
  const int cmin = 0;
  const int cmax = ncomps;

  // Periodic symmetries should have been translated to `none`
  static_assert(!any(symmetries == symmetry_t::periodic));

  if constexpr (all(symmetries == symmetry_t::none &&
                    boundaries == boundary_t::none)) {
    // If there are no boundary conditions to apply, then do
    // nothing.

    // do nothing

  } else {
    // This is the generic case for applying boundary conditions.

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
            dmin = dmin, dmax = dmax,
#endif
            xmin = xmin, dx = dx, layout = layout, destptr = destptr1,
            //
            cmin, cmax, dirichlet_values, neumann_source,
            linear_extrapolation_source, robin_source, robin_values,
            reflection_offset,
            reflection_parities] CCTK_DEVICE(const Arith::vect<int, dim> &dst)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
#ifdef CCTK_DEBUG
              assert(all(src >= dmin && src < dmax));
#endif

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
                           CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
