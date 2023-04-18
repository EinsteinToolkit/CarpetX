#include "boundaries.hxx"

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

// CarpetX::loop_region<CarpetX::BoundaryCondition::apply_on_face_symbcxyz<1, 1,
// 0, (CarpetX::symmetry_t)0, (CarpetX::boundary_t)0, (CarpetX::symmetry_t)0,
// (CarpetX::boundary_t)2, (CarpetX::symmetry_t)0,
// (CarpetX::boundary_t)0>(Arith::vect<int, 3> const&, Arith::vect<int, 3>
// const&) const::{lambda(Arith::vect<int, 3>
// const&)#1}>(CarpetX::BoundaryCondition::apply_on_face_symbcxyz<1, 1, 0,
// (CarpetX::symmetry_t)0, (CarpetX::boundary_t)0, (CarpetX::symmetry_t)0,
// (CarpetX::boundary_t)2, (CarpetX::symmetry_t)0,
// (CarpetX::boundary_t)0>(Arith::vect<int, 3> const&, Arith::vect<int, 3>
// const&) const::{lambda(Arith::vect<int, 3> const&)#1} const&,
// Arith::vect<int, 3> const&, Arith::vect<int, 3> const&)
//
// _ZN7CarpetX11loop_regionIZNKS_17BoundaryCondition22apply_on_face_symbcxyzILin1ELin1ELin1ELNS_10symmetry_tE0ELNS_10boundary_tE5ELS3_4ELS4_0ELS3_4ELS4_0EEEvRKN5Arith4vectIiLi3EEES9_EUlS9_E_EEvRKT_S9_S9_

////////////////////////////////////////////////////////////////////////////////

// TODO: Move this function to loop.hxx
template <typename F>
// inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_ATTRIBUTE_NOINLINE __attribute__((__flatten__, __hot__)) void
loop_region(const F &f, const Arith::vect<int, dim> &imin,
            const Arith::vect<int, dim> &imax) {
  assert(!any(imax <= imin));
  if (any(imax <= imin))
    return;

#if !defined __CUDACC__ && !defined __HIPCC__
  // CPU
  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      // #pragma omp simd
      for (int i = imin[0]; i < imax[0]; ++i) {
        const Arith::vect<int, dim> p{i, j, k};
        f(p);
      }
    }
  }
#else
  // GPU
  const amrex::Box box(amrex::IntVect(imin[0], imin[1], imin[2]),
                       amrex::IntVect(imax[0] - 1, imax[1] - 1, imax[2] - 1));
  amrex::ParallelFor(box, [=] CCTK_DEVICE(const int i, const int j, const int k)
                              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                const Arith::vect<int, dim> p{i, j, k};
                                f(p);
                              });
#endif
}

BoundaryCondition::BoundaryCondition(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const amrex::Box &box, amrex::FArrayBox &dest)
    : groupdata(groupdata), patchdata(ghext->patchdata.at(groupdata.patch)),
      geom(patchdata.amrcore->Geom(groupdata.level)),
      dest(dest), imin{geom.Domain().smallEnd(0), geom.Domain().smallEnd(1),
                       geom.Domain().smallEnd(2)},
      imax{geom.Domain().bigEnd(0) + 1 + !groupdata.indextype[0],
           geom.Domain().bigEnd(1) + 1 + !groupdata.indextype[1],
           geom.Domain().bigEnd(2) + 1 + !groupdata.indextype[2]},
      // Vertices are shifted outwards by 1/2 grid spacing
      xmin{geom.ProbLo(0) + CCTK_REAL(0.5) * groupdata.indextype[0] - 1,
           geom.ProbLo(1) + CCTK_REAL(0.5) * groupdata.indextype[1] - 1,
           geom.ProbLo(2) + CCTK_REAL(0.5) * groupdata.indextype[2] - 1},
      xmax{geom.ProbHi(0) - CCTK_REAL(0.5) * groupdata.indextype[0] + 1,
           geom.ProbHi(1) - CCTK_REAL(0.5) * groupdata.indextype[1] + 1,
           geom.ProbHi(2) - CCTK_REAL(0.5) * groupdata.indextype[2] + 1},
      dx{geom.CellSize(0), geom.CellSize(1), geom.CellSize(2)},
      amin{box.smallEnd(0), box.smallEnd(1), box.smallEnd(2)},
      amax{box.bigEnd(0) + 1, box.bigEnd(1) + 1, box.bigEnd(2) + 1},
      dmin{dest.box().smallEnd(0), dest.box().smallEnd(1),
           dest.box().smallEnd(2)},
      dmax{dest.box().bigEnd(0) + 1, dest.box().bigEnd(1) + 1,
           dest.box().bigEnd(2) + 1},
      layout(dmin, dmax), destptr(dest.dataPtr()) {
  // Check centering
  assert(box.ixType() == dest.box().ixType());
  for (int d = 0; d < dim; ++d)
    assert((box.type(d) == amrex::IndexType::CELL) == groupdata.indextype[d]);

  // Ensure we have enough ghost zones
  // TODO: Implement this
  // for (int d = 0; d < dim; ++d)
  //   assert(... <= groupdata.nghostzones[d]);

  // Ensure the destination array is large enough
  for (int d = 0; d < dim; ++d) {
    assert(dest.box().smallEnd(d) <= amin[d]);
    assert(dest.box().bigEnd(d) + 1 >= amax[d]);
  }
}

void BoundaryCondition::apply() const {
  constexpr int NEG = -1, INT = 0, POS = +1;

  apply_on_face<NEG, NEG, NEG>();
  apply_on_face<INT, NEG, NEG>();
  apply_on_face<POS, NEG, NEG>();
  apply_on_face<NEG, INT, NEG>();
  apply_on_face<INT, INT, NEG>();
  apply_on_face<POS, INT, NEG>();
  apply_on_face<NEG, POS, NEG>();
  apply_on_face<INT, POS, NEG>();
  apply_on_face<POS, POS, NEG>();

  apply_on_face<NEG, NEG, INT>();
  apply_on_face<INT, NEG, INT>();
  apply_on_face<POS, NEG, INT>();
  apply_on_face<NEG, INT, INT>();
  // apply_on_face<INT, INT, INT>();
  apply_on_face<POS, INT, INT>();
  apply_on_face<NEG, POS, INT>();
  apply_on_face<INT, POS, INT>();
  apply_on_face<POS, POS, INT>();

  apply_on_face<NEG, NEG, POS>();
  apply_on_face<INT, NEG, POS>();
  apply_on_face<POS, NEG, POS>();
  apply_on_face<NEG, INT, POS>();
  apply_on_face<INT, INT, POS>();
  apply_on_face<POS, INT, POS>();
  apply_on_face<NEG, POS, POS>();
  apply_on_face<INT, POS, POS>();
  apply_on_face<POS, POS, POS>();
}

template <int NI, int NJ, int NK>
void BoundaryCondition::apply_on_face() const {
  constexpr Arith::vect<int, dim> inormal{NI, NJ, NK};
  static_assert(!all(inormal == 0), "");

  using std::max, std::min;
  Arith::vect<int, dim> bmin, bmax;
  for (int d = 0; d < dim; ++d) {
    if (inormal[d] < 0) {
      // We have to fill the lower boundary
      bmin[d] = amin[d];
      bmax[d] = min(amax[d], imin[d]);
    } else if (inormal[d] > 0) {
      // We have to fill the upper boundary
      bmin[d] = max(amin[d], imax[d]);
      bmax[d] = amax[d];
    } else {
      // This direction is not a boundary
      bmin[d] = max(amin[d], imin[d]);
      bmax[d] = min(amax[d], imax[d]);
    }
  }

  // Do we actually own part of this boundary?
  bool npoints_is_zero = false;
  for (int d = 0; d < dim; ++d)
    npoints_is_zero |= bmax[d] <= bmin[d];
  if (npoints_is_zero)
    return;

  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners, multiple conditions will apply.
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
        (symmetry_x == symmetry_t::outer_boundary &&
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
      assert(0);
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI>
void BoundaryCondition::apply_on_face_symbcx(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners, multiple conditions will apply.
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
        (symmetry_y == symmetry_t::outer_boundary &&
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
      assert(0);
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
          symmetry_t SCJ, boundary_t BCJ>
void BoundaryCondition::apply_on_face_symbcxy(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  // Find which symmetry or boundary conditions apply to us. On edges
  // or in corners, multiple conditions will apply.
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
        (symmetry_z == symmetry_t::outer_boundary &&
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
      assert(0);
  }
}

template <int NI, int NJ, int NK, symmetry_t SCI, boundary_t BCI,
          symmetry_t SCJ, boundary_t BCJ, symmetry_t SCK, boundary_t BCK>
void BoundaryCondition::apply_on_face_symbcxyz(
    const Arith::vect<int, dim> &bmin,
    const Arith::vect<int, dim> &bmax) const {
  constexpr Arith::vect<int, dim> inormal{NI, NJ, NK};
  static_assert(!all(inormal == 0), "");
  constexpr Arith::vect<boundary_t, dim> boundaries{BCI, BCJ, BCK};
  constexpr Arith::vect<symmetry_t, dim> symmetries{SCI, SCJ, SCK};

  const int ncomps = dest.nComp();

  // Periodic symmetries should have been translated to `none`
  static_assert(!any(symmetries == symmetry_t::periodic), "");

  if constexpr (all(symmetries == symmetry_t::none &&
                    boundaries == boundary_t::none)) {
    // If there are no boundary conditions to apply, then do
    // nothing.

    // do nothing

  } else if constexpr (any(boundaries == boundary_t::dirichlet)) {
    // If we can apply Dirichlet boundary conditions, do so.
    // They are the easiest and cheapest.

    for (int comp = 0; comp < ncomps; ++comp) {
      const CCTK_REAL dirichlet_value = groupdata.dirichlet_values.at(comp);
      Loop::GF3D2<CCTK_REAL> var(layout, destptr + comp * layout.np);
      // #ifndef __CUDACC__
      // #pragma omp task final(true) untied
      // #endif
      loop_region(
          [
#ifdef CCTK_DEBUG
              amin = amin, amax = amax,
#endif
              dirichlet_value,
              var] CCTK_DEVICE(const Arith::vect<int, dim> &dst)
              CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
                for (int d = 0; d < dim; ++d)
                  assert(dst[d] >= amin[d] && dst[d] < amax[d]);
#endif
                var.store(dst, dirichlet_value);
              },
          bmin, bmax);
    }

  } else {
    // This is the generic case for applying boundary conditions.

    Arith::vect<int, dim> neumann_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::neumann) {
        if (inormal[d] != 0) {
          neumann_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
          if (inormal[d] < 0)
            assert(neumann_source[d] < amax[d]);
          else
            assert(neumann_source[d] >= amin[d]);
        }
      } else {
        neumann_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<int, dim> robin_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::robin) {
        assert(inormal[d] != 0);
        robin_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
        if (inormal[d] < 0)
          assert(robin_source[d] < amax[d]);
        else
          assert(robin_source[d] >= amin[d]);
      } else {
        robin_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<int, dim> linear_extrapolation_source;
    for (int d = 0; d < dim; ++d) {
      if (boundaries[d] == boundary_t::linear_extrapolation) {
        assert(inormal[d] != 0);
        linear_extrapolation_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
        if (inormal[d] < 0) {
          assert(linear_extrapolation_source[d] < amax[d]);
          assert(linear_extrapolation_source[d] - inormal[d] < amax[d]);
        } else {
          assert(linear_extrapolation_source[d] >= amin[d]);
          assert(linear_extrapolation_source[d] - inormal[d] >= amin[d]);
        }
      } else {
        linear_extrapolation_source[d] = 666666666; // go for a segfault
      }
    }

    Arith::vect<int, dim> reflection_offset;
    for (int d = 0; d < dim; ++d) {
      if (symmetries[d] == symmetry_t::reflection) {
        assert(inormal[d] != 0);
        reflection_offset[d] =
            inormal[d] < 0 ? 2 * imin[d] - groupdata.indextype.at(d)
                           : 2 * (imax[d] - 1) + groupdata.indextype.at(d);
        if (inormal[d] < 0)
          assert(reflection_offset[d] - bmin[d] < amax[d]);
        else
          assert(reflection_offset[d] - (bmax[d] - 1) >= amin[d]);
      } else {
        reflection_offset[d] = 666666666; // go for a segfault
      }
    }

    for (int comp = 0; comp < ncomps; ++comp) {
      const CCTK_REAL robin_value = groupdata.robin_values.at(comp);

      CCTK_REAL reflection_parity = +1;
      for (int d = 0; d < dim; ++d)
        if (symmetries[d] == symmetry_t::reflection)
          reflection_parity *= groupdata.parities.at(comp).at(d);
      using std::fabs;
      assert(fabs(reflection_parity) == 1);

      Loop::GF3D2<CCTK_REAL> var(layout, destptr + comp * layout.np);
      // #ifndef __CUDACC__
      // #pragma omp task final(true) untied
      // #endif
      loop_region(
          [
#ifdef CCTK_DEBUG
              amin = amin, amax = amax, dmin = dmin, dmax = dmax,
#endif
              neumann_source, linear_extrapolation_source, robin_source,
              robin_value, reflection_offset, reflection_parity, xmin = xmin,
              dx = dx, var] CCTK_DEVICE(const Arith::vect<int, dim>
                                            &dst) CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
              if (boundaries[d] == boundary_t::linear_extrapolation) {
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
            for (int d = 0; d < dim; ++d)
              assert(src[d] >= dmin[d] && src[d] < dmax[d]);
            assert(!all(src == dst));
#endif
            CCTK_REAL val = var(src);
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
            if constexpr (any(boundaries == boundary_t::linear_extrapolation)) {
              // Calculate gradient
              const CCTK_REAL grad = val - var(src + delta);
              using std::sqrt;
              val += sqrt(sum(pow2(dst - src)) / sum(pow2(delta))) * grad;
            }
#ifdef CCTK_DEBUG
            for (int d = 0; d < dim; ++d)
              assert(dst[d] >= amin[d] && dst[d] < amax[d]);
#endif
            if constexpr (any(symmetries == symmetry_t::reflection))
              val *= reflection_parity;
            var.store(dst, val);
          },
          bmin, bmax);
    }
  }
}

} // namespace CarpetX
