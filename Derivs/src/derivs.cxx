#include "derivs.hxx"

namespace Derivs {
using namespace Arith;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

// Tile-based multi-dimensional derivative operators

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const GF3D5<T> &gf, const GF3D5layout layout,
          const GridDescBaseDevice &grid, const GF3D2<const T> &gf0,
          const vect<T, dim> dx) {
  using vreal = simd<T>;
  using vbool = simdl<T>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<CI, CJ, CK, vsize>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index(layout, p.I);
        const auto val = gf0(mask, p.I);
        gf.store(mask, index, val);
      });
}

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
            const GF3D5layout layout, const GridDescBaseDevice &grid,
            const GF3D2<const T> &gf0, const vect<T, dim> dx,
            const int deriv_order) {
  using vreal = simd<T>;
  using vbool = simdl<T>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  switch (deriv_order) {

  case 2:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take account of ghost points
          const vbool mask1 =
              mask_for_loop_tail<vbool>(p.i, p.imax + deriv_order / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<2>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
        });
    break;

  case 4:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take account of ghost points
          const vbool mask1 =
              mask_for_loop_tail<vbool>(p.i, p.imax + deriv_order / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<4>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
        });
    break;

  default:
    CCTK_VERROR("Unsupported derivative order %d", deriv_order);
  }
}

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
             const smat<GF3D5<T>, dim> &ddgf, const GF3D5layout layout,
             const GridDescBaseDevice &grid, const GF3D2<const T> &gf0,
             const vect<T, dim> dx, const int deriv_order) {
  using vreal = simd<T>;
  using vbool = simdl<T>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  switch (deriv_order) {

  case 2:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take account of ghost points
          const vbool mask1 =
              mask_for_loop_tail<vbool>(p.i, p.imax + deriv_order / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<2>(gf0, mask1, p.I, dx);
          const auto ddval = calc_deriv2<2>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
          ddgf.store(mask, index, ddval);
        });
    break;

  case 4:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take account of ghost points
          const vbool mask1 =
              mask_for_loop_tail<vbool>(p.i, p.imax + deriv_order / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<4>(gf0, mask1, p.I, dx);
          const auto ddval = calc_deriv2<4>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
          ddgf.store(mask, index, ddval);
        });
    break;

  default:
    CCTK_VERROR("Unsupported derivative order %d", deriv_order);
  }
}

////////////////////////////////////////////////////////////////////////////////

// Template instantiations

using T = CCTK_REAL;

template CCTK_DEVICE CCTK_HOST Arith::vec<Arith::simd<T>, Loop::dim>
calc_deriv<2>(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
              const Arith::vect<int, Loop::dim> &I,
              const Arith::vect<T, Loop::dim> &dx);
template CCTK_DEVICE CCTK_HOST Arith::vec<Arith::simd<T>, Loop::dim>
calc_deriv<4>(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
              const Arith::vect<int, Loop::dim> &I,
              const Arith::vect<T, Loop::dim> &dx);

template CCTK_DEVICE CCTK_HOST Arith::vec<T, Loop::dim>
calc_deriv<2>(const Loop::GF3D2<const T> &gf,
              const Arith::vect<int, Loop::dim> &I,
              const Arith::vect<T, Loop::dim> &dx);
template CCTK_DEVICE CCTK_HOST Arith::vec<T, Loop::dim>
calc_deriv<4>(const Loop::GF3D2<const T> &gf,
              const Arith::vect<int, Loop::dim> &I,
              const Arith::vect<T, Loop::dim> &dx);

template CCTK_DEVICE CCTK_HOST Arith::smat<Arith::simd<T>, Loop::dim>
calc_deriv2<2>(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx);
template CCTK_DEVICE CCTK_HOST Arith::smat<Arith::simd<T>, Loop::dim>
calc_deriv2<4>(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx);

template CCTK_DEVICE CCTK_HOST Arith::smat<T, Loop::dim>
calc_deriv2<2>(const Loop::GF3D2<const T> &gf,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx);
template CCTK_DEVICE CCTK_HOST Arith::smat<T, Loop::dim>
calc_deriv2<4>(const Loop::GF3D2<const T> &gf,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx);

template void calc_copy<0, 0, 0>(const GF3D5<T> &gf, const GF3D5layout layout,
                                 const GridDescBaseDevice &grid,
                                 const GF3D2<const T> &gf0,
                                 const vect<T, dim> dx);

template void
calc_derivs<0, 0, 0>(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
                     const GF3D5layout layout, const GridDescBaseDevice &grid,
                     const GF3D2<const T> &gf0, const vect<T, dim> dx,
                     const int deriv_order);

template void
calc_derivs2<0, 0, 0>(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
                      const smat<GF3D5<T>, dim> &ddgf, const GF3D5layout layout,
                      const GridDescBaseDevice &grid, const GF3D2<const T> &gf0,
                      const vect<T, dim> dx, const int deriv_order);

} // namespace Derivs
