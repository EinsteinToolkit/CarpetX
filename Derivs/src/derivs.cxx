#include "derivs.hxx"

namespace Derivs {
using namespace Arith;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

// Tile-based multi-dimensional derivative operators

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const GF3D5<T> &gf, const GF3D5layout layout,
          const GridDescBaseDevice &grid, const GF3D2<const T> &gf0) {
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
calc_copy(const vec<GF3D5<T>, dim> &gf, const GF3D5layout layout,
          const GridDescBaseDevice &grid, const vec<GF3D2<const T>, dim> &gf0) {
  for (int a = 0; a < dim; ++a)
    calc_copy<CI, CJ, CK>(gf(a), layout, grid, gf0(a));
};

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_copy(const smat<GF3D5<T>, dim> &gf,
                                       const GF3D5layout layout,
                                       const GridDescBaseDevice &grid,
                                       const smat<GF3D2<const T>, dim> &gf0) {
  for (int a = 0; a < dim; ++a)
    for (int b = a; b < dim; ++b)
      calc_copy<CI, CJ, CK>(gf(a, b), layout, grid, gf0(a, b));
};

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
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 2 / 2);
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
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 4 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<4>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
        });
    break;

  case 6:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 6 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<6>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
        });
    break;

  case 8:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 8 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<8>(gf0, mask1, p.I, dx);
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
calc_derivs(const vec<GF3D5<T>, dim> &gf,
            const vec<vec<GF3D5<T>, dim>, dim> &dgf, const GF3D5layout layout,
            const GridDescBaseDevice &grid, const vec<GF3D2<const T>, dim> &gf0,
            const vect<T, dim> dx, const int deriv_order) {
  for (int a = 0; a < dim; ++a)
    calc_derivs<CI, CJ, CK>(gf(a), dgf(a), layout, grid, gf0(a), dx,
                            deriv_order);
}

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const smat<GF3D5<T>, dim> &gf,
            const smat<vec<GF3D5<T>, dim>, dim> &dgf, const GF3D5layout layout,
            const GridDescBaseDevice &grid,
            const smat<GF3D2<const T>, dim> &gf0, const vect<T, dim> dx,
            const int deriv_order) {
  for (int a = 0; a < dim; ++a)
    for (int b = a; b < dim; ++b)
      calc_derivs<CI, CJ, CK>(gf(a, b), dgf(a, b), layout, grid, gf0(a, b), dx,
                              deriv_order);
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
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 2 / 2);
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
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 4 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<4>(gf0, mask1, p.I, dx);
          const auto ddval = calc_deriv2<4>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
          ddgf.store(mask, index, ddval);
        });
    break;

  case 6:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 6 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<6>(gf0, mask1, p.I, dx);
          const auto ddval = calc_deriv2<6>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
          ddgf.store(mask, index, ddval);
        });
    break;

  case 8:
    grid.loop_int_device<CI, CJ, CK, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          // Take ghost points into account
          const vbool mask1 = mask_for_loop_tail<vbool>(p.i, p.imax + 8 / 2);
          const GF3D5index index(layout, p.I);
          const auto val = gf0(mask, p.I);
          const auto dval = calc_deriv<8>(gf0, mask1, p.I, dx);
          const auto ddval = calc_deriv2<8>(gf0, mask1, p.I, dx);
          gf.store(mask, index, val);
          dgf.store(mask, index, dval);
          ddgf.store(mask, index, ddval);
        });
    break;

  default:
    CCTK_VERROR("Unsupported derivative order %d", deriv_order);
  }
}

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const vec<GF3D5<T>, dim> &gf, const vec<vec<GF3D5<T>, dim>, dim> &dgf,
    const vec<smat<GF3D5<T>, dim>, dim> &ddgf, const GF3D5layout layout,
    const GridDescBaseDevice &grid, const vec<GF3D2<const T>, dim> &gf0,
    const vect<T, dim> dx, const int deriv_order) {
  for (int a = 0; a < dim; ++a)
    calc_derivs2<CI, CJ, CK>(gf(a), dgf(a), ddgf(a), layout, grid, gf0(a), dx,
                             deriv_order);
}

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const smat<GF3D5<T>, dim> &gf, const smat<vec<GF3D5<T>, dim>, dim> &dgf,
    const smat<smat<GF3D5<T>, dim>, dim> &ddgf, const GF3D5layout layout,
    const GridDescBaseDevice &grid, const smat<GF3D2<const T>, dim> &gf0,
    const vect<T, dim> dx, const int deriv_order) {
  for (int a = 0; a < dim; ++a)
    for (int b = a; b < dim; ++b)
      calc_derivs2<CI, CJ, CK>(gf(a, b), dgf(a, b), ddgf(a, b), layout, grid,
                               gf0(a, b), dx, deriv_order);
}

////////////////////////////////////////////////////////////////////////////////

// Template instantiations

using T = CCTK_REAL;

template void calc_copy<0, 0, 0>(const GF3D5<T> &gf, const GF3D5layout layout,
                                 const GridDescBaseDevice &grid,
                                 const GF3D2<const T> &gf0);

template void calc_copy<0, 0, 0>(const vec<GF3D5<T>, dim> &gf,
                                 const GF3D5layout layout,
                                 const GridDescBaseDevice &grid,
                                 const vec<GF3D2<const T>, dim> &gf0);

template void calc_copy<0, 0, 0>(const smat<GF3D5<T>, dim> &gf,
                                 const GF3D5layout layout,
                                 const GridDescBaseDevice &grid,
                                 const smat<GF3D2<const T>, dim> &gf0);
template void
calc_derivs<0, 0, 0>(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
                     const GF3D5layout layout, const GridDescBaseDevice &grid,
                     const GF3D2<const T> &gf0, const vect<T, dim> dx,
                     const int deriv_order);

template void calc_derivs<0, 0, 0>(const vec<GF3D5<T>, dim> &gf,
                                   const vec<vec<GF3D5<T>, dim>, dim> &dgf,
                                   const GF3D5layout layout,
                                   const GridDescBaseDevice &grid,
                                   const vec<GF3D2<const T>, dim> &gf0,
                                   const vect<T, dim> dx,
                                   const int deriv_order);

template void calc_derivs<0, 0, 0>(const smat<GF3D5<T>, dim> &gf,
                                   const smat<vec<GF3D5<T>, dim>, dim> &dgf,
                                   const GF3D5layout layout,
                                   const GridDescBaseDevice &grid,
                                   const smat<GF3D2<const T>, dim> &gf0,
                                   const vect<T, dim> dx,
                                   const int deriv_order);

template void
calc_derivs2<0, 0, 0>(const GF3D5<T> &gf, const vec<GF3D5<T>, dim> &dgf,
                      const smat<GF3D5<T>, dim> &ddgf, const GF3D5layout layout,
                      const GridDescBaseDevice &grid, const GF3D2<const T> &gf0,
                      const vect<T, dim> dx, const int deriv_order);

template void calc_derivs2<0, 0, 0>(
    const vec<GF3D5<T>, dim> &gf, const vec<vec<GF3D5<T>, dim>, dim> &dgf,
    const vec<smat<GF3D5<T>, dim>, dim> &ddgf, const GF3D5layout layout,
    const GridDescBaseDevice &grid, const vec<GF3D2<const T>, dim> &gf0,
    const vect<T, dim> dx, const int deriv_order);

template void calc_derivs2<0, 0, 0>(
    const smat<GF3D5<T>, dim> &gf, const smat<vec<GF3D5<T>, dim>, dim> &dgf,
    const smat<smat<GF3D5<T>, dim>, dim> &ddgf, const GF3D5layout layout,
    const GridDescBaseDevice &grid, const smat<GF3D2<const T>, dim> &gf0,
    const vect<T, dim> dx, const int deriv_order);

} // namespace Derivs
