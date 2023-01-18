#ifndef DERIVS_HXX
#define DERIVS_HXX

#include <div.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <ten3.hxx>
#include <vec.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Weyl {
using namespace Arith;
using namespace Loop;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

constexpr int deriv_order = 4;

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T pow3(const T &x) {
  return pow2(x) * x;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv1d(const simdl<T> &mask, const T *restrict const var, const ptrdiff_t di,
        const T dx) {
  const auto load = [&](const int n) {
    return maskz_loadu(mask, &var[n * di]) - maskz_loadu(mask, &var[-n * di]);
  };
  if constexpr (deriv_order == 2)
    return 1 / T(2) * load(1) / dx;
  if constexpr (deriv_order == 4)
    return (-1 / T(12) * load(2) + 2 / T(3) * load(1)) / dx;
  if constexpr (deriv_order == 6)
    return (1 / T(60) * load(3) - 3 / T(20) * load(2) + 3 / T(4) * load(1)) /
           dx;
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv2_1d(const simdl<T> &mask, const T *restrict const var, const ptrdiff_t di,
          const T dx) {
  const auto load = [&](const int n) {
    return maskz_loadu(mask, &var[n * di]) + maskz_loadu(mask, &var[-n * di]);
  };
  const auto load0 = [&]() { return maskz_loadu(mask, &var[0]); };
  if constexpr (deriv_order == 2)
    return (load(1) - 2 * load0()) / pow2(dx);
  if constexpr (deriv_order == 4)
    return (1 / T(12) * load(2) - 4 / T(3) * load(1) + 5 / T(2) * load0()) /
           pow2(dx);
  if constexpr (deriv_order == 6)
    return (1 / T(90) * load(3) - 3 / T(20) * load(2) + 3 / T(2) * load(1) -
            49 / T(18) * load0()) /
           pow2(dx);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv3_1d(const simdl<T> &mask, const T *restrict const var, const ptrdiff_t di,
          const T dx) {
  const auto load = [&](const int n) {
    return maskz_loadu(mask, &var[n * di]) - maskz_loadu(mask, &var[-n * di]);
  };
  if constexpr (deriv_order == 2)
    return (1 / T(2) * load(2) - load(1)) / pow3(dx);
  if constexpr (deriv_order == 4)
    return (-1 / T(8) * load(3) + load(2) - 13 / T(8) * load(1)) / pow3(dx);
  if constexpr (deriv_order == 6)
    return (7 / T(240) * load(4) - 3 / T(10) * load(3) +
            169 / T(129) * load(2) - 61 / T(30) * load(1)) /
           pow3(dx);
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv2_2d(const int vavail, const simdl<T> &mask, const T *restrict const var,
          const ptrdiff_t di, const ptrdiff_t dj, const T dx, const T dy) {
  constexpr size_t vsize = tuple_size_v<simd<T> >;
  if (di == 1) {
    assert(vavail > 0);
    constexpr int maxnpoints = deriv_order + 1 + vsize - 1;
    const int npoints = deriv_order + 1 + min(int(vsize), vavail) - 1;
    array<simd<T>, div_ceil(maxnpoints, int(vsize))> arrx;
    for (int i = 0; i < maxnpoints; i += vsize) {
      if (i < npoints) {
        const simdl<T> mask1 = mask_for_loop_tail<simdl<T> >(i, npoints);
        arrx[div_floor(i, int(vsize))] =
            deriv1d(mask1, &var[i - deriv_order / 2], dj, dy);
      }
    }
#ifdef CCTK_DEBUG
    for (int i = npoints; i < align_ceil(maxnpoints, int(vsize)); ++i)
      ((T *)&arrx[0])[i] = Arith::nan<T>()(); // unused
#endif
    const T *const varx = (T *)&arrx[0] + deriv_order / 2;
    return deriv1d(mask, varx, 1, dx);
  } else {
    assert(dj != 1);
    array<simd<T>, deriv_order + 1> arrx;
    for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
      if (j == 0) {
#ifdef CCTK_DEBUG
        arrx[deriv_order / 2 + j] = Arith::nan<simd<T> >()(); // unused
#endif
      } else {
        arrx[deriv_order / 2 + j] = deriv1d(mask, &var[j * dj], di, dx);
      }
    const T *const varx = (T *)(&arrx[deriv_order / 2]);
    return deriv1d(mask, varx, vsize, dy);
  }
}

////////////////////////////////////////////////////////////////////////////////

template <int dir, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST simd<T>
deriv(const simdl<T> &mask, const GF3D2<const T> &gf_, const vect<int, dim> &I,
      const vec<T, D> &dx) {
  static_assert(dir >= 0 && dir < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir));
  return deriv1d(mask, &gf_(I), di, dx(dir));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE
    ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 == dir2), simd<T> >
    deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  return deriv2_1d(mask, &gf_(I), di, dx(dir1));
}

template <int dir1, int dir2, typename T, int D>
inline ARITH_INLINE
    ARITH_DEVICE ARITH_HOST enable_if_t<(dir1 != dir2), simd<T> >
    deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  const ptrdiff_t dj = gf_.delta(DI(dir2));
  return deriv2_2d(vavail, mask, &gf_(I), di, dj, dx(dir1), dx(dir2));
}

template <int dir1, int dir2, int dir3, typename T, int D>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST
    enable_if_t<(dir1 == dir2 && dir1 == dir3), simd<T> >
    deriv3(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
           const vect<int, dim> &I, const vec<T, D> &dx) {
  static_assert(dir1 >= 0 && dir1 < D, "");
  static_assert(dir2 >= 0 && dir2 < D, "");
  static_assert(dir3 >= 0 && dir3 < D, "");
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.delta(DI(dir1));
  return deriv3_1d(mask, &gf_(I), di, dx(dir1));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<simd<T>, dim>
deriv(const simdl<T> &mask, const GF3D2<const T> &gf_, const vect<int, dim> &I,
      const vec<T, dim> &dx) {
  return {deriv<0>(mask, gf_, I, dx), deriv<1>(mask, gf_, I, dx),
          deriv<2>(mask, gf_, I, dx)};
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST smat<simd<T>, dim>
deriv2(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
       const vect<int, dim> &I, const vec<T, dim> &dx) {
  return {deriv2<0, 0>(vavail, mask, gf_, I, dx),
          deriv2<0, 1>(vavail, mask, gf_, I, dx),
          deriv2<0, 2>(vavail, mask, gf_, I, dx),
          deriv2<1, 1>(vavail, mask, gf_, I, dx),
          deriv2<1, 2>(vavail, mask, gf_, I, dx),
          deriv2<2, 2>(vavail, mask, gf_, I, dx)};
}

template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST sten3<simd<T>, dim>
deriv3(const int vavail, const simdl<T> &mask, const GF3D2<const T> &gf_,
       const vect<int, dim> &I, const vec<T, dim> &dx) {
  return {deriv3<0, 0, 0>(vavail, mask, gf_, I, dx),
          deriv3<0, 0, 1>(vavail, mask, gf_, I, dx),
          deriv3<0, 0, 2>(vavail, mask, gf_, I, dx),
          deriv3<0, 1, 1>(vavail, mask, gf_, I, dx),
          deriv3<0, 1, 2>(vavail, mask, gf_, I, dx),
          deriv3<0, 2, 2>(vavail, mask, gf_, I, dx),
          deriv3<1, 1, 1>(vavail, mask, gf_, I, dx),
          deriv3<1, 1, 2>(vavail, mask, gf_, I, dx),
          deriv3<1, 2, 2>(vavail, mask, gf_, I, dx),
          deriv3<2, 2, 2>(vavail, mask, gf_, I, dx)};
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
          const GF3D5<T> &gf0, const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const vec<CCTK_REAL, dim> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const auto val = gf1(mask, p.I);
        gf0.store(mask, index0, val);
      });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
            const GF3D5<T> &gf0, const vec<GF3D5<T>, dim> &dgf0,
            const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const vec<CCTK_REAL, dim> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const auto val = gf1(mask, p.I);
        gf0.store(mask, index0, val);
        const auto dval = deriv(mask, gf1, p.I, dx);
        dgf0.store(mask, index0, dval);
      });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
             const GF3D5<T> &gf0, const vec<GF3D5<T>, dim> &dgf0,
             const smat<GF3D5<T>, dim> &ddgf0, const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const vec<CCTK_REAL, dim> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const auto val = gf1(mask, p.I);
        gf0.store(mask, index0, val);
        const auto dval = deriv(mask, gf1, p.I, dx);
        dgf0.store(mask, index0, dval);
        const auto ddval = deriv2(p.imax - p.i, mask, gf1, p.I, dx);
        ddgf0.store(mask, index0, ddval);
      });
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_copy(const cGH *restrict const cctkGH,
                                       const vec<GF3D2<const T>, dim> &gf0_,
                                       const vec<GF3D5<T>, dim> &gf_,
                                       const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_copy(cctkGH, gf0_(a), gf_(a), layout);
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_derivs(
    const cGH *restrict const cctkGH, const vec<GF3D2<const T>, dim> &gf0_,
    const vec<GF3D5<T>, dim> &gf_, const vec<vec<GF3D5<T>, dim>, dim> &dgf_,
    const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs(cctkGH, gf0_(a), gf_(a), dgf_(a), layout);
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_derivs2(
    const cGH *restrict const cctkGH, const vec<GF3D2<const T>, dim> &gf0_,
    const vec<GF3D5<T>, dim> &gf_, const vec<vec<GF3D5<T>, dim>, dim> &dgf_,
    const vec<smat<GF3D5<T>, dim>, dim> &ddgf_, const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs2(cctkGH, gf0_(a), gf_(a), dgf_(a), ddgf_(a), layout);
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_copy(const cGH *restrict const cctkGH,
                                       const smat<GF3D2<const T>, dim> &gf0_,
                                       const smat<GF3D5<T>, dim> &gf_,
                                       const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_copy(cctkGH, gf0_(a, b), gf_(a, b), layout);
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_derivs(
    const cGH *restrict const cctkGH, const smat<GF3D2<const T>, dim> &gf0_,
    const smat<GF3D5<T>, dim> &gf_, const smat<vec<GF3D5<T>, dim>, dim> &dgf_,
    const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), layout);
}

template <typename T>
void CCTK_ATTRIBUTE_NOINLINE calc_derivs2(
    const cGH *restrict const cctkGH, const smat<GF3D2<const T>, dim> &gf0_,
    const smat<GF3D5<T>, dim> &gf_, const smat<vec<GF3D5<T>, dim>, dim> &dgf_,
    const smat<smat<GF3D5<T>, dim>, dim> &ddgf_, const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs2(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), ddgf_(a, b),
                   layout);
}

} // namespace Weyl

#endif // #ifndef DERIVS_HXX
