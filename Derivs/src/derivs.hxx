#ifndef CARPETX_DERIVS_DERIVS_HXX
#define CARPETX_DERIVS_DERIVS_HXX

#include <defs.hxx>
#include <div.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace Derivs {

////////////////////////////////////////////////////////////////////////////////

namespace stencils {
using namespace Arith;
using namespace Loop;

// Stencil coefficients

enum symmetry { none, symmetric, antisymmetric };

template <std::ptrdiff_t I0, std::ptrdiff_t I1, symmetry S> struct stencil {
  static constexpr std::ptrdiff_t N = I1 - I0 + 1;
  static_assert(N >= 0);
  static_assert(S == none || S == symmetric || S == antisymmetric);

  int divisor;
  std::array<int, N> coeffs;

  template <typename Array>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
      CCTK_DEVICE CCTK_HOST std::result_of_t<Array(std::ptrdiff_t)>
      apply(const Array &arr) const {
    using R = std::result_of_t<Array(std::ptrdiff_t)>;
    if constexpr (S == symmetric) {
      R r{0};
      for (std::ptrdiff_t n = 0; n < N / 2; ++n) {
        const std::ptrdiff_t n1 = N - 1 - n;
        r += coeffs[n] * (arr(n + I0) + arr(n1 + I0));
      }
      if (N % 2 != 0) {
        const std::ptrdiff_t n = N / 2 + 1;
        r += coeffs[n] * arr(n + I0);
      }
      r /= divisor;
      return r;
    }
    if constexpr (antisymmetric) {
      R r{0};
      for (std::ptrdiff_t n = 0; n < N / 2; ++n) {
        const std::ptrdiff_t n1 = N - 1 - n;
        r += coeffs[n] * (arr(n + I0) - arr(n1 + I0));
      }
      r /= divisor;
      return r;
    }
    R r{0};
    for (std::ptrdiff_t n = 0; n < N; ++n)
      r += coeffs[n] * arr(n + I0);
    return r;
  }
};

// Interpolate at i = 0
constexpr stencil<0, 0, symmetric> interp{1, {1}};

// Derivative at i = 0
constexpr stencil<-1, +1, antisymmetric> deriv1_o2{2, {-1, 0, +1}};

constexpr stencil<-2, +2, antisymmetric> deriv1_o4{12, {-1, +8, 0, -8, +1}};

constexpr stencil<-1, +1, symmetric> deriv2_o2{1, {-2, +1, -2}};

constexpr stencil<-2, +2, symmetric> deriv2_o4{12, {-1, +16, -30, +16, -1}};

// Interpolate at i = 1/2
constexpr stencil<-0, +1, symmetric> interp_c_o1{2, {1, 1}};

constexpr stencil<-1, +2, symmetric> interp_c_o3{16, {-1, +9, +9, -1}};

constexpr stencil<-2, +3, symmetric> interp_c_o5{
    256, {+3, -25, +150, +150, -25, +3}};

constexpr stencil<-3, +4, symmetric> interp_c_o7{
    2048, {-5, +49, -245, +1225, +1225, -245, +49, -5}};

// Derivative at i = 1/2
constexpr stencil<-0, +1, antisymmetric> deriv1_c_o2{1, {-1, +1}};

constexpr stencil<-1, +2, antisymmetric> deriv1_c_o4{12, {-1, +15, -15, +1}};

} // namespace stencils

namespace detail {
using namespace Arith;
using namespace Loop;

// Pointwise one-dimensional operators

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST simd<T>
interp1d(const simdl<T> &mask, const T *restrict const var) {
  return maskz_loadu(mask, var);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 1 / T(2);
  return (c1 * (var(1) - var(-1))) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 2 / T(3);
  constexpr T c2 = -1 / T(12);
  return (c2 * (var(2) - var(-2)) + c1 * (var(1) - var(-1))) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 6, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 3 / T(4);
  constexpr T c2 = -3 / T(20);
  constexpr T c3 = 1 / T(60);
  return (c3 * (var(3) - var(-3)) + c2 * (var(2) - var(-2)) +
          c1 * (var(1) - var(-1))) /
         dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 8, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 4 / T(5);
  constexpr T c2 = -1 / T(5);
  constexpr T c3 = 4 / T(105);
  constexpr T c4 = -1 / T(280);
  return (c4 * (var(4) - var(-4)) + c3 * (var(3) - var(-3)) +
          c2 * (var(2) - var(-2)) + c1 * (var(1) - var(-1))) /
         dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv1d_upwind(const TS var, const T dx, const R vel) {
  // arXiv:1111.2177 [gr-qc], (71)
  const T c1s = 1;
  const T c2s = -1 / T(4);
  const R symm = c2s * (var(2) - var(-2)) + c1s * (var(1) - var(-1));

  const T c0a = 3 / T(2);
  const T c1a = -1;
  const T c2a = 1 / T(4);
  const R anti =
      c2a * (var(2) + var(-2)) + c1a * (var(1) + var(-1)) + c0a * var(0);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv1d_upwind(const TS var, const T dx, const R vel) {
  // arXiv:1111.2177 [gr-qc], (71)
  const T c1s = 7 / T(8);
  const T c2s = -1 / T(4);
  const T c3s = 1 / T(24);
  const R symm = c3s * (var(3) - var(-3)) + c2s * (var(2) - var(-2)) +
                 c1s * (var(1) - var(-1));

  const T c0a = 5 / T(6);
  const T c1a = -5 / T(8);
  const T c2a = 1 / T(4);
  const T c3a = -1 / T(24);
  const R anti = c3a * (var(3) + var(-3)) + c2a * (var(2) + var(-2)) +
                 c1a * (var(1) + var(-1)) + c0a * var(0);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 6, R>
    deriv1d_upwind(const TS var, const T dx, const R vel) {
  // arXiv:1111.2177 [gr-qc], (71)
  const T c1s = 13 / T(15);
  const T c2s = -4 / T(15);
  const T c3s = 1 / T(15);
  const T c4s = -1 / T(120);
  const R symm = c4s * (var(4) - var(-4)) + c3s * (var(3) - var(-3)) +
                 c2s * (var(2) - var(-2)) + c1s * (var(1) - var(-1));

  const T c0a = 7 / T(12);
  const T c1a = -7 / T(15);
  const T c2a = 7 / T(30);
  const T c3a = -1 / T(15);
  const T c4a = 1 / T(120);
  const R anti = c4a * (var(4) + var(-4)) + c3a * (var(3) + var(-3)) +
                 c2a * (var(2) + var(-2)) + c1a * (var(1) + var(-1)) +
                 c0a * var(0);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 8, R>
    deriv1d_upwind(const TS var, const T dx, const R vel) {
  // arXiv:1111.2177 [gr-qc], (71)
  const T c1s = 7 / T(8);
  const T c2s = -2 / T(7);
  const T c3s = 29 / T(336);
  const T c4s = -1 / T(56);
  const T c5s = 1 / T(560);
  const R symm = c5s * (var(5) - var(-5)) + c4s * (var(4) - var(-4)) +
                 c3s * (var(3) - var(-3)) + c2s * (var(2) - var(-2)) +
                 c1s * (var(1) - var(-1));

  const T c0a = 9 / T(20);
  const T c1a = -3 / T(8);
  const T c2a = 3 / T(14);
  const T c3a = -9 / T(112);
  const T c4a = 1 / T(56);
  const T c5a = -1 / T(560);
  const R anti = c5a * (var(5) + var(-5)) + c4a * (var(4) + var(-4)) +
                 c3a * (var(3) + var(-3)) + c2a * (var(2) + var(-2)) +
                 c1a * (var(1) + var(-1)) + c0a * var(0);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv2_1d(const TS var, const T dx) {
  const T c1 = 1;
  return (c1 * ((var(1) - var(0)) - (var(0) - var(-1)))) / pow2(dx);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv2_1d(const TS var, const T dx) {
  constexpr T c1 = 4 / T(3);
  constexpr T c2 = -1 / T(12);
  return (c2 * ((var(2) - var(0)) - (var(0) - var(-2))) +
          c1 * ((var(1) - var(0)) - (var(0) - var(-1)))) /
         pow2(dx);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 6, R>
    deriv2_1d(const TS var, const T dx) {
  constexpr T c1 = 3 / T(2);
  constexpr T c2 = -3 / T(20);
  constexpr T c3 = 1 / T(90);
  return (c3 * ((var(3) - var(0)) - (var(0) - var(-3))) +
          c2 * ((var(2) - var(0)) - (var(0) - var(-2))) +
          c1 * ((var(1) - var(0)) - (var(0) - var(-1)))) /
         pow2(dx);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 8, R>
    deriv2_1d(const TS var, const T dx) {
  constexpr T c1 = 8 / T(5);
  constexpr T c2 = -1 / T(5);
  constexpr T c3 = 8 / T(315);
  constexpr T c4 = -1 / T(560);
  return (c4 * ((var(4) - var(0)) - (var(0) - var(-4))) +
          c3 * ((var(3) - var(0)) - (var(0) - var(-3))) +
          c2 * ((var(2) - var(0)) - (var(0) - var(-2))) +
          c1 * ((var(1) - var(0)) - (var(0) - var(-1)))) /
         pow2(dx);
}

template <int deriv_order, bool vectorize_di, typename T, typename TS,
          typename R = std::result_of_t<TS(int, int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
deriv2_2d(const TS var, const T dx, const T dy) {
  // We assume that the x-direction might be special since it might
  // be SIMD-vectorized. We assume that the y-direction is not
  // SIMD-vectorized.

  if constexpr (vectorize_di) {
    // Calculate y-derivative first
    static_assert(sizeof(R) % sizeof(T) == 0);
    static_assert(sizeof(R) / sizeof(T) > 0);
    constexpr std::ptrdiff_t vsize = sizeof(R) / sizeof(T);

    // We need fewer ndyvars than without vectorizon: Instead of `(2 *
    // deriv_order + 1) * vsize` scalars, we only need to calculate
    // `(2 * deriv_order + 1) + (vsize - 1)` scalars
    constexpr std::ptrdiff_t maxnpoints = deriv_order + 1 + vsize - 1;
    constexpr std::ptrdiff_t ndyvars = div_ceil(maxnpoints, vsize);
    std::array<R, ndyvars> dyvar;

    for (std::ptrdiff_t n = 0; n < maxnpoints; n += vsize) {
      const std::ptrdiff_t di = n - deriv_order / 2;
      // Skip the unused central point, but only if there is no vectorization
      if (vsize == 1 && di == 0)
        continue;
      dyvar[div_floor(n, vsize)] = deriv1d<deriv_order>(
          [&](int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
            assert(di >= -deriv_order / 2);
            assert(di <= +deriv_order / 2);
            assert(di >= -deriv_order / 2);
            assert(dj <= +deriv_order / 2);
#endif
            return var(di, dj);
          },
          dy);
    }

    // Calculate x-derivative next
    const T *const scalar_dyvar = (const T *)dyvar.data();
    return deriv1d<deriv_order>(
        [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
          assert(di >= -deriv_order / 2);
          assert(di <= +deriv_order / 2);
#endif
          if constexpr (vsize == 1)
            return scalar_dyvar[deriv_order / 2 + di];
          else
            return loadu<R>(&scalar_dyvar[deriv_order / 2 + di]);
        },
        dx);

  } else {
    // Calculate y-derivative first
    constexpr std::ptrdiff_t ndyvars = deriv_order + 1;
    std::array<R, ndyvars> dyvar;
#ifdef CCTK_DEBUG
    for (std::ptrdiff_t n = 0; n < ndyvars; ++n)
      dyvar[n] = Arith::nan<T>()();
#endif

    for (std::ptrdiff_t n = 0; n < ndyvars; ++n) {
      const std::ptrdiff_t di = n - deriv_order / 2;
      // Skip the unused central point
      if (di == 0)
        continue;
      dyvar[n] = deriv1d<deriv_order>(
          [&](int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
            assert(di >= -deriv_order / 2);
            assert(di <= +deriv_order / 2);
            assert(di >= -deriv_order / 2);
            assert(dj <= +deriv_order / 2);
#endif
            return var(di, dj);
          },
          dy);
    }

    // Calculate x-derivative next
    return deriv1d<deriv_order>(
        [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
#ifdef CCTK_DEBUG
          assert(di >= -deriv_order / 2);
          assert(di <= +deriv_order / 2);
#endif
          return dyvar[deriv_order / 2 + di];
        },
        dx);
  }
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    diss1d(const TS var, const T dx) {
  const T c0 = 6;
  const T c1 = -4;
  const T c2 = 1;
  return (c2 * (var(2) + var(-2)) + c1 * (var(1) + var(-1)) + c0 * var(0)) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    diss1d(const TS var, const T dx) {
  const T c0 = -20;
  const T c1 = 15;
  const T c2 = -6;
  const T c3 = 1;
  return (c3 * (var(3) + var(-3)) + c2 * (var(2) + var(-2)) +
          c1 * (var(1) + var(-1)) + c0 * var(0)) /
         dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 6, R>
    diss1d(const TS var, const T dx) {
  const T c0 = 70;
  const T c1 = -56;
  const T c2 = 28;
  const T c3 = -8;
  const T c4 = 1;
  return (c4 * (var(4) + var(-4)) + c3 * (var(3) + var(-3)) +
          c2 * (var(2) + var(-2)) + c1 * (var(1) + var(-1)) + c0 * var(0)) /
         dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 8, R>
    diss1d(const TS var, const T dx) {
  const T c0 = -252;
  const T c1 = 210;
  const T c2 = -120;
  const T c3 = 45;
  const T c4 = -10;
  const T c5 = 1;
  return (c5 * (var(5) + var(-5)) + c4 * (var(4) + var(-4)) +
          c3 * (var(3) + var(-3)) + c2 * (var(2) + var(-2)) +
          c1 * (var(1) + var(-1)) + c0 * var(0)) /
         dx;
}

} // namespace detail

////////////////////////////////////////////////////////////////////////////////

// Pointwise multi-dimensional derivative operators

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST Arith::simd<T>
calc_deriv_upwind(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
                  const Arith::vect<int, Loop::dim> &I,
                  const Arith::vect<T, Loop::dim> &dx,
                  const Arith::vec<Arith::simd<T>, Loop::dim> &vel) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return detail::deriv1d_upwind<deriv_order>(
             [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
               return maskz_loadu(mask, &ptr[di * offsets[0]]);
             },
             dx[0], vel(0)) +
         detail::deriv1d_upwind<deriv_order>(
             [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
               return maskz_loadu(mask, &ptr[di * offsets[1]]);
             },
             dx[1], vel(1)) +
         detail::deriv1d_upwind<deriv_order>(
             [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
               return maskz_loadu(mask, &ptr[di * offsets[2]]);
             },
             dx[2], vel(2));
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T calc_deriv_upwind(
    const Loop::GF3D2<const T> &gf, const Arith::vect<int, Loop::dim> &I,
    const Arith::vect<T, Loop::dim> &dx, const Arith::vec<T, Loop::dim> &vel) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return detail::deriv1d_upwind<deriv_order>(
             [&](int di)
                 CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
             dx[0], vel(0)) +
         detail::deriv1d_upwind<deriv_order>(
             [&](int di)
                 CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
             dx[1], vel(1)) +
         detail::deriv1d_upwind<deriv_order>(
             [&](int di)
                 CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
             dx[2], vel(2));
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST Arith::simd<T>
calc_diss(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
          const Arith::vect<int, Loop::dim> &I,
          const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  constexpr int diss_order = deriv_order + 2;
  constexpr int sign = diss_order % 4 == 0 ? -1 : +1;
  return sign / T(pown(2, deriv_order + 2)) *
         (detail::diss1d<deriv_order>(
              [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
                return maskz_loadu(mask, &ptr[di * offsets[0]]);
              },
              dx[0]) +
          detail::diss1d<deriv_order>(
              [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
                return maskz_loadu(mask, &ptr[di * offsets[1]]);
              },
              dx[1]) +
          detail::diss1d<deriv_order>(
              [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
                return maskz_loadu(mask, &ptr[di * offsets[2]]);
              },
              dx[2]));
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
calc_diss(const Loop::GF3D2<const T> &gf, const Arith::vect<int, Loop::dim> &I,
          const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  constexpr int diss_order = deriv_order + 2;
  constexpr int sign = diss_order % 4 == 0 ? -1 : +1;
  return sign / T(pown(2, deriv_order + 2)) *
         (detail::diss1d<deriv_order>(
              [&](int di)
                  CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
              dx[0]) +
          detail::diss1d<deriv_order>(
              [&](int di)
                  CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
              dx[1]) +
          detail::diss1d<deriv_order>(
              [&](int di)
                  CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
              dx[2]));
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::vec<Arith::simd<T>, Loop::dim>
    calc_deriv(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return {
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0]]);
          },
          dx[0]),
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1]]);
          },
          dx[1]),
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[2]]);
          },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::vec<T, Loop::dim>
    calc_deriv(const Loop::GF3D2<const T> &gf,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return {
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
          dx[0]),
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
          dx[1]),
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::smat<Arith::simd<T>, Loop::dim>
    calc_deriv2(const Loop::GF3D2<const T> &gf, const Arith::simdl<T> &mask,
                const Arith::vect<int, Loop::dim> &I,
                const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return {
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0]]);
          },
          dx[0]),
      detail::deriv2_2d<deriv_order, true>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0] + dj * offsets[1]]);
          },
          dx[0], dx[1]),
      detail::deriv2_2d<deriv_order, true>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0] + dj * offsets[2]]);
          },
          dx[0], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1]]);
          },
          dx[1]),
      detail::deriv2_2d<deriv_order, false>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1] + dj * offsets[2]]);
          },
          dx[1], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[2]]);
          },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::smat<T, Loop::dim>
    calc_deriv2(const Loop::GF3D2<const T> &gf,
                const Arith::vect<int, Loop::dim> &I,
                const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      gf.delta(1, 0, 0),
      gf.delta(0, 1, 0),
      gf.delta(0, 0, 1),
  };
  return {
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
          dx[0]),
      detail::deriv2_2d<deriv_order, true>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[0] + dj * offsets[1]];
          },
          dx[0], dx[1]),
      detail::deriv2_2d<deriv_order, true>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[0] + dj * offsets[2]];
          },
          dx[0], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
          dx[1]),
      detail::deriv2_2d<deriv_order, false>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[1] + dj * offsets[2]];
          },
          dx[1], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
          dx[2]),
  };
}

////////////////////////////////////////////////////////////////////////////////

// Tile-based multi-dimensional operators

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_copy(const Loop::GF3D5<T> &gf,
                                       const Loop::GF3D5layout layout,
                                       const Loop::GridDescBaseDevice &grid,
                                       const Loop::GF3D2<const T> &gf0);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const Arith::vec<Loop::GF3D5<T>, Loop::dim> &gf,
          const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
          const Arith::vec<Loop::GF3D2<const T>, Loop::dim> &gf0);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const Arith::smat<Loop::GF3D5<T>, Loop::dim> &gf,
          const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
          const Arith::smat<Loop::GF3D2<const T>, Loop::dim> &gf0);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs(
    const Loop::GF3D5<T> &gf, const Arith::vec<Loop::GF3D5<T>, Loop::dim> &dgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Loop::GF3D2<const T> &gf0, const Arith::vect<T, Loop::dim> dx,
    const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs(
    const Arith::vec<Loop::GF3D5<T>, Loop::dim> &gf,
    const Arith::vec<Arith::vec<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &dgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Arith::vec<Loop::GF3D2<const T>, Loop::dim> &gf0,
    const Arith::vect<T, Loop::dim> dx, const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs(
    const Arith::smat<Loop::GF3D5<T>, Loop::dim> &gf,
    const Arith::smat<Arith::vec<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &dgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Arith::smat<Loop::GF3D2<const T>, Loop::dim> &gf0,
    const Arith::vect<T, Loop::dim> dx, const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const Loop::GF3D5<T> &gf, const Arith::vec<Loop::GF3D5<T>, Loop::dim> &dgf,
    const Arith::smat<Loop::GF3D5<T>, Loop::dim> &ddgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Loop::GF3D2<const T> &gf0, const Arith::vect<T, Loop::dim> dx,
    const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const Arith::vec<Loop::GF3D5<T>, Loop::dim> &gf,
    const Arith::vec<Arith::vec<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &dgf,
    const Arith::vec<Arith::smat<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &ddgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Arith::vec<Loop::GF3D2<const T>, Loop::dim> &gf0,
    const Arith::vect<T, Loop::dim> dx, const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const Arith::smat<Loop::GF3D5<T>, Loop::dim> &gf,
    const Arith::smat<Arith::vec<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &dgf,
    const Arith::smat<Arith::smat<Loop::GF3D5<T>, Loop::dim>, Loop::dim> &ddgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Arith::smat<Loop::GF3D2<const T>, Loop::dim> &gf0,
    const Arith::vect<T, Loop::dim> dx, const int deriv_order);

} // namespace Derivs

#endif // #ifndef CARPETX_DERIVS_DERIVS_HXX
