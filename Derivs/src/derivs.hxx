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
      return std::move(r);
    }
    if constexpr (antisymmetric) {
      R r{0};
      for (std::ptrdiff_t n = 0; n < N / 2; ++n) {
        const std::ptrdiff_t n1 = N - 1 - n;
        r += coeffs[n] * (arr(n + I0) - arr(n1 + I0));
      }
      r /= divisor;
      return std::move(r);
    }
    R r{0};
    for (std::ptrdiff_t n = 0; n < N; ++n)
      r += coeffs[n] * arr(n + I0);
    return std::move(r);
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
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv1d(const TS var, const T dx) {
  const T c1 = 1 / (2 * dx);
  return c1 * (var(1) - var(-1));
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv1d(const TS var, const T dx) {
  const T c1 = 2 / (3 * dx);
  const T c2 = -1 / (12 * dx);
  return c2 * (var(2) - var(-2)) + c1 * (var(1) - var(-1));
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, simd<T> >
    deriv1d_upwind(const simdl<T> &mask, const T *restrict const var,
                   const std::ptrdiff_t di, const simd<T> &vel, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)

  // if (sign)
  //   // +     [ 0   -1   +1    0    0]
  //   // + 1/2 [+1   -2   +1    0    0]
  //   //       [+1/2 -2   +3/2  0    0]
  //   return (1 / T(2) * var[-2 * di] //
  //           - 2 * var[-di]          //
  //           + 3 / T(2) * var[0]) /
  //          dx;
  // else
  //   // +     [ 0    0   -1   +1    0  ]
  //   // - 1/2 [ 0    0   +1   -2   +1  ]
  //   //       [ 0    0   -3/2 +2   -1/2]
  //   return (-3 / T(2) * var[0] //
  //           + 2 * var[+di]     //
  //           - 1 / T(2) * var[+2 * di]) /
  //          dx;

  // + 1/2 [+1/2 -2   +3/2  0    0  ]
  // + 1/2 [ 0    0   -3/2 +2   -1/2]
  //       [+1/4 -1    0   +1   -1/4]
  constexpr T c1s = 1;
  constexpr T c2s = -1 / T(4);
  const simd<T> symm =
      c2s * (maskz_loadu(mask, &var[2 * di]) -
             maskz_loadu(mask, &var[-2 * di])) //
      + c1s(maskz_loadu(mask, &var[di]) - maskz_loadu(mask, &var[-di]));
  // + 1/2 [+1/2 -2   +3/2  0    0  ]
  // - 1/2 [ 0    0   -3/2 +2   -1/2]
  //       [+1/4 -1   +3/2 -1   +1/4]
  constexpr T c0a = 3 / T(2);
  constexpr T c1a = -1;
  constexpr T c2a = 1 / T(4);
  const simd<T> anti =
      c2a * (maskz_loadu(mask, &var[2 * di]) +
             maskz_loadu(mask, &var[-2 * di]))                          //
      + c1a(maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
      + c0a * maskz_loadu(mask, &var[0]);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, simd<T> >
    deriv1d_upwind(const simdl<T> &mask, const T *restrict const var,
                   const std::ptrdiff_t di, const simd<T> &vel, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)

  // A fourth order stencil for a first derivative, shifted by one grid point

  // if (sign)
  //   return (-1 / T(12) * var[-3 * di] //
  //           + 1 / T(2) * var[-2 * di] //
  //           - 3 / T(2) * var[-di]     //
  //           + 5 / T(6) * var[0]       //
  //           + 1 / T(4) * var[+di]) /
  //          dx;
  // else
  //   return (-1 / T(4) * var[-di]      //
  //           - 5 / T(6) * var[0]       //
  //           + 3 / T(2) * var[+di]     //
  //           - 1 / T(2) * var[+2 * di] //
  //           + 1 / T(12) * var[+3 * di]) /
  //          dx;

  constexpr T c1s = 7 / T(8);
  constexpr T c2s = -1 / T(4);
  constexpr T c3s = 1 / T(24);
  const simd<T> symm =
      c3s * (maskz_loadu(mask, &var[3 * di]) -
             maskz_loadu(mask, &var[-3 * di])) //
      + c2s * (maskz_loadu(mask, &var[2 * di]) -
               maskz_loadu(mask, &var[-2 * di])) //
      + c1s * (maskz_loadu(mask, &var[di]) - maskz_loadu(mask, &var[-di]));
  constexpr T c0a = 5 / T(6);
  constexpr T c1a = -5 / T(8);
  constexpr T c2a = 1 / T(4);
  constexpr T c3a = -1 / T(24);
  const simd<T> anti =
      c3a * (maskz_loadu(mask, &var[3 * di]) +
             maskz_loadu(mask, &var[-3 * di])) //
      + c2a * (maskz_loadu(mask, &var[2 * di]) +
               maskz_loadu(mask, &var[-2 * di]))                           //
      + c1a * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
      + c0a * maskz_loadu(mask, &var[0]);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv2_1d(const TS var, const T dx) {
  // constexpr T c0 = -2 / pow2(dx);
  // constexpr T c1 = 1 / pow2(dx);
  // return c1 * (var(-1) + var(1)) + c0 * var(0);
  const T c0 = 1 / pow2(dx);
  return c0 * ((var(1) - var(0)) - (var(0) - var(-1)));
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv2_1d(const TS var, const T dx) {
  // constexpr T c0 = -5 / T(2);
  // constexpr T c1 = 4 / T(3);
  // constexpr T c2 = -1 / T(12);
  // return (c2 * (var(-2) + var(2)) + c1 * (var(-1) + var(1)) + c0 * var(0)) /
  //        pow2(dx);
  const T c0 = 15 / (12 * pow2(dx));
  const T c1 = -1 / (12 * pow2(dx));
  return c1 * ((var(4) - var(1)) - (var(3) - var(0))) +
         c0 * ((var(3) - var(2)) - (var(2) - var(1)));
}

template <int deriv_order, bool vectorize_di, typename T, typename TS,
          typename R = std::result_of_t<TS(int, int)> >
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
    constexpr std::ptrdiff_t ndyvars =
        div_ceil(2 * deriv_order + 1 + vsize - 1, vsize);
    std::array<R, ndyvars> dyvar;

    for (std::ptrdiff_t n = 0; n < ndyvars; ++n) {
      std::ptrdiff_t di = n - deriv_order;
      // Skip the unused central point, but only if there is no vectorization
      if (vsize == 1 && di == 0)
        continue;
      dyvar[n] = deriv1d<deriv_order>(
          [&](int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE { return var(di, dj); }, dy);
    }

    // Calculate x-derivative next
    const T *const scalar_dyvar = (const T *)dyvar.data();
    return deriv1d<deriv_order>(
        [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if constexpr (vsize == 1)
            return scalar_dyvar[deriv_order + di];
          else
            return loadu<R>(&scalar_dyvar[deriv_order + di]);
        },
        dx);

  } else {
    // Calculate y-derivative first
    constexpr std::ptrdiff_t ndyvars = 2 * deriv_order + 1;
    std::array<R, ndyvars> dyvar;
#ifdef CCTK_DEBUG
    for (std::ptrdiff_t n = 0; n < ndyvars; ++n)
      dyvar[n] = Arith::nan<T>()();
#endif

    for (std::ptrdiff_t n = 0; n < ndyvars; ++n) {
      std::ptrdiff_t di = n - deriv_order;
      // Skip the unused central point
      if (di == 0)
        continue;
      dyvar[n] = deriv1d<deriv_order>(
          [&](int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE { return var(di, dj); }, dy);
    }

    // Calculate x-derivative next
    return deriv1d<deriv_order>(
        [&](int di)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { return dyvar[deriv_order + di]; },
        dx);
  }
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, simd<T> >
    diss1d(const simdl<T> &mask, const T *restrict const var,
           const std::ptrdiff_t di, const T dx) {
  constexpr T c0 = 6;
  constexpr T c1 = -4;
  constexpr T c2 = 1;
  return (c2 * (maskz_loadu(mask, &var[2 * di]) +
                maskz_loadu(mask, &var[-2 * di]))                             //
          + c1 * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
          + c0 * maskz_loadu(mask, &var[0])) /
         dx;
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, simd<T> >
    diss1d(const simdl<T> &mask, const T *restrict const var,
           const std::ptrdiff_t di, const T dx) {
  constexpr T c0 = -20;
  constexpr T c1 = 15;
  constexpr T c2 = -6;
  constexpr T c3 = 1;
  return (c3 * (maskz_loadu(mask, &var[3 * di]) +
                maskz_loadu(mask, &var[-3 * di])) //
          + c2 * *(maskz_loadu(mask, &var[2 * di]) +
                   maskz_loadu(mask, &var[-2 * di]))                          //
          + c1 * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
          + c0 * maskz_loadu(mask, &var[0])) /
         dx;
}

} // namespace detail

////////////////////////////////////////////////////////////////////////////////

// Pointwise multi-dimensional derivative operators

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
CCTK_ATTRIBUTE_NOINLINE void
calc_copy(const Loop::GF3D5<T> &gf, const Loop::GF3D5layout layout,
          const Loop::GridDescBaseDevice &grid, const Loop::GF3D2<const T> &gf0,
          const Arith::vect<T, Loop::dim> dx);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs(
    const Loop::GF3D5<T> &gf, const Arith::vec<Loop::GF3D5<T>, Loop::dim> &dgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Loop::GF3D2<const T> &gf0, const Arith::vect<T, Loop::dim> dx,
    const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const Loop::GF3D5<T> &gf, const Arith::vec<Loop::GF3D5<T>, Loop::dim> &dgf,
    const Arith::smat<Loop::GF3D5<T>, Loop::dim> &ddgf,
    const Loop::GF3D5layout layout, const Loop::GridDescBaseDevice &grid,
    const Loop::GF3D2<const T> &gf0, const Arith::vect<T, Loop::dim> dx,
    const int deriv_order);

} // namespace Derivs

#endif // #ifndef CARPETX_DERIVS_DERIVS_HXX
