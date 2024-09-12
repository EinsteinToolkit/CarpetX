#ifndef CARPETX_ARITH_SIMD_HXX
#define CARPETX_ARITH_SIMD_HXX

// This does not work -- when this file `simd.hxx` is `#include`d from
// other thorns, then `HAVE_CAPABILITY_NSIMD` will not be defined.
// This needs to be handled in a configuration script instead.

// // Disable SIMD when the `NSIMD` library is not available
// #ifndef HAVE_CAPABILITY_NSIMD
// #ifndef SIMD_DISABLE
// #define SIMD_DISABLE
// #endif
// #endif

// Accept `SIMD_CPU` as well for backward compatibility
#ifdef SIMD_CPU
#define SIMD_DISABLE
#endif

#include "defs.hxx"

#ifndef SIMD_DISABLE
#include <nsimd/nsimd-all.hpp>
#undef vec // This should arguably not be defined in C++
#else
#endif

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>
#include <type_traits>

namespace Arith {
using namespace std;

template <typename T> struct simd;
template <typename T> struct simdl;

#ifndef SIMD_DISABLE
namespace detail {
struct reinterpret32 {
  typedef i32 int_type;
  typedef u32 unsigned_type;
  typedef f32 float_type;
};
struct reinterpret64 {
  typedef i64 int_type;
  typedef u64 unsigned_type;
  typedef f64 float_type;
};

template <typename T> struct reinterpret;
template <> struct reinterpret<i32> : reinterpret32 {};
template <> struct reinterpret<u32> : reinterpret32 {};
template <> struct reinterpret<f32> : reinterpret32 {};
template <> struct reinterpret<i64> : reinterpret64 {};
template <> struct reinterpret<u64> : reinterpret64 {};
template <> struct reinterpret<f64> : reinterpret64 {};

template <typename T> using int_type = typename reinterpret<T>::int_type;
template <typename T>
using unsigned_type = typename reinterpret<T>::unsigned_type;
template <typename T> using float_type = typename reinterpret<T>::float_type;
} // namespace detail
#endif

////////////////////////////////////////////////////////////////////////////////

extern std::size_t flop_count, memop_count;
#pragma omp threadprivate(flop_count, memop_count)

void reset_counts();
std::size_t get_flop_count();
std::size_t get_memop_count();

////////////////////////////////////////////////////////////////////////////////

// A SIMD vector, holding elements of type `T`. `T` can be a
// floating-point or an integer type. To represent logical values
// (i.e. booleans), see `simdl` below.

// This class is based on NSIMD
// <https://github.com/agenium-scale/nsimd>. The NSIMD library is very
// close to what we need; with a few bugs corrected or features added,
// this wrapper class here might be avoided.

template <typename T> struct simd {
  using value_type = T;
#ifndef SIMD_DISABLE
  using storage_type = nsimd::pack<T>;
#else
  using storage_type = T;
#endif
  storage_type elts;
  static constexpr std::size_t storage_size = sizeof(storage_type) / sizeof(T);

  constexpr simd(const simd &) = default;
  constexpr simd(simd &&) = default;
  constexpr simd &operator=(const simd &) = default;
  constexpr simd &operator=(simd &&) = default;

  constexpr ARITH_DEVICE ARITH_HOST simd() {}
  constexpr ARITH_DEVICE ARITH_HOST simd(const T &a) : elts(a) {}
#ifndef SIMD_DISABLE
  template <typename U = T,
            enable_if_t<!is_same_v<nsimd::pack<T>, T> > * = nullptr>
  constexpr ARITH_DEVICE ARITH_HOST simd(const nsimd::pack<T> &elts)
      : elts(elts) {}
#endif

  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE ARITH_DEVICE ARITH_HOST void
  count_flop(const std::size_t n = 1) {
#ifdef ARITH_COUNT_OPS
    flop_count += n * storage_size;
#endif
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE ARITH_DEVICE ARITH_HOST void
  count_memop(const std::size_t n = 1) {
#ifdef ARITH_COUNT_OPS
    memop_count += n * storage_size;
#endif
  }

  constexpr ARITH_DEVICE ARITH_HOST std::size_t size() const {
#ifndef SIMD_DISABLE
    return sizeof(nsimd::pack<T>) / sizeof(T);
#else
    return 1;
#endif
  }

  constexpr ARITH_DEVICE ARITH_HOST T operator[](const std::ptrdiff_t n) const {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < int(storage_size));
#endif
#ifndef SIMD_DISABLE
    T xarr[storage_size];
    storeu(xarr, *this);
    return xarr[n];
#else
    return elts;
#endif
  }

  template <std::size_t I>
  friend constexpr
      ARITH_DEVICE ARITH_HOST std::enable_if_t<(I >= 0 && I < storage_size), T>
      get(const simd &x) {
    return x[I];
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator!(const simd &x) {
    count_flop();
    return !x.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator~(const simd &x) {
    count_flop();
    return ~x.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd operator+(const simd &x) {
    count_flop();
    return x.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator-(const simd &x) {
    count_flop();
    return -x.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd operator&(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts & y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator|(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts | y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator^(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts ^ y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator+(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts + y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator-(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts - y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator*(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts * y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator/(const simd &x,
                                                          const simd &y) {
    count_flop(10);
    return x.elts / y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator%(const simd &x,
                                                          const simd &y) {
    count_flop();
    return x.elts % y.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd operator&(const T &a,
                                                          const simd &y) {
    count_flop();
    return a & y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator|(const T &a,
                                                          const simd &y) {
    count_flop();
    return a | y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator^(const T &a,
                                                          const simd &y) {
    count_flop();
    return a ^ y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator+(const T &a,
                                                          const simd &y) {
    count_flop();
    return a + y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator-(const T &a,
                                                          const simd &y) {
    count_flop();
    return a - y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator*(const T &a,
                                                          const simd &y) {
    count_flop();
    return a * y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator/(const T &a,
                                                          const simd &y) {
    count_flop(10);
    return a / y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator%(const T &a,
                                                          const simd &y) {
    count_flop();
    return a % y.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd operator&(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts & b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator|(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts | b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator^(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts ^ b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator+(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts + b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator-(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts - b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator*(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts * b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator/(const simd &x,
                                                          const T &b) {
    count_flop(10);
    return x.elts / b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd operator%(const simd &x,
                                                          const T &b) {
    count_flop();
    return x.elts % b;
  }

  constexpr ARITH_DEVICE ARITH_HOST simd &operator&=(const simd &x) {
    return *this = *this & x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator|=(const simd &x) {
    return *this = *this | x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator^=(const simd &x) {
    return *this = *this ^ x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator+=(const simd &x) {
    return *this = *this + x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator-=(const simd &x) {
    return *this = *this - x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator*=(const simd &x) {
    return *this = *this * x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator/=(const simd &x) {
    return *this = *this / x;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator%=(const simd &x) {
    return *this = *this % x;
  }

  constexpr ARITH_DEVICE ARITH_HOST simd &operator&=(const T &a) {
    return *this = *this & a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator|=(const T &a) {
    return *this = *this | a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator^=(const T &a) {
    return *this = *this ^ a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator+=(const T &a) {
    return *this = *this + a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator-=(const T &a) {
    return *this = *this - a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator/=(const T &a) {
    return *this = *this / a;
  }
  constexpr ARITH_DEVICE ARITH_HOST simd &operator%=(const T &a) {
    return *this = *this % a;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator==(const simd &x,
                                                               const simd &y) {
    count_flop();
    return x.elts == y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator!=(const simd &x,
                                                               const simd &y) {
    count_flop();
    return x.elts != y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<(const simd &x,
                                                              const simd &y) {
    count_flop();
    return x.elts < y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>(const simd &x,
                                                              const simd &y) {
    count_flop();
    return x.elts > y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<=(const simd &x,
                                                               const simd &y) {
    count_flop();
    return x.elts <= y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>=(const simd &x,
                                                               const simd &y) {
    count_flop();
    return x.elts >= y.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator==(const T &a,
                                                               const simd &y) {
    count_flop();
    return a == y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator!=(const T &a,
                                                               const simd &y) {
    count_flop();
    return a != y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<(const T &a,
                                                              const simd &y) {
    count_flop();
    return a < y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>(const T &a,
                                                              const simd &y) {
    count_flop();
    return a > y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<=(const T &a,
                                                               const simd &y) {
    count_flop();
    return a <= y.elts;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>=(const T &a,
                                                               const simd &y) {
    count_flop();
    return a >= y.elts;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator==(const simd &x,
                                                               const T &b) {
    count_flop();
    return x.elts == b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator!=(const simd &x,
                                                               const T &b) {
    count_flop();
    return x.elts != b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<(const simd &x,
                                                              const T &b) {
    count_flop();
    return x.elts < b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>(const simd &x,
                                                              const T &b) {
    count_flop();
    return x.elts > b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator<=(const simd &x,
                                                               const T &b) {
    count_flop();
    return x.elts <= b;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> operator>=(const simd &x,
                                                               const T &b) {
    count_flop();
    return x.elts >= b;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd abs(const simd &x) {
    count_flop();
    using std::abs;
    return abs(x.elts);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd acos(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return acos_u10(x.elts);
#else
    using std::acos;
    return acos(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd acosh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return acosh_u10(x.elts);
#else
    using std::acosh;
    return acosh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST bool allisfinite(const simd &x) {
    using std::isfinite;
    return all(isfinite(x));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd andnot(const simd &x,
                                                       const simd &y) {
#ifndef SIMD_DISABLE
    count_flop();
    return andnotb(x.elts, y.elts);
#else
    return x & ~y;
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST bool anyisnan(const simd &x) {
    using std::isnan;
    return any(isnan(x));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd asin(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return asin_u10(x.elts);
#else
    using std::asin;
    return asin(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd asinh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return asinh_u10(x.elts);
#else
    using std::asinh;
    return asinh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd atan(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return atan_u10(x.elts);
#else
    using std::atan;
    return atan(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd atanh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return atanh_u10(x.elts);
#else
    using std::atanh;
    return atanh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd cbrt(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return cbrt_u10(x.elts);
#else
    using std::cbrt;
    return cbrt(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd copysign(const simd &x,
                                                         const simd &y) {
#ifndef SIMD_DISABLE
    typedef detail::unsigned_type<T> U;
    const T signmask =
        nsimd::scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return andnot(x, signmask) | (y & signmask);
#else
    count_flop();
    using std::copysign;
    return copysign(x.elts, y.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd cos(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return cos_u10(x.elts);
#else
    using std::cos;
    return cos(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd cosh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return cosh_u10(x.elts);
#else
    using std::cosh;
    return cosh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd cospi(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return cospi_u05(x.elts);
#else
    using std::acos, std::cos;
    const T pi = acos(T(-1));
    return cos(pi * x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd exp(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return exp_u10(x.elts);
#else
    using std::exp;
    return exp(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd exp10(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return exp10_u10(x.elts);
#else
    using std::exp2;
    constexpr T log2_10 = log2(T(10));
    return exp2(log2_10 * x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd exp2(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return exp2_u10(x.elts);
#else
    using std::exp2;
    return exp2(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd fabs(const simd &x) {
    count_flop();
    using std::abs;
    return abs(x.elts);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd flipsign(const simd &x,
                                                         const simd &y) {
#ifndef SIMD_DISABLE
    typedef detail::unsigned_type<T> U;
    const T signmask =
        nsimd::scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return x ^ (y & signmask);
#else
    count_flop();
    using std::copysign;
    return copysign(1, y) * x;
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd fmax(const simd &x,
                                                     const simd &y) {
    count_flop();
    using std::max;
    return max(x.elts, y.elts);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd fmax(const T &a,
                                                     const simd &y) {
    count_flop();
    using std::max;
    return max(simd(a), y);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd fmax(const simd &x,
                                                     const T &b) {
    count_flop();
    using std::max;
    return max(x, simd(b));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd fmin(const simd &x,
                                                     const simd &y) {
    count_flop();
    using std::min;
    return min(x.elts, y.elts);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd fmin(const T &a,
                                                     const simd &y) {
    count_flop();
    using std::min;
    return min(simd(a), y);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd fmin(const simd &x,
                                                     const T &b) {
    count_flop();
    using std::min;
    return min(x, simd(b));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd hypot(const simd &x,
                                                      const simd &y) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return hypot_u05(x.elts, y.elts);
#else
    using std::hypot;
    return hypot(x.elts, y.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd inv(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return rec(x.elts);
#else
    return inv(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> isfinite(const simd &x) {
    // using std::isfinite;
    // return isfinite(x.elts);
    return !isnan(x) && !isinf(x);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> isinf(const simd &x) {
    // using std::isinf;
    // return isinf(x.elts);
    return fabs(x) == numeric_limits<T>::infinity();
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> isnan(const simd &x) {
    // using std::isnan;
    // return isnan(x.elts);
    return x != x;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd log(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return log_u10(x.elts);
#else
    using std::log;
    return log(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd log10(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return log10_u10(x.elts);
#else
    using std::log10;
    return log10(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd log2(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return log2_u10(x.elts);
#else
    using std::log2;
    return log2(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd max(const simd &x,
                                                    const simd &y) {
    count_flop();
    using std::max;
    return max(x.elts, y.elts);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd max(const simd &x, const T &b) {
    count_flop();
    using std::max;
    return max(x, simd(b));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd max(const T &a, const simd &y) {
    count_flop();
    using std::max;
    return max(simd(a), y);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd
  max(std::initializer_list<simd> xs) {
    using std::max;
    simd r = simd(-std::numeric_limits<T>::infinity());
    for (const auto x : xs)
      r = max(r, x);
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd min(const simd &x,
                                                    const simd &y) {
    count_flop();
    using std::min;
    return min(x.elts, y.elts);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd min(const T &a, const simd &y) {
    count_flop();
    using std::min;
    return min(simd(a), y);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd min(const simd &x, const T &b) {
    count_flop();
    using std::min;
    return min(x, simd(b));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd
  min(std::initializer_list<simd> xs) {
    using std::min;
    simd r = simd(std::numeric_limits<T>::infinity());
    for (const auto x : xs)
      r = min(r, x);
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const simd &x,
                                                       const simd &y,
                                                       const simd &z) {
    count_flop(2);
#ifndef SIMD_DISABLE
    return nsimd::fma(x.elts, y.elts, z.elts);
#else
    return muladd(x.elts, y.elts, z.elts);
#endif
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const T &a,
                                                       const simd &y,
                                                       const simd &z) {
    return muladd(simd(a), y, z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const simd &x,
                                                       const T &b,
                                                       const simd &z) {
    return muladd(x, simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const T &a, const T &b,
                                                       const simd &z) {
    return muladd(simd(a), simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const simd &x,
                                                       const simd &y,
                                                       const T &c) {
    return muladd(x, y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const T &a,
                                                       const simd &y,
                                                       const T &c) {
    return muladd(simd(a), y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd muladd(const simd &x,
                                                       const T &b, const T &c) {
    return muladd(x, simd(b), simd(c));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const simd &x,
                                                       const simd &y,
                                                       const simd &z) {
    count_flop(2);
#ifndef SIMD_DISABLE
    return nsimd::fms(x.elts, y.elts, z.elts);
#else
    return mulsub(x.elts, y.elts, z.elts);
#endif
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const T &a,
                                                       const simd &y,
                                                       const simd &z) {
    return mulsub(simd(a), y, z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const simd &x,
                                                       const T &b,
                                                       const simd &z) {
    return mulsub(x, simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const T &a, const T &b,
                                                       const simd &z) {
    return mulsub(simd(a), simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const simd &x,
                                                       const simd &y,
                                                       const T &c) {
    return mulsub(x, y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const T &a,
                                                       const simd &y,
                                                       const T &c) {
    return mulsub(simd(a), y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd mulsub(const simd &x,
                                                       const T &b, const T &c) {
    return mulsub(x, simd(b), simd(c));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const simd &x,
                                                          const simd &y,
                                                          const simd &z) {
    count_flop(2);
#ifndef SIMD_DISABLE
    return nsimd::fnma(x.elts, y.elts, z.elts);
#else
    return negmuladd(x.elts, y.elts, z.elts);
#endif
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const T &a,
                                                          const simd &y,
                                                          const simd &z) {
    return negmuladd(simd(a), y, z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const simd &x,
                                                          const T &b,
                                                          const simd &z) {
    return negmuladd(x, simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const T &a,
                                                          const T &b,
                                                          const simd &z) {
    return negmuladd(simd(a), simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const simd &x,
                                                          const simd &y,
                                                          const T &c) {
    return negmuladd(x, y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const T &a,
                                                          const simd &y,
                                                          const T &c) {
    return negmuladd(simd(a), y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmuladd(const simd &x,
                                                          const T &b,
                                                          const T &c) {
    return negmuladd(x, simd(b), simd(c));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const simd &x,
                                                          const simd &y,
                                                          const simd &z) {
    count_flop(2);
#ifndef SIMD_DISABLE
    return nsimd::fnms(x.elts, y.elts, z.elts);
#else
    return negmulsub(x.elts, y.elts, z.elts);
#endif
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const T &a,
                                                          const simd &y,
                                                          const simd &z) {
    return negmulsub(simd(a), y, z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const simd &x,
                                                          const T &b,
                                                          const simd &z) {
    return negmulsub(x, simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const T &a,
                                                          const T &b,
                                                          const simd &z) {
    return negmulsub(simd(a), simd(b), z);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const simd &x,
                                                          const simd &y,
                                                          const T &c) {
    return negmulsub(x, y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const T &a,
                                                          const simd &y,
                                                          const T &c) {
    return negmulsub(simd(a), y, simd(c));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd negmulsub(const simd &x,
                                                          const T &b,
                                                          const T &c) {
    return negmulsub(x, simd(b), simd(c));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd pow(const simd &x,
                                                    const simd &y) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return pow_u10(x.elts, y.elts);
#else
    using std::pow;
    return pow(x.elts, y.elts);
#endif
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd pow(const simd &x, const T &b) {
    return pow(x, simd(b));
  }
  friend constexpr ARITH_DEVICE ARITH_HOST simd pow(const T &a, const simd &y) {
    return pow(simd(a), y);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> signbit(const simd &x) {
#ifndef SIMD_DISABLE
    typedef detail::unsigned_type<T> U;
    const T signmask =
        nsimd::scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return to_logical(x & signmask);
#else
    using std::signbit;
    return signbit(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd sin(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return sin_u10(x.elts);
#else
    using std::sin;
    return sin(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd sinh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return sinh_u10(x.elts);
#else
    using std::sinh;
    return sinh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd sinpi(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return sinpi_u05(x.elts);
#else
    using std::acos, std::sin;
    const T pi = acos(T(-1));
    return sin(pi * x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd sqrt(const simd &x) {
    count_flop(10);
    using std::sqrt;
    return sqrt(x.elts);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd tan(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return tan_u10(x.elts);
#else
    using std::tan;
    return tan(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simd tanh(const simd &x) {
    count_flop(10);
#ifndef SIMD_DISABLE
    return tanh_u10(x.elts);
#else
    using std::tanh;
    return tanh(x.elts);
#endif
  }

  friend constexpr ARITH_DEVICE ARITH_HOST simdl<T> to_logical(const simd &x) {
    return to_logical(x.elts);
  }

  friend ARITH_DEVICE ARITH_HOST void storea(T *ptr, const simd &x) {
    count_memop();
#ifndef SIMD_DISABLE
    storea(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST void storeu(T *ptr, const simd &x) {
    count_memop();
#ifndef SIMD_DISABLE
    storeu(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST void mask_storea(const simdl<T> &mask, T *ptr,
                                                  const simd &x) {
    count_memop();
#ifndef SIMD_DISABLE
    mask_storea(mask.elts, ptr, x.elts);
#else
    if (mask.elts)
      *ptr = x.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST void mask_storeu(const simdl<T> &mask, T *ptr,
                                                  const simd &x) {
    count_memop();
#ifndef SIMD_DISABLE
    mask_storeu(mask.elts, ptr, x.elts);
#else
    if (mask.elts) {
#if 0 && defined __CUDACC__
      // __stwb: Cache write-back all coherent levels
      // __stcg: Cache at global level (cache in L2 and below, not L1)
      // __stcs: Cache streaming, likely to be accessed once
      // __stwt: Cache write-through (to system memory)
      __stwt(ptr, x.elts);
#elif 0 && defined __HIPCC__
      __builtin_nontemporal_store(x.elts, ptr);
#else
      // CPU
#if 0 && defined __llvm__
      __builtin_nontemporal_store(x.elts, ptr);
#else
      *ptr = x.elts;
#endif
#endif
    }
#endif
  }

  friend ostream &operator<<(ostream &os, const simd &x) {
    os << "Ⓢ[";
    for (size_t n = 0; n < storage_size; ++n) {
      if (n != 0)
        os << ",";
      os << x[n];
    }
    os << "]";
    return os;
  }
};

} // namespace Arith
namespace std {
template <typename T>
struct tuple_size<Arith::simd<T> >
    : std::integral_constant<std::size_t, Arith::simd<T>::storage_size> {};
} // namespace std
namespace Arith {

template <typename T> struct zero<simd<T> > {
  typedef simd<T> value_type;
  // static constexpr value_type value = simd<T>(zero_v<T>);
  constexpr ARITH_INLINE operator value_type() const {
    return simd<T>(zero<T>());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return simd<T>(zero<T>());
  }
};

template <typename T> struct one<simd<T> > {
  typedef simd<T> value_type;
  // static constexpr value_type value = simd<T>(one_v<T>);
  constexpr ARITH_INLINE operator value_type() const {
    return simd<T>(one<T>());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return simd<T>(one<T>());
  }
};

template <typename T> struct nan<simd<T> > {
  typedef simd<T> value_type;
  // static constexpr value_type value = simd<T>(nan_v<T>);
  constexpr ARITH_INLINE operator value_type() const {
    return simd<T>(nan<T>());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return simd<T>(nan<T>());
  }
};

template <typename VT, typename T = typename VT::value_type>
ARITH_DEVICE ARITH_HOST inline simd<T> iota() {
#ifndef SIMD_DISABLE
  return nsimd::iota<nsimd::pack<T> >();
#else
  return 0;
#endif
}

template <typename VT, typename T = typename VT::value_type>
ARITH_DEVICE ARITH_HOST inline simdl<T> mask_for_loop_tail(const int i,
                                                           const int n) {
  simd<T>::count_flop();
#ifndef SIMD_DISABLE
  return nsimd::mask_for_loop_tail<nsimd::packl<T> >(i, n);
#else
  return i < n;
#endif
}

template <typename VT, typename T = typename VT::value_type>
ARITH_DEVICE ARITH_HOST inline simd<T> loada(const T *ptr) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return nsimd::loada<nsimd::pack<T> >(ptr);
#else
  return *ptr;
#endif
}

template <typename VT, typename T = typename VT::value_type>
ARITH_DEVICE ARITH_HOST inline simd<T> loadu(const T *ptr) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return nsimd::loadu<nsimd::pack<T> >(ptr);
#else
  return *ptr;
#endif
}

template <typename T>
ARITH_DEVICE ARITH_HOST inline simd<T> maskz_loada(const simdl<T> &mask,
                                                   const T *ptr) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return nsimd::maskz_loada(mask.elts, ptr);
#else
  return mask.elts ? *ptr : 0;
#endif
}

template <typename T>
ARITH_DEVICE ARITH_HOST inline simd<T> maskz_loadu(const simdl<T> &mask,
                                                   const T *ptr) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return nsimd::maskz_loadu(mask.elts, ptr);
#else
  return mask.elts ? *ptr : 0;
#endif
}

template <typename T>
ARITH_DEVICE ARITH_HOST inline simd<T>
masko_loada(const simdl<T> &mask, const T *ptr, const simd<T> &other) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return masko_loada(mask.elts, ptr, other.elts);
#else
  return mask.elts ? *ptr : other.elts;
#endif
}

template <typename T>
ARITH_DEVICE ARITH_HOST inline simd<T>
masko_loadu(const simdl<T> &mask, const T *ptr, const simd<T> &other) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return masko_loadu(mask.elts, ptr, other.elts);
#else
  return mask.elts ? *ptr : other.elts;
#endif
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
ARITH_DEVICE ARITH_HOST inline simd<T>
masko_loada(const simdl<T> &mask, const T *ptr, const U &other) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return masko_loada(mask, ptr, simd<T>(other));
#else
  return mask.elts ? *ptr : other;
#endif
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
ARITH_DEVICE ARITH_HOST inline simd<T>
masko_loadu(const simdl<T> &mask, const T *ptr, const U &other) {
  simd<T>::count_memop();
#ifndef SIMD_DISABLE
  return masko_loadu(mask, ptr, simd<T>(other));
#else
  return mask.elts ? *ptr : other;
#endif
}

////////////////////////////////////////////////////////////////////////////////

// A SIMD vector of booleans, usable with `simd<T>`.
//
// Booleans are special in SIMD operations. There is not just one
// boolean type, but rather one boolean type for each (floating-point
// or integer) type `T`. The reason is that SIMD operations do not
// like changing the number of bits, and hence there are 64-bit
// booleans, 32-bit booleans, 16-bit booleans, etc.
//
// This is described in the NSIMD documentation
// <https://github.com/agenium-scale/nsimd>.

template <typename T> struct simdl {
  typedef T value_type;
#ifndef SIMD_DISABLE
  using storage_type = nsimd::packl<T>;
#else
  using storage_type = bool;
#endif
  storage_type elts;
  static constexpr std::size_t storage_size = simd<T>::storage_size;

  constexpr simdl(const simdl &) = default;
  constexpr simdl(simdl &&) = default;
  constexpr simdl &operator=(const simdl &) = default;
  constexpr simdl &operator=(simdl &&) = default;

  constexpr ARITH_DEVICE ARITH_HOST simdl() {}
  constexpr ARITH_DEVICE ARITH_HOST simdl(bool a) : elts(a) {}
#ifndef SIMD_DISABLE
  constexpr ARITH_DEVICE ARITH_HOST simdl(const nsimd::packl<T> &elts)
      : elts(elts) {}
#endif

  constexpr ARITH_DEVICE ARITH_HOST std::size_t size() const {
#ifndef SIMD_DISABLE
    return sizeof(nsimd::packl<T>) / sizeof(T);
#else
    return 1;
#endif
  }

  constexpr ARITH_DEVICE ARITH_HOST bool
  operator[](const std::ptrdiff_t n) const {
#ifndef SIMD_DISABLE
    // TOOD: Introduce `to_mask` for simd/simdl
    return simd<T>(nsimd::to_mask(elts))[n];
#else
    return elts;
#endif
  }

  template <std::size_t I>
  friend constexpr ARITH_DEVICE ARITH_HOST
      std::enable_if_t<(I >= 0 && I < storage_size), bool>
      get(const simdl &x) {
    return x[I];
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator!(const simdl &x) {
    return !x.elts;
  }
  simdl operator~() const { return !elts; }

  friend ARITH_DEVICE ARITH_HOST simdl operator&&(const simdl &x,
                                                  const simdl &y) {
    return x.elts && y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator&(const simdl &x,
                                                 const simdl &y) {
    return x.elts && y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator||(const simdl &x,
                                                  const simdl &y) {
    return x.elts || y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator|(const simdl &x,
                                                 const simdl &y) {
    return x.elts || y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator^(const simdl &x,
                                                 const simdl &y) {
#ifndef SIMD_DISABLE
    return xorl(x.elts, y.elts);
#else
    return x.elts ^ y.elts;
#endif
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator&(const bool a, const simdl &y) {
    return a & y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator|(const bool a, const simdl &y) {
    return a | y.elts;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator^(const bool a, const simdl &y) {
    return a ^ y.elts;
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator&(const simdl &x, const bool b) {
    return x.elts & b;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator|(const simdl &x, const bool b) {
    return x.elts | b;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator^(const simdl &x, const bool b) {
    return x.elts ^ b;
  }

  ARITH_DEVICE ARITH_HOST simdl &operator&=(const simdl &x) {
    return *this = *this & x;
  }
  ARITH_DEVICE ARITH_HOST simdl &operator|=(const simdl &x) {
    return *this = *this | x;
  }
  ARITH_DEVICE ARITH_HOST simdl &operator^=(const simdl &x) {
    return *this = *this ^ x;
  }

  ARITH_DEVICE ARITH_HOST simdl &operator&=(const bool a) {
    return *this = *this & a;
  }
  ARITH_DEVICE ARITH_HOST simdl &operator|=(const bool a) {
    return *this = *this | a;
  }
  ARITH_DEVICE ARITH_HOST simdl &operator^=(const bool a) {
    return *this = *this ^ a;
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator==(const simdl &x,
                                                  const simdl &y) {
    return !(x != y);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator!=(const simdl &x,
                                                  const simdl &y) {
    return x ^ y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<(const simdl &x,
                                                 const simdl &y) {
    return !x & y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>(const simdl &x,
                                                 const simdl &y) {
    return y < x;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<=(const simdl &x,
                                                  const simdl &y) {
    return !(x > y);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>=(const simdl &x,
                                                  const simdl &y) {
    return !(x < y);
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator==(const bool a,
                                                  const simdl &y) {
    return simdl(a) == y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator!=(const bool a,
                                                  const simdl &y) {
    return simdl(a) != y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<(const bool a, const simdl &y) {
    return simdl(a) < y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>(const bool a, const simdl &y) {
    return simdl(a) > y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<=(const bool a,
                                                  const simdl &y) {
    return simdl(a) <= y;
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>=(const bool a,
                                                  const simdl &y) {
    return simdl(a) >= y;
  }

  friend ARITH_DEVICE ARITH_HOST simdl operator==(const simdl &x,
                                                  const bool b) {
    return x == simdl(b);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator!=(const simdl &x,
                                                  const bool b) {
    return x != simdl(b);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<(const simdl &x, const bool b) {
    return x < simdl(b);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>(const simdl &x, const bool b) {
    return x > simdl(b);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator<=(const simdl &x,
                                                  const bool b) {
    return x <= simdl(b);
  }
  friend ARITH_DEVICE ARITH_HOST simdl operator>=(const simdl &x,
                                                  const bool b) {
    return x >= simdl(b);
  }

  friend ARITH_DEVICE ARITH_HOST bool all(const simdl &x) {
#ifndef SIMD_DISABLE
    return all(x.elts);
#else
    return x.elts;
#endif
  }

  friend ARITH_DEVICE ARITH_HOST simdl andnot(const simdl &x, const simdl &y) {
#ifndef SIMD_DISABLE
    return andnotl(x.elts, y.elts);
#else
    return x && !y;
#endif
  }

  friend ARITH_DEVICE ARITH_HOST bool any(const simdl &x) {
#ifndef SIMD_DISABLE
    return any(x.elts);
#else
    return x.elts;
#endif
  }

  friend ARITH_DEVICE ARITH_HOST simd<T>
  if_else(const simdl &cond, const simd<T> &x, const simd<T> &y) {
#ifndef SIMD_DISABLE
    return if_else1(cond.elts, x.elts, y.elts);
#else
    return cond.elts ? x.elts : y.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST simd<T> if_else(const simdl &cond, const T &a,
                                                 const simd<T> &y) {
    return if_else(cond, simd(a), y);
  }
  friend ARITH_DEVICE ARITH_HOST simd<T> if_else(const simdl &cond,
                                                 const simd<T> &x, const T &b) {
    return if_else(cond, x, simd(b));
  }
  friend ARITH_DEVICE ARITH_HOST simd<T> if_else(const simdl &cond, const T &a,
                                                 const T &b) {
    return if_else(cond, simd(a), simd(b));
  }

  friend ARITH_DEVICE ARITH_HOST simdl<T>
  if_else(const simdl &cond, const simdl<T> &x, const simdl<T> &y) {
#ifndef SIMD_DISABLE
    return if_else1(cond.elts, x.elts, y.elts);
#else
    return cond.elts ? x.elts : y.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST simdl<T>
  if_else(const simdl &cond, const bool &a, const simdl<T> &y) {
    return if_else(cond, simdl(a), y);
  }
  friend ARITH_DEVICE ARITH_HOST simdl<T>
  if_else(const simdl &cond, const simdl<T> &x, const bool &b) {
    return if_else(cond, x, simdl(b));
  }
  friend ARITH_DEVICE ARITH_HOST simdl<T>
  if_else(const simdl &cond, const bool &a, const bool &b) {
    return if_else(cond, simdl(a), simdl(b));
  }

  friend ARITH_DEVICE ARITH_HOST void storela(T *ptr, const simdl &x) {
#ifndef SIMD_DISABLE
    storela(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend ARITH_DEVICE ARITH_HOST void storelu(T *ptr, const simdl &x) {
#ifndef SIMD_DISABLE
    storelu(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }

  friend ostream &operator<<(ostream &os, const simdl &x) {
    os << "Ⓑ[";
    for (size_t n = 0; n < storage_size; ++n) {
      if (n != 0)
        os << ",";
      os << x[n];
    }
    os << "]";
    return os;
  }
};

} // namespace Arith
namespace std {
template <typename T>
struct tuple_size<Arith::simdl<T> >
    : std::integral_constant<std::size_t, Arith::simdl<T>::storage_size> {};
} // namespace std
namespace Arith {

template <typename T> struct zero<simdl<T> > {
  typedef simdl<T> value_type;
  // static constexpr value_type value = simdl<T>(false);
  constexpr ARITH_INLINE operator value_type() const { return simdl<T>(false); }
  constexpr ARITH_INLINE value_type operator()() const {
    return simdl<T>(false);
  }
};

template <typename T> struct one<simdl<T> > {
  typedef simdl<T> value_type;
  // static constexpr value_type value = simdl<T>(true);
  constexpr ARITH_INLINE operator value_type() const { return simdl<T>(true); }
  constexpr ARITH_INLINE value_type operator()() const {
    return simdl<T>(true);
  }
};
} // namespace Arith

// Put these three things into the fake Cactus namespace
using CCTK_BOOLVEC = Arith::simdl<CCTK_REAL>;
using CCTK_REALVEC = Arith::simd<CCTK_REAL>;
constexpr std::size_t CCTK_VECSIZE = std::tuple_size_v<CCTK_REALVEC>;

#endif // #ifndef CARPETX_ARITH_SIMD_HXX
