#ifndef CARPETX_ARITH_DUAL_HXX
#define CARPETX_ARITH_DUAL_HXX

#include "defs.hxx"
#include "simd.hxx"
#include "vect.hxx"

#include <cassert>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <ostream>

namespace Arith {
using namespace std;

template <typename T, typename U = T> struct dual;

template <typename T, typename U> struct dual {
  // `T` is the value type.
  T val;
  // `U` is either `T`, or a vector (in the mathematical sense) of
  // `T`.
  U eps;

  constexpr ARITH_INLINE dual(const dual &) = default;
  constexpr ARITH_INLINE dual(dual &&) = default;
  constexpr ARITH_INLINE dual &operator=(const dual &) = default;
  constexpr ARITH_INLINE dual &operator=(dual &&) = default;

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual()
      : val(nan<T>()()), eps(nan<U>()()) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual(const T &x)
      : val(x), eps(zero<U>()()) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual(const T &x, const U &y)
      : val(x), eps(y) {}

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator+(const dual &x) {
    return {+x.val, +x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator-(const dual &x) {
    return {-x.val, -x.eps};
  }

  template <typename = T>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator+(const dual &x, const dual &y) {
    return {x.val + y.val, x.eps + y.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator-(const dual &x, const dual &y) {
    return {x.val - y.val, x.eps - y.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator+(const dual &x, const T &y) {
    return {x.val + y, x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator-(const dual &x, const T &y) {
    return {x.val - y, x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator*(const dual &x, const T &y) {
    return {x.val * y, x.eps * y};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator*(const T &x, const dual &y) {
    return {x * y.val, x * y.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator/(const dual &x, const T &y) {
    return {x.val / y, x.eps / y};
  }
  template <typename = T>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator*(const dual &x, const dual &y) {
    return {x.val * y.val, x.val * y.eps + x.eps * y.val};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  operator/(const dual &x, const dual &y) {
    return {x.val / y.val, (x.eps * y.val - x.val * y.eps) / (y.val * y.val)};
  }

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator+=(const dual &x) {
    return *this = *this + x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator-=(const dual &x) {
    return *this = *this - x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator*=(const T &x) {
    return *this = *this * x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator/=(const T &x) {
    return *this = *this / x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator*=(const dual &x) {
    return *this = *this * x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST dual &operator/=(const dual &x) {
    return *this = *this / x;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator==(const dual &x, const dual &y) {
    return x.val == y.val;
  };
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator<(const dual &x, const dual &y) {
    return x.val < y.val;
  };

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator!=(const dual &x, const dual &y) {
    return !(x == y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator>(const dual &x, const dual &y) {
    return y < x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator<=(const dual &x, const dual &y) {
    return !(x > y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator>=(const dual &x, const dual &y) {
    return !(x < y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  abs(const dual &x) {
    // return sqrt(pow2(x));
    using std::abs, std::copysign;
    return {abs(x.val), copysign(T{1}, x.val) * x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  allisfinite(const dual &x) {
    return allisfinite(x.val) &&
           all(fmap([](const auto &a) { return allisfinite(a); }, x.eps));
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const dual &x) {
    return anyisnan(x.val) ||
           any(fmap([](const auto &a) { return anyisnan(a); }, x.eps));
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  cbrt(const dual &x) {
    using std::cbrt;
    const T r = cbrt(x.val);
    return {r, r / (3 * x.val) * x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  cos(const dual &x) {
    using std::cos, std::sin;
    return {cos(x.val), -sin(x.val) * x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  exp(const dual &x) {
    using std::exp;
    const T r = exp(x.val);
    return {r, r * x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  fabs(const dual &x) {
    // return sqrt(pow2(x));
    using std::abs, std::copysign;
    return {abs(x.val), copysign(T{1}, x.val) * x.eps};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  isnan(const dual &x) {
    using std::isnan;
    return isnan(x.val) ||
           any(fmap([](const auto &a) { return isnan(a); }, x.eps));
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  max(const dual &x, const dual &y) {
    return if_else(x >= y, x, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  max(std::initializer_list<dual> xs) {
    using std::max;
    dual r = dual(-std::numeric_limits<T>::infinity());
    for (const auto x : xs)
      r = max(r, x);
    return r;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  min(const dual &x, const dual &y) {
    return if_else(x >= y, x, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  min(std::initializer_list<dual> xs) {
    using std::min;
    dual r = dual(std::numeric_limits<T>::infinity());
    for (const auto x : xs)
      r = min(r, x);
    return r;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual pow(const dual &x,
                                                                 const int n) {
    using std::pow;
    if (n == 0)
      return {1, U{}};
    return {pow(x.val, n), n * pow(x.val, n - 1) * x.eps};
  }
  template <typename = T>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  pow2(const dual &x) {
    return x * x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  sin(const dual &x) {
    using std::cos, std::sin;
    return {sin(x.val), cos(x.val) * x.eps};
  }
  template <typename = T>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST dual
  sqrt(const dual &x) {
    using std::sqrt;
    const T r = sqrt(x.val);
    return {r, x.eps / (2 * r)};
  }

  friend ostream &operator<<(ostream &os, const dual &x) {
    return os << x.val << "+ε*" << x.eps;
  }
};

template <typename T, typename U> struct zero<dual<T, U> > {
  typedef dual<T, U> value_type;
  // static constexpr value_type value = dual<T, U>(zero<T>(), zero<U>());
  constexpr ARITH_INLINE operator value_type() const {
    return dual<T, U>(zero<T>()(), zero<U>()());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return dual<T, U>(zero<T>()(), zero<U>()());
  }
};

template <typename T, typename U> struct one<dual<T, U> > {
  typedef dual<T, U> value_type;
  // static constexpr value_type value = dual<T, U>(one<T>(), zero<U>());
  constexpr ARITH_INLINE operator value_type() const {
    return dual<T, U>(one<T>()(), zero<U>()());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return dual<T, U>(one<T>()(), zero<U>()());
  }
};

template <typename T, typename U> struct nan<dual<T, U> > {
  typedef dual<T, U> value_type;
  // static constexpr value_type value = dual<T, U>(nan<T>(), nan<U>());
  constexpr ARITH_INLINE operator value_type() const {
    return dual<T, U>(nan<T>()(), nan<U>()());
  }
  constexpr ARITH_INLINE value_type operator()() const {
    return dual<T, U>(nan<T>()(), nan<U>()());
  }
};

} // namespace Arith
namespace std {
template <typename T, typename U> struct equal_to<Arith::dual<T, U> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::dual<T, U> &x, const Arith::dual<T, U> &y) const {
    return equal_to<T>()(x.val, y.val) && equal_to<U>()(x.eps, y.eps);
  }
};

template <typename T, typename U> struct less<Arith::dual<T, U> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::dual<T, U> &x, const Arith::dual<T, U> &y) const {
    if (less<T>(x.val, y.val))
      return true;
    if (less<T>(y.val, x.val))
      return false;
    if (less<U>(x.eps, y.eps))
      return true;
    if (less<U>(y.eps, x.eps))
      return false;
    return false;
  }
};
} // namespace std

namespace Arith {

////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
constexpr dual<simd<T>, U> if_else(const simdl<T> &cond,
                                   const dual<simd<T>, U> &x,
                                   const dual<simd<T>, U> &y) {
  return dual<simd<T>, U>(
      if_else(cond, x.val, y.val),
      fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
           x.eps, y.eps));
}

template <typename T, typename U>
constexpr dual<simd<T>, U> if_else(const simdl<T> &cond, const simd<T> &x,
                                   const dual<simd<T>, U> &y) {
  return dual<simd<T>, U>(
      if_else(cond, x, y.val),
      fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
           zero<U>()(), y.eps));
}

template <typename T, typename U>
constexpr dual<simd<T>, U>
if_else(const simdl<T> &cond, const dual<simd<T>, U> &x, const simd<T> &y) {
  return dual<simd<T>, U>(
      if_else(cond, x.val, y),
      fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
           x.eps, zero<U>()()));
}

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_DUAL_HXX
