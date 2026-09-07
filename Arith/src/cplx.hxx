#ifndef CARPETX_ARITH_CPLX_HXX
#define CARPETX_ARITH_CPLX_HXX

#include "defs.hxx"
#include "vect.hxx"

#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>

namespace Arith {

template <typename T> struct cplx {
  T real, imag;

  constexpr ARITH_INLINE cplx(const cplx &) = default;
  constexpr ARITH_INLINE cplx(cplx &&) = default;
  constexpr ARITH_INLINE cplx &operator=(const cplx &) = default;
  constexpr ARITH_INLINE cplx &operator=(cplx &&) = default;

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx() : real(), imag() {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx(const T &x)
      : real(x), imag() {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx(const T &x, const T &y)
      : real(x), imag(y) {}

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator+(const cplx &x) {
    return {+x.real, +x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator-(const cplx &x) {
    return {-x.real, -x.imag};
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator+(const cplx &x, const cplx &y) {
    return {x.real + y.real, x.imag + y.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator-(const cplx &x, const cplx &y) {
    return {x.real - y.real, x.imag - y.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator+(const cplx &x, const T &y) {
    return {x.real + y, x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator-(const cplx &x, const T &y) {
    return {x.real - y, x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator*(const cplx &x, const T &y) {
    return {x.real * y, x.imag * y};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator*(const T &x, const cplx &y) {
    return {x * y.real, x * y.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator/(const cplx &x, const T &y) {
    return {x.real / y, x.imag / y};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator*(const cplx &x, const cplx &y) {
    return {x.real * y.real - x.imag * y.imag,
            x.real * y.imag + x.imag * y.real};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  operator/(const cplx &x, const cplx &y) {
    return x * conj(y) / norm(y);
  }

  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator+=(const cplx &x) {
    return *this = *this + x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator-=(const cplx &x) {
    return *this = *this - x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator*=(const T &x) {
    return *this = *this * x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator/=(const T &x) {
    return *this = *this / x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator*=(const cplx &x) {
    return *this = *this * x;
  }
  inline ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx &operator/=(const cplx &x) {
    return *this = *this / x;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator==(const cplx &x, const cplx &y) {
    return x.real == y.real && x.imag == y.imag;
  };
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator<(const cplx &x, const cplx &y) {
    return x.real < y.real;
  };

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator!=(const cplx &x, const cplx &y) {
    return !(x == y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator>(const cplx &x, const cplx &y) {
    return y < x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator<=(const cplx &x, const cplx &y) {
    return !(x > y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator>=(const cplx &x, const cplx &y) {
    return !(x < y);
  }

  // Modulus |z|. Returns a real number, as `std::abs` does for
  // `std::complex`. `hypot` avoids the overflow of `sqrt(norm(x))`.
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T abs(const cplx &x) {
    using std::hypot;
    return hypot(x.real, x.imag);
  }
  // Argument arg(z), in (-pi, pi]
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T arg(const cplx &x) {
    using std::atan2;
    return atan2(x.imag, x.real);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  allisfinite(const cplx &x) {
    return allisfinite(x.real) && allisfinite(x.imag);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const cplx &x) {
    return anyisnan(x.real) || anyisnan(x.imag);
  }
  // Principal cube root, i.e. the root with the argument closest to zero.
  // (`std::complex` offers no `cbrt`; this is `exp(log(z)/3)`, evaluated in
  // polar form.)
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  cbrt(const cplx &x) {
    using std::cbrt, std::cos, std::sin;
    const T r = cbrt(abs(x));
    const T theta = arg(x) / T(3);
    return {r * cos(theta), r * sin(theta)};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  conj(const cplx &x) {
    return {x.real, -x.imag};
  }
  // cos(z) = cos(re) cosh(im) - i sin(re) sinh(im)
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  cos(const cplx &x) {
    using std::cos, std::cosh, std::sin, std::sinh;
    return {cos(x.real) * cosh(x.imag), -(sin(x.real) * sinh(x.imag))};
  }
  // exp(z) = exp(re) (cos(im) + i sin(im))
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  exp(const cplx &x) {
    using std::cos, std::exp, std::sin;
    const T r = exp(x.real);
    return {r * cos(x.imag), r * sin(x.imag)};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T imag(const cplx &x) {
    return x.imag;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  isnan(const cplx &x) {
    using std::isnan;
    return isnan(x.real) || isnan(x.imag);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T norm(const cplx &x) {
    return x.real * x.real + x.imag * x.imag;
  }
  // Exponentiation by squaring, spelled iteratively as `pown` in defs.hxx is.
  // It must not recurse: this function is ARITH_INLINE, and GCC rejects
  // `always_inline` on a recursive call it cannot unroll -- which is exactly
  // what happens once `n` is not a compile-time constant.
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx pow(cplx x,
                                                                 int n) {
    const bool invert = n < 0;
    // Negating `n` would overflow for INT_MIN, so take the magnitude unsigned
    unsigned m = invert ? -static_cast<unsigned>(n) : static_cast<unsigned>(n);
    cplx r = cplx(one<T>()());
    cplx y = x;
    // invariant: initial(x)^m == r * y^m
    while (m) {
      if (m & 1)
        r = r * y;
      y = y * y;
      m >>= 1;
    }
    return invert ? cplx(one<T>()()) / r : r;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  pow2(const cplx &x) {
    return x * x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T real(const cplx &x) {
    return x.real;
  }
  // sin(z) = sin(re) cosh(im) + i cos(re) sinh(im)
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  sin(const cplx &x) {
    using std::cos, std::cosh, std::sin, std::sinh;
    return {sin(x.real) * cosh(x.imag), cos(x.real) * sinh(x.imag)};
  }
  // Principal square root, i.e. the root with non-negative real part.
  // Evaluated via w = sqrt((|z| + |re|) / 2) so that neither component
  // suffers cancellation:
  //   re >= 0:  sqrt(z) = w + i im / (2 w)
  //   re <  0:  sqrt(z) = |im| / (2 w) + i copysign(w, im)
  // Written branch-free because `T` may be a SIMD vector.
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  sqrt(const cplx &x) {
    using std::copysign, std::fabs, std::sqrt;
    const T w = sqrt((abs(x) + fabs(x.real)) / T(2));
    const auto nonneg = x.real >= T(0);
    const T re = if_else(nonneg, w, fabs(x.imag) / (T(2) * w));
    const T im = if_else(nonneg, x.imag / (T(2) * w), copysign(w, x.imag));
    // sqrt(0) = 0, where the expressions above evaluate to 0/0
    const auto is_zero = w == T(0);
    return {if_else(is_zero, T(0), re), if_else(is_zero, T(0), im)};
  }

  friend std::ostream &operator<<(std::ostream &os, const cplx &x) {
    return os << x.real << "+" << x.imag << "*i";
  }
};

template <typename T> struct zero<cplx<T> > {
  typedef cplx<T> value_type;
  // static constexpr value_type value = cplx<T>(zero<T>::value,
  // zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return cplx<T>(zero<T>(), zero<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return cplx<T>(zero<T>(), zero<T>());
  }
};

template <typename T> struct one<cplx<T> > {
  typedef cplx<T> value_type;
  // static constexpr value_type value = cplx<T>(one<T>::value, zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return cplx<T>(one<T>(), zero<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return cplx<T>(one<T>(), zero<T>());
  }
};

template <typename T> struct nan<cplx<T> > {
  typedef cplx<T> value_type;
  // static constexpr value_type value = cplx<T>(nan<T>::value, nan<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return cplx<T>(nan<T>(), nan<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return cplx<T>(nan<T>(), nan<T>());
  }
};

} // namespace Arith
namespace std {
// These are the comparisons the standard containers use: `equal_to` and
// `hash` consider both components, `less` orders lexicographically. Note that
// `cplx::operator<` and friends compare only the real part -- complex numbers
// have no natural order, and that operator exists only so that a `cplx` can be
// branched on like a number. `std::less<cplx<T>>` is therefore deliberately
// *not* the same relation as `operator<`; use it, not the operator, whenever an
// ordering has to be a strict weak ordering (`std::map`, `std::sort`).
template <typename T> struct equal_to<Arith::cplx<T> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    return equal_to<T>()(x.real, y.real) && equal_to<T>()(x.imag, y.imag);
  }
};

template <typename T> struct less<Arith::cplx<T> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    if (less<T>()(x.real, y.real))
      return true;
    if (less<T>()(y.real, x.real))
      return false;
    if (less<T>()(x.imag, y.imag))
      return true;
    if (less<T>()(y.imag, x.imag))
      return false;
    return false;
  }
};

// Consistent with `equal_to` above, i.e. it hashes both components
template <typename T> struct hash<Arith::cplx<T> > {
  std::size_t operator()(const Arith::cplx<T> &x) const {
    const std::size_t h = hash<T>()(x.real);
    return h ^ (hash<T>()(x.imag) + 0x9e3779b9 + (h << 6) + (h >> 2));
  }
};
} // namespace std

#endif // #ifndef CARPETX_ARITH_CPLX_HXX
