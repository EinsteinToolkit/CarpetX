#ifndef CARPETX_ARITH_CPLX_HXX
#define CARPETX_ARITH_CPLX_HXX

#include "defs.hxx"
#include "vect.hxx"

#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>

namespace Arith {
using namespace std;

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
    return x.real == y.real;
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

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  abs(const cplx &x) {
    return sqrt(norm(x));
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  allisfinite(const cplx &x) {
    return allisfinite(x.real) && allisfinite(x.imag);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const cplx &x) {
    return anyisnan(x.real) || anyisnan(x.imag);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  cbrt(const cplx &x) {
    using std::cbrt;
    const T r = cbrt(x.real);
    return {r, r / (3 * x.real) * x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  conj(const cplx &x) {
    return {x.real, -x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  cos(const cplx &x) {
    using std::cos, std::sin;
    return {cos(x.real), -sin(x.real) * x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  exp(const cplx &x) {
    using std::exp;
    const T r = exp(x.real);
    return {r, r * x.imag};
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
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx pow(cplx x,
                                                                 int n) {
    if (n == 0)
      return 1;
    if (n == 1)
      return x;
    // Older compilers cannot handle variable definitions in constexpr functions
    // const auto y = pow(x * x, n / 2);
    // if (n % 2 == 0)
    //   return y;
    // return x * y;
    if (n % 2 == 0)
      return pow(x * x, n / 2);
    return x * pow(x * x, n / 2);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  pow2(const cplx &x) {
    return x * x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T real(const cplx &x) {
    return x.real;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  sin(const cplx &x) {
    using std::cos, std::sin;
    return {sin(x.real), cos(x.real) * x.imag};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST cplx
  sqrt(const cplx &x) {
    using std::sqrt;
    const T r = sqrt(x.real);
    return {r, x.imag / (2 * r)};
  }

  friend ostream &operator<<(ostream &os, const cplx &x) {
    return os << x.real << "+" << x.real << "*i";
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
template <typename T> struct equal_to<Arith::cplx<T> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    return equal_to<T>()(x.real, y.real) && equal_to<T>()(x.imag, y.imag);
  }
};

template <typename T> struct less<Arith::cplx<T> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::cplx<T> &x, const Arith::cplx<T> &y) const {
    if (less<T>(x.real, y.real))
      return true;
    if (less<T>(y.real, x.real))
      return false;
    if (less<T>(x.imag, y.imag))
      return true;
    if (less<T>(y.imag, x.imag))
      return false;
    return false;
  }
};
} // namespace std

#endif // #ifndef CARPETX_ARITH_CPLX_HXX
