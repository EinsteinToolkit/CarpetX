#ifndef CARPETX_ARITH_RATIONAL_HXX
#define CARPETX_ARITH_RATIONAL_HXX

#include <cctk.h>

#ifdef HAVE_CAPABILITY_yaml_cpp
#include <yaml-cpp/yaml.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <type_traits>

namespace Arith {
using namespace std;

// Rational numbers

template <typename I> struct rational {
  I num, den;

  struct no_normalize {};

  constexpr void normalize() {
    const I x = gcd(num, den);
    num /= x;
    den /= x;
  }

  constexpr rational() : num(0), den(1) {}

  constexpr rational(const rational &) = default;
  constexpr rational(rational &&) = default;
  constexpr rational &operator=(const rational &) = default;
  constexpr rational &operator=(rational &&) = default;

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational(const J &i) : num(i), den(1) {}
  template <typename J, typename K, enable_if_t<is_integral_v<J> > * = nullptr,
            enable_if_t<is_integral_v<K> > * = nullptr>
  constexpr rational(const J &num, const K &den, no_normalize)
      : num(num), den(den) {}
  template <typename J, typename K, enable_if_t<is_integral_v<J> > * = nullptr,
            enable_if_t<is_integral_v<K> > * = nullptr>
  constexpr rational(const J &num, const K &den) : num(num), den(den) {
    if (this->den < 0) {
      this->num = -this->num;
      this->den = -this->den;
    }
    normalize();
  }

  template <typename F, enable_if_t<is_floating_point_v<F> > * = nullptr>
  constexpr rational(const F &f) : num(0), den(1) {
    const F f1 = nextafter(f);
    const F df = f1 - f;
    const F df1 = 1 / df;
    assert(df1 == rint(df1));
    den = llrint(df1);
    num = llrint(f * df1);
    normalize();
  }

  template <typename F, enable_if_t<is_floating_point_v<F> > * = nullptr>
  constexpr operator F() const {
    return F(num) / F(den);
  }

  friend constexpr rational operator+(const rational &x) {
    return rational(+x.num, x.den, no_normalize());
  }
  friend constexpr rational operator-(const rational &x) {
    return rational(-x.num, x.den, no_normalize());
  }

  friend constexpr rational operator+(const rational &x, const rational &y) {
    return rational(x.num * y.den + x.den * y.num, x.den * y.den);
  }
  friend constexpr rational operator-(const rational &x, const rational &y) {
    return rational(x.num * y.den - x.den * y.num, x.den * y.den);
  }
  friend constexpr rational operator*(const rational &x, const rational &y) {
    return rational(x.num * y.num, x.den * y.den);
  }
  friend constexpr rational operator/(const rational &x, const rational &y) {
    return rational(x.num * y.den, x.den * y.num);
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator+(const rational &x, const J &y) {
    return x + rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator-(const rational &x, const J &y) {
    return x - rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator*(const rational &x, const J &y) {
    return x * rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator/(const rational &x, const J &y) {
    return x / rational(y);
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator+(const J &x, const rational &y) {
    return rational(x) + y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator-(const J &x, const rational &y) {
    return rational(x) - y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator*(const J &x, const rational &y) {
    return rational(x) * y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational operator/(const J &x, const rational &y) {
    return rational(x) / y;
  }

  constexpr rational &operator+=(const rational &x) {
    return *this = *this + x;
  }
  constexpr rational &operator-=(const rational &x) {
    return *this = *this - x;
  }
  constexpr rational &operator*=(const rational &x) {
    return *this = *this * x;
  }
  constexpr rational &operator/=(const rational &x) {
    return *this = *this / x;
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational &operator+=(const J &x) {
    return *this = *this + rational(x);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational &operator-=(const J &x) {
    return *this = *this - rational(x);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational &operator*=(const J &x) {
    return *this = *this * rational(x);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  constexpr rational &operator/=(const J &x) {
    return *this = *this / rational(x);
  }

  friend constexpr rational abs(const rational &x) {
    // std::abs is not constexpr
    return rational(x.num >= 0 ? x.num : -x.num, x.den);
  }

  friend constexpr rational max(const rational &x, const rational &y) {
    return x >= y ? x : y;
  }
  friend constexpr rational min(const rational &x, const rational &y) {
    return x <= y ? x : y;
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr rational pown(rational x, J a) {
    if (a < 0) {
      x = 1 / x;
      a = -a;
    }
    rational r = 1;
    while (a > 0) {
      if (a % 2)
        r *= x;
      x *= x;
      a /= 2;
    }
    return r;
  }
  friend constexpr double pow(const rational &x, const rational &y) {
    using std::pow;
    return pow(double(x), double(y));
  }

  friend constexpr bool operator==(const rational &x, const rational &y) {
    return x.num * y.den == x.den * y.num;
  }
  friend constexpr bool operator!=(const rational &x, const rational &y) {
    return !(x == y);
  }
  friend constexpr bool operator<(const rational &x, const rational &y) {
    return x.num * y.den < x.den * y.num;
  }
  friend constexpr bool operator>(const rational &x, const rational &y) {
    return y < x;
  }
  friend constexpr bool operator<=(const rational &x, const rational &y) {
    return !(x > y);
  }
  friend constexpr bool operator>=(const rational &x, const rational &y) {
    return y <= x;
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator==(const rational &x, const J &y) {
    return x == rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator!=(const rational &x, const J &y) {
    return x != rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator<(const rational &x, const J &y) {
    return x < rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator>(const rational &x, const J &y) {
    return x > rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator<=(const rational &x, const J &y) {
    return x <= rational(y);
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator>=(const rational &x, const J &y) {
    return x >= rational(y);
  }

  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator==(const J &x, const rational &y) {
    return rational(x) == y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator!=(const J &x, const rational &y) {
    return rational(x) != y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator<(const J &x, const rational &y) {
    return rational(x) < y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator>(const J &x, const rational &y) {
    return rational(x) > y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator<=(const J &x, const rational &y) {
    return rational(x) <= y;
  }
  template <typename J, enable_if_t<is_integral_v<J> > * = nullptr>
  friend constexpr bool operator>=(const J &x, const rational &y) {
    return rational(x) >= y;
  }

  friend ostream &operator<<(ostream &os, const rational &x) {
    return os << "(" << x.num << "/" << x.den << ")";
  }

#ifdef HAVE_CAPABILITY_yaml_cpp
  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const rational &x) {
    yaml << YAML::LocalTag("rational-1.0.0");
    yaml << YAML::Flow << YAML::BeginSeq << x.num << x.den << YAML::EndSeq;
    return yaml;
  }
#endif
};

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_RATIONAL_HXX
