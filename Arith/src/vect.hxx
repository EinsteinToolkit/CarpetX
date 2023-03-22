#ifndef VECT_HXX
#define VECT_HXX

#include "defs.hxx"
#include "div.hxx"
#include "simd.hxx"
#include "tuple.hxx"

#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace Arith {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

#if 0
// This does not work with thrust::tuple

#if 0
namespace detail {
template <typename T, std::size_t... I>
decltype(auto) make_tuple_type(const T &t, std::index_sequence<I...>) {
  // return std_make_tuple((I, t)...);
  return std_make_tuple(std::enable_if_t<(I || true), T>(t)...);
}
} // namespace detail

template <typename T, size_t N> struct ntuple {
  using type = decltype(detail::make_tuple_type(declval<T>(),
                                                std::make_index_sequence<N>{}));
};
#endif
template <typename T, size_t N> struct ntuple {
  using type = decltype(std_tuple_cat(
      declval<std_tuple<T> >(), declval<typename ntuple<T, N - 1>::type>()));
};
template <typename T> struct ntuple<T, 0> { using type = std_tuple<>; };
template <typename T, size_t N> using ntuple_t = typename ntuple<T, N>::type;

#endif

////////////////////////////////////////////////////////////////////////////////

namespace detail {
template <typename T, class Tuple, std::size_t... Is>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_from_tuple(const Tuple &t, std::index_sequence<Is...>) {
  return std::array<T, sizeof...(Is)>{std::get<Is>(t)...};
}
template <typename T, class Tuple, std::size_t... Is>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_from_tuple(Tuple &&t, std::index_sequence<Is...>) {
  return std::array<T, sizeof...(Is)>{std::move(std::get<Is>(t))...};
}

template <typename T, class Tuple, std::size_t... Is, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_push(const Tuple &t, std::index_sequence<Is...>, const U &x) {
  return std::array<T, sizeof...(Is) + 1>{std::get<Is>(t)..., T(x)};
}
} // namespace detail

template <typename T, std::size_t N, typename Tuple>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST std::array<T, N>
array_from_tuple(const Tuple &t) {
  return detail::array_from_tuple<T>(t, std::make_index_sequence<N>());
}
template <typename T, std::size_t N, typename Tuple>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST std::array<T, N>
array_from_tuple(Tuple &&t) {
  return detail::array_from_tuple<T>(std::move(t),
                                     std::make_index_sequence<N>());
}

template <class T, std::size_t N, typename U>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto
array_push(const std::array<T, N> &a, const U &e) {
  return detail::array_push<T>(a, std::make_index_sequence<N>(), e);
}

template <typename T, std::size_t N, typename F>
constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST std::array<T, N>
construct_array(const F &f) {
  if constexpr (N == 0)
    return std::array<T, N>();
  if constexpr (N > 0)
    return array_push<T>(construct_array<T, N - 1>(f), f(N - 1));
}

////////////////////////////////////////////////////////////////////////////////

// A small vector with a length that is known at compile time, similar
// to `std::array`. The main difference is that `vect` supports
// arithmetic operations, which is most useful for multi-dimensional
// array indices.

template <typename T, int D> struct vect {
  array<T, D> elts;

  typedef T value_type;
  typedef int size_type;
  static constexpr int size() { return D; }

  // (no it doesn't) initializes all elements to zero
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect() : elts() {}

  constexpr ARITH_INLINE vect(const vect &) = default;
  constexpr ARITH_INLINE vect(vect &&) = default;
  constexpr ARITH_INLINE vect &operator=(const vect &) = default;
  constexpr ARITH_INLINE vect &operator=(vect &&) = default;

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(vect<U, D> x)
      : elts(construct_array<T, D>(
            [&](int d) ARITH_INLINE { return std::move(x[d]); })) {}

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(array<T, D> arr)
      : elts(std::move(arr)) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(const T (&arr)[D])
      : elts(construct_array<T, D>([&](int d)
                                       ARITH_INLINE { return arr[d]; })) {}
#if 0
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(ntuple_t<T, D> tup)
      : elts(array_from_tuple<T, D>(std::move(tup))) {}
#endif
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(initializer_list<T> lst)
      : elts(construct_array<T, D>([&](size_t d) ARITH_INLINE {
#ifdef CCTK_DEBUG
          assert(d < lst.size());
#endif
          return lst.begin()[d];
        })) {
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(const vector<T> &vec)
      : elts(construct_array<T, D>([&](size_t d) ARITH_INLINE {
#ifdef CCTK_DEBUG
          return vec.at(d);
#else
          return vec[d];
#endif
        })) {
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect(vector<T> &&vec)
      : elts(construct_array<T, D>([&](size_t d) ARITH_INLINE {
#ifdef CCTK_DEBUG
          return std::move(vec.at(d));
#else
          return std::move(vec[d]);
#endif
        })) {
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator std::array<T, D>() const {
    return elts;
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
  operator[](int d) const {
#ifdef CCTK_DEBUG
    assert(d >= 0 && d < D);
#endif
    return elts[d];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &operator[](int d) {
#ifdef CCTK_DEBUG
    assert(d >= 0 && d < D);
#endif
    return elts[d];
  }

  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void loop(const F &f) {
    for (int d = 0; d < D; ++d)
      f(d);
  }
  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect make(const F &f) {
    return construct_array<T, D>(f);
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect pure(T a) {
    return make([&](int) ARITH_INLINE { return a; });
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect unit(int dir) {
    return make([&](int d) ARITH_INLINE { return d == dir; });
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<int, D> iota() {
    return vect<int, D>::make([&](int d) ARITH_INLINE { return d; });
  }

  template <typename F, typename... Args,
            typename R =
                remove_cv_t<remove_reference_t<result_of_t<F(T, Args...)> > > >
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  fmap(const F &f, const vect &x, const vect<Args, D> &...args) {
    return vect<R, D>::make(
        [&](int d) ARITH_INLINE { return f(x.elts[d], args.elts[d]...); });
  }
  template <typename F, typename... Args>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  fmap_(const F &f, const vect &x, const vect<Args, D> &...args) {
    loop([&](int d) ARITH_INLINE { f(x.elts[d], args.elts[d]...); });
  }

  template <typename Op, typename R, typename... Args,
            typename = remove_cv_t<
                remove_reference_t<result_of_t<Op(R, T, Args...)> > > >
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST R
  fold(const Op &op, R r, const vect &x, const vect<Args, D> &...args) {
    loop([&](int d) ARITH_INLINE { r = op(r, x[d], args[d]...); });
    return r;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  reversed(const vect &x) {
    return vect::make([&](int d) ARITH_INLINE { return x.elts[D - 1 - d]; });
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator+(const vect &x) {
    return fmap([](const T &a) ARITH_INLINE { return +a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator-(const vect &x) {
    return fmap([](const T &a) ARITH_INLINE { return -a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator~(const vect &x) {
    return fmap([](const T &a) ARITH_INLINE { return ~a; }, x);
  }

  template <typename U,
            typename R = decltype(std::declval<T>() + std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator+(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a + b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() - std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator-(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a - b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() * std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator*(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a * b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() / std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator/(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a / b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() % std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator%(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a % b; }, x,
                y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  div_floor(const vect &x, const vect &y) {
    return fmap([](const T &a, const T &b)
                    ARITH_INLINE { return div_floor(a, b); },
                x, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  mod_floor(const vect &x, const vect &y) {
    return fmap([](const T &a, const T &b)
                    ARITH_INLINE { return mod_floor(a, b); },
                x, y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() & std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator&(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a & b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() | std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator|(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a | b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() ^ std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator^(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a ^ b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() << std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator<<(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a << b; }, x,
                y);
  }
  template <typename U,
            typename R = decltype(std::declval<T>() >> std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<R, D>
  operator>>(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const U &b) ARITH_INLINE { return a >> b; }, x,
                y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator+(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a + b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator-(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a - b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator*(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a * b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator/(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a / b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator%(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a % b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  div_floor(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return div_floor(a, b); }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  mod_floor(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return mod_floor(a, b); }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator&(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a & b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator|(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a | b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator^(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a ^ b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator<<(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a << b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator>>(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a >> b; }, x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator+(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b + a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator-(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b - a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator*(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b * a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator/(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b / a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator%(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b % a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  div_floor(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return div_floor(b, a); }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  mod_floor(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return mod_floor(b, a); }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator&(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b & a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator|(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b | a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator^(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b ^ a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator<<(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b << a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator>>(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b >> a; }, x);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator+=(const vect &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator-=(const vect &x) {
    return *this = *this - x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator+=(const vect<U, D> &x) {
    return *this = *this + x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator-=(const vect<U, D> &x) {
    return *this = *this - x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator*=(const vect<U, D> &x) {
    return *this = *this * x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator/=(const vect<U, D> &x) {
    return *this = *this / x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator%=(const vect<U, D> &x) {
    return *this = *this % x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator&=(const vect &x) {
    return *this = *this & x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator|=(const vect &x) {
    return *this = *this | x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator^=(const vect<U, D> &x) {
    return *this = *this ^ x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator<<=(const vect<U, D> &x) {
    return *this = *this << x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  operator>>=(const vect<U, D> &x) {
    return *this = *this >> x;
  }

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator+=(const U &a) {
    return *this = *this + a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator-=(const U &a) {
    return *this = *this - a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator*=(const U &a) {
    return *this = *this * a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator/=(const U &a) {
    return *this = *this / a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator%=(const U &a) {
    return *this = *this % a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator&=(const U &a) {
    return *this = *this & a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator|=(const U &a) {
    return *this = *this | a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator^=(const U &a) {
    return *this = *this ^ a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator<<=(const U &a) {
    return *this = *this << a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect operator>>=(const U &a) {
    return *this = *this >> a;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator!(const vect &x) {
    return fmap([](const T &a) ARITH_INLINE { return !a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator&&(const vect &x, const vect &y) {
    return fmap([](const T &a, const T &b) ARITH_INLINE { return a && b; }, x,
                y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator||(const vect &x, const vect &y) {
    return fmap([](const T &a, const T &b) ARITH_INLINE { return a || b; }, x,
                y);
  }

  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator==(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const T &b) ARITH_INLINE { return a == b; }, x,
                y);
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator!=(const vect &x, const vect<U, D> &y) {
    return !(x == y);
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<(const vect &x, const vect<U, D> &y) {
    return fmap([](const T &a, const T &b) ARITH_INLINE { return a < b; }, x,
                y);
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>(const vect &x, const vect<U, D> &y) {
    return y < x;
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<=(const vect &x, const vect<U, D> &y) {
    return !(x > y);
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>=(const vect &x, const vect<U, D> &y) {
    return !(x < y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator==(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a == b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator!=(const T &a, const vect &x) {
    return !(x == a);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<(const T &a, const vect &x) {
    return fmap([&](const T &b) ARITH_INLINE { return a < b; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>(const T &a, const vect &x) {
    return a < x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<=(const T &a, const vect &x) {
    return !(x > a);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>=(const T &a, const vect &x) {
    return !(x < a);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator==(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b == a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator!=(const vect &x, const T &a) {
    return !(x == a);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<(const vect &x, const T &a) {
    return fmap([&](const T &b) ARITH_INLINE { return b < a; }, x);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>(const vect &x, const T &a) {
    return a < x;
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator<=(const vect &x, const T &a) {
    return !(x > a);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect<bool, D>
  operator>=(const vect &x, const T &a) {
    return !(x < a);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  abs(const vect &x) {
    using std::abs;
    return fmap([](const T &a) ARITH_INLINE { return abs(a); }, x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  all(const vect &x) {
    return fold([](const T &a, const T &b) ARITH_INLINE { return a && b; },
                one<T>()(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  any(const vect &x) {
    return fold([](const T &a, const T &b) ARITH_INLINE { return a || b; },
                zero<T>()(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  allisfinite(const vect &x) {
    return all(
        fmap([](const auto &a) ARITH_INLINE { return allisfinite(a); }, x));
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  anyisnan(const vect &x) {
    return any(fmap([](const auto &a) ARITH_INLINE { return anyisnan(a); }, x));
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T dot(const vect &x,
                                                              const vect &y) {
    return sum(x * y);
  }

  template <typename C>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  if_else(const vect<C, D> &cond, const vect &x, const vect &y) {
    return fmap([](const C &c, const T &x, const T &y)
                    ARITH_INLINE { return if_else(c, x, y); },
                cond, x, y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*vect<bool, D>*/
  isnan(const vect &x) {
    using std::isnan;
    return fmap([](const T &a) ARITH_INLINE { return isnan(a); }, x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  max(const vect &x, const vect &y) {
    using std::max;
    return fmap([](const T &x, const T &y) ARITH_INLINE { return max(x, y); },
                x, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  max(const T &a, const vect &y) {
    using std::max;
    return fmap([&](const T &y) ARITH_INLINE { return max(a, y); }, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect max(const vect &x,
                                                                 const T &a) {
    using std::max;
    return fmap([&](const T &x) ARITH_INLINE { return max(x, a); }, x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  max1(const vect &x, const vect &y) {
    return fmap([](const T &x, const T &y) ARITH_INLINE { return max1(x, y); },
                x, y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  maxabs(const vect &x) {
    using std::abs;
    return fold([](const T &r, const T &a)
                    ARITH_INLINE { return max1(r, abs(a)); },
                zero<T>()(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  maximum(const vect &x) {
    return fold([&](const T &a, const T &b) ARITH_INLINE { return max1(a, b); },
                -numeric_limits<T>::infinity(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  min(const vect &x, const vect &y) {
    using std::min;
    return fmap([](const T &x, const T &y) ARITH_INLINE { return min(x, y); },
                x, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  min(const T &a, const vect &y) {
    using std::min;
    return fmap([&](const T &y) ARITH_INLINE { return min(a, y); }, y);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect min(const vect &x,
                                                                 const T &a) {
    using std::min;
    return fmap([&](const T &x) ARITH_INLINE { return min(x, a); }, x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vect
  min1(const vect &x, const vect &y) {
    return fmap([](const T &x, const T &y) ARITH_INLINE { return min1(x, y); },
                x, y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  minimum(const vect &x) {
    return fold(min1, numeric_limits<T>::infinity(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T prod(const vect &x) {
    return fold([](const T &a, const T &b) ARITH_INLINE { return a * b; },
                one<T>()(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T sum(const vect &x) {
    return fold([](const T &a, const T &b) ARITH_INLINE { return a + b; },
                zero<T>()(), x);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  sumabs(const vect &x) {
    using std::abs;
    return fold([](const T &a, const T &b) ARITH_INLINE { return a + abs(b); },
                zero<T>()(), x);
  }

  friend ostream &operator<<(ostream &os, const vect &x) {
    os << "[";
    for (int d = 0; d < D; ++d) {
      if (d > 0)
        os << ",";
      os << x.elts[d];
    }
    os << "]";
    return os;
  }
};

template <typename T, int D> struct zero<vect<T, D> > {
  typedef vect<T, D> value_type;
  // static constexpr value_type value = vect<T, D>::pure(zero_v<T>);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vect<T, D>::pure(zero<T>()());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vect<T, D>::pure(zero<T>()());
  }
};

template <typename T, int D> struct nan<vect<T, D> > {
  typedef vect<T, D> value_type;
  // static constexpr value_type value = vect<T, D>::pure(nan_v<T>);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vect<T, D>::pure(nan<T>()());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

} // namespace Arith
namespace std {
template <typename T, int D> struct equal_to<Arith::vect<T, D> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::vect<T, D> &lhs, const Arith::vect<T, D> &rhs) const {
    // This is not yet constexpr in C++17
    // return equal_to<array<T, D> >()(lhs.elts, rhs.elts);
    for (int d = 0; d < D; ++d)
      if (!equal_to<T>()(lhs.elts[d], rhs.elts[d]))
        return false;
    return true;
  }
};

template <typename T, int D> struct less<Arith::vect<T, D> > {
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  operator()(const Arith::vect<T, D> &lhs, const Arith::vect<T, D> &rhs) const {
    // return less()(lhs.elts, rhs.elts);
    for (int d = 0; d < D; ++d) {
      if (less<T>()(lhs.elts[d], rhs.elts[d]))
        return true;
      if (less<T>()(rhs.elts[d], lhs.elts[d]))
        return false;
    }
    return false;
  }
};

template <typename T, int D>
struct tuple_size<Arith::vect<T, D> >
    : integral_constant<typename Arith::vect<T, D>::size_type, D> {};

} // namespace std
namespace Arith {

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D>
constexpr vect<simd<T>, D> if_else(const simdl<T> &cond,
                                   const vect<simd<T>, D> &x,
                                   const vect<simd<T>, D> &y) {
  return fmap([&](const auto &x, const auto &y)
                  ARITH_INLINE { return if_else(cond, x, y); },
              x, y);
}

} // namespace Arith

#endif // #ifndef VECT_HXX
