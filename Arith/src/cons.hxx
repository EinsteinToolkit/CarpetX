#ifndef CARPETX_ARITH_CONS_HXX
#define CARPETX_ARITH_CONS_HXX

#include "defs.hxx"

#include <type_traits>
#include <utility>

namespace Arith {

template <typename T1, typename T2> struct cons {
  using first_type = T1;
  using second_type = T2;

  T1 first;
  T2 second;

  constexpr cons(const cons &) = default;
  constexpr cons(cons &&) = default;
  constexpr cons &operator=(const cons &) = default;
  constexpr cons &operator=(cons &&) = default;

  ARITH_DEVICE ARITH_HOST void swap(cons &c1, cons &c2) {
    using std::swap;
    swap(c1.first, c2.first);
    swap(c1.second, c2.second);
  }

  constexpr ARITH_DEVICE ARITH_HOST cons() : first{}, second{} {}

  constexpr ARITH_DEVICE ARITH_HOST cons(T1 first_, T2 second_)
      : first(std::move(first_)), second(std::move(second_)) {}

  friend constexpr ARITH_DEVICE ARITH_HOST cons make_cons(T1 first, T2 second) {
    return cons(std::move(first), std::move(second));
  }

  friend constexpr ARITH_DEVICE ARITH_HOST bool operator==(const cons &c1,
                                                           const cons &c2) {
    return c1.first == c2.first && c1.second == c2.second;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool operator!=(const cons &c1,
                                                           const cons &c2) {
    return !(c1 == c2);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool operator<(const cons &c1,
                                                          const cons &c2) {
    if (c1.first < c2.first)
      return true;
    if (c2.first < c1.first)
      return false;
    if (c1.second < c2.second)
      return true;
    return false;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool operator>(const cons &c1,
                                                          const cons &c2) {
    return c2 < c1;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool operator<=(const cons &c1,
                                                           const cons &c2) {
    return !(c1 > c2);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool operator>=(const cons &c1,
                                                           const cons &c2) {
    return !(c1 < c2);
  }
};

} // namespace Arith

namespace std {

template <size_t I, typename T1, typename T2>
constexpr ARITH_DEVICE ARITH_HOST enable_if_t<I == 0, const T1 &>
get(const Arith::cons<T1, T2> &c) {
  return c.first;
}
template <size_t I, typename T1, typename T2>
constexpr ARITH_DEVICE ARITH_HOST enable_if_t<I == 1, const T2 &>
get(const Arith::cons<T1, T2> &c) {
  return c.second;
}

template <size_t I, typename T1, typename T2>
constexpr ARITH_DEVICE ARITH_HOST enable_if_t<I == 0, T1 &>
get(Arith::cons<T1, T2> &c) {
  return c.first;
}
template <size_t I, typename T1, typename T2>
constexpr ARITH_DEVICE ARITH_HOST enable_if_t<I == 1, T2 &>
get(Arith::cons<T1, T2> &c) {
  return c.second;
}

template <typename T1, typename T2>
struct tuple_size<Arith::cons<T1, T2> > : integral_constant<size_t, 2> {};

template <typename T1, typename T2>
struct tuple_element<0, Arith::cons<T1, T2> > {
  using type = T1;
};
template <typename T1, typename T2>
struct tuple_element<1, Arith::cons<T1, T2> > {
  using type = T2;
};

} // namespace std

#endif // #ifndef CARPETX_ARITH_CONS_HXX
