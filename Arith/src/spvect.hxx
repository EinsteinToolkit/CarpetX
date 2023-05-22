#ifndef CARPETX_ARITH_SPVECT_HXX
#define CARPETX_ARITH_SPVECT_HXX

#include "cons.hxx"
#include "defs.hxx"

#ifdef __CUDACC__
#include <thrust/sort.h>
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace Arith {

template <typename T> using vector1 = std::vector<T>;

template <typename I, typename T, template <typename> typename V = vector1>
class spvect {
  mutable V<cons<I, T> > elts;
  mutable bool sorted;

  template <typename F, typename R = std::invoke_result_t<F, T> >
  friend constexpr ARITH_DEVICE ARITH_HOST spvect<I, R, V>
  fmap(const F &f, const spvect &x) {
    spvect<I, R, V> r;
    for (const auto &elt : x.elts) {
      const auto &i = elt.first;
      const auto &a = elt.second;
      r.elts.emplace_back(i, f(a));
    }
    return r;
  }

  template <typename F, typename U, typename R = std::invoke_result_t<F, T, U> >
  friend constexpr ARITH_DEVICE ARITH_HOST spvect<I, R, V>
  fmap(const F &f, const spvect &x, const spvect<I, U, V> &y) {
    spvect<I, R, V> r;
    x.make_sorted();
    y.make_sorted();
    auto xi = x.elts.begin();
    auto yi = y.elts.begin();
    while (xi != x.elts.end() || yi != y.elts.end()) {
      const bool xvalid = xi != x.elts.end();
      const bool yvalid = yi != y.elts.end();
      if (xvalid && yvalid && std::equal_to<I>()(xi->first, yi->first)) {
        r.emplace_back(xi->first, f(xi->second, yi->second));
        ++xi;
        ++yi;
      } else if (xvalid && (!yvalid || std::less<I>()(xi->first, yi->first))) {
        r.emplace_back(xi->first, f(xi->second, U(0)));
        ++xi;
      } else {
        assert(yvalid && (!xvalid || std::less<I>()(yi->first, xi->first)));
        r.emplace_back(yi->first, f(T(0), yi->second));
        ++yi;
      }
    }
    return r;
  }

  template <typename Op, typename R,
            typename = std::enable_if_t<
                std::is_same_v<R, std::invoke_result_t<Op, R, T> > > >
  friend constexpr ARITH_DEVICE ARITH_HOST R fold(const Op &op, R r,
                                                  const spvect &x) {
    for (const auto &elt : x.elts) {
      const auto &a = elt.second;
      r = op(r, a);
    }
    return r;
  }

  template <typename F, typename Op, typename R,
            typename S = std::invoke_result_t<F, T>,
            typename = std::enable_if_t<
                std::is_same_v<R, std::invoke_result_t<Op, R, S> > > >
  friend constexpr ARITH_DEVICE ARITH_HOST R foldmap(const F &f, const Op &op,
                                                     R r, const spvect &x) {
    for (const auto &elt : x.elts) {
      const auto &a = elt.second;
      r = op(r, f(a));
    }
    return r;
  }

  template <typename F, typename Op, typename R, typename U,
            typename S = std::invoke_result_t<F, T, U>,
            typename = std::enable_if_t<
                std::is_same_v<R, std::invoke_result_t<Op, R, S> > > >
  friend constexpr ARITH_DEVICE ARITH_HOST R foldmap(const F &f, const Op &op,
                                                     R r, const spvect &x,
                                                     const spvect<I, U, V> &y) {
    x.make_sorted();
    y.make_sorted();
    auto xi = x.elts.begin();
    auto yi = y.elts.begin();
    while (xi != x.elts.end() || yi != y.elts.end()) {
      const bool xvalid = xi != x.elts.end();
      const bool yvalid = yi != y.elts.end();
      if (xvalid && yvalid && xi->first == yi->first) {
        r = op(r, f(xi->second, yi->second));
        ++xi;
        ++yi;
      } else if (xvalid && (!yvalid || xi->first < yi->first)) {
        r = op(r, f(xi->second, U(0)));
        ++xi;
      } else {
        assert(yvalid && (!xvalid || xi->first > yi->first));
        r = op(r, f(T(0), yi->second));
        ++yi;
      }
    }
    return r;
  }

public:
  constexpr spvect(const spvect &) = default;
  constexpr spvect(spvect &&) = default;
  constexpr spvect &operator=(const spvect &) = default;
  constexpr spvect &operator=(spvect &&) = default;

  constexpr ARITH_DEVICE ARITH_HOST spvect() : sorted(true) {}

  inline ARITH_DEVICE ARITH_HOST void make_sorted() const {
    if (!sorted)
      do_make_sorted();
  }

private:
  CCTK_ATTRIBUTE_NOINLINE ARITH_DEVICE ARITH_HOST void do_make_sorted() const {
    // sort
#ifndef __CUDACC__
    std::sort(elts.begin(), elts.end(), [](const auto &elt1, const auto &elt2) {
      return std::less<I>()(elt1.first, elt2.first);
    });
#else
    thrust::sort(elts.begin(), elts.end(),
                 [](const auto &elt1, const auto &elt2) {
                   return std::less<I>()(elt1.first, elt2.first);
                 });
#endif
    // combine elements
    std::size_t iwrite = 0, iread = 0;
    while (iread < elts.size()) {
      if (iread > iwrite)
        elts[iwrite] = std::move(elts[iread++]);
      while (iread < elts.size() &&
             std::equal_to<I>()(elts[iread].first, elts[iwrite].first))
        elts[iwrite].second += std::move(elts[iread++]).second;
      ++iwrite;
    }
    elts.resize(iwrite);
    sorted = true;
  }

public:
  ARITH_DEVICE ARITH_HOST void emplace_back(const I &idx, const T &val) {
    if (empty()) {
      elts.emplace_back(idx, val);
    } else {
      const auto &lastidx = elts.back().first;
      if (std::equal_to<I>()(idx, lastidx)) {
        // add
        elts.back().second += val;
      } else {
        // emplace
        sorted &= std::less<I>()(lastidx, idx);
        elts.emplace_back(idx, val);
      }
    }
  }

  constexpr ARITH_DEVICE ARITH_HOST bool empty() const { return elts.empty(); }
  ARITH_DEVICE ARITH_HOST std::size_t size() const {
    make_sorted();
    return elts.size();
  }

  ARITH_DEVICE ARITH_HOST std::size_t count(const I &idx) const {
    make_sorted();
    const auto &pos = std::lower_bound(
        elts.begin(), elts.end(), idx,
        [](const auto &elt, const auto &val) { return elt.first < val; });
    return pos != elts.end() && pos->first == idx;
  }

  ARITH_DEVICE ARITH_HOST const T &operator[](const I &idx) const {
    make_sorted();
    const auto &pos = std::lower_bound(
        elts.begin(), elts.end(), idx,
        [](const auto &elt, const auto &val) { return elt.first < val; });
    assert(pos != elts.end() && pos->first == idx);
    return pos->second;
  }
  ARITH_DEVICE ARITH_HOST T &operator[](const I &idx) {
    make_sorted();
    const auto &pos = std::lower_bound(
        elts.begin(), elts.end(), idx,
        [](const auto &elt, const auto &val) { return elt.first < val; });
    assert(pos != elts.end() && pos->first == idx);
    return pos->second;
  }

  ARITH_DEVICE ARITH_HOST auto begin() const {
    make_sorted();
    return elts.begin();
  }
  ARITH_DEVICE ARITH_HOST auto end() const {
    make_sorted();
    return elts.end();
  }

  template <typename = T>
  friend constexpr ARITH_DEVICE ARITH_HOST spvect operator+(const spvect &x) {
    return fmap([](const auto &a) { return +a; }, x);
  }
  template <typename = T>
  friend constexpr ARITH_DEVICE ARITH_HOST spvect operator-(const spvect &x) {
    return fmap([](const auto &a) { return -a; }, x);
  }

  friend ARITH_DEVICE ARITH_HOST spvect operator+(const spvect &x,
                                                  const spvect &y) {
    return fmap([](const auto &a, const auto &b) { return a + b; }, x, y);
  }
  friend ARITH_DEVICE ARITH_HOST spvect operator-(const spvect &x,
                                                  const spvect &y) {
    return fmap([](const auto &a, const auto &b) { return a - b; }, x, y);
  }

  friend ARITH_DEVICE ARITH_HOST spvect operator*(const T &a, const spvect &y) {
    return fmap([&a](const auto &b) { return a * b; }, y);
  }
  friend ARITH_DEVICE ARITH_HOST spvect operator*(const spvect &x, const T &b) {
    return fmap([&b](const auto &a) { return a * b; }, x);
  }
  friend ARITH_DEVICE ARITH_HOST spvect operator/(const spvect &x, const T &b) {
    return fmap([&b](const auto &a) { return a / b; }, x);
  }
  template <typename T1 = T,
            typename = std::enable_if_t<std::is_integral_v<T1> > >
  friend ARITH_DEVICE ARITH_HOST spvect operator%(const spvect &x, const T &b) {
    return fmap([&b](const auto &a) { return a % b; }, x);
  }

  ARITH_DEVICE ARITH_HOST spvect &operator+=(const spvect &x) {
    return *this = *this + x;
  }
  ARITH_DEVICE ARITH_HOST spvect &operator-=(const spvect &x) {
    return *this = *this - x;
  }
  ARITH_DEVICE ARITH_HOST spvect &operator*=(const T &a) {
    return *this = *this * a;
  }
  ARITH_DEVICE ARITH_HOST spvect &operator/=(const T &a) {
    return *this = *this / a;
  }
  template <typename T1 = T,
            typename = std::enable_if_t<std::is_integral_v<T1> > >
  ARITH_DEVICE ARITH_HOST spvect &operator%=(const T &a) {
    return *this = *this % a;
  }

  friend ARITH_DEVICE ARITH_HOST spvect abs(const spvect &x) {
    using std::abs;
    return fmap([](const auto &a) { return abs(a); }, x);
  }

  friend ARITH_DEVICE ARITH_HOST bool iszero(const spvect &x) {
    return foldmap([](const auto &a) { return a == 0; },
                   [](const bool r, const bool a) { return r && a; }, true, x);
  }
  // nvcc @11.2.2 does not find this function when it is a friend
  // friend ARITH_DEVICE ARITH_HOST T maximum(const spvect &x) {
  //   return fold(
  //       [] ARITH_DEVICE ARITH_HOST(const T r, const T a) { return max1(r, a);
  //       }, -std::numeric_limits<T>::infinity(), x);
  // }
  friend ARITH_DEVICE ARITH_HOST T minimum(const spvect &x) {
    return fold([](const T r, const T a) { return min1(r, a); },
                std::numeric_limits<T>::infinity(), x);
  }
  friend ARITH_DEVICE ARITH_HOST T prod(const spvect &x) {
    return fold([](const T r, const T a) { return r * a; }, T(1), x);
  }
  friend ARITH_DEVICE ARITH_HOST T sum(const spvect &x) {
    return fold([](const T r, const T a) { return r + a; }, T(0), x);
  }

  friend ARITH_DEVICE ARITH_HOST bool operator==(const spvect &x,
                                                 const spvect &y) {
    return foldmap([](const auto &a, const auto &b) { return a == b; },
                   [](const bool r, const bool a) { return r && a; }, true, x,
                   y);
  }
  friend ARITH_DEVICE ARITH_HOST bool operator!=(const spvect &x,
                                                 const spvect &y) {
    return !(x == y);
  }

  friend std::ostream &operator<<(std::ostream &os, const spvect &x) {
    x.make_sorted();
    os << "⠎⠧[";
    for (const auto &elt : x.elts) {
      os << elt.first << ":" << elt.second << ",";
    }
    os << "]";
    return os;
  }
};

template <typename I, typename T, template <typename> typename V = vector1>
ARITH_DEVICE ARITH_HOST T maximum(const spvect<I, T, V> &x) {
  return fold(
      [] ARITH_DEVICE ARITH_HOST(const T r, const T a) { return max1(r, a); },
      -std::numeric_limits<T>::infinity(), x);
}

template <typename I, typename T, template <typename> typename V>
struct zero<spvect<I, T, V> > {
  typedef spvect<I, T, V> value_type;
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return spvect<I, T, V>();
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return spvect<I, T, V>();
  }
};

template <typename I, typename T, template <typename> typename V>
struct nan<spvect<I, T, V> > {
  typedef spvect<I, T, V> value_type;
  // static constexpr value_type value = spvect<I,T,V>::pure(nan_v<T>);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return spvect<I, T, V>::pure(nan<T>()());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return spvect<I, T, V>::pure(nan<T>()());
  }
};

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_SPVECT_HXX
