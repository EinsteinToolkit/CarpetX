#ifndef CARPETX_ARITH_SMALLVECTOR_HXX
#define CARPETX_ARITH_SMALLVECTOR_HXX

#include "arr.hxx"
#include "defs.hxx"
#include "vect.hxx"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <initializer_list>
#include <iostream>
#include <limits>

#ifdef CCTK_DEBUG
#define debug_assert(arg) assert(arg)
#else
#define debug_assert(arg) ((void)0)
#endif

namespace Arith {

template <typename T, std::size_t N> class smallvector {
  arr<T, N> elts;
  std::size_t sz;

public:
  using value_type = typename arr<T, N>::value_type;
  using size_type = typename arr<T, N>::size_type;
  using difference_type = typename arr<T, N>::difference_type;
  using reference = typename arr<T, N>::reference;
  using const_reference = typename arr<T, N>::const_reference;
  using pointer = typename arr<T, N>::pointer;
  using const_pointer = typename arr<T, N>::const_pointer;
  using iterator = typename arr<T, N>::iterator;
  using const_iterator = typename arr<T, N>::const_iterator;
  using reverse_iterator = typename arr<T, N>::reverse_iterator;
  using const_reverse_iterator = typename arr<T, N>::const_reverse_iterator;

  constexpr smallvector(const smallvector &) = default;
  constexpr smallvector(smallvector &&) = default;
  constexpr smallvector &operator=(const smallvector &) = default;
  constexpr smallvector &operator=(smallvector &&) = default;

  constexpr ARITH_DEVICE ARITH_HOST smallvector() : elts{}, sz(0) {}

  // constexpr smallvector(const std::initializer_list<T> &lst)
  //     : elts(arr<T, N>(vect<T, N>(lst))), sz(lst.size()) {}
  constexpr ARITH_DEVICE ARITH_HOST
  smallvector(const std::initializer_list<T> &lst)
      : elts(construct_array<T, N>([&lst](const size_type i) {
          return i < lst.size() ? lst.begin()[i] : T{};
        })),
        sz(lst.size()) {}
  template <size_type N1>
  constexpr ARITH_DEVICE ARITH_HOST smallvector(const arr<T, N1> &elts_)
      : elts(construct_array<T, N>(
            [&elts_](const size_type i) { return i < N1 ? elts_[i] : T{}; })),
        sz(N1) {}
  template <size_type N1>
  constexpr ARITH_DEVICE ARITH_HOST smallvector(const std::array<T, N1> &elts_)
      : elts(construct_array<T, N>(
            [&elts_](const size_type i) { return i < N1 ? elts_[i] : T{}; })),
        sz(N1) {}

  constexpr ARITH_DEVICE ARITH_HOST bool empty() const { return sz == 0; }
  constexpr ARITH_DEVICE ARITH_HOST size_type size() const { return sz; }
  constexpr ARITH_DEVICE ARITH_HOST size_type max_size() const { return N; }

  ARITH_DEVICE ARITH_HOST void resize(const size_type newsz) {
    debug_assert(newsz <= N);
    sz = newsz;
  }

  ARITH_DEVICE ARITH_HOST const_reference at(const size_type i) const {
    if (CCTK_BUILTIN_EXPECT(i >= sz, false))
      throw std::out_of_range("smallvector::at");
    return elts[i];
  }
  ARITH_DEVICE ARITH_HOST reference at(const size_type i) {
    if (CCTK_BUILTIN_EXPECT(i >= sz, false))
      throw std::out_of_range("smallvector::at");
    return elts[i];
  }
  constexpr ARITH_DEVICE ARITH_HOST const_reference
  operator[](const size_type i) const {
    debug_assert(i < sz);
    return elts[i];
  }
  constexpr ARITH_DEVICE ARITH_HOST reference operator[](const size_type i) {
    debug_assert(i < sz);
    return elts[i];
  }

  constexpr ARITH_DEVICE ARITH_HOST reference front() {
    debug_assert(sz > 0);
    return elts.front();
  }
  constexpr ARITH_DEVICE ARITH_HOST const_reference front() const {
    debug_assert(sz > 0);
    return elts.front();
  }
  constexpr ARITH_DEVICE ARITH_HOST reference back() {
    debug_assert(sz > 0);
    return elts.back();
  }
  constexpr ARITH_DEVICE ARITH_HOST const_reference back() const {
    debug_assert(sz > 0);
    return elts.back();
  }

  ARITH_DEVICE ARITH_HOST const T *data() const { return elts.data(); }
  ARITH_DEVICE ARITH_HOST T *data() { return elts.data(); }

  ARITH_DEVICE ARITH_HOST iterator begin() { return elts.begin(); }
  ARITH_DEVICE ARITH_HOST const_iterator begin() const { return elts.begin(); }
  ARITH_DEVICE ARITH_HOST const_iterator cbegin() const {
    return elts.cbegin();
  }
  ARITH_DEVICE ARITH_HOST iterator end() { return elts.end(); }
  ARITH_DEVICE ARITH_HOST const_iterator end() const { return elts.end(); }
  ARITH_DEVICE ARITH_HOST const_iterator cend() const { return elts.cend(); }

  // constexpr const T *rbegin() const { return &elts[0]; }
  // constexpr T *rbegin() { return &elts[0]; }
  // constexpr const T *rend() const { return &elts[sz]; }
  // constexpr T *rend() { return &elts[sz]; }

  ARITH_DEVICE ARITH_HOST void clear() { sz = 0; }

  ARITH_DEVICE ARITH_HOST void push_back(const T &x) {
    debug_assert(sz < N);
    elts[sz++] = x;
  }
  ARITH_DEVICE ARITH_HOST void push_back(T &&x) {
    debug_assert(sz < N);
    elts[sz++] = x;
  }
  ARITH_DEVICE ARITH_HOST void pop_back() {
    debug_assert(sz > 0);
    --sz;
  }

  template <typename... Args>
  ARITH_DEVICE ARITH_HOST reference emplace_back(Args &&...args) {
    debug_assert(sz < N);
    return elts[sz++] = T(std::forward<Args>(args)...);
  }

  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator+(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = +x.elts[i];
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator-(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = -x.elts[i];
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator+(const smallvector &x, const smallvector &y) {
    smallvector r;
    debug_assert(x.sz == y.sz);
    using std::min;
    r.sz = min(x.sz, y.sz);
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] + y.elts[i];
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator-(const smallvector &x, const smallvector &y) {
    smallvector r;
    debug_assert(x.sz == y.sz);
    using std::min;
    r.sz = min(x.sz, y.sz);
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] - y.elts[i];
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator*(const T &a, const smallvector &y) {
    smallvector r;
    r.sz = y.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = a * y.elts[i];
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator*(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] * b;
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator/(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] / b;
    return r;
  }
  template <typename T1 = T,
            typename = std::enable_if_t<std::is_integral_v<T1> > >
  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  operator%(const smallvector &x, const T &b) {
    smallvector r;
    r.sz = x.sz;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = x.elts[i] % b;
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST smallvector
  abs(const smallvector &x) {
    smallvector r;
    r.sz = x.sz;
    using std::abs;
    for (size_type i = 0; i < r.sz; ++i)
      r.elts[i] = abs(x.elts[i]);
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST T maximum(const smallvector &x) {
    T r = -std::numeric_limits<T>::infinity();
    for (size_type i = 0; i < x.sz; ++i)
      r = max1(r, x.elts[i]);
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST T minimum(const smallvector &x) {
    T r = std::numeric_limits<T>::infinity();
    for (size_type i = 0; i < x.sz; ++i)
      r = min1(r, x.elts[i]);
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST T prod(const smallvector &x) {
    T r = 1;
    for (size_type i = 0; i < x.sz; ++i)
      r *= x.elts[i];
    return r;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST T sum(const smallvector &x) {
    T r = 0;
    for (size_type i = 0; i < x.sz; ++i)
      r += x.elts[i];
    return r;
  }

  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator==(const smallvector &x, const smallvector &y) {
    if (x.sz != y.sz)
      return false;
    for (size_type i = 0; i < x.sz; ++i)
      if (x.elts[i] != y.elts[i])
        return false;
    return true;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator!=(const smallvector &x, const smallvector &y) {
    return !(x == y);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator<(const smallvector &x, const smallvector &y) {
    using std::min;
    for (size_type i = 0; i < min(x.sz, y.sz); ++i) {
      if (x.elts[i] < y.elts[i])
        return true;
      if (y.elts[i] < x.elts[i])
        return false;
    }
    return x.sz < y.sz;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator>(const smallvector &x, const smallvector &y) {
    return y < x;
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator<=(const smallvector &x, const smallvector &y) {
    return !(y < x);
  }
  friend constexpr ARITH_DEVICE ARITH_HOST bool
  operator>=(const smallvector &x, const smallvector &y) {
    return !(x < y);
  }

  friend std::ostream &operator<<(std::ostream &os, const smallvector &x) {
    os << "[";
    for (const auto &a : x)
      os << a << ",";
    os << "]";
    return os;
  }
};

template <typename T, std::size_t N> struct zero<smallvector<T, N> > {
  typedef smallvector<T, N> value_type;
  // static constexpr value_type value = smallvector<T, N>::pure(zero_v<T>);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return smallvector<T, N>::pure(zero<T>()());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return smallvector<T, N>::pure(zero<T>()());
  }
};

template <typename T, std::size_t N> struct nan<smallvector<T, N> > {
  typedef smallvector<T, N> value_type;
  // static constexpr value_type value = smallvector<T, N>::pure(nan_v<T>);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return smallvector<T, N>::pure(nan<T>()());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return smallvector<T, N>::pure(nan<T>()());
  }
};

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_SMALLVECTOR_HXX
