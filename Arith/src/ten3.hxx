#ifndef CARPETX_ARITH_TEN3_HXX
#define CARPETX_ARITH_TEN3_HXX

#include "arr.hxx"
#include "defs.hxx"
#include "simd.hxx"
#include "tuple.hxx"
#include "vect.hxx"

#include "vec.hxx" // for symm_t

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

namespace Arith {
using namespace std;

// A rank-3 tensor with various symmetries
template <typename T, int D, symm_t symm> struct gten3 {
  constexpr static int Nfull = D * D * D;
  constexpr static int Nsymm = D * (D + 1) * (D + 2) / 6;
  constexpr static int Nanti = D * (D - 1) * (D - 2) / 6;
  constexpr static int N = if_symm<symm>(Nfull, Nsymm, Nanti);
  vect<T, N> elts;

public:
  typedef T value_type;
  typedef int size_type;
  static constexpr int size_value = N;

private:
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST arr<int, 3>
  sorted(const int i, const int j, const int k) {
    using std::max, std::min;
    if (i <= j && i <= k)
      return arr<int, 3>{i, min(j, k), max(j, k)};
    else if (j <= i && j <= k)
      return arr<int, 3>{j, min(i, k), max(i, k)};
    else
      return arr<int, 3>{k, min(i, j), max(i, j)};
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_full(const int i, const int j, const int k) {
    return i * D * D + j * D + k;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_symm(const int i0, const int j0, const int k0) {
    const int i = get<0>(sorted(i0, j0, k0));
    const int j = get<1>(sorted(i0, j0, k0));
    const int k = get<2>(sorted(i0, j0, k0));
    return (i * (i * i - 3 * D * i + 3 * D * D - 1) / 6) -
           (j * (j - 2 * D + 1) / 2) + k;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_anti(const int i0, const int j0, const int k0) {
    if (i0 == j0 || i0 == k0 || j0 == k0)
      return 0;
    const int i = get<0>(sorted(i0, j0, k0));
    const int j = get<1>(sorted(i0, j0, k0));
    const int k = get<2>(sorted(i0, j0, k0));
    return (i * (i * i + (-3 * D + 6) * i + 3 * D * D - 12 * D + 11) / 6) -
           (j * (j - 2 * D + 3) / 2) + k - D;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind(const int i, const int j, const int k) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < D && j >= 0 && j < D && k >= 0 && k < D);
#endif
    int n{};
    if constexpr (symm == FULL)
      n = ind_full(i, j, k);
    if constexpr (symm == SYMM)
      n = ind_symm(i, j, k);
    if constexpr (symm == ANTI)
      n = ind_anti(i, j, k);
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < N);
#endif
    return n;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  sign(const int i, const int j, const int k) {
    if constexpr (symm != ANTI)
      return 1;
    if (i == j || i == k || j == k)
      return 0;
    return i < j && j < k   ? +1
           : i < k && k < j ? -1
           : j < i && i < k ? -1
           : j < k && k < i ? +1
           : k < i && i < j ? +1
                            : -1;
  }

public:
  // initializes all elements to nan
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3()
      : elts(nan<vect<T, N> >()) {}

  constexpr ARITH_INLINE gten3(const gten3 &) = default;
  constexpr ARITH_INLINE gten3(gten3 &&) = default;
  constexpr ARITH_INLINE gten3 &operator=(const gten3 &) = default;
  constexpr ARITH_INLINE gten3 &operator=(gten3 &&) = default;

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(gten3<U, D, symm> x)
      : elts(std::move(x.elts)) {}

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(initializer_list<T> x)
      : elts(x) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(vect<T, N> elts)
      : elts(std::move(elts)) {}
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(arr<T, N> x)
      : elts(std::move(x)) {}
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  gten3(std::array<T, N> x)
      : elts(std::move(x)) {}
  // explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(const
  // vector<T> &x) : elts(x) {} explicit constexpr ARITH_INLINE ARITH_DEVICE
  // ARITH_HOST gten3(vector<T> &&x) : elts(std::move(x)) {}

  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator vect<T, N>() const {
    return elts;
  }

  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void loop(const F &f) {
    if constexpr (symm == symm_t::full)
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          for (int k = 0; k < D; ++k)
            f(i, j, k);
    if constexpr (symm == symm_t::symm)
      for (int i = 0; i < D; ++i)
        for (int j = i; j < D; ++j)
          for (int k = j; k < D; ++k)
            f(i, j, k);
    if constexpr (symm == symm_t::anti)
      for (int i = 0; i < D; ++i)
        for (int j = i + 1; j < D; ++j)
          for (int k = j + 1; k < D; ++k)
            f(i, j, k);
  }
  template <typename F, typename = result_of_t<F(int, int, int)> >
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3(const F &f)
      : gten3(fmap(f, iota1(), iota2(), iota3())) {}

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3 pure(const T &a) {
    return {vect<T, N>::pure(a)};
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  unit(const int i, const int j, const int k) {
    gten3 r = zero<gten3>();
    r(i, j, k) = one<T>();
    return r;
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<int, D, symm>
  iota() {
    gten3<int, D, symm> r;
    gten3<int, D, symm>::loop(
        [&](int i, int j, int k) { r(i, j, k) = {i, j, k}; });
    return r;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<int, D, symm>
  iota1() {
    gten3<int, D, symm> r;
    gten3<int, D, symm>::loop([&](int i, int j, int k) { r.set(i, j, k, i); });
    return r;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<int, D, symm>
  iota2() {
    gten3<int, D, symm> r;
    gten3<int, D, symm>::loop([&](int i, int j, int k) { r.set(i, j, k, j); });
    return r;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<int, D, symm>
  iota3() {
    gten3<int, D, symm> r;
    gten3<int, D, symm>::loop([&](int i, int j, int k) { r.set(i, j, k, k); });
    return r;
  }

  template <typename F, typename... Args,
            typename R =
                remove_cv_t<remove_reference_t<result_of_t<F(T, Args...)> > > >
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<R, D, symm>
  fmap(const F &f, const gten3 &x, const gten3<Args, D, symm> &...args) {
    return fmap(f, x.elts, args.elts...);
  }
  template <typename F, typename... Args>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  fmap_(const F &f, const gten3 &x, const gten3<Args, D, symm> &...args) {
    fmap_(f, x.elts, args.elts...);
  }

  template <
      typename... Args,
      typename R = remove_cv_t<remove_reference_t<result_of_t<T(Args...)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<R, D, symm>
  operator()(const Args &...args) const {
    return fmap([&](const T &var) { return var(args...); }, *this);
  }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename R = remove_cv_t<
                remove_reference_t<result_of_t<T(Arg1, Arg2, U)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3<R, D, symm>
  operator()(const Arg1 &arg1, const Arg2 &arg2,
             const gten3<U, D, symm> &val) const {
    return fmap([&](const T &var, const U &x) { return var(arg1, arg2, x); },
                *this, val);
  }
  // template <typename... Args, typename U, typename T1 = T,
  //           typename = decltype(declval<T1>().store(declval<Args>()...,
  //                                                   declval<U>()))>
  // ARITH_INLINE ARITH_DEVICE ARITH_HOST void store(const Args &...args,
  //                         const gten3<U, D, symm> &val)
  //                         const {
  //   fmap_([&](const auto &var, const auto &x) { return var.store(args..., x);
  //   },
  //        val);
  // }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename = decltype(declval<T1>().store(
                declval<Arg1>(), declval<Arg2>(), declval<U>()))>
  ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  store(const Arg1 &arg1, const Arg2 &arg2,
        const gten3<U, D, symm> &val) const {
    fmap_([&](const T &var, const U &x) { return var.store(arg1, arg2, x); },
          *this, val);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type size() const {
    return N;
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 != ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
  operator()(int i, int j, int k) const {
    return elts[ind(i, j, k)];
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 == ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T operator()(int i, int j,
                                                              int k) const {
    return sign(i, j, k) * elts[ind(i, j, k)];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void set(int i, int j, int k,
                                                          const T &x) {
#ifdef CCTK_DEBUG
    assert(sign(i, j, k) == 1);
#endif
    elts[ind(i, j, k)] = x;
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 != ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &operator()(int i, int j,
                                                               int k) {
    return elts[ind(i, j, k)];
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator+(const gten3 &x) {
    return {+x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator-(const gten3 &x) {
    return {-x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator+(const gten3 &x, const gten3 &y) {
    return {x.elts + y.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator-(const gten3 &x, const gten3 &y) {
    return {x.elts - y.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator*(const T &a, const gten3 &x) {
    return {a * x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator*(const gten3 &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator/(const gten3 &x, const T &a) {
    return {x.elts / a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator%(const gten3 &x, const T &a) {
    return {x.elts % a};
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator+=(const gten3 &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3
  operator-=(const gten3 &x) {
    return *this = *this - x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3 operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3 operator/=(const T &a) {
    return *this = *this / a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gten3 operator%=(const T &a) {
    return *this = *this % a;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator==(const gten3 & x, const gten3 & y) {
    return all(x.elts == y.elts);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator!=(const gten3 & x, const gten3 & y) {
    return !(x == y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  all(const gten3 &x) {
    return all(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  any(const gten3 &x) {
    return any(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  allisfinite(const gten3 &x) {
    return allisfinite(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const gten3 &x) {
    return anyisnan(x.elts);
  }

  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*gten3<bool, D, symm>*/
      isnan(const gten3 &x) {
    return isnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  maxabs(const gten3 &x) {
    return maxabs(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  sumabs(const gten3 &x) {
    return sumabs(x.elts);
  }

  friend ostream &operator<<(ostream &os, const gten3 &x) {
    os << "[";
    for (int k = 0; k < D; ++k) {
      if (k > 0)
        os << ",";
      os << "[";
      for (int j = 0; j < D; ++j) {
        if (j > 0)
          os << ",";
        os << "[";
        for (int i = 0; i < D; ++i) {
          if (i > 0)
            os << ",";
          os << x(i, j, k);
        }
        os << "]";
      }
      os << "]";
    }
    os << "]";
    return os;
  }
};

template <typename T, int D, symm_t symm> struct zero<gten3<T, D, symm> > {
  typedef gten3<T, D, symm> value_type;
  static constexpr value_type value = gten3<T, D, symm>::pure(zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return value;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return value;
  }
};

template <typename T, int D, symm_t symm> struct nan<gten3<T, D, symm> > {
  typedef gten3<T, D, symm> value_type;
  // static constexpr value_type value =
  //     gten3<T, D, symm>::pure(nan<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return gten3<T, D, symm>::pure(nan<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return gten3<T, D, symm>::pure(nan<T>());
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D> using ten3 = gten3<T, D, symm_t::full>;
template <typename T, int D> using sten3 = gten3<T, D, symm_t::symm>;
template <typename T, int D> using aten3 = gten3<T, D, symm_t::anti>;

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D, symm_t symm>
constexpr gten3<simd<T>, D, symm> if_else(const simdl<T> &cond,
                                          const gten3<simd<T>, D, symm> &x,
                                          const gten3<simd<T>, D, symm> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_TEN3_HXX
