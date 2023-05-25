#ifndef CARPETX_ARITH_VEC_HXX
#define CARPETX_ARITH_VEC_HXX

#include "defs.hxx"
#include "dual.hxx"
#include "simd.hxx"
#include "tuple.hxx"
#include "vect.hxx"

#include <cctk.h>

#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

namespace Arith {
using namespace std;

// Symmetries

enum class symm_t : int { full, symm, anti };
constexpr symm_t FULL = symm_t::full;
constexpr symm_t SYMM = symm_t::symm;
constexpr symm_t ANTI = symm_t::anti;
template <symm_t symm, typename T1, typename T2, typename T3>
constexpr auto if_symm(const T1 &if_full, const T2 &if_symm,
                       const T3 &if_anti) {
  if constexpr (symm == FULL)
    return if_full;
  if constexpr (symm == SYMM)
    return if_symm;
  if constexpr (symm == ANTI)
    return if_anti;
}
template <typename T1, typename T2, typename T3>
constexpr auto if_symm(symm_t symm, const T1 &if_full, const T2 &if_symm,
                       const T3 &if_anti) {
  if (symm == FULL)
    return if_full;
  if (symm == SYMM)
    return if_symm;
  if (symm == ANTI)
    return if_anti;
  abort();
}
inline ostream &operator<<(ostream &os, const symm_t symm) {
  return os << if_symm(symm, "", "s", "a");
}

////////////////////////////////////////////////////////////////////////////////

// A small vector.

// There are two small vector classes, `vect` and `vec`. The former
// (`vect`) is a drop-in replacement for `std::array` which supports
// arithmetic operations. The latter (`vec`, defined here) is a rank-1
// tensor compatible with `mat`, `ten3`, and `rten`. While similar,
// these classes have different APIs because they are used
// differently.
//
// The rule of thumb is: If you are looking for a type similar to a
// C++ array, use `vect`. If you are looking for a type for a physics
// quantity, use `vec`.
//
// It might be possible to unify these two types. In practice, there
// is little confusion, so this is not urgent.

template <typename T, int D> struct vec {

  constexpr static int N = D;
  vect<T, N> elts;

public:
  typedef T value_type;
  typedef int size_type;
  static constexpr int size_value = N;

private:
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < D);
#endif
    return n;
  }

public:
  // initializes all elements to nan
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec()
      : elts(nan<vect<T, N> >()) {}

  constexpr ARITH_INLINE vec(const vec &) = default;
  constexpr ARITH_INLINE vec(vec &&) = default;
  constexpr ARITH_INLINE vec &operator=(const vec &) = default;
  constexpr ARITH_INLINE vec &operator=(vec &&) = default;

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(vec<U, D> x)
      : elts(std::move(x.elts)) {}

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(initializer_list<T> x)
      : elts(x) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(vect<T, N> elts)
      : elts(std::move(elts)) {}
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(array<T, N> x)
      : elts(std::move(x)) {}
  // explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(const vector<T>
  // &x) : elts(x) {} explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  // vec(vector<T> &&x) : elts(std::move(x)) {}

  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator vect<T, N>() const {
    return elts;
  }

  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void loop(const F &f) {
    vect<T, D>::loop(f);
  }
  template <typename F, typename = result_of_t<F(int)> >
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(const F &f)
      : elts(vect<T, D>::make(f)) {}

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec pure(const T &a) {
    return vect<T, N>::pure(a);
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec unit(const int i) {
    return vect<T, N>::unit(i);
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<int, D> iota() {
    return {vect<int, N>::iota()};
  }

  template <typename F, typename... Args,
            typename R =
                remove_cv_t<remove_reference_t<result_of_t<F(T, Args...)> > > >
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D>
  fmap(const F &f, const vec &x, const vec<Args, D> &...args) {
    return fmap(f, x.elts, args.elts...);
  }
  template <typename F, typename... Args>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  fmap_(const F &f, const vec &x, const vec<Args, D> &...args) {
    fmap_(f, x.elts, args.elts...);
  }

  template <
      typename... Args,
      typename R = remove_cv_t<remove_reference_t<result_of_t<T(Args...)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D>
  operator()(const Args &...args) const {
    return fmap([&](const T &var) { return var(args...); }, *this);
  }
  template <typename Arg1, typename Arg2, typename U,
            typename R = remove_cv_t<
                remove_reference_t<result_of_t<T(Arg1, Arg2, U)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D>
  operator()(const Arg1 &arg1, const Arg2 &arg2, const vec<U, D> &val) const {
    return fmap([&](const T &var, const U &x) { return var(arg1, arg2, x); },
                *this, val);
  }
  // template <typename... Args, typename U, typename T1 = T,
  //           typename = decltype(declval<T1>().store(declval<Args>()...,
  //                                                   declval<U>()))>
  // ARITH_INLINE ARITH_DEVICE ARITH_HOST void store(const Args &...args,
  //                         const vec<U, D> &val) const {
  //   fmap_([&](const auto &var, const auto &x) { return var.store(args..., x);
  //   },
  //        val);
  // }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename = decltype(declval<T1>().store(
                declval<Arg1>(), declval<Arg2>(), declval<U>()))>
  ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  store(const Arg1 &arg1, const Arg2 &arg2, const vec<U, D> &val) const {
    fmap_([&](const T &var, const U &x) { return var.store(arg1, arg2, x); },
          *this, val);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type size() const {
    return N;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
  operator()(int i) const {
    return elts[ind(i)];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &operator()(int i) {
    return elts[ind(i)];
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator+(const vec<T, D> &x) {
    return {+x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator-(const vec<T, D> &x) {
    return {-x.elts};
  }
  template <typename U,
            typename R = decltype(std::declval<T>() + std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D>
  operator+(const vec<T, D> &x, const vec<U, D> &y) {
    return {x.elts + y.elts};
  }
  template <typename U,
            typename R = decltype(std::declval<T>() - std::declval<U>())>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D>
  operator-(const vec<T, D> &x, const vec<U, D> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator*(const T &a, const vec<T, D> &x) {
    return {a * x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator*(const vec<T, D> &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator/(const vec<T, D> &x, const T &a) {
    return {x.elts / a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D>
  operator%(const vec<T, D> &x, const T &a) {
    return {x.elts % a};
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator+=(const vec &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator-=(const vec &x) {
    return *this = *this - x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec
  operator+=(const vec<U, D> &x) {
    return *this = *this + x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec
  operator-=(const vec<U, D> &x) {
    return *this = *this - x;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator*=(const U &a) {
    return *this = *this * a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator/=(const U &a) {
    return *this = *this / a;
  }
  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator%=(const U &a) {
    return *this = *this % a;
  }

  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator==(const vec<T, D> &x, const vec<U, D> &y) {
    return all(x.elts == y.elts);
  }
  template <typename U>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator!=(const vec<T, D> &x, const vec<U, D> &y) {
    return !(x == y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  all(const vec &x) {
    return all(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  any(const vec &x) {
    return any(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  allisfinite(const vec &x) {
    return allisfinite(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const vec &x) {
    return anyisnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*vec<bool, D>*/
  isnan(const vec &x) {
    return isnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T maxabs(const vec &x) {
    return maxabs(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T sumabs(const vec &x) {
    return sumabs(x.elts);
  }

  friend ostream &operator<<(ostream &os, const vec<T, D> &v) {
    os << "[";
    for (int i = 0; i < D; ++i) {
      if (i > 0)
        os << ",";
      os << v(i);
    }
    os << "]";
    return os;
  }
};

template <typename T, int D> struct zero<vec<T, D> > {
  typedef vec<T, D> value_type;
  // static constexpr value_type value = vec<T, D>::pure(zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vec<T, D>::pure(zero<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vec<T, D>::pure(zero<T>());
  }
};

template <typename T, int D> struct nan<vec<T, D> > {
  typedef vec<T, D> value_type;
  // static constexpr value_type value = vec<T, D>::pure(nan<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vec<T, D>::pure(nan<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vec<T, D>::pure(nan<T>());
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D>
constexpr vec<simd<T>, D> if_else(const simdl<T> &cond,
                                  const vec<simd<T>, D> &x,
                                  const vec<simd<T>, D> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

template <typename T, typename U, int D>
constexpr vec<dual<simd<T>, U>, D> if_else(const simdl<T> &cond,
                                           const vec<dual<simd<T>, U>, D> &x,
                                           const vec<dual<simd<T>, U>, D> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

} // namespace Arith

#endif // #ifndef CARPETX_ARITH_VEC_HXX
