#ifndef ALGO_HXX
#define ALGO_HXX

#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <boost/math/tools/roots.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>

namespace Algo {

template <typename F, typename T>
std::pair<T, T> bisect(F &&f, T min, T max, int min_bits, int max_iters,
                       int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::bisect(
      std::forward<F>(f), min, max,
      boost::math::tools::eps_tolerance<double>(min_bits), max_iter);
  iters = max_iter;
  return res;
}

template <typename F, typename T>
std::pair<T, T> bracket_and_solve_root(F &&f, T guess, T factor, bool rising,
                                       int min_bits, int max_iters,
                                       int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::bracket_and_solve_root(
      std::forward<F>(f), guess, factor, rising,
      boost::math::tools::eps_tolerance<double>(min_bits), max_iter);
  iters = max_iter;
  return res;
}

// See <https://en.wikipedia.org/wiki/Brent%27s_method>
template <typename F, typename T>
std::pair<T, T> brent(F f, T a, T b, int min_bits, int max_iters, int &iters) {
  using std::abs, std::min, std::max, std::swap;

  auto tol = boost::math::tools::eps_tolerance<T>(min_bits);

  iters = 0;
  auto fa = f(a);
  auto fb = f(b);
  if (abs(fa) < abs(fb)) {
    swap(a, b);
    swap(fa, fb);
  }
  if (fb == 0)
    return {b, b};
  if (fa * fb >= 0) {
    // Root is not bracketed
    iters = max_iters;
    return {min(a, b), max(a, b)};
  }
  T c = a;
  auto fc = fa;
  bool mflag = true;
  T d{};

  while (fb != 0 && !tol(a, b) && iters < max_iters) {
    T s;
    if (fa != fc && fb != fc)
      // inverse quadratic interpolation
      s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
          (b * fa * fc) / ((fb - fa) * (fb - fc)) +
          (c * fa * fb) / ((fc - fa) * (fc - fb));
    else
      // secant method
      s = (a + b) / 2 - (fa + fb) / 2 * (b - a) / (fb - fa);
    T u = (3 * a + b) / 4;
    T v = b;
    if (u > v)
      swap(u, v);
    bool cond1 = !(u <= s && s <= v);
    bool cond2 = mflag && abs(s - b) >= abs(b - c) / 2;
    bool cond3 = !mflag && abs(s - b) >= abs(c - d) / 2;
    bool cond4 = mflag && tol(c, b);
    bool cond5 = !mflag && tol(c, d);
    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      // bisection
      s = (a + b) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    auto fs = f(s);
    // `d` is assigned for the first time here; it won't be used above on the
    // first iteration because `mflag` is set
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    // CCTK_VINFO("iters=%d mflag=%d   a=%.17g b=%.17g c=%.17g d=%.17g fa=%.17g"
    //            "fb=%.17g fc=%.17g",
    //            iters, int(mflag), double(a), double(b), double(c), double(d),
    //            double(fa), double(fb), double(fc));
    assert(fa * fb <= 0);
    if (abs(fa) < abs(fb)) {
      swap(a, b);
      swap(fa, fb);
    }
    ++iters;
  }

  if (fb == 0)
    return {b, b};
  return {min(a, b), max(a, b)};
}

// Requires function and its derivative
template <typename F, typename T>
T newton_raphson(F f, T guess, T min, T max, int min_bits, int max_iters,
                 int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::newton_raphson_iterate(
      std::forward<F>(f), guess, min, max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

// Requires function and first two derivatives
template <typename F, typename T>
T halley(F f, T guess, T min, T max, int min_bits, int max_iters, int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::halley_iterate(std::forward<F>(f), guess, min,
                                                max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

// Requires function and first two derivatives
template <typename F, typename T>
T schroder(F f, T guess, T min, T max, int min_bits, int max_iters,
           int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::halley_iterate(std::forward<F>(f), guess, min,
                                                max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

template <typename F, typename T, int N>
Arith::vec<T, N> newton_raphson_nd(F f, const Arith::vec<T, N> &guess,
                                   const Arith::vec<T, N> &min,
                                   const Arith::vec<T, N> &max, int min_bits,
                                   int max_iters, int &iters, bool &failed) {
  using vec = Arith::vec<T, N>;
  using mat = Arith::mat<T, N>;
  failed = false;
  auto tolfx = boost::math::tools::eps_tolerance<T>(min_bits);
  vec x = guess;
  for (iters = 1; iters <= max_iters; ++iters) {
    const auto [fx0, jac0] = f(x);
    const vec fx = fx0;
    const mat jac = jac0;
    const T errfx = sumabs(fx);
    if (tolfx(1 + errfx, 1))
      return x;
    const T det_jac = calc_det(jac);
    const mat inv_jac = calc_inv(jac, det_jac);
    const vec dx([&](int i) {
      return -Arith::sum<2>([&](int j) { return inv_jac(i, j) * fx(j); });
    });
    x += dx;
  }
  failed = true;
  return x;
}

} // namespace Algo

#endif // #ifndef ALGO_HXX
