#ifndef CARPETX_ALGO_ROOTS_HXX
#define CARPETX_ALGO_ROOTS_HXX

#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <cctk.h>

#include <boost/math/tools/roots.hpp>

#ifdef __HIPCC__
#include <hip/hip_runtime.h>
#endif

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <limits>
#include <utility>

// For accelerators: Declare whether a function should live on the device or the
// host (or both)
#if defined __CUDACC__
#define ALGO_DEVICE __device__
#define ALGO_HOST __host__
#elif defined __HIPCC__
#define ALGO_DEVICE __device__
#define ALGO_HOST __host__
#else
#define ALGO_DEVICE
#define ALGO_HOST
#endif

namespace Algo {

// These helpers exist so that the `using std::...` declaration can live in a
// function body. That matters for `clamp1`, which is called with arguments
// named `min` and `max` that would otherwise shadow `std::min` and `std::max`.

namespace {
template <typename T>
constexpr ALGO_HOST ALGO_DEVICE void swap1(T &x, T &y) {
  T z = std::move(x);
  x = std::move(y);
  y = std::move(z);
}
} // namespace

namespace {
template <typename T>
constexpr ALGO_HOST ALGO_DEVICE T ldexp1(const T &x, const int n) {
  using std::ldexp;
  return ldexp(x, n);
}
} // namespace

namespace {
template <typename T>
constexpr ALGO_HOST ALGO_DEVICE T max1(const T &x, const T &y) {
  using std::max;
  return max(x, y);
}
} // namespace

namespace {
template <typename T>
constexpr ALGO_HOST ALGO_DEVICE T clamp1(const T &x, const T &lo, const T &hi) {
  using std::max, std::min;
  return min(max(x, lo), hi);
}
} // namespace

// A relative convergence criterion, equivalent to
// `boost::math::tools::eps_tolerance`: two values are considered equal once
// they agree to `min_bits` binary digits.
//
// Note that this is a purely relative criterion, as it is in Boost. A bracket
// that straddles zero therefore never satisfies it, because `min(|x|, |y|)`
// tends to zero along with the bracket; such a search terminates on its
// iteration count instead. Callers looking for a root at zero need an absolute
// criterion, which this class deliberately does not provide.
template <typename T> class eps_tolerance {
  T eps;

public:
  constexpr ALGO_HOST ALGO_DEVICE eps_tolerance()
      : eps(4 * std::numeric_limits<T>::epsilon()) {}
  // `eps` is clamped from below at `4 * epsilon`, as Boost does. Without the
  // clamp, any `min_bits >= digits` asks for a precision that no two distinct
  // floating-point numbers can satisfy, and the search silently runs until it
  // exhausts its iteration count.
  explicit constexpr ALGO_HOST ALGO_DEVICE eps_tolerance(const int min_bits)
      : eps(max1(ldexp1(T(1), 1 - min_bits),
                 T(4 * std::numeric_limits<T>::epsilon()))) {}
  constexpr ALGO_HOST ALGO_DEVICE bool operator()(const T &x,
                                                  const T &y) const {
    using std::abs, std::min;
    return abs(x - y) <= eps * min(abs(x), abs(y));
  }
};

template <typename F, typename T>
std::pair<T, T> bisect(F &&f, T min, T max, int min_bits, int max_iters,
                       int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::bisect(
      std::forward<F>(f), min, max,
      boost::math::tools::eps_tolerance<T>(min_bits), max_iter);
  iters = max_iter;
  return res;
}

template <typename F, typename T>
std::pair<T, T> bracket_and_solve_root(F &&f, T guess, T factor, bool rising,
                                       int min_bits, int max_iters,
                                       int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::bracket_and_solve_root(
      std::forward<F>(f), guess, factor, rising, eps_tolerance<T>(min_bits),
      max_iter);
  iters = max_iter;
  return res;
}

// See <https://en.wikipedia.org/wiki/Brent%27s_method>
//
template <typename F, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE ALGO_HOST ALGO_DEVICE std::pair<T, T>
brent(F f, T a, T b, int min_bits, int max_iters, int &iters) {
  using std::abs, std::min, std::max;

  const auto tol = eps_tolerance<T>(min_bits);

  iters = 0;
  auto fa = f(a);
  auto fb = f(b);
  if (abs(fa) < abs(fb)) {
    swap1(a, b);
    swap1(fa, fb);
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
      swap1(u, v);
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
    assert(fa * fb <= 0);
    if (abs(fa) < abs(fb)) {
      swap1(a, b);
      swap1(fa, fb);
    }
    ++iters;
  }

  if (fb == 0)
    return {b, b};
  return {min(a, b), max(a, b)};
}

// Requires function and its derivative
template <typename F, typename T>
T newton_raphson(F &&f, T guess, T min, T max, int min_bits, int max_iters,
                 int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::newton_raphson_iterate(
      std::forward<F>(f), guess, min, max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

// Requires function and first two derivatives
template <typename F, typename T>
T halley(F &&f, T guess, T min, T max, int min_bits, int max_iters,
         int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::halley_iterate(std::forward<F>(f), guess, min,
                                                max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

// Requires function and first two derivatives
template <typename F, typename T>
T schroder(F &&f, T guess, T min, T max, int min_bits, int max_iters,
           int &iters) {
  std::uintmax_t max_iter = max_iters;
  auto res = boost::math::tools::schroder_iterate(
      std::forward<F>(f), guess, min, max, min_bits, max_iter);
  iters = max_iter;
  return res;
}

// Multi-dimensional Newton-Raphson iteration. `f` returns both the residual
// and its Jacobian. Iterates are confined to the box `[min, max]`.
//
// `iters` is the number of Newton steps taken, at most `max_iters`. `failed`
// is set if no solution was found, i.e. if the Jacobian became singular, if the
// iteration stalled on a bound, or if `max_iters` was reached.
template <typename F, typename T, int N>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE ALGO_HOST ALGO_DEVICE Arith::vec<T, N>
newton_raphson_nd(F f, const Arith::vec<T, N> &guess,
                  const Arith::vec<T, N> &min, const Arith::vec<T, N> &max,
                  int min_bits, int max_iters, int &iters, bool &failed) {
  using vec = Arith::vec<T, N>;
  using mat = Arith::mat<T, N>;
  using std::isfinite;

  const auto tolfx = eps_tolerance<T>(min_bits);

  failed = true;
  vec x([&](int i) { return clamp1(guess(i), min(i), max(i)); });

  for (iters = 0; iters <= max_iters; ++iters) {
    const auto [fx0, jac0] = f(x);
    const vec fx = fx0;
    const mat jac = jac0;

    // Comparing `1 + errfx` against `1` reduces to `errfx <= eps`, i.e. this is
    // an absolute criterion on the residual, and hence sensitive to how `f` is
    // scaled.
    const T errfx = sumabs(fx);
    if (tolfx(1 + errfx, 1)) {
      failed = false;
      return x;
    }
    if (iters == max_iters)
      break;

    const T det_jac = calc_det(jac);
    if (!isfinite(det_jac) || det_jac == 0)
      // Singular Jacobian; there is no Newton step to take
      break;
    const mat inv_jac = calc_inv(jac, det_jac);
    const vec dx([&](int i) {
      return -Arith::sum<N>([&](int j) { return inv_jac(i, j) * fx(j); });
    });
    const vec xnew([&](int i) { return clamp1(x(i) + dx(i), min(i), max(i)); });

    bool moved = false;
    for (int i = 0; i < N; ++i)
      moved |= xnew(i) != x(i);
    if (!moved)
      // The iteration is stuck, most likely against a bound
      break;
    x = xnew;
  }

  return x;
}

} // namespace Algo

#endif // #ifndef CARPETX_ALGO_ROOTS_HXX
