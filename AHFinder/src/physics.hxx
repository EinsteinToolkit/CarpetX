#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include <defs.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <cmath>

namespace AHFinder {
using namespace Arith;

template <typename T> constexpr T det(const smat<T, 3> &g) {
  return 2 * g(0, 1) * g(0, 2) * g(1, 2) - g(2, 2) * pow2(g(0, 1)) -
         g(1, 1) * pow2(g(0, 2)) +
         g(0, 0) * (g(1, 1) * g(2, 2) - pow2(g(1, 2)));
}

template <typename T>
constexpr smat<T, 3> inv(const smat<T, 3> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 3>{
                     g(1, 1) * g(2, 2) - pow2(g(1, 2)),
                     g(0, 2) * g(1, 2) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),
                     g(0, 0) * g(2, 2) - pow2(g(0, 2)),
                     g(0, 1) * g(0, 2) - g(0, 0) * g(1, 2),
                     g(0, 0) * g(1, 1) - pow2(g(0, 1)),
                 };
}
template <typename T> constexpr smat<T, 3> inv(const smat<T, 3> &g) {
  return inv(g, det(g));
}

template <typename T> constexpr T det(const mat<T, 3> &g) {
  return g(0, 2) * (-(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1)) +
         g(0, 1) * (g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2)) +
         g(0, 0) * (-(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2));
}

template <typename T>
constexpr mat<T, 3> inv(const mat<T, 3> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 3>{
                     -(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2),
                     g(0, 2) * g(2, 1) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),

                     g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2),
                     -(g(0, 2) * g(2, 0)) + g(0, 0) * g(2, 2),
                     g(0, 2) * g(1, 0) - g(0, 0) * g(1, 2),

                     -(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1),
                     g(0, 1) * g(2, 0) - g(0, 0) * g(2, 1),
                     -(g(0, 1) * g(1, 0)) + g(0, 0) * g(1, 1),
                 };
}
template <typename T> constexpr mat<T, 3> inv(const mat<T, 3> &g) {
  return inv(g, det(g));
}

template <typename T1, typename T2, int D, symm_t symm,
          typename R = decltype(std::declval<T1>() * std::declval<T2>())>
constexpr vec<R, D> matmul(const gmat<T1, D, symm> &A, const vec<T2, D> &x) {
  return vec<R, D>(
      [&](int a) { return sum<D>([&](int b) { return A(a, b) * x(b); }); });
}

template <typename T, int D, symm_t symm>
void gram_schmidt(const gmat<T, D, symm> g, vec<T, D> &x1, vec<T, D> &x2,
                  vec<T, D> &x3) {
  using std::sqrt;

  // First vector
  T x1x1 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x1(b); });
  x1 /= sqrt(x1x1);

  // Second vector
  T x1x2 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x2(b); });
  x2 -= x1x2 * x1;
  T x2x2 = sum<D>([&](int a, int b) { return g(a, b) * x2(a) * x2(b); });
  x2 /= sqrt(x2x2);

  // Third vector
  T x1x3 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x3(b); });
  x3 -= x1x3 * x1;
  T x2x3 = sum<D>([&](int a, int b) { return g(a, b) * x2(a) * x3(b); });
  x3 -= x2x3 * x2;
  T x3x3 = sum<D>([&](int a, int b) { return g(a, b) * x3(a) * x3(b); });
  x3 /= sqrt(x3x3);

#if 0
  x1x1 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x1(b); });
  x1x2 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x2(b); });
  x2x2 = sum<D>([&](int a, int b) { return g(a, b) * x2(a) * x2(b); });
  x1x3 = sum<D>([&](int a, int b) { return g(a, b) * x1(a) * x3(b); });
  x2x3 = sum<D>([&](int a, int b) { return g(a, b) * x2(a) * x3(b); });
  x3x3 = sum<D>([&](int a, int b) { return g(a, b) * x3(a) * x3(b); });
  using std::abs;
  assert(abs(x1x1 - 1) <= 1.0e-12);
  assert(abs(x1x2 - 0) <= 1.0e-12);
  assert(abs(x2x2 - 1) <= 1.0e-12);
  assert(abs(x1x3 - 0) <= 1.0e-12);
  assert(abs(x2x3 - 0) <= 1.0e-12);
  assert(abs(x3x3 - 1) <= 1.0e-12);
#endif
}

} // namespace AHFinder

#endif // #ifndef PHYSICS_HXX
