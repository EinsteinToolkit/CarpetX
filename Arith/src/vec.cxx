#include "vec.hxx"

#include <cctk.h>

#include <cmath>
#include <functional>

namespace Arith {
using namespace std;

template <typename T, typename U> constexpr bool eq(const T &x, const U &y) {
  return x == y || (x != x && y != y);
}

template <typename T, int D>
constexpr bool eqv(const vec<T, D> &x, const vec<T, D> &y) {
  using std::isnan;
  return all(fmap([](const auto &a, const auto &b) { return eq(a, b); }, x, y));
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestVec() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#if !defined __CUDACC__ && !defined __HIPCC__
  using V3D = vec<CCTK_REAL, 3>;

  constexpr CCTK_REAL N = nan<CCTK_REAL>();
  static_assert(eq(vec<CCTK_REAL, 3>()(0), N));
  static_assert(eq(vec<CCTK_REAL, 3>()(1), N));
  static_assert(eq(vec<CCTK_REAL, 3>()(2), N));
  static_assert(eq(vec<CCTK_REAL, 3>()(0), N));
  static_assert(eq(vec<CCTK_REAL, 3>()(1), N));
  static_assert(eq(vec<CCTK_REAL, 3>()(2), N));

  static_assert(eq(vec<CCTK_REAL, 3>{1, 2, 3}(0), 1));
  static_assert(eq(vec<CCTK_REAL, 3>{1, 2, 3}(1), 2));
  static_assert(eq(vec<CCTK_REAL, 3>{1, 2, 3}(2), 3));

  static_assert(eqv(vec<CCTK_REAL, 3>::iota(), {0, 1, 2}));

  static_assert(eqv(vec<CCTK_REAL, 3>::unit(0), {1, 0, 0}));
  static_assert(eqv(vec<CCTK_REAL, 3>::unit(1), {0, 1, 0}));
  static_assert(eqv(vec<CCTK_REAL, 3>::unit(2), {0, 0, 1}));

  constexpr V3D x = {5, 6, 8};
  static_assert(eqv(+x, {5, 6, 8}));
  static_assert(eqv(-x, {-5, -6, -8}));

  constexpr V3D y = {10, 4, 6};
  static_assert(eqv(x + y, {15, 10, 14}));
  static_assert(eqv(x - y, {-5, 2, 2}));

  static_assert(eqv(3 * x, {15, 18, 24}));
  static_assert(eqv(x * 3, {15, 18, 24}));
#endif
}

} // namespace Arith
