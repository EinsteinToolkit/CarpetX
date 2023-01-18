#include "mat.hxx"

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

template <typename T, typename U> constexpr bool eq(const T &x, const U &y) {
  using std::isnan;
  return x == y || (isnan(x) && isnan(y));
}

template <typename T, int D, symm_t symm>
constexpr bool eqm(const gmat<T, D, symm> &x, const gmat<T, D, symm> &y) {
  return all(fmap([](const auto &a, const auto &b) { return eq(a, b); }, x, y));
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestMat() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
  using M3D = smat<CCTK_REAL, 3>;

  constexpr CCTK_REAL N = nan<CCTK_REAL>();
  static_assert(eq(smat<CCTK_REAL, 3>()(0, 0), N));
  static_assert(eq(smat<CCTK_REAL, 3>()(0, 1), N));
  static_assert(eq(smat<CCTK_REAL, 3>()(0, 2), N));
  static_assert(eq(smat<CCTK_REAL, 3>()(1, 1), N));
  static_assert(eq(smat<CCTK_REAL, 3>()(1, 2), N));
  static_assert(eq(smat<CCTK_REAL, 3>()(2, 2), N));

  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(0, 0), 1));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(0, 1), 2));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(0, 2), 3));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(1, 0), 2));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(1, 1), 4));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(1, 2), 5));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(2, 0), 3));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(2, 1), 5));
  static_assert(eq(smat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6})(2, 2), 6));

  static_assert(eqm(smat<CCTK_REAL, 3>::iota1(), {0, 0, 0, 1, 1, 2}));
  static_assert(eqm(smat<CCTK_REAL, 3>::iota2(), {0, 1, 2, 1, 2, 2}));

  static_assert(eqm(smat<CCTK_REAL, 3>::unit(0, 0), {1, 0, 0, 0, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(0, 1), {0, 1, 0, 0, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(0, 2), {0, 0, 1, 0, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(1, 0), {0, 1, 0, 0, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(1, 1), {0, 0, 0, 1, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(1, 2), {0, 0, 0, 0, 1, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(2, 0), {0, 0, 1, 0, 0, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(2, 1), {0, 0, 0, 0, 1, 0}));
  static_assert(eqm(smat<CCTK_REAL, 3>::unit(2, 2), {0, 0, 0, 0, 0, 1}));

  constexpr M3D x = {71, 90, 14, 50, 41, 65};
  static_assert(eqm(+x, {71, 90, 14, 50, 41, 65}));
  static_assert(eqm(-x, {-71, -90, -14, -50, -41, -65}));

  constexpr M3D y = {77, 32, 67, 5, 81, 10};
  static_assert(eqm(x + y, {148, 122, 81, 55, 122, 75}));
  static_assert(eqm(x - y, {-6, 58, -53, 45, -40, 55}));
  static_assert(eqm(3 * x, {213, 270, 42, 150, 123, 195}));
  static_assert(eqm(x * 3, {213, 270, 42, 150, 123, 195}));

  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(0, 0), 1));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(0, 1), 2));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(0, 2), 3));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(1, 0), 4));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(1, 1), 5));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(1, 2), 6));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(2, 0), 7));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(2, 1), 8));
  static_assert(eq(mat<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9})(2, 2), 9));

  static_assert(eqm(mat<CCTK_REAL, 3>::iota1(), {0, 0, 0, 1, 1, 1, 2, 2, 2}));
  static_assert(eqm(mat<CCTK_REAL, 3>::iota2(), {0, 1, 2, 0, 1, 2, 0, 1, 2}));

  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(0, 0), 0));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(0, 1), 1));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(0, 2), 2));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(1, 0), -1));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(1, 1), 0));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(1, 2), 3));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(2, 0), -2));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(2, 1), -3));
  static_assert(eq(amat<CCTK_REAL, 3>({1, 2, 3})(2, 2), 0));

  static_assert(eqm(amat<CCTK_REAL, 3>::iota1(), {0, 0, 1}));
  static_assert(eqm(amat<CCTK_REAL, 3>::iota2(), {1, 2, 2}));
#endif
}

} // namespace Arith
