#include "ten3.hxx"

#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

template <typename T, typename U> constexpr bool eq(const T &x, const U &y) {
  using std::isnan;
  return x == y || (isnan(x) && isnan(y));
}

template <typename T, int D, symm_t symm>
constexpr bool eqm(const gten3<T, D, symm> &x, const gten3<T, D, symm> &y) {
  return all(fmap([](const auto &a, const auto &b) { return eq(a, b); }, x, y));
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestTen3() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#if !defined __CUDACC__ && !defined __HIPCC__
  using T3D = sten3<CCTK_REAL, 3>;

  constexpr CCTK_REAL N = nan<CCTK_REAL>();
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 0, 0), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 0, 1), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 0, 2), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 1, 1), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 1, 2), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(0, 2, 2), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(1, 1, 1), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(1, 1, 2), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(1, 2, 2), N));
  static_assert(eq(sten3<CCTK_REAL, 3>()(2, 2, 2), N));

  constexpr auto t3 = sten3<CCTK_REAL, 3>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
  static_assert(eq(t3(0, 0, 0), 1));
  static_assert(eq(t3(0, 0, 1), 2));
  static_assert(eq(t3(0, 0, 2), 3));
  static_assert(eq(t3(0, 1, 0), 2));
  static_assert(eq(t3(0, 1, 1), 4));
  static_assert(eq(t3(0, 1, 2), 5));
  static_assert(eq(t3(0, 2, 0), 3));
  static_assert(eq(t3(0, 2, 1), 5));
  static_assert(eq(t3(0, 2, 2), 6));
  static_assert(eq(t3(1, 0, 0), 2));
  static_assert(eq(t3(1, 0, 1), 4));
  static_assert(eq(t3(1, 0, 2), 5));
  static_assert(eq(t3(1, 1, 0), 4));
  static_assert(eq(t3(1, 1, 1), 7));
  static_assert(eq(t3(1, 1, 2), 8));
  static_assert(eq(t3(1, 2, 0), 5));
  static_assert(eq(t3(1, 2, 1), 8));
  static_assert(eq(t3(1, 2, 2), 9));
  static_assert(eq(t3(2, 0, 0), 3));
  static_assert(eq(t3(2, 0, 1), 5));
  static_assert(eq(t3(2, 0, 2), 6));
  static_assert(eq(t3(2, 1, 0), 5));
  static_assert(eq(t3(2, 1, 1), 8));
  static_assert(eq(t3(2, 1, 2), 9));
  static_assert(eq(t3(2, 2, 0), 6));
  static_assert(eq(t3(2, 2, 1), 9));
  static_assert(eq(t3(2, 2, 2), 10));

  static_assert(eqm(T3D::iota1(), {0, 0, 0, 0, 0, 0, 1, 1, 1, 2}));
  static_assert(eqm(T3D::iota2(), {0, 0, 0, 1, 1, 2, 1, 1, 2, 2}));
  static_assert(eqm(T3D::iota3(), {0, 1, 2, 1, 2, 2, 1, 2, 2, 2}));

  static_assert(eqm(T3D::unit(0, 0, 0), {1, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 0, 1), {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 0, 2), {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 1, 0), {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 1, 1), {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 1, 2), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 2, 0), {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 2, 1), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(0, 2, 2), {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 0, 0), {0, 1, 0, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 0, 1), {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 0, 2), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 1, 0), {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 1, 1), {0, 0, 0, 0, 0, 0, 1, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 1, 2), {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}));
  static_assert(eqm(T3D::unit(1, 2, 0), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(1, 2, 1), {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}));
  static_assert(eqm(T3D::unit(1, 2, 2), {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}));
  static_assert(eqm(T3D::unit(2, 0, 0), {0, 0, 1, 0, 0, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(2, 0, 1), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(2, 0, 2), {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(2, 1, 0), {0, 0, 0, 0, 1, 0, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(2, 1, 1), {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}));
  static_assert(eqm(T3D::unit(2, 1, 2), {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}));
  static_assert(eqm(T3D::unit(2, 2, 0), {0, 0, 0, 0, 0, 1, 0, 0, 0, 0}));
  static_assert(eqm(T3D::unit(2, 2, 1), {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}));
  static_assert(eqm(T3D::unit(2, 2, 2), {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}));

  constexpr T3D x = {71, 90, 14, 50, 41, 65, 81, 18, 61, 22};
  static_assert(eqm(+x, {71, 90, 14, 50, 41, 65, 81, 18, 61, 22}));
  static_assert(eqm(-x, {-71, -90, -14, -50, -41, -65, -81, -18, -61, -22}));

  constexpr T3D y = {77, 32, 67, 5, 81, 10, 80, 57, 92, 92};
  static_assert(eqm(x + y, {148, 122, 81, 55, 122, 75, 161, 75, 153, 114}));
  static_assert(eqm(x - y, {-6, 58, -53, 45, -40, 55, 1, -39, -31, -70}));
  static_assert(eqm(3 * x, {213, 270, 42, 150, 123, 195, 243, 54, 183, 66}));
  static_assert(eqm(x * 3, {213, 270, 42, 150, 123, 195, 243, 54, 183, 66}));

  // TODO: Complete this test
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 0, 0),
         1));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 0, 1),
         2));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 0, 2),
         3));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 1, 0),
         4));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 1, 1),
         5));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 1, 2),
         6));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 2, 0),
         7));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 2, 1),
         8));
  static_assert(
      eq(ten3<CCTK_REAL, 3>({1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27})(0, 2, 2),
         9));

  static_assert(eqm(ten3<CCTK_REAL, 3>::iota1(),
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
                     1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2}));
  static_assert(eqm(ten3<CCTK_REAL, 3>::iota2(),
                    {0, 0, 0, 1, 1, 1, 2, 2, 2, 0, 0, 0, 1, 1,
                     1, 2, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2}));
  static_assert(eqm(ten3<CCTK_REAL, 3>::iota3(),
                    {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1,
                     2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}));

  // TODO: Convert this test to 4D
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 0, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 0, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 0, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 1, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 1, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 1, 2), 1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 2, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 2, 1), -1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(0, 2, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 0, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 0, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 0, 2), -1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 1, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 1, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 1, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 2, 0), 1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 2, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(1, 2, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 0, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 0, 1), 1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 0, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 1, 0), -1));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 1, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 1, 2), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 2, 0), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 2, 1), 0));
  static_assert(eq(aten3<CCTK_REAL, 3>({1})(2, 2, 2), 0));

  static_assert(eqm(aten3<CCTK_REAL, 3>::iota1(), {0}));
  static_assert(eqm(aten3<CCTK_REAL, 3>::iota2(), {1}));
  static_assert(eqm(aten3<CCTK_REAL, 3>::iota3(), {2}));
#endif

  static_assert(sten3<CCTK_REAL, 4>().size() == 20);
  static_assert(ten3<CCTK_REAL, 4>().size() == 64);
  static_assert(aten3<CCTK_REAL, 4>().size() == 4);
}

} // namespace Arith
