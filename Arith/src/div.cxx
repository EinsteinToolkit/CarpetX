#include "div.hxx"

namespace Arith {

template <typename T> constexpr bool test_div_floor(T x, T y) {
  T d = div_floor(x, y);
  T m = mod_floor(x, y);
  T absy = y >= 0 ? y : -y;
  return d * y + m == x && m >= 0 && m < absy;
}

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestDiv() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
  // clang in HIP mode does not have constexpr for isnan
  // https://sep5.readthedocs.io/en/latest/Programming_Guides/HIP-porting-guide.html
#if(!defined(__CUDACC__) && !defined(__clang__))
  static_assert(test_div_floor(0, 1));
  static_assert(test_div_floor(0, 2));
  static_assert(test_div_floor(0, 3));

  static_assert(test_div_floor(1, 1));
  static_assert(test_div_floor(1, 2));
  static_assert(test_div_floor(1, 3));
  static_assert(test_div_floor(2, 1));
  static_assert(test_div_floor(2, 2));
  static_assert(test_div_floor(2, 3));
  static_assert(test_div_floor(3, 1));
  static_assert(test_div_floor(3, 2));
  static_assert(test_div_floor(3, 3));

  static_assert(test_div_floor(-1, 1));
  static_assert(test_div_floor(-1, 2));
  static_assert(test_div_floor(-1, 3));
  static_assert(test_div_floor(-2, 1));
  static_assert(test_div_floor(-2, 2));
  static_assert(test_div_floor(-2, 3));
  static_assert(test_div_floor(-3, 1));
  static_assert(test_div_floor(-3, 2));
  static_assert(test_div_floor(-3, 3));
#endif
}

} // namespace Arith
