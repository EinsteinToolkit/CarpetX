#include "cplx.hxx"

#include <cctk.h>

#include <functional>

namespace Arith {
using namespace std;

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestCplx() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
  // clang in HIP mode does not have constexpr for isnan
  // https://sep5.readthedocs.io/en/latest/Programming_Guides/HIP-porting-guide.html
#if(!defined(__CUDACC__) && !defined(__clang__))
  typedef cplx<CCTK_REAL> CREAL;
  constexpr equal_to<CCTK_REAL> eq;
  constexpr equal_to<CREAL> eqc;

  static_assert(eq(CREAL().real, 0));
  static_assert(eq(CREAL().imag, 0));

  static_assert(eq(CREAL(1).real, 1));
  static_assert(eq(CREAL(1).imag, 0));

  static_assert(eq(CREAL(1, 2).real, 1));
  static_assert(eq(CREAL(1, 2).imag, 2));

  static_assert(eqc(CREAL(1, 2), CREAL(1, 2)));
  static_assert(!eqc(CREAL(1, 2), CREAL(2, 3)));

  static_assert(eqc(+CREAL(1, 2), CREAL(1, 2)));
  static_assert(eqc(-CREAL(1, 2), CREAL(-1, -2)));

  static_assert(eqc(CREAL(1, 2) + CREAL(3, 4), CREAL(4, 6)));
  static_assert(eqc(CREAL(1, 2) - CREAL(3, 4), CREAL(-2, -2)));
  static_assert(eqc(2 * CREAL(3, 4), CREAL(6, 8)));
  static_assert(eqc(CREAL(3, 4) * 2, CREAL(6, 8)));
  static_assert(eqc(CREAL(3, 4) / 2, CREAL(1.5, 2)));

  static_assert(eqc(CREAL(2, 3) * CREAL(4, 5), CREAL(-7, 22)));
  static_assert(eqc(CREAL(4, 5) / CREAL(3, 4), CREAL(1.28, -0.04)));

  // static_assert(eqc(sqrt(CREAL(4, 3)), CREAL(2, 0.75)));

  static_assert(eqc(pow2(CREAL(1, 2)), CREAL(1, 2) * CREAL(1, 2)));

  static_assert(eqc(pow(CREAL(1, 2), 0), CREAL(1, 0)));
  static_assert(eqc(pow(CREAL(1, 2), 1), CREAL(1, 2)));
  static_assert(eqc(pow(CREAL(1, 2), 2), pow2(CREAL(1, 2))));
  static_assert(eqc(pow(CREAL(1, 2), 3), CREAL(1, 2) * pow2(CREAL(1, 2))));
  static_assert(eqc(pow(CREAL(1, 2), 4), pow2(pow2(CREAL(1, 2)))));
#endif
}

} // namespace Arith
