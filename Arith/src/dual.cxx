#include "dual.hxx"

#include <cctk.h>

#include <functional>

using namespace std;

namespace Arith {

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestDual() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
  typedef dual<CCTK_REAL> DREAL;
  constexpr equal_to<CCTK_REAL> eq;
  constexpr equal_to<DREAL> eqd;

  static_assert(eq(DREAL(0).val, 0));
  static_assert(eq(DREAL(0).eps, 0));

  static_assert(eq(DREAL(1).val, 1));
  static_assert(eq(DREAL(1).eps, 0));

  static_assert(eq(DREAL(1, 2).val, 1));
  static_assert(eq(DREAL(1, 2).eps, 2));

  static_assert(eqd(DREAL(1, 2), DREAL(1, 2)));
  static_assert(!eqd(DREAL(1, 2), DREAL(2, 3)));

  static_assert(eqd(+DREAL(1, 2), DREAL(1, 2)));
  static_assert(eqd(-DREAL(1, 2), DREAL(-1, -2)));

  static_assert(eqd(DREAL(1, 2) + DREAL(3, 4), DREAL(4, 6)));
  static_assert(eqd(DREAL(1, 2) - DREAL(3, 4), DREAL(-2, -2)));
  static_assert(eqd(2 * DREAL(3, 4), DREAL(6, 8)));
  static_assert(eqd(DREAL(3, 4) * 2, DREAL(6, 8)));
  static_assert(eqd(DREAL(3, 4) / 2, DREAL(1.5, 2)));

  static_assert(eqd(DREAL(2, 3) * DREAL(4, 5), DREAL(8, 22)));
  static_assert(eqd(DREAL(4, 5) / DREAL(2, 3), DREAL(2, -0.5)));

  // static_assert(eqd(sqrt(DREAL(4, 3)), DREAL(2, 0.75)));
#endif
}

} // namespace Arith
