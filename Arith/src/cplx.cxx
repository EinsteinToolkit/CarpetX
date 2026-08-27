#include "cplx.hxx"
#include "simd.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>

using namespace std;

namespace Arith {

// This function is compiled, but not executed. The tests are "run" at
// compile time. If this function compiles, the tests pass.
void TestCplx() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
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
  // `operator==` must consider both components too, not just the real part.
  // (`eqc` is `std::equal_to`, which is specialised separately and already
  // compared both, so it would not have caught this.)
  static_assert(CREAL(1, 2) == CREAL(1, 2));
  static_assert(CREAL(1, 2) != CREAL(1, 3));
  static_assert(CREAL(1, 2) != CREAL(2, 2));

  static_assert(eqc(+CREAL(1, 2), CREAL(1, 2)));
  static_assert(eqc(-CREAL(1, 2), CREAL(-1, -2)));

  static_assert(eqc(CREAL(1, 2) + CREAL(3, 4), CREAL(4, 6)));
  static_assert(eqc(CREAL(1, 2) - CREAL(3, 4), CREAL(-2, -2)));
  static_assert(eqc(2 * CREAL(3, 4), CREAL(6, 8)));
  static_assert(eqc(CREAL(3, 4) * 2, CREAL(6, 8)));
  static_assert(eqc(CREAL(3, 4) / 2, CREAL(1.5, 2)));

  static_assert(eqc(CREAL(2, 3) * CREAL(4, 5), CREAL(-7, 22)));
  static_assert(eqc(CREAL(4, 5) / CREAL(3, 4), CREAL(1.28, -0.04)));

  // The transcendental functions cannot be checked here because `std::sqrt`
  // and friends are not constexpr; see `test_cplx` below.

  static_assert(eqc(pow2(CREAL(1, 2)), CREAL(1, 2) * CREAL(1, 2)));

  static_assert(eqc(pow(CREAL(1, 2), 0), CREAL(1, 0)));
  static_assert(eqc(pow(CREAL(1, 2), 1), CREAL(1, 2)));
  static_assert(eqc(pow(CREAL(1, 2), 2), pow2(CREAL(1, 2))));
  static_assert(eqc(pow(CREAL(1, 2), 3), CREAL(1, 2) * pow2(CREAL(1, 2))));
  static_assert(eqc(pow(CREAL(1, 2), 4), pow2(pow2(CREAL(1, 2)))));

  // Negative exponents invert
  static_assert(eqc(pow(CREAL(0, 2), -1), CREAL(0, -0.5)));
  static_assert(eqc(pow(CREAL(0, 2), -2), CREAL(-0.25, 0)));
  static_assert(eqc(pow(CREAL(1, 0), -3), CREAL(1, 0)));
#endif
}

namespace {

// Compiled, but not executed. Instantiates every `cplx` function for both a
// scalar and a SIMD element type. Comparisons on SIMD types yield masks rather
// than booleans, so all of these have to be branch-free; if this compiles,
// they are.
template <typename T> void instantiate_cplx(const cplx<T> &z, const int n) {
  const cplx<T> r = exp(z) + sqrt(z) + cos(z) + sin(z) + cbrt(z) + conj(z) +
                    pow2(z) + pow(z, n);
  const cplx<T> s(abs(z) + arg(z) + norm(z) + real(z) + imag(z));
  // Consume both results
  if (!allisfinite(r) && anyisnan(s))
    CCTK_ERROR("never reached");
}

using std::complex;

// Relative error, absolute for small values
CCTK_REAL relerr(const complex<CCTK_REAL> &got,
                 const complex<CCTK_REAL> &want) {
  return std::abs(got - want) / std::max(CCTK_REAL(1), std::abs(want));
}
complex<CCTK_REAL> conv(const cplx<CCTK_REAL> &z) { return {z.real, z.imag}; }

// Check the transcendental functions against `std::complex`, which defines the
// branch cuts we want. `std::complex` is the reference, not the
// implementation: `cplx` has to work for SIMD element types as well.
void test_cplx() {
  typedef CCTK_REAL R;
  typedef cplx<R> C;
  typedef complex<R> S;
  const R eps = std::numeric_limits<R>::epsilon();
  const R pi = std::acos(R(-1));
  const R vals[] = {0, 1, -1, R(0.5), R(-0.5), 2, -2, R(3.7), R(-3.7)};

  for (const R re : vals) {
    for (const R im : vals) {
      const C z(re, im);
      const S w(re, im);

      assert(std::abs(abs(z) - std::abs(w)) <=
             16 * eps * std::max(R(1), std::abs(w)));
      assert(std::abs(arg(z) - std::arg(w)) <= 16 * eps);
      assert(relerr(conv(exp(z)), std::exp(w)) <= 16 * eps);
      assert(relerr(conv(sqrt(z)), std::sqrt(w)) <= 16 * eps);
      assert(relerr(conv(cos(z)), std::cos(w)) <= 16 * eps);
      assert(relerr(conv(sin(z)), std::sin(w)) <= 16 * eps);

      // `std::complex` has no `cbrt`, so check the defining property and the
      // choice of branch instead
      const C c = cbrt(z);
      assert(relerr(conv(c * c * c), w) <= 64 * eps);
      assert(std::abs(std::arg(conv(c))) <= pi / 3 + 16 * eps);

      for (int n = -4; n <= 4; ++n) {
        if (re == 0 && im == 0 && n <= 0)
          continue;
        assert(relerr(conv(pow(z, n)), std::pow(w, n)) <= 512 * eps);
      }
    }
  }

  // Exact values, including the principal branch of `sqrt` and its cut along
  // the negative real axis
  assert(sqrt(C(0, 0)) == C(0, 0));
  assert(sqrt(C(4, 0)) == C(2, 0));
  assert(sqrt(C(-1, 0)) == C(0, 1));
  assert(sqrt(C(-4, 0)) == C(0, 2));
  assert(sqrt(C(-1, -R(0.0))) == C(0, -1));
  assert(abs(C(3, 4)) == 5);
  assert(relerr(conv(exp(C(0, pi))), S(-1, 0)) <= 16 * eps);
  assert(relerr(conv(cbrt(C(-8, 0))), S(1, std::sqrt(R(3)))) <= 16 * eps);
}

} // namespace

extern "C" void Test_cplx(CCTK_ARGUMENTS) {
#ifndef __CUDACC__
  CCTK_INFO("Test_cplx");

  test_cplx();
#endif
}

// Referenced so that `instantiate_cplx` is compiled for both element types
void InstantiateCplx(const cplx<CCTK_REAL> &z, const cplx<simd<CCTK_REAL> > &vz,
                     const int n) {
  instantiate_cplx(z, n);
  instantiate_cplx(vz, n);
}

} // namespace Arith
