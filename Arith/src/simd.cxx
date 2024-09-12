#include "simd.hxx"

#include <cctk.h>

#include <cassert>
#include <cmath>

namespace Arith {

std::size_t flop_count, memop_count;

void reset_counts() {
#pragma omp parallel
  flop_count = memop_count = 0;
}
std::size_t get_flop_count() {
  std::size_t count = 0;
#pragma omp parallel reduction(+ : count)
  count += flop_count;
  return count;
}
std::size_t get_memop_count() {
  std::size_t count = 0;
#pragma omp parallel reduction(+ : count)
  count += memop_count;
  return count;
}

template <typename T> constexpr bool isequal(const simd<T> x, const T a) {
  return all(x == a);
}
template <typename T> constexpr bool isequal(const simdl<T> x, const bool a) {
  return all(x == a);
}
template <typename T> constexpr bool isapprox(const simd<T> x, const T a) {
  return all(fabs(x - a) < 1.0e-14);
}

void check(const bool isgood) {
  if (isgood)
    return;
  CCTK_ERROR("Test failure");
}

void TestSIMD() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
  using real = CCTK_REAL;
  using realv = simd<real>;

  realv x;
  realv y = 0;
  realv z = zero<realv>();

  check(all(y == 0));
  check(all(y == z));

  realv a = 2;
  realv b = 3;
  realv c = 4;

  realv r0 = muladd(a, b, c);
  realv r1 = mulsub(a, b, c);
  realv r2 = negmuladd(a, b, c);
  realv r3 = negmulsub(a, b, c);
  check(all(r0 == muladd(2, 3, 4)));
  check(all(r1 == mulsub(2, 3, 4)));
  check(all(r2 == negmuladd(2, 3, 4)));
  check(all(r3 == negmulsub(2, 3, 4)));

  real s = 2;
  real t = 2;
  real u = 4;

  check(isequal(+a, +s));
  check(isequal(-a, -s));

  check(isequal(a + b, s + t));
  check(isequal(a - b, s - t));
  check(isequal(a * b, s * t));
  check(isequal(a / b, s / t));
  check(isequal(s + b, s + t));
  check(isequal(s - b, s - t));
  check(isequal(s * b, s * t));
  check(isequal(s / b, s / t));
  check(isequal(a + t, s + t));
  check(isequal(a - t, s - t));
  check(isequal(a * t, s * t));
  check(isequal(a / t, s / t));

  check(isequal(a == b, s == t));
  check(isequal(a != b, s != t));
  check(isequal(a < b, s < t));
  check(isequal(a > b, s > t));
  check(isequal(a <= b, s <= t));
  check(isequal(a >= b, s >= t));
  check(isequal(s == b, s == t));
  check(isequal(s != b, s != t));
  check(isequal(s < b, s < t));
  check(isequal(s > b, s > t));
  check(isequal(s <= b, s <= t));
  check(isequal(s >= b, s >= t));
  check(isequal(a == t, s == t));
  check(isequal(a != t, s != t));
  check(isequal(a < t, s < t));
  check(isequal(a > t, s > t));
  check(isequal(a <= t, s <= t));
  check(isequal(a >= t, s >= t));

  check(isapprox(abs(a), abs(s)));
  check(isapprox(acos(a), acos(s)));
  check(isapprox(acosh(a), acosh(s)));
  check(allisfinite(a) == allisfinite(s));
  check(anyisnan(a) == anyisnan(s));
  check(isapprox(asin(a), asin(s)));
  check(isapprox(asinh(a), asinh(s)));
  check(isapprox(atan(a), atan(s)));
  check(isapprox(atanh(a), atanh(s)));
  check(isapprox(cbrt(a), cbrt(s)));
  // check(isapprox(cis(a), cis(s)));
  // check(isapprox(cispi(a), cispi(s)));
  check(isapprox(copysign(a, b), copysign(s, t)));
  check(isapprox(cos(a), cos(s)));
  check(isapprox(cosh(a), cosh(s)));
  check(isapprox(cospi(a), cos(CCTK_REAL(M_PI) * s)));
  check(isapprox(exp(a), exp(s)));
  check(isapprox(exp10(a), exp10(s)));
  check(isapprox(exp2(a), exp2(s)));
  check(isapprox(fabs(a), fabs(s)));
  check(isapprox(flipsign(a, b), flipsign(s, t)));
  check(isapprox(fmax(a, b), fmax(s, t)));
  check(isapprox(fmax(a, t), fmax(s, t)));
  check(isapprox(fmax(s, b), fmax(s, t)));
  check(isapprox(fmin(a, b), fmin(s, t)));
  check(isapprox(fmin(a, t), fmin(s, t)));
  check(isapprox(fmin(s, b), fmin(s, t)));
  check(isapprox(hypot(a, b), hypot(s, t)));
  check(isapprox(inv(a), inv(s)));
  check(isequal(isfinite(a), isfinite(s)));
  check(isequal(isinf(a), isinf(s)));
  check(isequal(isnan(a), isnan(s)));
  check(isapprox(log(a), log(s)));
  check(isapprox(log10(a), log10(s)));
  check(isapprox(log2(a), log2(s)));
  check(isapprox(max(a, b), max(s, t)));
  check(isapprox(max(a, t), max(s, t)));
  check(isapprox(max(s, b), max(s, t)));
  check(isapprox(min(a, b), min(s, t)));
  check(isapprox(min(a, t), min(s, t)));
  check(isapprox(min(s, b), min(s, t)));
  check(isapprox(muladd(a, b, c), muladd(s, t, u)));
  check(isapprox(muladd(a, b, u), muladd(s, t, u)));
  check(isapprox(muladd(a, t, c), muladd(s, t, u)));
  check(isapprox(muladd(a, t, u), muladd(s, t, u)));
  check(isapprox(muladd(s, b, c), muladd(s, t, u)));
  check(isapprox(muladd(s, b, u), muladd(s, t, u)));
  check(isapprox(muladd(s, t, c), muladd(s, t, u)));
  check(isapprox(mulsub(a, b, c), mulsub(s, t, u)));
  check(isapprox(mulsub(a, b, u), mulsub(s, t, u)));
  check(isapprox(mulsub(a, t, c), mulsub(s, t, u)));
  check(isapprox(mulsub(a, t, u), mulsub(s, t, u)));
  check(isapprox(mulsub(s, b, c), mulsub(s, t, u)));
  check(isapprox(mulsub(s, b, u), mulsub(s, t, u)));
  check(isapprox(mulsub(s, t, c), mulsub(s, t, u)));
  check(isapprox(negmuladd(a, b, c), negmuladd(s, t, u)));
  check(isapprox(negmuladd(a, b, u), negmuladd(s, t, u)));
  check(isapprox(negmuladd(a, t, c), negmuladd(s, t, u)));
  check(isapprox(negmuladd(a, t, u), negmuladd(s, t, u)));
  check(isapprox(negmuladd(s, b, c), negmuladd(s, t, u)));
  check(isapprox(negmuladd(s, b, u), negmuladd(s, t, u)));
  check(isapprox(negmuladd(s, t, c), negmuladd(s, t, u)));
  check(isapprox(negmulsub(a, b, c), negmulsub(s, t, u)));
  check(isapprox(negmulsub(a, b, u), negmulsub(s, t, u)));
  check(isapprox(negmulsub(a, t, c), negmulsub(s, t, u)));
  check(isapprox(negmulsub(a, t, u), negmulsub(s, t, u)));
  check(isapprox(negmulsub(s, b, c), negmulsub(s, t, u)));
  check(isapprox(negmulsub(s, b, u), negmulsub(s, t, u)));
  check(isapprox(negmulsub(s, t, c), negmulsub(s, t, u)));
  check(isapprox(pow(a, b), pow(s, t)));
  check(isapprox(pow(a, t), pow(s, t)));
  check(isapprox(pow(s, b), pow(s, t)));
  check(isequal(signbit(a), signbit(s)));
  check(isapprox(sin(a), sin(s)));
  check(isapprox(sinh(a), sinh(s)));
  check(isapprox(sinpi(a), sin(CCTK_REAL(M_PI) * s)));
  check(isapprox(sqrt(a), sqrt(s)));
  check(isapprox(tan(a), tan(s)));
  check(isapprox(tanh(a), tanh(s)));
  // check(isapprox(tanpi(a), tanpi(s)));

#endif
}

} // namespace Arith
