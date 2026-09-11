namespace HybridMethods {

// Frozen coefficients for the 3-step, 2-stage fourth order method
// b3 and a32 are not listed here: solve.cxx derives them as b3 = 1 - (b0+b1+b2)
// and a32 = c3 - (a30+a31)

template <typename T> static inline auto hrk432_sol_1_c3() -> T {
  return T(9) / T(25);
}
template <typename T> static inline auto hrk432_sol_1_b0() -> T {
  return T(-85) / T(1416);
}
template <typename T> static inline auto hrk432_sol_1_b1() -> T {
  return T(131) / T(408);
}
template <typename T> static inline auto hrk432_sol_1_b2() -> T {
  return T(-29) / T(24);
}
template <typename T> static inline auto hrk432_sol_1_a30() -> T {
  return T(2511) / T(62500);
}
template <typename T> static inline auto hrk432_sol_1_a31() -> T {
  return T(-2268) / T(15625);
}

} // namespace HybridMethods
