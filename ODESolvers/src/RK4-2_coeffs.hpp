namespace MultiStepRungeKutta {

// Frozen coefficients for the two published solutions of the 2-step, 3-stage
// fourth order methods. b3, a21 and a32 are not listed here: solve.cxx derives
// them as b3 = 1 - (b0+b1+b2), a21 = c2 - a20, a32 = c3 - (a30+a31)

template <typename T> static inline auto rk4_dash_2_sol_1_c2() -> T {
  return T(7) / T(25);
}
template <typename T> static inline auto rk4_dash_2_sol_1_c3() -> T {
  return T(-13) / T(25);
}
template <typename T> static inline auto rk4_dash_2_sol_1_b0() -> T {
  return T(-643) / T(1536);
}
template <typename T> static inline auto rk4_dash_2_sol_1_b1() -> T {
  return T(-4237) / T(1092);
}
template <typename T> static inline auto rk4_dash_2_sol_1_b2() -> T {
  return T(38125) / T(10752);
}
template <typename T> static inline auto rk4_dash_2_sol_1_a20() -> T {
  return T(-49) / T(1250);
}
template <typename T> static inline auto rk4_dash_2_sol_1_a30() -> T {
  return T(7033) / T(960000);
}
template <typename T> static inline auto rk4_dash_2_sol_1_a31() -> T {
  return T(-217633) / T(210000);
}

template <typename T> static inline auto rk4_dash_2_sol_2_c2() -> T {
  return T(-99) / T(50);
}
template <typename T> static inline auto rk4_dash_2_sol_2_c3() -> T {
  return T(101) / T(100);
}
template <typename T> static inline auto rk4_dash_2_sol_2_b0() -> T {
  return T(-191) / T(882);
}
template <typename T> static inline auto rk4_dash_2_sol_2_b1() -> T {
  return T(48241) / T(59994);
}
template <typename T> static inline auto rk4_dash_2_sol_2_b2() -> T {
  return T(193750) / T(4351347);
}
template <typename T> static inline auto rk4_dash_2_sol_2_a20() -> T {
  return T(1309) / T(15500);
}
template <typename T> static inline auto rk4_dash_2_sol_2_a30() -> T {
  return T(-241289) / T(5880000);
}
template <typename T> static inline auto rk4_dash_2_sol_2_a31() -> T {
  return T(22846301) / T(16170000);
}

} // namespace MultiStepRungeKutta
