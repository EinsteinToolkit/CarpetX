#include <cmath>

namespace HybridMethods {

template <typename T> static inline auto rk423_sol_1_b0(T c2, T c3) -> T {
  using std::pow;
  return (-3 + c2 * (4 - 6 * c3) + 4 * c3) / (12. * (1 + c2) * (1 + c3));
}

template <typename T> static inline auto rk423_sol_1_b1(T c2, T c3) -> T {
  using std::pow;
  return (7 - 10 * c3 + 2 * c2 * (-5 + 9 * c3)) / (12. * c2 * c3);
}

template <typename T> static inline auto rk423_sol_1_b2(T c2, T c3) -> T {
  using std::pow;
  return (7 - 10 * c3) / (12. * c2 * (1 + c2) * (c2 - c3));
}

template <typename T> static inline auto rk423_sol_1_a20(T c2, T c3) -> T {
  using std::pow;
  return -0.5 * pow(c2, 2);
}

template <typename T> static inline auto rk423_sol_1_a30(T c2, T c3) -> T {
  using std::pow;
  return (c3 *
          (7 * c2 * (3 + 2 * c2) - 3 * c2 * (-4 + 5 * c2 * (1 + 2 * c2)) * c3 -
           2 * (7 + 12 * c2) * pow(c3, 2))) /
         (6. * pow(1 + c2, 2) * (-7 + 10 * c2));
}

template <typename T> static inline auto rk423_sol_1_a31(T c2, T c3) -> T {
  using std::pow;
  return (c3 * (pow(c2, 2) * (4 - 15 * c3) + 30 * pow(c2, 3) * (2 + c3) +
                7 * c3 * (3 + 2 * c3) + 3 * c2 * (-21 + c3 * (-7 + 8 * c3)))) /
         (6. * c2 * (1 + c2) * (-7 + 10 * c2));
}

template <typename T> static inline auto rk423_sol_2_b0(T c2, T c3) -> T {
  using std::pow;
  return (-3 + c2 * (4 - 6 * c3) + 4 * c3) / (12. * (1 + c2) * (1 + c3));
}

template <typename T> static inline auto rk423_sol_2_b1(T c2, T c3) -> T {
  using std::pow;
  return (7 - 10 * c3 + 2 * c2 * (-5 + 9 * c3)) / (12. * c2 * c3);
}

template <typename T> static inline auto rk423_sol_2_b2(T c2, T c3) -> T {
  using std::pow;
  return (7 - 10 * c3) / (12. * c2 * (1 + c2) * (c2 - c3));
}

template <typename T> static inline auto rk423_sol_2_a20(T c2, T c3) -> T {
  using std::pow;
  return (c2 * (-21 + 2 * c2 * (7 + 12 * c3) + 4 * c3 * (8 + 15 * c3))) /
         (12. * (1 + c3) * (-7 + 10 * c3));
}

template <typename T> static inline auto rk423_sol_2_a30(T c2, T c3) -> T {
  using std::pow;
  return (c3 * (7 * (-3 + 8 * c2) - 2 * (5 + 6 * c2 * (1 + 5 * c2)) * c3 +
                12 * (8 - 5 * c2) * pow(c3, 2) + 120 * pow(c3, 3))) /
         (12. * (1 + c2) * (-7 + 10 * c2));
}

template <typename T> static inline auto rk423_sol_2_a31(T c2, T c3) -> T {
  using std::pow;
  return (c3 * (c2 * (-147 + 20 * c2 * (-1 + 6 * c2)) +
                2 * (42 + c2 * (23 + 6 * c2 * (1 + 5 * c2))) * c3 +
                12 * (1 + c2) * (-3 + 5 * c2) * pow(c3, 2) -
                120 * (1 + c2) * pow(c3, 3))) /
         (12. * c2 * (1 + c2) * (-7 + 10 * c2));
}

} // namespace HybridMethods