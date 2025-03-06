#include <cmath>

namespace HybridMethods {

template <typename T> static inline auto hrk432_sol_1_b0(T c3) -> T {
  using std::pow;
  return 0.4166666666666667 - 9 / (8. * (2 + c3));
}

template <typename T> static inline auto hrk432_sol_1_b1(T c3) -> T {
  using std::pow;
  return -1.3333333333333333 + 9 / (4. * (1 + c3));
}

template <typename T> static inline auto hrk432_sol_1_b2(T c3) -> T {
  using std::pow;
  return 1.9166666666666667 - 9 / (8. * c3);
}

template <typename T> static inline auto hrk432_sol_1_a30(T c3) -> T {
  using std::pow;
  return (pow(c3, 2) * (3 + 2 * c3)) / 12.;
}

template <typename T> static inline auto hrk432_sol_1_a31(T c3) -> T {
  using std::pow;
  return -0.3333333333333333 * (pow(c3, 2) * (3 + c3));
}

} // namespace HybridMethods