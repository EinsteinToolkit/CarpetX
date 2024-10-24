#ifndef TEST_RKAB_STANDING_WAVE_HXX
#define TEST_RKAB_STANDING_WAVE_HXX

#include <cctk.h>

#include <cmath>

namespace TestRKAB::sw {

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE phi(T A, T kx, T ky, T kz, T t, T x,
                                             T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return A * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Pi(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * omega * pi * cos(2 * kx * pi * x) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * omega * pi * t);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dx(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * kx * pi * cos(2 * omega * pi * t) * cos(2 * ky * pi * y) *
         cos(2 * kz * pi * z) * sin(2 * kx * pi * x);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dy(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * ky * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * kz * pi * z) * sin(2 * ky * pi * y);
}

template <typename T>
static inline auto CCTK_HOST CCTK_DEVICE Dz(T A, T kx, T ky, T kz, T t, T x,
                                            T y, T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos, std::acos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * kz * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

} // namespace TestRKAB::sw

#endif // TEST_RKAB_STANDING_WAVE_HXX