#include <cstddef>
#include <iterator>
#include <loop_device.hxx>

#include <sum.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace TestRKAB {
using namespace Arith;

constexpr int dim = 3;

// Standing wave functions
template <typename T>
static inline auto sw_phi(T A, T kx, T ky, T kz, T t, T x, T y,
                          T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return A * cos(2 * omega * pi * t) * sin(2 * kx * pi * x) *
         sin(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto sw_Pi(T A, T kx, T ky, T kz, T t, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return -2 * A * omega * pi * sin(2 * omega * pi * t) * sin(2 * kx * pi * x) *
         sin(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto sw_Dx(T A, T kx, T ky, T kz, T t, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 2 * A * kx * pi * cos(2 * omega * pi * t) * cos(2 * kx * pi * x) *
         sin(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto sw_Dy(T A, T kx, T ky, T kz, T t, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 2 * A * ky * pi * cos(2 * omega * pi * t) * sin(2 * kx * pi * x) *
         cos(2 * ky * pi * y) * sin(2 * kz * pi * z);
}

template <typename T>
static inline auto sw_Dz(T A, T kx, T ky, T kz, T t, T x, T y,
                         T z) noexcept -> T {
  using std::sqrt, std::sin, std::cos;
  const auto pi{acos(T{-1})};
  const auto omega{sqrt(kx * kx + ky * ky + kz * kz)};
  return 2 * A * kz * pi * cos(2 * omega * pi * t) * sin(2 * kx * pi * x) *
         sin(2 * ky * pi * y) * cos(2 * kz * pi * z);
}

// Finite difference helpers
enum class fd_dir : std::size_t { x = 0, y = 1, z = 2 };

template <fd_dir direction, typename T>
static inline auto fd_c_1_4(const Loop::PointDesc &p,
                            const Loop::GF3D2<T> &gf) noexcept -> T {
  // (1*f[i-2]-8*f[i-1]+8*f[i+1]-1*f[i+2])
  constexpr auto d{static_cast<size_t>(direction)};
  const auto num{gf(p.I - 2 * p.DI[d]) - 8.0 * gf(p.I - 1 * p.DI[d]) +
                 8.0 * gf(p.I + 1 * p.DI[d]) - 1.0 * gf(p.I + 2 * p.DI[d])};
  const auto den{1.0 / (12.0 * p.DX[d])};
  return den * num;
}

extern "C" void TestRKAB_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};
          const auto t_pre{t - cctk_delta_time};

          const auto A{amplitude};
          const auto kx{standing_wave_kx};
          const auto ky{standing_wave_ky};
          const auto kz{standing_wave_kz};

          phi(p.I) = sw_phi(A, kx, ky, kz, t, p.x, p.y, p.z);
          Pi(p.I) = sw_Pi(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dx(p.I) = sw_Dx(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dy(p.I) = sw_Dy(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dz(p.I) = sw_Dz(A, kx, ky, kz, t, p.x, p.y, p.z);

          phi_pre(p.I) = sw_phi(A, kx, ky, kz, t_pre, p.x, p.y, p.z);
          Pi_pre(p.I) = sw_Pi(A, kx, ky, kz, t_pre, p.x, p.y, p.z);
          Dx_pre(p.I) = sw_Dx(A, kx, ky, kz, t_pre, p.x, p.y, p.z);
          Dy_pre(p.I) = sw_Dy(A, kx, ky, kz, t_pre, p.x, p.y, p.z);
          Dz_pre(p.I) = sw_Dz(A, kx, ky, kz, t_pre, p.x, p.y, p.z);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    CCTK_ERROR("Gaussian initial data not implemented yet");
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void TestRKAB_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        phi_rhs(p.I) = Pi(p.I);
        Pi_rhs(p.I) = fd_c_1_4<fd_dir::x>(p, Dx) + fd_c_1_4<fd_dir::y>(p, Dx) +
                      fd_c_1_4<fd_dir::z>(p, Dz);
        Dx_rhs(p.I) = fd_c_1_4<fd_dir::x>(p, Pi);
        Dy_rhs(p.I) = fd_c_1_4<fd_dir::y>(p, Pi);
        Dz_rhs(p.I) = fd_c_1_4<fd_dir::z>(p, Pi);
      });
}

extern "C" void TestRKAB_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

extern "C" void TestRKAB_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};

          const auto A{amplitude};
          const auto kx{standing_wave_kx};
          const auto ky{standing_wave_ky};
          const auto kz{standing_wave_kz};

          const auto expected_phi{sw_phi(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Pi{sw_Pi(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dx{sw_Dx(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dy{sw_Dy(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dz{sw_Dz(A, kx, ky, kz, t, p.x, p.y, p.z)};

          const auto actual_phi{phi(p.I)};
          const auto actual_Pi{Pi(p.I)};
          const auto actual_Dx{Dx(p.I)};
          const auto actual_Dy{Dy(p.I)};
          const auto actual_Dz{Dz(p.I)};

          phi_err(p.I) = fabs(expected_phi - actual_phi);
          Pi_err(p.I) = fabs(expected_Pi - actual_Pi);
          Dx_err(p.I) = fabs(expected_Dx - actual_Dx);
          Dy_err(p.I) = fabs(expected_Dy - actual_Dy);
          Dz_err(p.I) = fabs(expected_Dz - actual_Dz);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    CCTK_VERROR("Gaussian error not implemented yet");
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace TestRKAB
