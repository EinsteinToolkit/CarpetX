#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace WaveToyX {

constexpr int dim = 3;

// u(t,x,y,z) =
//   A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
template <typename T>
constexpr void standing_wave(const T t, const T x, const T y, const T z, T &u,
                             T &rho) {
  DECLARE_CCTK_PARAMETERS;

  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T A = amplitude;
  const T kx = standing_wave_kx;
  const T ky = standing_wave_ky;
  const T kz = standing_wave_kz;
  const T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  rho = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
        cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

extern "C" void WaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          standing_wave(cctk_time, p.x, p.y, p.z, u(p.I), rho(p.I));
        });
  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void WaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_RHS;
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::pow;
        CCTK_REAL ddu = 0;
        for (int d = 0; d < dim; ++d)
          ddu += (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
                 pow(p.DX[d], 2);

        udot(p.I) = rho(p.I);
        rhodot(p.I) = ddu;
      });
}

extern "C" void WaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, rho0;
          standing_wave(cctk_time, p.x, p.y, p.z, u0, rho0);
          uerr(p.I) = u(p.I) - u0;
          rhoerr(p.I) = rho(p.I) - rho0;
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void WaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Energy;
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::pow;
        CCTK_REAL du2 = 0;
        for (int d = 0; d < dim; ++d)
          du2 += pow((u(p.I + p.DI[d]) - u(p.I + p.DI[d])) / (2 * p.DX[d]), 2);

        eps(p.I) = (pow(rho(p.I), 2) + du2) / 2;
      });
}

} // namespace WaveToyX
