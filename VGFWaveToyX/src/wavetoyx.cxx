#include <loop_device.hxx>

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>

namespace VGFWaveToyX {

constexpr int dim = 3;

// u(t,x,y,z) =
//   A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
template <typename T>
constexpr void standing_wave(const T A, const T kx, const T ky, const T kz,
                             const T t, const T x, const T y, const T z, T &u,
                             T &rho) {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  rho = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
        cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t, const T x, const T y,
                        const T z, T &u, T &rho) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    rho = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
  } else {
    u = (f(t - r) - f(t + r)) / r;
    rho = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
  }
}

extern "C" void VGFWaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VGFWaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          for (int n = 0; n < nvars; ++n)
            standing_wave(amplitude[n], standing_wave_kx[n],
                          standing_wave_ky[n], standing_wave_kz[n], cctk_time,
                          p.x, p.y, p.z, u(p.I, n), rho(p.I, n));
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          for (int n = 0; n < nvars; ++n)
            gaussian(amplitude[n], gaussian_width[n], cctk_time, p.x, p.y, p.z,
                     u(p.I, n), rho(p.I, n));
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void VGFWaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VGFWaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < dim; ++d)
    if (cctk_nghostzones[d] < 1)
      CCTK_ERROR("Too few ghost zones");

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int n = 0; n < nvars; ++n) {
          using std::pow;
          CCTK_REAL ddu = 0;
          for (int d = 0; d < dim; ++d)
            ddu += (u(p.I - p.DI[d], n) - 2 * u(p.I, n) + u(p.I + p.DI[d], n)) /
                   pow(p.DX[d], 2);

          u_rhs(p.I, n) = rho(p.I, n);
          rho_rhs(p.I, n) = ddu;
        }
      });
}

extern "C" void VGFWaveToyX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VGFWaveToyX_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void VGFWaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VGFWaveToyX_Energy;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < dim; ++d)
    if (cctk_nghostzones[d] < 1)
      CCTK_ERROR("Too few ghost zones");

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL eps1 = 0;
        for (int n = 0; n < nvars; ++n) {
          using std::pow;
          CCTK_REAL du2 = 0;
          for (int d = 0; d < dim; ++d)
            du2 += pow(
                (u(p.I + p.DI[d], n) - u(p.I - p.DI[d], n)) / (2 * p.DX[d]), 2);

          eps1 += (pow(rho(p.I, n), 2) + du2) / 2;
        }
        eps(p.I) = eps1;
      });
}

extern "C" void VGFWaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_VGFWaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          for (int n = 0; n < nvars; ++n) {
            CCTK_REAL u0, rho0;
            standing_wave(amplitude[n], standing_wave_kx[n],
                          standing_wave_ky[n], standing_wave_kz[n], cctk_time,
                          p.x, p.y, p.z, u0, rho0);
            u_err(p.I, n) = u(p.I, n) - u0;
            rho_err(p.I, n) = rho(p.I, n) - rho0;
          }
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          for (int n = 0; n < nvars; ++n) {
            CCTK_REAL u0, rho0;
            gaussian(amplitude[n], gaussian_width[n], cctk_time, p.x, p.y, p.z,
                     u0, rho0);
            u_err(p.I, n) = u(p.I, n) - u0;
            rho_err(p.I, n) = rho(p.I, n) - rho0;
          }
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace VGFWaveToyX
