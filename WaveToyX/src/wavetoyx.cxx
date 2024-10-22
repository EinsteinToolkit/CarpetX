#include <loop_device.hxx>

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>

namespace WaveToyX {

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

extern "C" void WaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, p.x, p.y, p.z, u(p.I),
                        rho(p.I));
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u(p.I),
                   rho(p.I));
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

#if 0
extern "C" void WaveToyX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::abs, std::max;
        CCTK_REAL maxabs_du = 0;
        for (int d = 0; d < dim; ++d) {
          const int ni = d == 0 ? 1 : 2;
          const int nj = d == 1 ? 1 : 2;
          const int nk = d == 2 ? 1 : 2;
          for (int dk = 0; dk < nk; ++dk) {
            for (int dj = 0; dj < nj; ++dj) {
              for (int di = 0; di < ni; ++di) {
                const auto I = p.I + di * p.DI[0] + dj * p.DI[1] + dk * p.DI[2];
                maxabs_du = max(maxabs_du, abs(u(I + p.DI[d]) - u(I)));
              }
            }
          }
        }
        regrid_error(p.I) = maxabs_du;
      });
}
#endif

extern "C" void WaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < dim; ++d)
    if (cctk_nghostzones[d] < fd_order / 2)
      CCTK_ERROR("Too few ghost zones");

  if (CCTK_EQUALS(boundary_condition, "CarpetX")) {

    switch (fd_order) {

    case 2:
      grid.loop_int_device<0, 0, 0>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p)
              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                using std::pow;
                CCTK_REAL ddu = 0;
                for (int d = 0; d < dim; ++d)
                  ddu += (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
                         pow(p.DX[d], 2);

                u_rhs(p.I) = rho(p.I);
                rho_rhs(p.I) = ddu;
              });
      break;

    case 4:
      grid.loop_int_device<0, 0, 0>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p)
              CCTK_ATTRIBUTE_ALWAYS_INLINE {
                using std::pow;
                CCTK_REAL ddu = 0;
                for (int d = 0; d < dim; ++d)
                  ddu += (-1 / CCTK_REAL(12) * u(p.I - 2 * p.DI[d]) +
                          4 / CCTK_REAL(3) * u(p.I - p.DI[d]) -
                          5 / CCTK_REAL(2) * u(p.I) +
                          4 / CCTK_REAL(3) * u(p.I + p.DI[d]) -
                          1 / CCTK_REAL(2) * u(p.I + 2 * p.DI[d])) /
                         pow(p.DX[d], 2);

                u_rhs(p.I) = rho(p.I);
                rho_rhs(p.I) = ddu;
              });
      break;

    default:
      CCTK_ERROR("Unsupported finite differencing order");
      break;
    }

  } else if (CCTK_EQUALS(boundary_condition, "reflecting")) {

    if (fd_order != 2)
      CCTK_ERROR("Unsupported finite differencing order");

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::pow;

          Arith::vect<CCTK_REAL, dim> ddu;
          for (int d = 0; d < dim; ++d)
            if (p.BI[d] < 0)
              ddu[d] = (u(p.I) - 2 * u(p.I + p.DI[d]) + u(p.I + 2 * p.DI[d])) /
                       pow(p.DX[d], 2);
            else if (p.BI[d] > 0)
              ddu[d] = (u(p.I - 2 * p.DI[d]) - 2 * u(p.I - p.DI[d]) + u(p.I)) /
                       pow(p.DX[d], 2);
            else
              ddu[d] = (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
                       pow(p.DX[d], 2);

          u_rhs(p.I) = rho(p.I);
          rho_rhs(p.I) = ddu[0] + ddu[1] + ddu[2];
        });

  } else if (CCTK_EQUALS(boundary_condition, "radiative")) {

    if (fd_order != 2)
      CCTK_ERROR("Unsupported finite differencing order");

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // Scalar wave equation in a first-order formulation:
          //     rho = dt u
          //     v   = grad u

          //     dot u   = rho
          //     dot rho = div v
          //     dot v   = grad rho

          // Non-trivial characteristics (with outward normal n):
          //     a = rho + n v   (incoming)
          //     b = rho - n v   (outgoing)

          //     rho = (a + b) / 2
          //     n v = (a - b) / 2

          // Set incoming mode to zero ("radiative"):
          //     a = 0

          // Calculate new state variables:

          //     rho_new = (a + b) / 2
          //     n v_new = (a - b) / 2

          //     rho_new = + b / 2
          //     n v_new = - b / 2

          //     rho_new = + (rho - n v) / 2
          //     n v_new = - (rho - n v) / 2

          // We're actually using a second order formulation.
          // Insert the definitions for rho and v:

          //     dot u_new    = + (dot u - n grad u) / 2
          //     n grad u_new = - (dot u - n grad u) / 2

          // Take a time derivative:

          //     dot rho_new = + (dot rho - n grad rho) / 2
          //     n grad rho_new = - (dot rho - n grad rho) / 2

          // Extract the definitions for dot rho and dot u from this.

          using std::pow;

          assert(all(p.loop_min <= p.I && p.I < p.loop_max));
          assert(all(if_else(p.BI < 0, p.I == p.bnd_min,
                             if_else(p.BI > 0, p.I == p.bnd_max - 1,
                                     p.bnd_min < p.I && p.I < p.bnd_max - 1))));

          Arith::vect<CCTK_REAL, dim> ddu;
          for (int d = 0; d < dim; ++d)
            if (p.BI[d] < 0)
              ddu[d] = (u(p.I) - 2 * u(p.I + p.DI[d]) + u(p.I + 2 * p.DI[d])) /
                       pow(p.DX[d], 2);
            else if (p.BI[d] > 0)
              ddu[d] = (u(p.I - 2 * p.DI[d]) - 2 * u(p.I - p.DI[d]) + u(p.I)) /
                       pow(p.DX[d], 2);
            else
              ddu[d] = (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
                       pow(p.DX[d], 2);

          Arith::vect<CCTK_REAL, dim> du;
          for (int d = 0; d < dim; ++d)
            if (p.BI[d] < 0)
              du[d] = (u(p.I + p.DI[d]) - u(p.I)) / p.DX[d];
            else if (p.BI[d] > 0)
              du[d] = (u(p.I) - u(p.I - p.DI[d])) / p.DX[d];
            else
              du[d] = (u(p.I + p.DI[d]) - u(p.I - p.DI[d])) / (2 * p.DX[d]);

          Arith::vect<CCTK_REAL, dim> drho;
          for (int d = 0; d < dim; ++d)
            if (p.BI[d] < 0)
              drho[d] = (rho(p.I + p.DI[d]) - rho(p.I)) / p.DX[d];
            else if (p.BI[d] > 0)
              drho[d] = (rho(p.I) - rho(p.I - p.DI[d])) / p.DX[d];
            else
              drho[d] =
                  (rho(p.I + p.DI[d]) - rho(p.I - p.DI[d])) / (2 * p.DX[d]);

          u_rhs(p.I) = rho(p.I);
          rho_rhs(p.I) = ddu[0] + ddu[1] + ddu[2];

          if (any(p.BI != 0)) {
            // We are on the boundary
            u_rhs(p.I) = (u_rhs(p.I) - dot(p.BI, du)) / 2;
            rho_rhs(p.I) = (rho_rhs(p.I) - dot(p.BI, drho)) / 2;
          }
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

extern "C" void WaveToyX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Energy;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < dim; ++d)
    if (cctk_nghostzones[d] < fd_order / 2)
      CCTK_ERROR("Too few ghost zones");

  switch (fd_order) {

  case 2:
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::pow;
          CCTK_REAL du2 = 0;
          for (int d = 0; d < dim; ++d)
            du2 +=
                pow((u(p.I + p.DI[d]) - u(p.I - p.DI[d])) / (2 * p.DX[d]), 2);

          eps(p.I) = (pow(rho(p.I), 2) + du2) / 2;
        });
    break;

  case 4:
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::pow;
          CCTK_REAL du2 = 0;
          for (int d = 0; d < dim; ++d)
            du2 += pow((1 / CCTK_REAL(12) * u(p.I - 2 * p.DI[d]) -
                        2 / CCTK_REAL(3) * u(p.I - p.DI[d]) +
                        2 / CCTK_REAL(3) * u(p.I + p.DI[d]) -
                        1 / CCTK_REAL(12) * u(p.I + 2 * p.DI[d])) /
                           p.DX[d],
                       2);

          eps(p.I) = (pow(rho(p.I), 2) + du2) / 2;
        });
    break;

  default:
    CCTK_ERROR("Unsupported finite differencing order");
    break;
  }
}

extern "C" void WaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, rho0;
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, p.x, p.y, p.z, u0, rho0);
          u_err(p.I) = u(p.I) - u0;
          rho_err(p.I) = rho(p.I) - rho0;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, rho0;
          gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u0,
                   rho0);
          u_err(p.I) = u(p.I) - u0;
          rho_err(p.I) = rho(p.I) - rho0;
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace WaveToyX
