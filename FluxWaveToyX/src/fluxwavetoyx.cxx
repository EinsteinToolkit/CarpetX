#include <loop_device.hxx>

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace FluxWaveToyX {

constexpr int dim = 3;

// Average from cells to faces ("reconstruction")
template <typename T>
CCTK_DEVICE CCTK_HOST inline CCTK_ATTRIBUTE_ALWAYS_INLINE auto
average(const Loop::GF3D2<T> &var, const Loop::PointDesc &p, const int dir) {
  return (var(p.I - p.DI[dir]) + var(p.I)) / 2;
}

// Flux divergence
// Expects fluxes on faces, calculates divergence in cell
template <typename T>
CCTK_DEVICE CCTK_HOST inline CCTK_ATTRIBUTE_ALWAYS_INLINE auto
flux_div(const Loop::GF3D2<T> &flux_x, const Loop::GF3D2<T> &flux_y,
         const Loop::GF3D2<T> &flux_z, const Loop::PointDesc &p) {
  return (flux_x(p.I + p.DI[0]) - flux_x(p.I)) / p.DX[0] +
         (flux_y(p.I + p.DI[1]) - flux_y(p.I)) / p.DX[1] +
         (flux_z(p.I + p.DI[2]) - flux_z(p.I)) / p.DX[2];
}

// Derivative
template <typename T>
CCTK_DEVICE CCTK_HOST inline CCTK_ATTRIBUTE_ALWAYS_INLINE auto
deriv(const Loop::GF3D2<T> &var, const Loop::PointDesc &p, const int dir) {
  return (var(p.I + p.DI[dir]) - var(p.I - p.DI[dir])) / (2 * p.DX[dir]);
}

// u(t,x,y,z) =
//   A cos(2 pi omega t) cos(2 pi kx x) cos(2 pi ky y) cos(2 pi kz z)
template <typename T>
constexpr void standing_wave(const T A, const T kx, const T ky, const T kz,
                             const T t, const T x, const T y, const T z, T &u,
                             T &ft, T &fx, T &fy, T &fz) {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  ft = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
       cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  fx = A * (-2 * pi * kx) * cos(2 * pi * omega * t) * sin(2 * pi * kx * x) *
       cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  fy = A * (-2 * pi * ky) * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
       sin(2 * pi * ky * y) * cos(2 * pi * kz * z);
  fz = A * (-2 * pi * kz) * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
       cos(2 * pi * ky * y) * sin(2 * pi * kz * z);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t, const T x, const T y,
                        const T z, T &u, T &ft, T &fx, T &fy, T &fz) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    ft = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
    fx = 0;
    fy = 0;
    fz = 0;
  } else {
    u = (f(t - r) - f(t + r)) / r;
    ft = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
    const T fr = -(f(t - r) - f(t + r)) / pow(r, 2) +
                 (f(t - r) * (t - r) + f(t + r) * (t + r)) / (pow(W, 2) * r);
    fx = fr * x / r;
    fy = fr * y / r;
    fz = fr * z / r;
  }
}

extern "C" void FluxWaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, p.x, p.y, p.z, u(p.I),
                        ft(p.I), fx(p.I), fy(p.I), fz(p.I));
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u(p.I),
                   ft(p.I), fx(p.I), fy(p.I), fz(p.I));
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

enum class bc_t { CarpetX, reflecting, radiative };

bc_t get_bc() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(boundary_condition, "CarpetX"))
    return bc_t::CarpetX;
  if (CCTK_EQUALS(boundary_condition, "reflecting"))
    return bc_t::reflecting;
  if (CCTK_EQUALS(boundary_condition, "radiative"))
    return bc_t::radiative;
  assert(0);
}

extern "C" void FluxWaveToyX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_Fluxes;

  const bc_t bc = get_bc();

  // "Reconstructing" at the cell interface is just averaging here, and flux
  // limiting is not necessary since the solution is smooth

  // Calculate x-flux
  grid.loop_int_device<0, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_flux_x(p.I) = 0;
        ft_flux_x(p.I) = average(fx, p, 0);
        fx_flux_x(p.I) = average(ft, p, 0);
        fy_flux_x(p.I) = 0;
        fz_flux_x(p.I) = 0;

        if (bc != bc_t::CarpetX && p.BI[0] != 0) {
          auto chi_m = ft_flux_x(p.I) + p.BI[0] * fx_flux_x(p.I);
          auto chi_p = ft_flux_x(p.I) - p.BI[0] * fx_flux_x(p.I);

          switch (bc) {
          case bc_t::reflecting:
            chi_m = -chi_p;
            break;
          case bc_t::radiative:
            chi_m = 0;
            break;
          default:
            assert(0);
          }

          ft_flux_x(p.I) = (chi_m + chi_p) / 2;
          fx_flux_x(p.I) = p.BI[0] * (chi_m - chi_p) / 2;
        }
      });

  // Calculate y-flux
  grid.loop_int_device<1, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_flux_y(p.I) = 0;
        ft_flux_y(p.I) = average(fy, p, 1);
        fx_flux_y(p.I) = 0;
        fy_flux_y(p.I) = average(ft, p, 1);
        fz_flux_y(p.I) = 0;

        if (bc != bc_t::CarpetX && p.BI[1] != 0) {
          auto chi_m = ft_flux_y(p.I) + p.BI[1] * fy_flux_y(p.I);
          auto chi_p = ft_flux_y(p.I) - p.BI[1] * fy_flux_y(p.I);

          switch (bc) {
          case bc_t::reflecting:
            chi_m = -chi_p;
            break;
          case bc_t::radiative:
            chi_m = 0;
            break;
          default:
            assert(0);
          }

          ft_flux_y(p.I) = (chi_m + chi_p) / 2;
          fy_flux_y(p.I) = p.BI[1] * (chi_m - chi_p) / 2;
        }
      });

  // Calculate z-flux
  grid.loop_int_device<1, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_flux_z(p.I) = 0;
        ft_flux_z(p.I) = average(fz, p, 2);
        fx_flux_z(p.I) = 0;
        fy_flux_z(p.I) = 0;
        fz_flux_z(p.I) = average(ft, p, 2);

        if (bc != bc_t::CarpetX && p.BI[2] != 0) {
          auto chi_m = ft_flux_z(p.I) + p.BI[2] * fz_flux_z(p.I);
          auto chi_p = ft_flux_z(p.I) - p.BI[2] * fz_flux_z(p.I);

          switch (bc) {
          case bc_t::reflecting:
            chi_m = -chi_p;
            break;
          case bc_t::radiative:
            chi_m = 0;
            break;
          default:
            assert(0);
          }

          ft_flux_z(p.I) = (chi_m + chi_p) / 2;
          fz_flux_z(p.I) = p.BI[2] * (chi_m - chi_p) / 2;
        }
      });
}

extern "C" void FluxWaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_rhs(p.I) = flux_div(u_flux_x, u_flux_y, u_flux_z, p) + ft(p.I);
        ft_rhs(p.I) = flux_div(ft_flux_x, ft_flux_y, ft_flux_z, p);
        fx_rhs(p.I) = flux_div(fx_flux_x, fx_flux_y, fx_flux_z, p);
        fy_rhs(p.I) = flux_div(fy_flux_x, fy_flux_y, fy_flux_z, p);
        fz_rhs(p.I) = flux_div(fz_flux_x, fz_flux_y, fz_flux_z, p);
      });
}

extern "C" void FluxWaveToyX_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_Constraints;

  const std::array<Loop::GF3D2<const CCTK_REAL>, dim> f = {fx, fy, fz};
  const std::array<Loop::GF3D2<CCTK_REAL>, dim> curlf = {curlfx, curlfy,
                                                         curlfz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int d = 0; d < dim; ++d) {
          const int a = (d + 1) % dim;
          const int b = (d + 2) % dim;
          curlf[d](p.I) = deriv(f[a], p, b) - deriv(f[b], p, a);
        }
      });
}

extern "C" void FluxWaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_Energy;

  using std::pow;

  const std::array<Loop::GF3D2<const CCTK_REAL>, dim> f = {fx, fy, fz};

  grid.loop_int_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      CCTK_REAL f2 = 0;
                                      for (int d = 0; d < dim; ++d)
                                        f2 += pow(f[d](p.I), 2);

                                      eps(p.I) = (pow(fx(p.I), 2) + f2) / 2;
                                    });
}

extern "C" void FluxWaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_FluxWaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0, fx0, fy0, fz0;
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, p.x, p.y, p.z, u0, ft0,
                        fx0, fy0, fz0);
          u_error(p.I) = u(p.I) - u0;
          ft_error(p.I) = ft(p.I) - ft0;
          fx_error(p.I) = fx(p.I) - fx0;
          fy_error(p.I) = fy(p.I) - fy0;
          fz_error(p.I) = fz(p.I) - fz0;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0, fx0, fy0, fz0;
          gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u0, ft0,
                   fx0, fy0, fz0);
          u_error(p.I) = u(p.I) - u0;
          ft_error(p.I) = ft(p.I) - ft0;
          fx_error(p.I) = fx(p.I) - fx0;
          fy_error(p.I) = fy(p.I) - fy0;
          fz_error(p.I) = fz(p.I) - fz0;
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace FluxWaveToyX
