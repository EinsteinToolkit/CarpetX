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

extern "C" void WaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_WaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(boundary_condition, "Carpetx")) {

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

  } else if (CCTK_EQUALS(boundary_condition, "radiative")) {

    const Arith::vect<bool, dim> bboxlo = {
        (bool)cctk_bbox[0],
        (bool)cctk_bbox[2],
        (bool)cctk_bbox[4],
    };
    const Arith::vect<bool, dim> bboxhi = {
        (bool)cctk_bbox[1],
        (bool)cctk_bbox[3],
        (bool)cctk_bbox[5],
    };
    const Arith::vect<int, dim> domain_imin = {
        cctk_nghostzones[0],
        cctk_nghostzones[1],
        cctk_nghostzones[2],
    };
    const Arith::vect<int, dim> domain_imax = {
        cctk_lsh[0] - cctk_nghostzones[0],
        cctk_lsh[1] - cctk_nghostzones[1],
        cctk_lsh[2] - cctk_nghostzones[2],
    };

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const Arith::vect<bool, dim> isbndlo = bboxlo && p.I == domain_imin;
          const Arith::vect<bool, dim> isbndhi =
              bboxhi && p.I == domain_imax - 1;

          if (any(isbndlo || isbndhi)) {
            // boundary

            // Split potential into its two characteristics:
            // (We use the notation for the upper x boundary here.)
            //     u(t,x) = f(t-x) + g(t+x)
            // "Radiative" means the incoming characteristic is zero:
            //     g(v) = 0
            // This is equivalent to:
            //     d_t u(t,x) + d_x u(t,x) = 0
            // The boundary condition then follows as:
            //     d_t u(t,x) = - d_x u(t,x)

            // Take another time derivative for a first-order-in-time system:
            //     rho(t,x) = - d_x rho(t,x)

            // In general, with the outward boundary normal n^i:
            //     d_t u = - n^i d_i u

            const Arith::vect<int, dim> n = isbndhi - isbndlo; // outward normal

            Arith::vect<CCTK_REAL, dim> drho;
            for (int d = 0; d < dim; ++d)
              if (isbndlo[d])
                drho[d] = (rho(p.I + p.DI[d]) - rho(p.I)) / p.DX[d];
              else if (isbndhi[d])
                drho[d] = (rho(p.I) - rho(p.I - p.DI[d])) / p.DX[d];
              else
                drho[d] = 0;

            udot(p.I) = rho(p.I);
            rhodot(p.I) = -dot(n, drho);

          } else {
            // interior

            using std::pow;
            CCTK_REAL ddu = 0;
            for (int d = 0; d < dim; ++d)
              ddu += (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
                     pow(p.DX[d], 2);

            udot(p.I) = rho(p.I);
            rhodot(p.I) = ddu;
          }
        });

  } else {
    CCTK_ERROR("Internal error");
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
          uerr(p.I) = u(p.I) - u0;
          rhoerr(p.I) = rho(p.I) - rho0;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, rho0;
          gaussian(amplitude, gaussian_width, cctk_time, p.x, p.y, p.z, u0,
                   rho0);
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
