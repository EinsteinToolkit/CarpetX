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

namespace StaggeredWaveToyX {

constexpr int dim = 3;

// u(t,x,y,z) =
//   A cos(2 π omega t) cos(2 π kx x) cos(2 π ky y) cos(2 π kz z)
template <typename T>
constexpr void standing_wave(const T A, const Arith::vect<T, dim> &k, const T t,
                             const Arith::vect<T, dim> &x, T &u, T &ft) {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(dot(k, k));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * k[0] * x[0]) *
      cos(2 * pi * k[1] * x[1]) * cos(2 * pi * k[2] * x[2]);
  ft = A * (-2 * pi * omega) * sin(2 * pi * omega * t) *
       cos(2 * pi * k[0] * x[0]) * cos(2 * pi * k[1] * x[1]) *
       cos(2 * pi * k[2] * x[2]);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t,
                        const Arith::vect<T, dim> &x, T &u, T &ft) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(dot(x, x));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
    ft = -2 / pow(W, 4) * f(t) * (pow(t, 2) - pow(W, 2));
  } else {
    u = (f(t - r) - f(t + r)) / r;
    ft = -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow(W, 2) * r);
  }
}

extern "C" void StaggeredWaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, p.X, u0, ft0);
          u(p.I) = u0;
          ft(p.I) = ft0;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fx(p.I) = (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fy(p.I) = (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fz(p.I) = (up - um) / p.DX[2];
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0;
          gaussian(amplitude, gaussian_width, cctk_time, p.X, u0, ft0);
          u(p.I) = u0;
          ft(p.I) = ft0;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fx(p.I) = (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fy(p.I) = (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fz(p.I) = (up - um) / p.DX[2];
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void StaggeredWaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_RHS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_rhs(p.I) = ft(p.I);
        ft_rhs(p.I) = (fx(p.I) - fx(p.I - p.DI[0])) / p.DX[0] +
                      (fy(p.I) - fy(p.I - p.DI[1])) / p.DX[1] +
                      (fz(p.I) - fz(p.I - p.DI[2])) / p.DX[2];
      });

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fx_rhs(p.I) = (ft(p.I + p.DI[0]) - ft(p.I)) / p.DX[0];
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fy_rhs(p.I) = (ft(p.I + p.DI[1]) - ft(p.I)) / p.DX[1];
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fz_rhs(p.I) = (ft(p.I + p.DI[2]) - ft(p.I)) / p.DX[2];
      });
}

extern "C" void StaggeredWaveToyX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_Boundaries;

  // Do nothing
}

extern "C" void StaggeredWaveToyX_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_Constraints;

  grid.loop_int_device<0, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        curlfx(p.I) = (fz(p.I + p.DI[1]) - fz(p.I)) / p.DX[1] -
                      (fy(p.I + p.DI[2]) - fy(p.I)) / p.DX[2];
      });

  grid.loop_int_device<1, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        curlfy(p.I) = (fx(p.I + p.DI[2]) - fx(p.I)) / p.DX[2] -
                      (fz(p.I + p.DI[0]) - fz(p.I)) / p.DX[0];
      });

  grid.loop_int_device<1, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        curlfz(p.I) = (fy(p.I + p.DI[0]) - fy(p.I)) / p.DX[0] -
                      (fx(p.I + p.DI[1]) - fx(p.I)) / p.DX[1];
      });
}

extern "C" void StaggeredWaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_Energy;

  using std::pow;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        eps(p.I) =
            (pow(ft(p.I), 2) + pow((fx(p.I) + fx(p.I + p.DI[0])) / 2, 2) +
             pow((fy(p.I) + fy(p.I + p.DI[1])) / 2, 2) +
             pow((fz(p.I) + fz(p.I + p.DI[2])) / 2, 2)) /
            2;
      });
}

extern "C" void StaggeredWaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_StaggeredWaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, p.X, u0, ft0);
          u_err(p.I) = u(p.I) - u0;
          ft_err(p.I) = ft(p.I) - ft0;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fx_err(p.I) = fx(p.I) - (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fy_err(p.I) = fy(p.I) - (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um, ftm;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up, ftp;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up, ftp);
          fz_err(p.I) = fz(p.I) - (up - um) / p.DX[2];
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0, ft0;
          gaussian(amplitude, gaussian_width, cctk_time, p.X, u0, ft0);
          u_err(p.I) = u(p.I) - u0;
          ft_err(p.I) = ft(p.I) - ft0;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fx_err(p.I) = fx(p.I) - (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fy_err(p.I) = fy(p.I) - (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um, ftm;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um, ftm);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up, ftp;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up, ftp);
          fz_err(p.I) = fz(p.I) - (up - um) / p.DX[2];
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace StaggeredWaveToyX
