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

namespace TestSubcyclingMC {

constexpr int dim = 3;

// u(t,x,y,z) =
//   A cos(2 π omega t) cos(2 π kx x) cos(2 π ky y) cos(2 π kz z)
template <typename T>
constexpr void standing_wave(const T A, const Arith::vect<T, dim> &k, const T t,
                             const Arith::vect<T, dim> &x, T &u) {
  using std::acos, std::cos, std::pow, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(dot(k, k));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * k[0] * x[0]) *
      cos(2 * pi * k[1] * x[1]) * cos(2 * pi * k[2] * x[2]);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T>
constexpr void gaussian(const T A, const T W, const T t,
                        const Arith::vect<T, dim> &x, T &u) {
  using std::exp, std::pow, std::sqrt;

  const T r = sqrt(dot(x, x));

  const auto f = [&](const T v) {
    return A * exp(-pow(v, 2) / (2 * pow(W, 2)));
  };

  if (r < sqrt(std::numeric_limits<T>::epsilon())) {
    // L'Hôpital
    u = 2 / pow(W, 2) * f(t) * t;
  } else {
    u = (f(t - r) - f(t + r)) / r;
  }
}

extern "C" void TestSubcyclingMC_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, p.X, u0);
          u(p.I) = u0;
        });

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto tm = cctk_time - CCTK_DELTA_TIME / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        tm, p.X, um);
          const auto tp = cctk_time + CCTK_DELTA_TIME / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        tp, p.X, up);
          ft(p.I) = (up - um) / CCTK_DELTA_TIME;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fx(p.I) = (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fy(p.I) = (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fz(p.I) = (up - um) / p.DX[2];
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0;
          gaussian(amplitude, gaussian_width, cctk_time, p.X, u0);
          u(p.I) = u0;
        });

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto tm = cctk_time - CCTK_DELTA_TIME / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, tm, p.X, um);
          const auto tp = cctk_time + CCTK_DELTA_TIME / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, tp, p.X, up);
          ft(p.I) = (up - um) / CCTK_DELTA_TIME;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fx(p.I) = (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fy(p.I) = (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fz(p.I) = (up - um) / p.DX[2];
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void TestSubcyclingMC_Evol1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Evol1;

  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u(p.I) += CCTK_DELTA_TIME * ft(p.I);
                                    });

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fx(p.I) += CCTK_DELTA_TIME * (ft(p.I + p.DI[0]) - ft(p.I)) / p.DX[0];
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fy(p.I) += CCTK_DELTA_TIME * (ft(p.I + p.DI[1]) - ft(p.I)) / p.DX[1];
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        fz(p.I) += CCTK_DELTA_TIME * (ft(p.I + p.DI[2]) - ft(p.I)) / p.DX[2];
      });
}

extern "C" void TestSubcyclingMC_Evol2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Evol2;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        ft(p.I) += CCTK_DELTA_TIME * ((fx(p.I) - fx(p.I - p.DI[0])) / p.DX[0] +
                                      (fy(p.I) - fy(p.I - p.DI[1])) / p.DX[1] +
                                      (fz(p.I) - fz(p.I - p.DI[2])) / p.DX[2]);
      });
}

extern "C" void TestSubcyclingMC_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Energy;

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

extern "C" void TestSubcyclingMC_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, p.X, u0);
          u_err(p.I) = u(p.I) - u0;
        });

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto tm = cctk_time - CCTK_DELTA_TIME / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        tm, p.X, um);
          const auto tp = cctk_time + CCTK_DELTA_TIME / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        tp, p.X, up);
          ft_err(p.I) = ft(p.I) - (up - um) / CCTK_DELTA_TIME;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fx_err(p.I) = fx(p.I) - (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fy_err(p.I) = fy(p.I) - (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up;
          standing_wave(amplitude,
                        {standing_wave_kx, standing_wave_ky, standing_wave_kz},
                        cctk_time, Xp, up);
          fz_err(p.I) = fz(p.I) - (up - um) / p.DX[2];
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0;
          gaussian(amplitude, gaussian_width, cctk_time, p.X, u0);
          u_err(p.I) = u(p.I) - u0;
        });

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto tm = cctk_time - CCTK_DELTA_TIME / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, tm, p.X, um);
          const auto tp = cctk_time + CCTK_DELTA_TIME / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, tp, p.X, up);
          ft_err(p.I) = ft(p.I) - (up - um) / CCTK_DELTA_TIME;
        });

    grid.loop_int_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[0] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[0] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fx_err(p.I) = fx(p.I) - (up - um) / p.DX[0];
        });

    grid.loop_int_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[1] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[1] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fy_err(p.I) = fy(p.I) - (up - um) / p.DX[1];
        });

    grid.loop_int_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto Xm = p.X - p.DI[2] * p.DX / 2;
          CCTK_REAL um;
          gaussian(amplitude, gaussian_width, cctk_time, Xm, um);
          const auto Xp = p.X + p.DI[2] * p.DX / 2;
          CCTK_REAL up;
          gaussian(amplitude, gaussian_width, cctk_time, Xp, up);
          fz_err(p.I) = fz(p.I) - (up - um) / p.DX[2];
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace TestSubcyclingMC
