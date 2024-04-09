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
using namespace Arith;

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
          rho(p.I) = 0;
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL u0;
          gaussian(amplitude, gaussian_width, cctk_time, p.X, u0);
          u(p.I) = u0;
          rho(p.I) = 0;
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

void calcRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY2;
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::pow;
        Arith::vect<CCTK_REAL, dim> ddu;
        for (int d = 0; d < dim; ++d) {
          ddu[d] = (u_w(p.I - p.DI[d]) - 2 * u_w(p.I) + u_w(p.I + p.DI[d])) /
                   pow(p.DX[d], 2);
        }
        u_rhs(p.I) = rho_w(p.I);
        rho_rhs(p.I) = ddu[0] + ddu[1] + ddu[2];
      });
}

void updateU(CCTK_ARGUMENTS, CCTK_REAL dt) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY2;
  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u(p.I) += u_rhs(p.I) * dt;
                                      rho(p.I) += rho_rhs(p.I) * dt;
                                    });
}

void calcYs(CCTK_ARGUMENTS, CCTK_REAL dt) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY2;
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_w(p.I) = u_p(p.I) + u_rhs(p.I) * dt;
        rho_w(p.I) = rho_p(p.I) + rho_rhs(p.I) * dt;
      });
}

extern "C" void TestSubcyclingMC_CalcY1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY1;
  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u_p(p.I) = u(p.I);
                                      rho_p(p.I) = rho(p.I);
                                      u_w(p.I) = u(p.I);
                                      rho_w(p.I) = rho(p.I);
                                    });
}

extern "C" void TestSubcyclingMC_CalcY2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY2;
  CCTK_REAL dt = CCTK_DELTA_TIME;
  calcRHS(CCTK_PASS_CTOC); // k1
  updateU(CCTK_PASS_CTOC, dt / CCTK_REAL(6.));
  calcYs(CCTK_PASS_CTOC, dt * CCTK_REAL(0.5)); // Y2
}

extern "C" void TestSubcyclingMC_CalcY3(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY3;
  CCTK_REAL dt = CCTK_DELTA_TIME;
  calcRHS(CCTK_PASS_CTOC); // k2
  updateU(CCTK_PASS_CTOC, dt / CCTK_REAL(3.));
  calcYs(CCTK_PASS_CTOC, dt * CCTK_REAL(0.5)); // Y3
}

extern "C" void TestSubcyclingMC_CalcY4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY4;
  CCTK_REAL dt = CCTK_DELTA_TIME;
  calcRHS(CCTK_PASS_CTOC); // k3
  updateU(CCTK_PASS_CTOC, dt / CCTK_REAL(3.));
  calcYs(CCTK_PASS_CTOC, dt); // Y4
}

extern "C" void TestSubcyclingMC_UpdateU(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_UpdateU;
  CCTK_REAL dt = CCTK_DELTA_TIME;
  calcRHS(CCTK_PASS_CTOC); // k4
  updateU(CCTK_PASS_CTOC, dt / CCTK_REAL(6.));
}

} // namespace TestSubcyclingMC
