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

extern "C" void TestSubcyclingMC_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_Initial;
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

void CalcRhsAndUpdateU(const Loop::GridDescBaseDevice &grid,
                       // output 1
                       const Loop::GF3D2<CCTK_REAL> &u_rhs,
                       const Loop::GF3D2<CCTK_REAL> &rho_rhs,
                       // input 1
                       const Loop::GF3D2<CCTK_REAL> &u_w,
                       const Loop::GF3D2<CCTK_REAL> &rho_w,
                       // output 2
                       const Loop::GF3D2<CCTK_REAL> &u,
                       const Loop::GF3D2<CCTK_REAL> &rho,
                       // input
                       CCTK_REAL dt) {
  // calculate rhs
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
  // update ustate
  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u(p.I) += u_rhs(p.I) * dt;
                                      rho(p.I) += rho_rhs(p.I) * dt;
                                    });
}

void CalcYs(const Loop::GridDescBaseDevice &grid,
            // output
            const Loop::GF3D2<CCTK_REAL> &u_w,
            const Loop::GF3D2<CCTK_REAL> &rho_w,
            // input
            const Loop::GF3D2<const CCTK_REAL> &u_p,
            const Loop::GF3D2<const CCTK_REAL> &rho_p,
            const Loop::GF3D2<CCTK_REAL> &u_rhs,
            const Loop::GF3D2<CCTK_REAL> &rho_rhs, CCTK_REAL dt) {
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        u_w(p.I) = u_p(p.I) + u_rhs(p.I) * dt;
        rho_w(p.I) = rho_p(p.I) + rho_rhs(p.I) * dt;
      });
}

/**
 * \brief Calculate Ys ghost points for fine grid using Ks on coarse grid
 *
 * \param Yf        RK substage Ys on the fine side to be interperated into
 *                  the ghost zones
 * \param kcs       RK ks on the coarset side
 * \param u0        u at t0
 * \param dtc       Time step size on coarse side
 * \param xsi       which substep on fine level during a coarse time
 *                  step.  For an AMR simulation with subcycling and a
 *                  refinement ratio of 2, the number is either 0 or 0.5,
 *                  denoting the first and second substep, respectively.
 * \param stage     RK stage number starting from 1
 */
void CalcYfFromKcs(const Loop::GridDescBaseDevice &grid,
                   // output
                   const Loop::GF3D2<CCTK_REAL> &Yf,
                   // input
                   const array<const Loop::GF3D2<CCTK_REAL>, 4> &kcs,
                   CCTK_REAL u0, CCTK_REAL dtc, CCTK_REAL xsi, CCTK_INT stage) {
  assert(stage > 0 && stage <= 4);

  CCTK_REAL r = 0.5; // ratio between coarse and fine cell size (2 to 1 MR case)
  CCTK_REAL xsi2 = xsi * xsi;
  CCTK_REAL xsi3 = xsi2 * xsi;
  // coefficients for U
  CCTK_REAL b1 = xsi - CCTK_REAL(1.5) * xsi2 + CCTK_REAL(2. / 3.) * xsi3;
  CCTK_REAL b2 = xsi2 - CCTK_REAL(2. / 3.) * xsi3;
  CCTK_REAL b3 = b2;
  CCTK_REAL b4 = CCTK_REAL(-0.5) * xsi2 + CCTK_REAL(2. / 3.) * xsi3;
  // coefficients for Ut
  CCTK_REAL c1 = CCTK_REAL(1.) - CCTK_REAL(3.) * xsi + CCTK_REAL(2.) * xsi2;
  CCTK_REAL c2 = CCTK_REAL(2.) * xsi - CCTK_REAL(2.) * xsi2;
  CCTK_REAL c3 = c2;
  CCTK_REAL c4 = -xsi + CCTK_REAL(2.) * xsi2;
  // coefficients for Utt
  CCTK_REAL d1 = CCTK_REAL(-3.) + CCTK_REAL(4.) * xsi;
  CCTK_REAL d2 = CCTK_REAL(2.) - CCTK_REAL(4.) * xsi;
  CCTK_REAL d3 = d2;
  CCTK_REAL d4 = CCTK_REAL(-1.) + CCTK_REAL(4.) * xsi;
  // coefficients for Uttt
  constexpr CCTK_REAL e1 = CCTK_REAL(4.);
  constexpr CCTK_REAL e2 = CCTK_REAL(-4.);
  constexpr CCTK_REAL e3 = CCTK_REAL(-4.);
  constexpr CCTK_REAL e4 = CCTK_REAL(4.);

  if (stage == 1) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL k1 = kcs[0](p.I);
          CCTK_REAL k2 = kcs[1](p.I);
          CCTK_REAL k3 = kcs[2](p.I);
          CCTK_REAL k4 = kcs[3](p.I);
          CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
          Yf(p.I) = u0 + dtc * uu;
        });
  } else if (stage == 2) {

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL k1 = kcs[0](p.I);
          CCTK_REAL k2 = kcs[1](p.I);
          CCTK_REAL k3 = kcs[2](p.I);
          CCTK_REAL k4 = kcs[3](p.I);
          CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
          CCTK_REAL ut = c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4;
          Yf(p.I) = u0 + dtc * (uu + CCTK_REAL(0.5) * r * ut);
        });
  } else if (stage == 3 || stage == 4) {
    CCTK_REAL r2 = r * r;
    CCTK_REAL r3 = r2 * r;
    CCTK_REAL at = (stage == 3) ? CCTK_REAL(0.5) * r : r;
    CCTK_REAL att = (stage == 3) ? CCTK_REAL(0.25) * r2 : CCTK_REAL(0.5) * r2;
    CCTK_REAL attt =
        (stage == 3) ? CCTK_REAL(0.0625) * r3 : CCTK_REAL(0.125) * r3;
    CCTK_REAL ak = (stage == 3) ? CCTK_REAL(-4.) : CCTK_REAL(4.);

    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL k1 = kcs[0](p.I);
          CCTK_REAL k2 = kcs[1](p.I);
          CCTK_REAL k3 = kcs[2](p.I);
          CCTK_REAL k4 = kcs[3](p.I);
          CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
          CCTK_REAL ut = c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4;
          CCTK_REAL utt = d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4;
          CCTK_REAL uttt = e1 * k1 + e2 * k2 + e3 * k3 + e4 * k4;
          Yf(p.I) = u0 + dtc * (uu + at * ut + att * utt +
                                attt * (uttt + ak * (k3 - k2)));
        });
  }
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
  CalcRhsAndUpdateU(grid, u_k1, rho_k1, u_w, rho_w, u, rho,
                    CCTK_DELTA_TIME / CCTK_REAL(6.)); // k1
  CalcYs(grid, u_w, rho_w, u_p, rho_p, u_k1, rho_k1,
         CCTK_DELTA_TIME * CCTK_REAL(0.5)); // Y2
}

extern "C" void TestSubcyclingMC_CalcY3(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY3;
  CalcRhsAndUpdateU(grid, u_k2, rho_k2, u_w, rho_w, u, rho,
                    CCTK_DELTA_TIME / CCTK_REAL(3.)); // k2
  CalcYs(grid, u_w, rho_w, u_p, rho_p, u_k2, rho_k2,
         CCTK_DELTA_TIME * CCTK_REAL(0.5)); // Y3
}

extern "C" void TestSubcyclingMC_CalcY4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_CalcY4;
  CalcRhsAndUpdateU(grid, u_k3, rho_k3, u_w, rho_w, u, rho,
                    CCTK_DELTA_TIME / CCTK_REAL(3.)); // k3
  CalcYs(grid, u_w, rho_w, u_p, rho_p, u_k3, rho_k3,
         CCTK_DELTA_TIME); // Y4
}

extern "C" void TestSubcyclingMC_UpdateU(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingMC_UpdateU;
  CalcRhsAndUpdateU(grid, u_k4, rho_k4, u_w, rho_w, u, rho,
                    CCTK_DELTA_TIME / CCTK_REAL(6.)); // k4
}

} // namespace TestSubcyclingMC
