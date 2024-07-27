#include <loop_device.hxx>
#include <subcycling.hxx>

#include <sum.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace TestSubcycling {
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

extern "C" void TestSubcycling_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_Initial;
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

  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u_k1(p.I) = 0.0;
                                      u_k2(p.I) = 0.0;
                                      u_k3(p.I) = 0.0;
                                      u_k4(p.I) = 0.0;
                                      rho_k1(p.I) = 0.0;
                                      rho_k2(p.I) = 0.0;
                                      rho_k3(p.I) = 0.0;
                                      rho_k4(p.I) = 0.0;
                                    });
}

/**
 * \brief Calculate RHSs from Us
 *          rhs = RHS(u),
 *
 * \param vlr       RHSs of state vector Us
 * \param vlu       state vector Us
 */
template <int D, typename tVarOut, typename tVarIn>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcRhs(const Loop::GridDescBaseDevice &grid, const array<tVarOut, D> &vlr,
        const array<tVarIn, D> &vlu) {
  tVarOut &u_rhs = vlr[0];
  tVarOut &rho_rhs = vlr[1];
  tVarIn &u = vlu[0];
  tVarIn &rho = vlu[1];

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::pow;
        Arith::vect<CCTK_REAL, dim> ddu;
        for (int d = 0; d < dim; ++d) {
          // ddu[d] = (u(p.I - p.DI[d]) - 2 * u(p.I) + u(p.I + p.DI[d])) /
          //          pow(p.DX[d], 2);
          ddu[d] = (-(u(p.I - 2 * p.DI[d]) + u(p.I + 2 * p.DI[d])) +
                    16 * (u(p.I - p.DI[d]) + u(p.I + p.DI[d])) - 30 * u(p.I)) /
                   (12 * pow(p.DX[d], 2));
        }
        u_rhs(p.I) = rho(p.I);
        rho_rhs(p.I) = ddu[0] + ddu[1] + ddu[2];
      });
}

template <int D, typename tVarOut, typename tVarIn>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcU(const Loop::GridDescBaseDevice &grid, const array<tVarOut, D> &vlu,
      const array<tVarIn, D> &vlp, const array<tVarIn, D> &k1,
      const array<tVarIn, D> &k2, const array<tVarIn, D> &k3,
      const array<tVarIn, D> &k4, const CCTK_REAL dt) {
  for (size_t v = 0; v < D; ++v) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          vlu[v](p.I) =
              vlp[v](p.I) +
              (k1[v](p.I) + (k2[v](p.I) + k3[v](p.I)) * 2 + k4[v](p.I)) * dt /
                  6;
        });
  }
}

template <int D, typename tVarOut, typename tVarIn>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
SetK(const Loop::GridDescBaseDevice &grid, const array<tVarOut, D> &vlk,
     const array<tVarIn, D> &vlr) {
  for (size_t v = 0; v < D; ++v) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { vlk[v](p.I) = vlr[v](p.I); });
  }
}

/**
 * \brief Calculate Ys from U0 and rhs
 *          w = u_p + rhs * dt
 *
 * \param vlw       RK substage Ys
 * \param vlp       state vector U0
 * \param vlr       RK Ks
 * \param dt        time factor of each RK substage
 */
template <int D, typename tVarOut, typename tVarIn>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcYs(const Loop::GridDescBaseDevice &grid, const array<tVarOut, D> &vlw,
       const array<tVarIn, D> &vlp, const array<tVarIn, D> &vlr,
       const CCTK_REAL dt) {
  for (size_t v = 0; v < D; ++v) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          vlw[v](p.I) = vlp[v](p.I) + vlr[v](p.I) * dt;
        });
  }
}

extern "C" void TestSubcycling_CalcYfFromKcs(CCTK_ARGUMENTS, const CCTK_REAL dt,
                                             const CCTK_INT stage) {
  vector<int> u_groups, p_groups;
  array<vector<int>, 4> ks_groups;
  u_groups.push_back(CCTK_GroupIndex("TestSubcycling::ustate"));
  p_groups.push_back(CCTK_GroupIndex("TestSubcycling::pstate"));
  ks_groups[0].push_back(CCTK_GroupIndex("TestSubcycling::k1"));
  ks_groups[1].push_back(CCTK_GroupIndex("TestSubcycling::k2"));
  ks_groups[2].push_back(CCTK_GroupIndex("TestSubcycling::k3"));
  ks_groups[3].push_back(CCTK_GroupIndex("TestSubcycling::k4"));
  const CCTK_REAL xsi = (cctkGH->cctk_iteration % 2) ? 0.0 : 0.5;
  Subcycling::CalcYfFromKcs<4>(CCTK_PASS_CTOC, u_groups, p_groups, ks_groups,
                               dt * 2, xsi, stage);
}

extern "C" void TestSubcycling_SetP(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_SetP;

  CCTK_VINFO("Updating grid function at iteration %d level %d time %g",
             cctk_iteration, cctk_level, cctk_time);
  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      u_p(p.I) = u(p.I);
                                      rho_p(p.I) = rho(p.I);
                                    });
}

extern "C" void TestSubcycling_CalcK1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcK1;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> k1{u_k1, rho_k1};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlu{u, rho};
  constexpr size_t nvars = vlu.size();
  if (use_subcycling_wip)
    TestSubcycling_CalcYfFromKcs(CCTK_PASS_CTOC, dt, 1);
  // k1 = f(Y1)
  CalcRhs<nvars>(grid, k1, vlu);
}

extern "C" void TestSubcycling_CalcY2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcY2;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> vlu{u, rho};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlp{u_p, rho_p};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k1{u_k1, rho_k1};
  constexpr size_t nvars = vlu.size();
  // Y2 = y0 + h/2 k1
  CalcYs<nvars>(grid, vlu, vlp, k1, dt * CCTK_REAL(0.5));
}

extern "C" void TestSubcycling_CalcK2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcK2;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> k2{u_k2, rho_k2};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlu{u, rho};
  constexpr size_t nvars = vlu.size();
  if (use_subcycling_wip)
    TestSubcycling_CalcYfFromKcs(CCTK_PASS_CTOC, dt, 2);
  // k2 = f(Y2)
  CalcRhs<nvars>(grid, k2, vlu);
}

extern "C" void TestSubcycling_CalcY3(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcY3;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> vlu{u, rho};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlp{u_p, rho_p};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k2{u_k2, rho_k2};
  constexpr size_t nvars = vlu.size();
  // Y3 = y0 + h/2 k2
  CalcYs<nvars>(grid, vlu, vlp, k2, dt * CCTK_REAL(0.5));
}

extern "C" void TestSubcycling_CalcK3(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcK3;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> k3{u_k3, rho_k3};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlu{u, rho};
  constexpr size_t nvars = vlu.size();
  if (use_subcycling_wip)
    TestSubcycling_CalcYfFromKcs(CCTK_PASS_CTOC, dt, 3);
  // k3 = f(Y3)
  CalcRhs<nvars>(grid, k3, vlu);
}

extern "C" void TestSubcycling_CalcY4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcY4;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> vlu{u, rho};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlp{u_p, rho_p};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k3{u_k3, rho_k3};
  constexpr size_t nvars = vlu.size();
  // Y4 = y0 + h k3
  CalcYs<nvars>(grid, vlu, vlp, k3, dt);
}

extern "C" void TestSubcycling_CalcK4(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_CalcK4;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> k4{u_k4, rho_k4};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlu{u, rho};
  constexpr size_t nvars = vlu.size();
  if (use_subcycling_wip)
    TestSubcycling_CalcYfFromKcs(CCTK_PASS_CTOC, dt, 4);
  // k4 = f(Y4)
  CalcRhs<nvars>(grid, k4, vlu);
}

extern "C" void TestSubcycling_UpdateU(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_UpdateU;
  DECLARE_CCTK_PARAMETERS;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const array<const Loop::GF3D2<CCTK_REAL>, 2> vlu{u, rho};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> vlp{u_p, rho_p};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k1{u_k1, rho_k1};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k2{u_k2, rho_k2};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k3{u_k3, rho_k3};
  const array<const Loop::GF3D2<const CCTK_REAL>, 2> k4{u_k4, rho_k4};
  constexpr size_t nvars = vlu.size();
  // Update U
  CalcU<nvars>(grid, vlu, vlp, k1, k2, k3, k4, dt);
}

extern "C" void TestSubcycling_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

extern "C" void TestSubcycling_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSubcycling_Error;
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

} // namespace TestSubcycling
