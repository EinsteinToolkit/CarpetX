#include <loop_device.hxx>

#include <defs.hxx>
#include <simd.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>
#include <type_traits>

namespace SIMDWaveToyX {

constexpr int dim = 3;

// u(t,x,y,z) =
//   A cos(2 pi omega t) sin(2 pi kx x) sin(2 pi ky y) sin(2 pi kz z)
template <typename T, typename VT>
constexpr void standing_wave(const T A, const T kx, const T ky, const T kz,
                             const T t, const VT x, const T y, const T z, VT &u,
                             VT &rho) {
  using Arith::pow2;
  using std::acos, std::cos, std::sin, std::sqrt;

  const T pi = acos(-T(1));
  const T omega = sqrt(pow2(kx) + pow2(ky) + pow2(kz));

  u = A * cos(2 * pi * omega * t) * cos(2 * pi * kx * x) *
      cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
  rho = A * (-2 * pi * omega) * sin(2 * pi * omega * t) * cos(2 * pi * kx * x) *
        cos(2 * pi * ky * y) * cos(2 * pi * kz * z);
}

// u(t,r) = (f(t-r) - f(t+r)) / r
// f(v) = A exp(-1/2 (r/W)^2)
template <typename T, typename VT>
constexpr void gaussian(const T A, const T W, const T t, const VT x, const T y,
                        const T z, VT &u, VT &rho) {
  using Arith::if_else, Arith::pow2, Arith::pown;
  using std::exp, std::sqrt;

  const VT r = sqrt(pow2(x) + pow2(y) + pow2(z));

  const auto f = [&](const VT v) { return A * exp(-pow2(v) / (2 * pow2(W))); };

  // L'Hôpital
  u = if_else(r < sqrt(std::numeric_limits<T>::epsilon()),
              2 / pow2(W) * f(t) * t, (f(t - r) - f(t + r)) / r);
  rho = if_else(r < sqrt(std::numeric_limits<T>::epsilon()),
                -2 / pown(W, 4) * f(t) * (pow2(t) - pow2(W)),
                -(f(t - r) * (t - r) - f(t + r) * (t + r)) / (pow2(W) * r));
}

extern "C" void SIMDWaveToyX_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SIMDWaveToyX_Initial;
  DECLARE_CCTK_PARAMETERS;

  using real = CCTK_REAL;
  using vreal = Arith::simd<CCTK_REAL>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vreal x0 = p.x + Arith::iota<vreal>() * p.dx;
          const real y0 = p.y;
          const real z0 = p.z;
          vreal u0, rho0;
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, x0, y0, z0, u0, rho0);
          u.store(p.mask, p.I, u0);
          rho.store(p.mask, p.I, rho0);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vreal x0 = p.x + Arith::iota<vreal>() * p.dx;
          const real y0 = p.y;
          const real z0 = p.z;
          vreal u0, rho0;
          gaussian(amplitude, gaussian_width, cctk_time, x0, y0, z0, u0, rho0);
          u.store(p.mask, p.I, u0);
          rho.store(p.mask, p.I, rho0);
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

extern "C" void SIMDWaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SIMDWaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  using vreal = Arith::simd<CCTK_REAL>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using Arith::pow2;

        vreal ddu = 0;
        for (int d = 0; d < dim; ++d)
          ddu += (u(p.mask, p.I - p.DI[d]) - 2 * u(p.mask, p.I) +
                  u(p.mask, p.I + p.DI[d])) /
                 pow2(p.DX[d]);

        const vreal udot = rho(p.mask, p.I);
        const vreal rhodot = ddu;

        u_rhs.store(p.mask, p.I, udot);
        rho_rhs.store(p.mask, p.I, rhodot);
      });
}

extern "C" void SIMDWaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SIMDWaveToyX_Energy;

  using vreal = Arith::simd<CCTK_REAL>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using Arith::pow2;

        vreal du2 = 0;
        for (int d = 0; d < dim; ++d)
          du2 += pow2((u(p.mask, p.I + p.DI[d]) - u(p.mask, p.I + p.DI[d])) /
                      (2 * p.DX[d]));

        const vreal eps0 = (pow2(rho(p.mask, p.I)) + du2) / 2;

        eps.store(p.mask, p.I, eps0);
      });
}

extern "C" void SIMDWaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_SIMDWaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  using real = CCTK_REAL;
  using vreal = Arith::simd<CCTK_REAL>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vreal x0 = p.x + Arith::iota<vreal>() * p.dx;
          const real y0 = p.y;
          const real z0 = p.z;
          vreal u0, rho0;
          standing_wave(amplitude, standing_wave_kx, standing_wave_ky,
                        standing_wave_kz, cctk_time, x0, y0, z0, u0, rho0);
          u_err.store(p.mask, p.I, u(p.mask, p.I) - u0);
          rho_err.store(p.mask, p.I, rho(p.mask, p.I) - rho0);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const vreal x0 = p.x + Arith::iota<vreal>() * p.dx;
          const real y0 = p.y;
          const real z0 = p.z;
          vreal u0, rho0;
          gaussian(amplitude, gaussian_width, cctk_time, x0, y0, z0, u0, rho0);
          u_err.store(p.mask, p.I, u(p.mask, p.I) - u0);
          rho_err.store(p.mask, p.I, rho(p.mask, p.I) - rho0);
        });

  } else {
    CCTK_ERROR("Unknown initial condition");
  }
}

} // namespace SIMDWaveToyX
