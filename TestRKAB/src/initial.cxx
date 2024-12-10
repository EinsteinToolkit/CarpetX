#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "gaussian.hxx"

#include <random>

namespace TestRKAB {
using namespace Arith;

extern "C" void TestRKAB_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_Initial;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};

          const auto A{amplitude};
          const auto kx{standing_wave_kx};
          const auto ky{standing_wave_ky};
          const auto kz{standing_wave_kz};

          phi(p.I) = sw::phi(A, kx, ky, kz, t, p.x, p.y, p.z);
          Pi(p.I) = sw::Pi(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dx(p.I) = sw::Dx(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dy(p.I) = sw::Dy(A, kx, ky, kz, t, p.x, p.y, p.z);
          Dz(p.I) = sw::Dz(A, kx, ky, kz, t, p.x, p.y, p.z);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto t{cctk_time};

          phi(p.I) = gauss::phi(amplitude, gaussian_width, t, p.x, p.y, p.z);
          Pi(p.I) = gauss::Pi(amplitude, gaussian_width, t, p.x, p.y, p.z);
          Dx(p.I) = gauss::Dx(amplitude, gaussian_width, t, p.x, p.y, p.z);
          Dy(p.I) = gauss::Dy(amplitude, gaussian_width, t, p.x, p.y, p.z);
          Dz(p.I) = gauss::Dz(amplitude, gaussian_width, t, p.x, p.y, p.z);
        });

  } else if (CCTK_EQUALS(initial_condition, "noise")) {
    grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [&] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::cos;

          const auto noise_value{noise_boundary * cos(10.0 * p.x / p.dx) *
                                 cos(10.0 * p.y / p.dy) *
                                 cos(10.0 * p.z / p.dz)};
          phi(p.I) = noise_value;
          Pi(p.I) = noise_value;
          Dx(p.I) = noise_value;
          Dy(p.I) = noise_value;
          Dz(p.I) = noise_value;
        });
  } else {
    CCTK_VERROR("Unknown initial condition \"%s\"", initial_condition);
  }

  // Regardless of ID choice, we initialize the previous RHS with zeros
  grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      phi_p_rhs(p.I) = 0.0;
                                      Pi_p_rhs(p.I) = 0.0;
                                      Dx_p_rhs(p.I) = 0.0;
                                      Dy_p_rhs(p.I) = 0.0;
                                      Dz_p_rhs(p.I) = 0.0;

                                      phi_pp_rhs(p.I) = 0.0;
                                      Pi_pp_rhs(p.I) = 0.0;
                                      Dx_pp_rhs(p.I) = 0.0;
                                      Dy_pp_rhs(p.I) = 0.0;
                                      Dz_pp_rhs(p.I) = 0.0;
                                    });
}

} // namespace TestRKAB