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
    std::uniform_real_distribution<CCTK_REAL> noise_distrib{-noise_boundary,
                                                            noise_boundary};
    std::mt19937 noise_engine{noise_seed};

    grid.loop_int_device<0, 0, 0>(grid.nghostzones,
                                  [&] CCTK_DEVICE(const Loop::PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        const auto t{cctk_time};

                                        phi(p.I) = noise_distrib(noise_engine);
                                        Pi(p.I) = noise_distrib(noise_engine);
                                        Dx(p.I) = noise_distrib(noise_engine);
                                        Dy(p.I) = noise_distrib(noise_engine);
                                        Dz(p.I) = noise_distrib(noise_engine);
                                      });
  } else {
    CCTK_VERROR("Unknown initial condition \"%s\"", initial_condition);
  }
}

} // namespace TestRKAB