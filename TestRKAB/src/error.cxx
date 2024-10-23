#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include "standing_wave.hxx"
#include "gaussian.hxx"

namespace TestRKAB {
using namespace Arith;

extern "C" void TestRKAB_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_Error;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_condition, "standing wave")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};

          const auto A{amplitude};
          const auto kx{standing_wave_kx};
          const auto ky{standing_wave_ky};
          const auto kz{standing_wave_kz};

          const auto expected_phi{sw::phi(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Pi{sw::Pi(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dx{sw::Dx(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dy{sw::Dy(A, kx, ky, kz, t, p.x, p.y, p.z)};
          const auto expected_Dz{sw::Dz(A, kx, ky, kz, t, p.x, p.y, p.z)};

          const auto actual_phi{phi(p.I)};
          const auto actual_Pi{Pi(p.I)};
          const auto actual_Dx{Dx(p.I)};
          const auto actual_Dy{Dy(p.I)};
          const auto actual_Dz{Dz(p.I)};

          phi_err(p.I) = fabs(expected_phi - actual_phi);
          Pi_err(p.I) = fabs(expected_Pi - actual_Pi);
          Dx_err(p.I) = fabs(expected_Dx - actual_Dx);
          Dy_err(p.I) = fabs(expected_Dy - actual_Dy);
          Dz_err(p.I) = fabs(expected_Dz - actual_Dz);
        });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          using std::fabs;

          const auto t{cctk_time};

          const auto expected_phi{
              gauss::phi(amplitude, gaussian_width, t, p.x, p.y, p.z)};
          const auto expected_Pi{
              gauss::Pi(amplitude, gaussian_width, t, p.x, p.y, p.z)};
          const auto expected_Dx{
              gauss::Dx(amplitude, gaussian_width, t, p.x, p.y, p.z)};
          const auto expected_Dy{
              gauss::Dy(amplitude, gaussian_width, t, p.x, p.y, p.z)};
          const auto expected_Dz{
              gauss::Dz(amplitude, gaussian_width, t, p.x, p.y, p.z)};

          const auto actual_phi{phi(p.I)};
          const auto actual_Pi{Pi(p.I)};
          const auto actual_Dx{Dx(p.I)};
          const auto actual_Dy{Dy(p.I)};
          const auto actual_Dz{Dz(p.I)};

          phi_err(p.I) = fabs(expected_phi - actual_phi);
          Pi_err(p.I) = fabs(expected_Pi - actual_Pi);
          Dx_err(p.I) = fabs(expected_Dx - actual_Dx);
          Dy_err(p.I) = fabs(expected_Dy - actual_Dy);
          Dz_err(p.I) = fabs(expected_Dz - actual_Dz);
        });
  } else {
    CCTK_VERROR("Unknown initial condition \"%s\"", initial_condition);
  }
}

} // namespace TestRKAB