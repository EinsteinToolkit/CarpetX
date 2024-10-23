#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

namespace TestRKAB {
using namespace Arith;

extern "C" void TestRKAB_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestRKAB_Energy;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto local_Pi{Pi(p.I)};
        const auto local_Dx{Dx(p.I)};
        const auto local_Dy{Dy(p.I)};
        const auto local_Dz{Dz(p.I)};

        energy_density(p.I) = 0.5 * (local_Pi * local_Pi + local_Dx * local_Dx +
                                     local_Dy * local_Dy + local_Dz * local_Dz);
      });
}

} // namespace TestRKAB