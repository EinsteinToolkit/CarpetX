#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop.hxx>
#include <loop_device.hxx>

namespace TestLoop {
using namespace Loop;

extern "C" void TestLoop_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoop_Init;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) = 0.0;
      });
}

extern "C" void TestLoop_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoop_Sync;
}

extern "C" void TestLoop_OutermostInterior(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoop_OutermostInterior;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_bnd<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 10.0;
      });

  grid.loop_int_bnd_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 1.0;
      });
}

} // namespace TestLoop
