#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop.hxx>
#include <loop_device.hxx>

namespace TestLoopX {
using namespace Loop;

extern "C" void TestLoopX_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoopX_Init;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) = 0.0;
      });
}

extern "C" void TestLoopX_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoopX_Sync;
}

extern "C" void TestLoopX_OutermostInterior(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestLoopX_OutermostInterior;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_outermost_int<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 10.0;
      });

  grid.loop_outermost_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 1.0;
      });
}

} // namespace TestLoopX
