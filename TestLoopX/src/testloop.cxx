#include <cctk.h>
#include "CarpetX/CarpetX/src/driver.hxx"
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

  const auto symmetries = CarpetX::ghext->patchdata.at(cctk_patch).symmetries;
  const vect<vect<bool, Loop::dim>, 2> is_sym_bnd {
    {
      symmetries[0][0] != CarpetX::symmetry_t::none,
      symmetries[0][1] != CarpetX::symmetry_t::none,
      symmetries[0][2] != CarpetX::symmetry_t::none
    },
    {
      symmetries[1][0] != CarpetX::symmetry_t::none,
      symmetries[1][1] != CarpetX::symmetry_t::none,
      symmetries[1][2] != CarpetX::symmetry_t::none
    }
  };
  grid.loop_outermost_int<0, 0, 0>(
      grid.nghostzones, is_sym_bnd,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 10.0;
      });

  grid.loop_outermost_int_device<0, 0, 0>(
      grid.nghostzones, is_sym_bnd,
      [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
          testloop_gf(p.I) += 1.0;
      });
}

} // namespace TestLoopX
