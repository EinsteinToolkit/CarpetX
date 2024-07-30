#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace CarpetXRegrid {
using namespace std;

extern "C" void CarpetXRegrid_InitError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CarpetXRegrid_InitError;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { regrid_error(p.I) = 0; });
}

} // namespace CarpetXRegrid
