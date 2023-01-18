#include "schedule.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace CarpetX {
using namespace std;

extern "C" void CarpetX_InitError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CarpetX_InitError;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_device<1, 1, 1, where_t::everywhere>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { regrid_error(p.I) = 0; });
}

} // namespace CarpetX
