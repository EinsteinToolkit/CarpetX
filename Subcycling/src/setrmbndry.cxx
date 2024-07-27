#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Subcycling {

extern "C" void Subcycling_SetLevelNeighbor(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Subcycling_SetLevelNeighbor;
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { level_neighbor(p.I) = cctk_level; });
}

extern "C" void Subcycling_SetIsRMBndry(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Subcycling_SetIsRMBndry;
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { isrmbndry(p.I) = 0; });
  grid.loop_ghosts_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        isrmbndry(p.I) = (level_neighbor(p.I) == cctk_level) ? 0 : 1;
      });
}

} // namespace Subcycling
