#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Coordinates {

extern "C" void Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { vcoordx(p.I) = p.x; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { vcoordy(p.I) = p.y; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { vcoordz(p.I) = p.z; });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { ccoordx(p.I) = p.x; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { ccoordy(p.I) = p.y; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { ccoordz(p.I) = p.z; });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { cvol(p.I) = p.dx * p.dy * p.dz; });
}

} // namespace Coordinates
