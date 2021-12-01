#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <array>

namespace HydroBase {
using namespace Loop;

extern "C" void HydroBase_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroBase_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> rho_(layout, rho);

  const GF3D2<CCTK_REAL> velx_(layout, velx);
  const GF3D2<CCTK_REAL> vely_(layout, vely);
  const GF3D2<CCTK_REAL> velz_(layout, velz);

  const GF3D2<CCTK_REAL> eps_(layout, eps);

  const GF3D2<CCTK_REAL> press_(layout, press);

  const GF3D2<CCTK_REAL> Bvecx_(layout, Bvecx);
  const GF3D2<CCTK_REAL> Bvecy_(layout, Bvecy);
  const GF3D2<CCTK_REAL> Bvecz_(layout, Bvecz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { rho_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { velx_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { vely_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { velz_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eps_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { press_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Bvecx_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Bvecy_(p.I) = 0; });
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Bvecz_(p.I) = 0; });
}

} // namespace HydroBase
