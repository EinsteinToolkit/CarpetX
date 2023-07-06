#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

namespace TmunuBaseX {
using namespace std;
using namespace Loop;

extern "C" void TmunuBaseX_ZeroTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TmunuBaseX_ZeroTmunu;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> eTtt_(layout, eTtt);
  const GF3D2<CCTK_REAL> eTtx_(layout, eTtx);
  const GF3D2<CCTK_REAL> eTty_(layout, eTty);
  const GF3D2<CCTK_REAL> eTtz_(layout, eTtz);
  const GF3D2<CCTK_REAL> eTxx_(layout, eTxx);
  const GF3D2<CCTK_REAL> eTxy_(layout, eTxy);
  const GF3D2<CCTK_REAL> eTxz_(layout, eTxz);
  const GF3D2<CCTK_REAL> eTyy_(layout, eTyy);
  const GF3D2<CCTK_REAL> eTyz_(layout, eTyz);
  const GF3D2<CCTK_REAL> eTzz_(layout, eTzz);

  const GridDescBaseDevice grid(cctkGH);
#if 0
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTtt_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTtx_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTty_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTtz_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTxx_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTxy_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTxz_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTyy_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTyz_(p.I) = 0; });
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { eTzz_(p.I) = 0; });
#endif

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      eTtt_(p.I) = 0;
                                      eTtx_(p.I) = 0;
                                      eTty_(p.I) = 0;
                                      eTtz_(p.I) = 0;
                                      eTxx_(p.I) = 0;
                                      eTxy_(p.I) = 0;
                                      eTxz_(p.I) = 0;
                                      eTyy_(p.I) = 0;
                                      eTyz_(p.I) = 0;
                                      eTzz_(p.I) = 0;
                                    });
}

} // namespace TmunuBaseX
