#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

namespace ADMBaseX {
using namespace std;
using namespace Loop;

extern "C" void ADMBaseX_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBaseX_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gxx_(layout, gxx);
  const GF3D2<CCTK_REAL> gxy_(layout, gxy);
  const GF3D2<CCTK_REAL> gxz_(layout, gxz);
  const GF3D2<CCTK_REAL> gyy_(layout, gyy);
  const GF3D2<CCTK_REAL> gyz_(layout, gyz);
  const GF3D2<CCTK_REAL> gzz_(layout, gzz);

  const GF3D2<CCTK_REAL> kxx_(layout, kxx);
  const GF3D2<CCTK_REAL> kxy_(layout, kxy);
  const GF3D2<CCTK_REAL> kxz_(layout, kxz);
  const GF3D2<CCTK_REAL> kyy_(layout, kyy);
  const GF3D2<CCTK_REAL> kyz_(layout, kyz);
  const GF3D2<CCTK_REAL> kzz_(layout, kzz);

  const GridDescBaseDevice grid(cctkGH);

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      gxx_(p.I) = 1;
                                      gxy_(p.I) = 0;
                                      gxz_(p.I) = 0;
                                      gyy_(p.I) = 1;
                                      gyz_(p.I) = 0;
                                      gzz_(p.I) = 1;

                                      kxx_(p.I) = 0;
                                      kxy_(p.I) = 0;
                                      kxz_(p.I) = 0;
                                      kyy_(p.I) = 0;
                                      kyz_(p.I) = 0;
                                      kzz_(p.I) = 0;
                                    });
}

extern "C" void ADMBaseX_initial_lapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBaseX_initial_lapse;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> alp_(layout, alp);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { alp_(p.I) = 1; });
}

extern "C" void ADMBaseX_initial_dtlapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBaseX_initial_dtlapse;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> dtalp_(layout, dtalp);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { dtalp_(p.I) = 0; });
}

extern "C" void ADMBaseX_initial_shift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBaseX_initial_shift;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      betax_(p.I) = 0;
                                      betay_(p.I) = 0;
                                      betaz_(p.I) = 0;
                                    });
}

extern "C" void ADMBaseX_initial_dtshift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBaseX_initial_dtshift;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> dtbetax_(layout, dtbetax);
  const GF3D2<CCTK_REAL> dtbetay_(layout, dtbetay);
  const GF3D2<CCTK_REAL> dtbetaz_(layout, dtbetaz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      dtbetax_(p.I) = 0;
                                      dtbetay_(p.I) = 0;
                                      dtbetaz_(p.I) = 0;
                                    });
}

} // namespace ADMBaseX
