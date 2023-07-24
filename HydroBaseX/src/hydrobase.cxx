#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace HydroBaseX {
using namespace Loop;

extern "C" void HydroBaseX_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_initial_data;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      rho(p.I) = 0;
                                      velx(p.I) = 0;
                                      vely(p.I) = 0;
                                      velz(p.I) = 0;
                                      eps(p.I) = 0;
                                      press(p.I) = 0;
                                      temperature(p.I) = 0;
                                      entropy(p.I) = 0;
                                      Ye(p.I) = 0;
                                      Bvecx(p.I) = 0;
                                      Bvecy(p.I) = 0;
                                      Bvecz(p.I) = 0;
                                    });

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecx(p.I) = 0; });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecy(p.I) = 0; });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecz(p.I) = 0; });
}

} // namespace HydroBaseX
