#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace HydroBaseX {
using namespace Loop;

extern "C" void HydroBaseX_Zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_Zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      rho(p.I) = 0.0;
                                      velx(p.I) = 0.0;
                                      vely(p.I) = 0.0;
                                      velz(p.I) = 0.0;
                                      eps(p.I) = 0.0;
                                      press(p.I) = 0.0;
                                    });
}

extern "C" void HydroBaseX_temperature_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_temperature_zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              temperature(p.I) = 0.0;
                            });
}

extern "C" void HydroBaseX_entropy_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_entropy_zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              entropy(p.I) = 0.0;
                            });
}

extern "C" void HydroBaseX_Y_e_one(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_Y_e_one;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Ye(p.I) = 1.0; });
}

extern "C" void HydroBaseX_Bvec_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_Bvec_zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      Bvecx(p.I) = 0.0;
                                      Bvecy(p.I) = 0.0;
                                      Bvecz(p.I) = 0.0;
                                    });
}

extern "C" void HydroBaseX_Avec_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_Avec_zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecx(p.I) = 0.0; });
  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecy(p.I) = 0.0; });
  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avecz(p.I) = 0.0; });
}

extern "C" void HydroBaseX_Aphi_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_Aphi_zero;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Aphi(p.I) = 0.0; });
}

} // namespace HydroBaseX
