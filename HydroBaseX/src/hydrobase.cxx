#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace HydroBaseX {
using namespace Loop;

extern "C" void HydroBaseX_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroBaseX_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const bool vacuum = CCTK_EQUALS(initial_hydro, "vacuum");
  const bool temperature_zero = CCTK_EQUALS(initial_temperature, "zero");
  const bool entropy_zero = CCTK_EQUALS(initial_entropy, "zero");
  const bool Ye_one = CCTK_EQUALS(initial_Y_e, "one");
  const bool Bvec_zero = CCTK_EQUALS(initial_Bvec, "zero");

  if (vacuum || temperature_zero || entropy_zero || Ye_one || Bvec_zero)
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (vacuum) {
            rho(p.I) = 0;
            velx(p.I) = 0;
            vely(p.I) = 0;
            velz(p.I) = 0;
            eps(p.I) = 0;
            press(p.I) = 0;
          }
          if (temperature_zero)
            temperature(p.I) = 0;
          if (entropy_zero)
            entropy(p.I) = 0;
          if (Ye_one)
            Ye(p.I) = 1;
          if (Bvec_zero) {
            Bvecx(p.I) = 0;
            Bvecy(p.I) = 0;
            Bvecz(p.I) = 0;
          }
        });

  if (CCTK_EQUALS(initial_Avec, "zero")) {
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

  if (CCTK_EQUALS(initial_Aphi, "zero"))
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                              CCTK_ATTRIBUTE_ALWAYS_INLINE { Aphi(p.I) = 0; });
}

} // namespace HydroBaseX
