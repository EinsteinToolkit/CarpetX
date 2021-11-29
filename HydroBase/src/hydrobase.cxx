#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace HydroBase {
using namespace Loop;

extern "C" void HydroBase_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroBase_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);

  const GF3D<CCTK_REAL, 1, 1, 1> velx_(cctkGH, velx);
  const GF3D<CCTK_REAL, 1, 1, 1> vely_(cctkGH, vely);
  const GF3D<CCTK_REAL, 1, 1, 1> velz_(cctkGH, velz);

  const GF3D<CCTK_REAL, 1, 1, 1> eps_(cctkGH, eps);

  const GF3D<CCTK_REAL, 1, 1, 1> press_(cctkGH, press);

  const GF3D<CCTK_REAL, 1, 1, 1> Bvecx_(cctkGH, Bvecx);
  const GF3D<CCTK_REAL, 1, 1, 1> Bvecy_(cctkGH, Bvecy);
  const GF3D<CCTK_REAL, 1, 1, 1> Bvecz_(cctkGH, Bvecz);

  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { rho_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { velx_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { vely_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { velz_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { eps_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { press_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { Bvecx_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { Bvecy_(p.I) = 0; });
  loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) { Bvecz_(p.I) = 0; });
}

} // namespace HydroBase
