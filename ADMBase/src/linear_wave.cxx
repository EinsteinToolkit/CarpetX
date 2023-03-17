#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace ADMBase {
using namespace Loop;
using namespace std;

extern "C" void ADMBase_linear_wave(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_linear_wave;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

  const GF3D<CCTK_REAL, 0, 0, 0> kxx_(cctkGH, kxx);
  const GF3D<CCTK_REAL, 0, 0, 0> kxy_(cctkGH, kxy);
  const GF3D<CCTK_REAL, 0, 0, 0> kxz_(cctkGH, kxz);
  const GF3D<CCTK_REAL, 0, 0, 0> kyy_(cctkGH, kyy);
  const GF3D<CCTK_REAL, 0, 0, 0> kyz_(cctkGH, kyz);
  const GF3D<CCTK_REAL, 0, 0, 0> kzz_(cctkGH, kzz);

  // See arXiv:1111.2177 [gr-qc], (74-75)

  const auto b = [&](const PointDesc &p) {
    return linear_wave_amplitude *
           sin(2 * CCTK_REAL(M_PI) * (p.x - t) / linear_wave_wavelength);
  };
  const auto bt = [&](const PointDesc &p) {
    return -2 * CCTK_REAL(M_PI) * linear_wave_amplitude /
           linear_wave_wavelength *
           cos(2 * CCTK_REAL(M_PI) * (p.x - t) / linear_wave_wavelength);
  };

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxx_(p.I) = 1; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyy_(p.I) = 1 + b(p); });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gzz_(p.I) = 1 - b(p); });

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxx_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyy_(p.I) = bt(p) / 2; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { kzz_(p.I) = -bt(p) / 2; });
}

} // namespace ADMBase
