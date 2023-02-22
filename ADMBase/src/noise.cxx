#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <random>

namespace ADMBase {
using namespace Loop;
using namespace std;

extern "C" void ADMBase_add_noise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_add_noise;
  DECLARE_CCTK_PARAMETERS;

  // Hardware random device
  random_device device;
  // Create and seed software random number engine from hardware random number
  default_random_engine engine(device());
  // Random number distribution
  uniform_real_distribution<CCTK_REAL> distribution(-noise_amplitude,
                                                    noise_amplitude);
  const auto add_noise = [&](CCTK_REAL &restrict var) {
    var += distribution(engine);
  };

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

  const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  const GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

  const GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
  const GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
  const GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

  const GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gxx_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gxy_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gxz_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gyy_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gyz_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(gzz_(p.I)); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kxx_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kxy_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kxz_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kyy_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kyz_(p.I)); });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(kzz_(p.I)); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { add_noise(alp_(p.I)); });

  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(dtalp_(p.I)); });

  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(betax_(p.I)); });
  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(betay_(p.I)); });
  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(betaz_(p.I)); });

  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(dtbetax_(p.I)); });
  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(dtbetay_(p.I)); });
  loop_int<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { add_noise(dtbetaz_(p.I)); });
}

} // namespace ADMBase
