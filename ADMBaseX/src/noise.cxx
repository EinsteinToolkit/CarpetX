#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <noisex.hxx>

namespace ADMBaseX {
using namespace Loop;

extern "C" void ADMBaseX_add_noise(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_ADMBaseX_add_noise;

  NoiseX::add(cctkGH, gxx);
  NoiseX::add(cctkGH, gxx);
  NoiseX::add(cctkGH, gxy);
  NoiseX::add(cctkGH, gxz);
  NoiseX::add(cctkGH, gyy);
  NoiseX::add(cctkGH, gyz);
  NoiseX::add(cctkGH, gzz);

  NoiseX::add(cctkGH, kxx);
  NoiseX::add(cctkGH, kxy);
  NoiseX::add(cctkGH, kxz);
  NoiseX::add(cctkGH, kyy);
  NoiseX::add(cctkGH, kyz);
  NoiseX::add(cctkGH, kzz);

  NoiseX::add(cctkGH, alp);

  NoiseX::add(cctkGH, dtalp);

  NoiseX::add(cctkGH, betax);
  NoiseX::add(cctkGH, betay);
  NoiseX::add(cctkGH, betaz);

  NoiseX::add(cctkGH, dtbetax);
  NoiseX::add(cctkGH, dtbetay);
  NoiseX::add(cctkGH, dtbetaz);
}

} // namespace ADMBaseX
