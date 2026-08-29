#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace DGWaveToyX {

extern "C" void DGWaveToyX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGWaveToyX_Boundaries;
  DECLARE_CCTK_PARAMETERS;
}

} // namespace DGWaveToyX
