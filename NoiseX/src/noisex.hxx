#include <cctk.h>

#include <loop_device.hxx>

namespace NoiseX {

void add(const cGH *restrict const cctkGH, Loop::GF3D2<CCTK_REAL> var);

} // namespace NoiseX