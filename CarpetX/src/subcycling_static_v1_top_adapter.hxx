#ifndef CARPETX_SUBCYCLING_STATIC_V1_TOP_ADAPTER_HXX
#define CARPETX_SUBCYCLING_STATIC_V1_TOP_ADAPTER_HXX

#include "subcycling_static_v1_contract.hxx"

#include <cctk.h>

namespace CarpetX {

int EvolveStaticV1(tFleshConfig *config);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_STATIC_V1_TOP_ADAPTER_HXX
