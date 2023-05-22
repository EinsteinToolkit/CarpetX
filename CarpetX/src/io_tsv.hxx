#ifndef CARPETX_CARPETX_IO_TSV_HXX
#define CARPETX_CARPETX_IO_TSV_HXX

#include <cctk.h>

namespace CarpetX {

void OutputTSVold(const cGH *restrict cctkGH);
void OutputTSV(const cGH *restrict cctkGH);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_IO_TSV_HXX
