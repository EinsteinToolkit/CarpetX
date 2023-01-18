#ifndef IO_TSV_HXX
#define IO_TSV_HXX

#include <cctk.h>

namespace CarpetX {

void OutputTSVold(const cGH *restrict cctkGH);
void OutputTSV(const cGH *restrict cctkGH);

} // namespace CarpetX

#endif // #ifndef IO_TSV_HXX
