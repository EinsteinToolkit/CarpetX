#ifndef CARPETX_CARPETX_SYNC_RESTRICT_INTERNAL_HXX
#define CARPETX_CARPETX_SYNC_RESTRICT_INTERNAL_HXX

#include <cctk.h>

namespace CarpetX {

// Intra-CarpetX symbols defined in sync_restrict.cxx and called from
// schedule.cxx. Not part of the public schedule.hxx surface.

void Reflux(const cGH *cctkGH, int level);
void Restrict(const cGH *cctkGH, int level);
void SyncRestrictedGFs(const cGH *cctkGH);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_SYNC_RESTRICT_INTERNAL_HXX
