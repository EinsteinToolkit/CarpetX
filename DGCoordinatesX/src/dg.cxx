#include <dg.hxx>

#include <cctk.h>

namespace DG {

template struct dg_basis<CCTK_REAL, 2>;
template struct dg_basis<CCTK_REAL, 3>;
template struct dg_basis<CCTK_REAL, 4>;
template struct dg_basis<CCTK_REAL, 8>;

} // namespace
