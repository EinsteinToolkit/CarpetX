#ifndef CARPETX_CF_MASK_HXX
#define CARPETX_CF_MASK_HXX

namespace CarpetX {

// Truthy value written by BuildMask at coarse-fine ghost cells in the mask
// `iMultiFab` produced by LevelData::get_cf_mask. All other cells are 0.
inline constexpr int cf_ghost = 1;

} // namespace CarpetX

#endif // CARPETX_CF_MASK_HXX
