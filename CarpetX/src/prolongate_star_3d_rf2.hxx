#ifndef CARPETX_CARPETX_PROLONGATE_STAR_3D_RF2_HXX
#define CARPETX_CARPETX_PROLONGATE_STAR_3D_RF2_HXX

#include <AMReX_Interpolater.H>

#include <array>
#include <map>

namespace CarpetX {

using interpolator_family_t =
    std::map<int, std::array<amrex::Interpolater *, 8> >;

extern const interpolator_family_t prolongate_star_poly_3d_rf2;
extern const interpolator_family_t prolongate_star_cons_3d_rf2;

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_PROLONGATE_STAR_3D_RF2_HXX
