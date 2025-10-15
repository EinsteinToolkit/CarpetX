#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Interpolate polynomially in vertex centred directions and conserve
// with 3rd order accuracy in cell centred directions

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_eno_star_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO_STAR, 5, 5, 2, FB_NONE>
    prolongate_eno_star_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO_STAR, POLY, 5, 2, 5, FB_NONE>
    prolongate_eno_star_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO_STAR, ENO_STAR, 5, 2, 2, FB_NONE>
    prolongate_eno_star_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, ENO_STAR, POLY, POLY, 2, 5, 5, FB_NONE>
    prolongate_eno_star_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, ENO_STAR, POLY, ENO_STAR, 2, 5, 2, FB_NONE>
    prolongate_eno_star_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, ENO_STAR, ENO_STAR, POLY, 2, 2, 5, FB_NONE>
    prolongate_eno_star_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, ENO_STAR, ENO_STAR, ENO_STAR, 2, 2, 2, FB_NONE>
    prolongate_eno_star_3d_rf2_c111_o5;

} // namespace CarpetX
