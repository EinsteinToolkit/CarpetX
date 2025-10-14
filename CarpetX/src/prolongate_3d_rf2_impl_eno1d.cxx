#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Interpolate polynomially in vertex centred directions and conserve
// with 3rd order accuracy in cell centred directions

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_eno1d_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO1D, 5, 5, 2, FB_NONE>
    prolongate_eno1d_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO1D, POLY, 5, 2, 5, FB_NONE>
    prolongate_eno1d_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO1D, ENO1D, 5, 2, 2, FB_NONE>
    prolongate_eno1d_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, ENO1D, POLY, POLY, 2, 5, 5, FB_NONE>
    prolongate_eno1d_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, ENO1D, POLY, ENO1D, 2, 5, 2, FB_NONE>
    prolongate_eno1d_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, ENO1D, ENO1D, POLY, 2, 2, 5, FB_NONE>
    prolongate_eno1d_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, ENO1D, ENO1D, ENO1D, 2, 2, 2, FB_NONE>
    prolongate_eno1d_3d_rf2_c111_o5;

} // namespace CarpetX
