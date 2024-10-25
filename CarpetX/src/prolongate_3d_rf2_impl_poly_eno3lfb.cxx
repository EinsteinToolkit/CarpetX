#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Interpolate polynomially in vertex centred directions and conserve
// with 3rd order accuracy and a linear fallback in cell centred
// directions

// prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c000_o1;
// prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c001_o1;
// prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c010_o1;
// prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c011_o1;
// prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c100_o1;
// prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c101_o1;
// prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c110_o1;
// prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 1, 1, 1, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c111_o1;
//
// prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c000_o3;
// prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c001_o3;
// prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c010_o3;
// prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c011_o3;
// prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c100_o3;
// prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c101_o3;
// prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c110_o3;
// prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c111_o3;
//
// prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c000_o5;
// prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 5, 5, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c001_o5;
// prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 5, 3, 5, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c010_o5;
// prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 5, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c011_o5;
// prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 3, 5, 5, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c100_o5;
// prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 3, 5, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c101_o5;
// prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 3, 3, 5, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c110_o5;
// prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 3, 3, 3, FB_LINEAR>
//     prolongate_poly_eno3lfb_3d_rf2_c111_o5;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 0, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 0, 1, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 0, 0, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 0, 1, 1, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 0, 1, 0, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 0, 0, 1, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 3, 3, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 3, 2, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 3, 2, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 2, 3, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 2, 3, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 2, 2, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 2, 2, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c111_o3;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 5, 5, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 5, 2, 5, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 5, 2, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 2, 5, 5, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 2, 5, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 2, 2, 5, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 2, 2, 2, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c111_o5;

#if 0
prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c000_o7;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 7, 7, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c001_o7;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 7, 3, 7, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c010_o7;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 7, 3, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c011_o7;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 3, 7, 7, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c100_o7;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 3, 7, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c101_o7;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 3, 3, 7, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c110_o7;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 3, 3, 3, FB_LINEAR>
    prolongate_poly_eno3lfb_3d_rf2_c111_o7;
#endif

} // namespace CarpetX
