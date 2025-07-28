#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Hermite interpolation

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o3;

prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o5;

} // namespace CarpetX
