#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Order 5 operators of prolongate_3d_rf2_impl_ddf.cxx. Each order is compiled
// separately to keep compile times short.

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 5, 5, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 5, 4, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 5, 4, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 4, 5, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 4, 5, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 4, 4, 5, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 4, 4, 4, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o5;

} // namespace CarpetX
