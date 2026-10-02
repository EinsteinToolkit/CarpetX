#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Order 7 operators of prolongate_3d_rf2_impl_ddf.cxx. Each order is compiled
// separately to keep compile times short.

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 7, 7, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o7;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 7, 7, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o7;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 7, 6, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o7;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 7, 6, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o7;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 6, 7, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o7;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 6, 7, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o7;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 6, 6, 7, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o7;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 6, 6, 6, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o7;

} // namespace CarpetX
