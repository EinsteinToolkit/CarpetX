#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Order 3 operators of prolongate_3d_rf2_impl_ddf.cxx. Each order is compiled
// separately to keep compile times short.

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 2, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 2, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 2, 3, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 2, 3, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 2, 2, 3, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 2, 2, 2, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o3;

} // namespace CarpetX
