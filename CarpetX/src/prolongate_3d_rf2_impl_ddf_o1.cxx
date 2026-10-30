#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Order 1 operators of prolongate_3d_rf2_impl_ddf.cxx. Each order is compiled
// separately to keep compile times short.

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 0, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 0, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 0, 1, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 0, 1, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 0, 0, 1, FB_NONE>
    prolongate_ddf_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 0, 0, 0, FB_NONE>
    prolongate_ddf_3d_rf2_c111_o1;

} // namespace CarpetX
