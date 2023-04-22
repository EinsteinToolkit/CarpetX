#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// DDF ENO (tensor product) interpolation

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c000_o1;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 1, 1, 0, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c001_o1;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 1, 0, 1, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c010_o1;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 1, 0, 0, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c011_o1;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 0, 1, 1, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c100_o1;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 0, 1, 0, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c101_o1;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 0, 0, 1, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c110_o1;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 0, 0, 0, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c111_o1;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c000_o3;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 3, 3, 2, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c001_o3;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 3, 2, 3, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c010_o3;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 3, 2, 2, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c011_o3;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 2, 3, 3, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c100_o3;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 2, 3, 2, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c101_o3;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 2, 2, 3, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c110_o3;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 2, 2, 2, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c111_o3;

prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c000_o5;
prolongate_3d_rf2<VC, VC, CC, POLY, POLY, ENO, 5, 5, 4, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c001_o5;
prolongate_3d_rf2<VC, CC, VC, POLY, ENO, POLY, 5, 4, 5, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c010_o5;
prolongate_3d_rf2<VC, CC, CC, POLY, ENO, ENO, 5, 4, 4, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c011_o5;
prolongate_3d_rf2<CC, VC, VC, ENO, POLY, POLY, 4, 5, 5, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c100_o5;
prolongate_3d_rf2<CC, VC, CC, ENO, POLY, ENO, 4, 5, 4, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c101_o5;
prolongate_3d_rf2<CC, CC, VC, ENO, ENO, POLY, 4, 4, 5, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c110_o5;
prolongate_3d_rf2<CC, CC, CC, ENO, ENO, ENO, 4, 4, 4, FB_NONE>
    prolongate_ddf_eno_3d_rf2_c111_o5;

} // namespace CarpetX
