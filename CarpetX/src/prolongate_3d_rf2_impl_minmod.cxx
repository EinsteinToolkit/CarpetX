#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Minmod (tensor product) interpolation

static prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c000_o1;
static prolongate_3d_rf2<VC, VC, CC, POLY, POLY, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c001_o1;
static prolongate_3d_rf2<VC, CC, VC, POLY, MINMOD, POLY, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c010_o1;
static prolongate_3d_rf2<VC, CC, CC, POLY, MINMOD, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c011_o1;
static prolongate_3d_rf2<CC, VC, VC, MINMOD, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c100_o1;
static prolongate_3d_rf2<CC, VC, CC, MINMOD, POLY, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c101_o1;
static prolongate_3d_rf2<CC, CC, VC, MINMOD, MINMOD, POLY, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c110_o1;
static prolongate_3d_rf2<CC, CC, CC, MINMOD, MINMOD, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c111_o1;

static prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_minmod_3d_rf2_c000_o3;
static prolongate_3d_rf2<VC, VC, CC, POLY, POLY, MINMOD, 3, 3, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c001_o3;
static prolongate_3d_rf2<VC, CC, VC, POLY, MINMOD, POLY, 3, 1, 3, FB_NONE>
    prolongate_minmod_3d_rf2_c010_o3;
static prolongate_3d_rf2<VC, CC, CC, POLY, MINMOD, MINMOD, 3, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c011_o3;
static prolongate_3d_rf2<CC, VC, VC, MINMOD, POLY, POLY, 1, 3, 3, FB_NONE>
    prolongate_minmod_3d_rf2_c100_o3;
static prolongate_3d_rf2<CC, VC, CC, MINMOD, POLY, MINMOD, 1, 3, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c101_o3;
static prolongate_3d_rf2<CC, CC, VC, MINMOD, MINMOD, POLY, 1, 1, 3, FB_NONE>
    prolongate_minmod_3d_rf2_c110_o3;
static prolongate_3d_rf2<CC, CC, CC, MINMOD, MINMOD, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c111_o3;

static prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 5, 5, 5, FB_NONE>
    prolongate_minmod_3d_rf2_c000_o5;
static prolongate_3d_rf2<VC, VC, CC, POLY, POLY, MINMOD, 5, 5, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c001_o5;
static prolongate_3d_rf2<VC, CC, VC, POLY, MINMOD, POLY, 5, 1, 5, FB_NONE>
    prolongate_minmod_3d_rf2_c010_o5;
static prolongate_3d_rf2<VC, CC, CC, POLY, MINMOD, MINMOD, 5, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c011_o5;
static prolongate_3d_rf2<CC, VC, VC, MINMOD, POLY, POLY, 1, 5, 5, FB_NONE>
    prolongate_minmod_3d_rf2_c100_o5;
static prolongate_3d_rf2<CC, VC, CC, MINMOD, POLY, MINMOD, 1, 5, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c101_o5;
static prolongate_3d_rf2<CC, CC, VC, MINMOD, MINMOD, POLY, 1, 1, 5, FB_NONE>
    prolongate_minmod_3d_rf2_c110_o5;
static prolongate_3d_rf2<CC, CC, CC, MINMOD, MINMOD, MINMOD, 1, 1, 1, FB_NONE>
    prolongate_minmod_3d_rf2_c111_o5;

const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_minmod_3d_rf2{
        {1,
         {
             &prolongate_minmod_3d_rf2_c000_o1,
             &prolongate_minmod_3d_rf2_c001_o1,
             &prolongate_minmod_3d_rf2_c010_o1,
             &prolongate_minmod_3d_rf2_c011_o1,
             &prolongate_minmod_3d_rf2_c100_o1,
             &prolongate_minmod_3d_rf2_c101_o1,
             &prolongate_minmod_3d_rf2_c110_o1,
             &prolongate_minmod_3d_rf2_c111_o1,
         }},
        {3,
         {
             &prolongate_minmod_3d_rf2_c000_o3,
             &prolongate_minmod_3d_rf2_c001_o3,
             &prolongate_minmod_3d_rf2_c010_o3,
             &prolongate_minmod_3d_rf2_c011_o3,
             &prolongate_minmod_3d_rf2_c100_o3,
             &prolongate_minmod_3d_rf2_c101_o3,
             &prolongate_minmod_3d_rf2_c110_o3,
             &prolongate_minmod_3d_rf2_c111_o3,
         }},
        {5,
         {
             &prolongate_minmod_3d_rf2_c000_o5,
             &prolongate_minmod_3d_rf2_c001_o5,
             &prolongate_minmod_3d_rf2_c010_o5,
             &prolongate_minmod_3d_rf2_c011_o5,
             &prolongate_minmod_3d_rf2_c100_o5,
             &prolongate_minmod_3d_rf2_c101_o5,
             &prolongate_minmod_3d_rf2_c110_o5,
             &prolongate_minmod_3d_rf2_c111_o5,
         }},
    };

} // namespace CarpetX
