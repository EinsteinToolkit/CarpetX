#include "prolongate_3d_rf2_impl.hxx"

namespace CarpetX {

// Hermite interpolation

static prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c000_o1;
static prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o1;
static prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o1;
static prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o1;
static prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o1;
static prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o1;
static prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o1;
static prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 1, 1, 1, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o1;

static prolongate_3d_rf2<VC, VC, VC, POLY, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c000_o3;
static prolongate_3d_rf2<VC, VC, CC, POLY, POLY, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o3;
static prolongate_3d_rf2<VC, CC, VC, POLY, CONS, POLY, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o3;
static prolongate_3d_rf2<VC, CC, CC, POLY, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o3;
static prolongate_3d_rf2<CC, VC, VC, CONS, POLY, POLY, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o3;
static prolongate_3d_rf2<CC, VC, CC, CONS, POLY, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o3;
static prolongate_3d_rf2<CC, CC, VC, CONS, CONS, POLY, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o3;
static prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 3, 3, 3, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o3;

static prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 5, 5, 5,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o5;
static prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o5;
static prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o5;
static prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o5;
static prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o5;
static prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o5;
static prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o5;
static prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 5, 5, 5, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o5;

static prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 7, 7, 7,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o7;
static prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o7;
static prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o7;
static prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o7;
static prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o7;
static prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o7;
static prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o7;
static prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 7, 7, 7, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o7;

static prolongate_3d_rf2<VC, VC, VC, HERMITE, HERMITE, HERMITE, 9, 9, 9,
                         FB_NONE>
    prolongate_hermite_3d_rf2_c000_o9;
static prolongate_3d_rf2<VC, VC, CC, HERMITE, HERMITE, CONS, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c001_o9;
static prolongate_3d_rf2<VC, CC, VC, HERMITE, CONS, HERMITE, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c010_o9;
static prolongate_3d_rf2<VC, CC, CC, HERMITE, CONS, CONS, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c011_o9;
static prolongate_3d_rf2<CC, VC, VC, CONS, HERMITE, HERMITE, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c100_o9;
static prolongate_3d_rf2<CC, VC, CC, CONS, HERMITE, CONS, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c101_o9;
static prolongate_3d_rf2<CC, CC, VC, CONS, CONS, HERMITE, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c110_o9;
static prolongate_3d_rf2<CC, CC, CC, CONS, CONS, CONS, 9, 9, 9, FB_NONE>
    prolongate_hermite_3d_rf2_c111_o9;

const std::map<int, std::array<amrex::Interpolater *, 8> >
    prolongate_hermite_3d_rf2{
        {1,
         {
             &prolongate_hermite_3d_rf2_c000_o1,
             &prolongate_hermite_3d_rf2_c001_o1,
             &prolongate_hermite_3d_rf2_c010_o1,
             &prolongate_hermite_3d_rf2_c011_o1,
             &prolongate_hermite_3d_rf2_c100_o1,
             &prolongate_hermite_3d_rf2_c101_o1,
             &prolongate_hermite_3d_rf2_c110_o1,
             &prolongate_hermite_3d_rf2_c111_o1,
         }},
        {3,
         {
             &prolongate_hermite_3d_rf2_c000_o3,
             &prolongate_hermite_3d_rf2_c001_o3,
             &prolongate_hermite_3d_rf2_c010_o3,
             &prolongate_hermite_3d_rf2_c011_o3,
             &prolongate_hermite_3d_rf2_c100_o3,
             &prolongate_hermite_3d_rf2_c101_o3,
             &prolongate_hermite_3d_rf2_c110_o3,
             &prolongate_hermite_3d_rf2_c111_o3,
         }},
        {5,
         {
             &prolongate_hermite_3d_rf2_c000_o5,
             &prolongate_hermite_3d_rf2_c001_o5,
             &prolongate_hermite_3d_rf2_c010_o5,
             &prolongate_hermite_3d_rf2_c011_o5,
             &prolongate_hermite_3d_rf2_c100_o5,
             &prolongate_hermite_3d_rf2_c101_o5,
             &prolongate_hermite_3d_rf2_c110_o5,
             &prolongate_hermite_3d_rf2_c111_o5,
         }},
        {7,
         {
             &prolongate_hermite_3d_rf2_c000_o7,
             &prolongate_hermite_3d_rf2_c001_o7,
             &prolongate_hermite_3d_rf2_c010_o7,
             &prolongate_hermite_3d_rf2_c011_o7,
             &prolongate_hermite_3d_rf2_c100_o7,
             &prolongate_hermite_3d_rf2_c101_o7,
             &prolongate_hermite_3d_rf2_c110_o7,
             &prolongate_hermite_3d_rf2_c111_o7,
         }},
        {9,
         {
             &prolongate_hermite_3d_rf2_c000_o9,
             &prolongate_hermite_3d_rf2_c001_o9,
             &prolongate_hermite_3d_rf2_c010_o9,
             &prolongate_hermite_3d_rf2_c011_o9,
             &prolongate_hermite_3d_rf2_c100_o9,
             &prolongate_hermite_3d_rf2_c101_o9,
             &prolongate_hermite_3d_rf2_c110_o9,
             &prolongate_hermite_3d_rf2_c111_o9,
         }},
    };

} // namespace CarpetX
