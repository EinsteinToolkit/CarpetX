#include "vect.hxx"

namespace Arith {

namespace TestVect {
// constexpr vect<int, 3> vect1(make_tuple(0, 1, 2));
constexpr vect<int, 3> vect1(std::array<int, 3>{0, 1, 2});
static_assert(vect1[0] == 0, "");
static_assert(vect1[1] == 1, "");
static_assert(vect1[2] == 2, "");

// constexpr vect<vect<int, 3>, 3>
//     vect2(make_tuple(vect<int, 3>(make_tuple(100, 101, 102)),
//                      vect<int, 3>(make_tuple(110, 111, 112)),
//                      vect<int, 3>(make_tuple(120, 121, 122))));
constexpr vect<vect<int, 3>, 3> vect2(std::array<vect<int, 3>, 3>{
    vect<int, 3>(std::array<int, 3>{100, 101, 102}),
    vect<int, 3>(std::array<int, 3>{110, 111, 112}),
    vect<int, 3>(std::array<int, 3>{120, 121, 122})});
static_assert(vect2[0][0] == 100, "");
static_assert(vect2[0][1] == 101, "");
static_assert(vect2[0][2] == 102, "");
static_assert(vect2[1][0] == 110, "");
static_assert(vect2[1][1] == 111, "");
static_assert(vect2[1][2] == 112, "");
static_assert(vect2[2][0] == 120, "");
static_assert(vect2[2][1] == 121, "");
static_assert(vect2[2][2] == 122, "");

// Comparisons against a scalar, with the scalar on either side. These must
// agree component by component with the underlying scalar comparison.
constexpr vect<int, 3> vect3(std::array<int, 3>{1, 2, 3});
static_assert(all((vect3 < 2) == vect<bool, 3>{true, false, false}), "");
static_assert(all((vect3 > 2) == vect<bool, 3>{false, false, true}), "");
static_assert(all((vect3 <= 2) == vect<bool, 3>{true, true, false}), "");
static_assert(all((vect3 >= 2) == vect<bool, 3>{false, true, true}), "");
static_assert(all((2 < vect3) == vect<bool, 3>{false, false, true}), "");
static_assert(all((2 > vect3) == vect<bool, 3>{true, false, false}), "");
static_assert(all((2 <= vect3) == vect<bool, 3>{false, true, true}), "");
static_assert(all((2 >= vect3) == vect<bool, 3>{true, true, false}), "");

} // namespace TestVect

} // namespace Arith
