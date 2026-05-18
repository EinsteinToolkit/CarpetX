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

} // namespace TestVect

} // namespace Arith
