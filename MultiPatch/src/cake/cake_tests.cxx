#include "cake.hxx"
#include "tests.hxx"

#include <cctk_Parameters.h>

#include <cstddef>
#include <string>
#include <sstream>
#include <utility>
#include <cmath>

namespace MultiPatch {
namespace CakeTests {

/**
 * Tests the get_owner_patch function.
 *
 * @param pt The patch transformations structure
 * @param x A global point to test.
 * @param expected The expected result of the get_owner_patch execution.
 * @return A string indicating the test status.
 */
std::string patch_owner_test(const PatchTransformations &pt,
                             const MultiPatch::Cake::svec &x,
                             MultiPatch::Cake::patch_piece expected) {

  using namespace MultiPatch::Cake;
  using namespace MultiPatchTests;

  std::ostringstream msg;
  msg << "has ";

  const auto owner_patch = get_owner_patch(pt, x);

  if (owner_patch == expected) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: Expected to get patch "
        << piece_name(expected) << " and got " << piece_name(owner_patch);
  }

  msg << ".";
  return msg.str();
}

/**
 * Tests if local2global(global2local(global)) == global
 *
 * @param pt The patch transformations structure
 * @param global_vars A global point to test.
 * @return A string indicating the test status.
 */
std::string global_identity_test(const PatchTransformations &pt,
                                 const MultiPatch::Cake::svec &global_vars) {

  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "has ";

  const auto g2l = pt.global2local(pt, global_vars);

  const auto owner_patch_idx = std::get<0>(g2l);
  const auto local_vars = std::get<1>(g2l);

  const auto l2g = pt.local2global(pt, owner_patch_idx, local_vars);

  const auto test1 = isapprox(l2g(0), global_vars(0));
  const auto test2 = isapprox(l2g(1), global_vars(1));
  const auto test3 = isapprox(l2g(2), global_vars(2));

  if (test1 && test2 && test3) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << l2g(0) << " =/= " << global_vars(0) << ". ";
    }

    if (!test2) {
      msg << l2g(1) << " =/= " << global_vars(1) << ". ";
    }

    if (!test3) {
      msg << l2g(2) << " =/= " << global_vars(2) << ". ";
    }
  }

  return msg.str();
}

/**
 * Tests if global2local(local2global(local, patch)) == (local, patch)
 *
 * @param pt The patch transformations structure
 * @param patch The patch index to test.
 * @param local_point A local point to test.
 * @return A string indicating the test status.
 */
std::string local_identity_test(const PatchTransformations &pt, int patch,
                                const MultiPatch::Cake::svec &local_point) {

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatchTests::colored;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "has ";

  const auto l2g = pt.local2global(pt, patch, local_point);
  const auto g2l = pt.global2local(pt, l2g);

  const auto computed_patch_idx = std::get<0>(g2l);
  const auto computed_local_point = std::get<1>(g2l);

  const bool test1 = computed_patch_idx == patch;
  const bool test2 = isapprox(computed_local_point(0), local_point(0));
  const bool test3 = isapprox(computed_local_point(1), local_point(1));
  const bool test4 = isapprox(computed_local_point(2), local_point(2));

  if (test1 && test2 && test3 && test4) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << "the computed patch is "
          << piece_name(static_cast<patch_piece>(computed_patch_idx)) << ". ";
    }

    if (!test2) {
      msg << "the computed local coordinate a is " << computed_local_point(0)
          << ". ";
    }

    if (!test3) {
      msg << "the computed local coordinate a is " << computed_local_point(1)
          << ". ";
    }

    if (!test4) {
      msg << "the computed local coordinate a is " << computed_local_point(2)
          << ". ";
    }
  }

  return msg.str();
}

/**
 * Tests the cake jacobian implementations by approximating da^i/dx^i with
 * fourth order finite differences. The "step size" and comparison tolerance is
 * customizable at compile time.
 *
 * @param pt The patch transformations structure
 * @param patch The patch index to test.
 * @param global_point A global point to test.
 * @return A string indicating the test status.
 */
std::string jacobian_test(const PatchTransformations &pt, int patch,
                          const MultiPatch::Cake::svec &global_point) {

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatch::Cake::svec;
  using MultiPatchTests::colored;
  using MultiPatchTests::fd_4;
  using MultiPatchTests::fd_comp_tol;
  using MultiPatchTests::fd_delta;
  using MultiPatchTests::fd_direction;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "patch ";

  // Compute a local point and patch number from the global poiint
  const auto local_data = pt.global2local(pt, global_point);
  msg << piece_name(static_cast<patch_piece>(std::get<0>(local_data)))
      << " has ";

  // From local point and patch number, we can compute the jacobian
  const auto J_data =
      pt.dlocal_dglobal(pt, std::get<0>(local_data), std::get<1>(local_data));
  const auto &J = std::get<1>(J_data);

  // Compute the second derivative of local2global by finite differencing
  const auto g2l_wrapper = [&](const svec &point) -> svec {
    return std::get<1>(pt.global2local(pt, point));
  };

  const auto expected_dx =
      fd_4<fd_direction::x, svec>(g2l_wrapper, global_point);
  const auto expected_dy =
      fd_4<fd_direction::y, svec>(g2l_wrapper, global_point);
  const auto expected_dz =
      fd_4<fd_direction::z, svec>(g2l_wrapper, global_point);

  const bool test1 = isapprox(expected_dx(0), J(0)(0), fd_comp_tol);
  const bool test2 = isapprox(expected_dx(1), J(1)(0), fd_comp_tol);
  const bool test3 = isapprox(expected_dx(2), J(2)(0), fd_comp_tol);

  const bool test4 = isapprox(expected_dy(0), J(0)(1), fd_comp_tol);
  const bool test5 = isapprox(expected_dy(1), J(1)(1), fd_comp_tol);
  const bool test6 = isapprox(expected_dy(2), J(2)(1), fd_comp_tol);

  const bool test7 = isapprox(expected_dz(0), J(0)(2), fd_comp_tol);
  const bool test8 = isapprox(expected_dz(1), J(1)(2), fd_comp_tol);
  const bool test9 = isapprox(expected_dz(2), J(2)(2), fd_comp_tol);

  const bool all_tests = test1 && test2 && test3 && test4 && test5 && test6;

  if (all_tests) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << "computed J(0)(0) value is " << (J(0)(0))
          << " and the expected value is " << expected_dx(0) << ". ";
    }

    if (!test2) {
      msg << "computed J(1)(0) value is " << (J(1)(0))
          << " and the expected value is " << expected_dx(1) << ". ";
    }

    if (!test3) {
      msg << "computed J(2)(0) value is " << (J(2)(0))
          << " and the expected value is " << expected_dx(2) << ". ";
    }

    if (!test4) {
      msg << "computed  J(0)(1) value is " << (J(0)(1))
          << " and the expected value is " << expected_dy(0) << ". ";
    }

    if (!test5) {
      msg << "computed  J(1)(1) value is " << (J(1)(1))
          << " and the expected value is " << expected_dy(1) << ". ";
    }

    if (!test6) {
      msg << "computed  J(2)(1) value is " << (J(2)(1))
          << " and the expected value is " << expected_dy(2) << ". ";
    }

    if (!test7) {
      msg << "computed  J(0)(2) value is " << (J(0)(2))
          << " and the expected value is " << expected_dz(0) << ". ";
    }

    if (!test8) {
      msg << "computed  J(1)(2) value is " << (J(1)(2))
          << " and the expected value is " << expected_dz(1) << ". ";
    }

    if (!test9) {
      msg << "computed  J(2)(2) value is " << (J(2)(2))
          << " and the expected value is " << expected_dz(2) << ". ";
    }
  }

  return msg.str();
}

std::string djacobian_test(const PatchTransformations &pt, int patch,
                           const MultiPatch::Cake::svec &global_point) {

  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatch::Cake::svec;
  using MultiPatchTests::colored;
  using MultiPatchTests::fd2_4;
  using MultiPatchTests::fd_4;
  using MultiPatchTests::fd_comp_tol;
  using MultiPatchTests::fd_delta;
  using MultiPatchTests::fd_direction;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::string_color;

  std::ostringstream msg;
  msg << "patch ";

  // Compute a local point and patch number from the global poiint
  const auto local_data = pt.global2local(pt, global_point);
  msg << piece_name(static_cast<patch_piece>(std::get<0>(local_data)))
      << " has ";

  // From local point and patch number, we can compute the jacobian derivative
  const auto J_data =
      pt.d2local_dglobal2(pt, std::get<0>(local_data), std::get<1>(local_data));
  const auto &dJ = std::get<2>(J_data);

  // Compute the derivative of the jacobian by finite differencing
  const auto g2l_wrapper = [&](const svec &point) -> svec {
    return std::get<1>(pt.global2local(pt, point));
  };

  const auto expected_dxdx =
      fd2_4<fd_direction::x, fd_direction::x, svec>(g2l_wrapper, global_point);
  const auto expected_dxdy =
      fd2_4<fd_direction::x, fd_direction::y, svec>(g2l_wrapper, global_point);
  const auto expected_dxdz =
      fd2_4<fd_direction::x, fd_direction::z, svec>(g2l_wrapper, global_point);
  const auto expected_dydy =
      fd2_4<fd_direction::y, fd_direction::y, svec>(g2l_wrapper, global_point);
  const auto expected_dydz =
      fd2_4<fd_direction::y, fd_direction::z, svec>(g2l_wrapper, global_point);
  const auto expected_dzdz =
      fd2_4<fd_direction::z, fd_direction::z, svec>(g2l_wrapper, global_point);

  const bool test1 = isapprox(expected_dxdx(0), dJ(0)(0, 0), fd_comp_tol);
  const bool test2 = isapprox(expected_dxdx(1), dJ(1)(0, 0), fd_comp_tol);
  const bool test3 = isapprox(expected_dxdx(2), dJ(2)(0, 0), fd_comp_tol);

  const bool test4 = isapprox(expected_dxdy(0), dJ(0)(0, 1), fd_comp_tol);
  const bool test5 = isapprox(expected_dxdy(1), dJ(1)(0, 1), fd_comp_tol);
  const bool test6 = isapprox(expected_dxdy(2), dJ(2)(0, 1), fd_comp_tol);

  const bool test7 = isapprox(expected_dxdz(0), dJ(0)(0, 2), fd_comp_tol);
  const bool test8 = isapprox(expected_dxdz(1), dJ(1)(0, 2), fd_comp_tol);
  const bool test9 = isapprox(expected_dxdz(2), dJ(2)(0, 2), fd_comp_tol);

  const bool test10 = isapprox(expected_dydy(0), dJ(0)(1, 1), fd_comp_tol);
  const bool test11 = isapprox(expected_dydy(1), dJ(1)(1, 1), fd_comp_tol);
  const bool test12 = isapprox(expected_dydy(2), dJ(2)(1, 1), fd_comp_tol);

  const bool test13 = isapprox(expected_dydz(0), dJ(0)(1, 2), fd_comp_tol);
  const bool test14 = isapprox(expected_dydz(1), dJ(1)(1, 2), fd_comp_tol);
  const bool test15 = isapprox(expected_dydz(2), dJ(2)(1, 2), fd_comp_tol);

  const bool test16 = isapprox(expected_dzdz(0), dJ(0)(2, 2), fd_comp_tol);
  const bool test17 = isapprox(expected_dzdz(1), dJ(1)(2, 2), fd_comp_tol);
  const bool test18 = isapprox(expected_dzdz(2), dJ(2)(2, 2), fd_comp_tol);

  const bool all_tests = test1 && test2 && test3 && test4 && test5 && test6 &&
                         test7 && test8 && test9 && test10 && test11 &&
                         test12 && test13 && test14 && test15 && test16 &&
                         test17 && test18;

  if (all_tests) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason: ";

    if (!test1) {
      msg << "computed dJ(0)(0, 0) value is " << (dJ(0)(0, 0))
          << " and the expected value is " << expected_dxdx(0) << ". ";
    }

    if (!test2) {
      msg << "computed dJ(1)(0, 0) value is " << (dJ(1)(0, 0))
          << " and the expected value is " << expected_dxdx(1) << ". ";
    }

    if (!test3) {
      msg << "computed dJ(2)(0, 0) value is " << (dJ(2)(0, 0))
          << " and the expected value is " << expected_dxdx(2) << ". ";
    }

    if (!test4) {
      msg << "computed dJ(0)(0, 1) value is " << (dJ(0)(0, 1))
          << " and the expected value is " << expected_dxdy(0) << ". ";
    }

    if (!test5) {
      msg << "computed dJ(1)(0, 1) value is " << (dJ(1)(0, 1))
          << " and the expected value is " << expected_dxdy(1) << ". ";
    }

    if (!test6) {
      msg << "computed dJ(2)(0, 1) value is " << (dJ(2)(0, 1))
          << " and the expected value is " << expected_dxdy(2) << ". ";
    }

    if (!test7) {
      msg << "computed dJ(0)(0, 2) value is " << (dJ(0)(0, 2))
          << " and the expected value is " << expected_dxdz(0) << ". ";
    }

    if (!test8) {
      msg << "computed dJ(1)(0, 2) value is " << (dJ(1)(0, 2))
          << " and the expected value is " << expected_dxdz(1) << ". ";
    }

    if (!test9) {
      msg << "computed dJ(2)(0, 2) value is " << (dJ(2)(0, 2))
          << " and the expected value is " << expected_dxdz(2) << ". ";
    }

    if (!test10) {
      msg << "computed dJ(0)(1, 1) value is " << (dJ(0)(1, 1))
          << " and the expected value is " << expected_dydy(0) << ". ";
    }

    if (!test11) {
      msg << "computed dJ(1)(1, 1) value is " << (dJ(1)(1, 1))
          << " and the expected value is " << expected_dydy(1) << ". ";
    }

    if (!test12) {
      msg << "computed dJ(2)(1, 1) value is " << (dJ(2)(1, 1))
          << " and the expected value is " << expected_dydy(2) << ". ";
    }

    if (!test13) {
      msg << "computed dJ(0)(1, 2) value is " << (dJ(0)(1, 2))
          << " and the expected value is " << expected_dydz(0) << ". ";
    }

    if (!test14) {
      msg << "computed dJ(1)(1, 2) value is " << (dJ(1)(1, 2))
          << " and the expected value is " << expected_dydz(1) << ". ";
    }

    if (!test15) {
      msg << "computed dJ(2)(1, 2) value is " << (dJ(2)(1, 2))
          << " and the expected value is " << expected_dydz(2) << ". ";
    }

    if (!test16) {
      msg << "computed dJ(0)(2, 2) value is " << (dJ(0)(2, 2))
          << " and the expected value is " << expected_dzdz(0) << ". ";
    }

    if (!test17) {
      msg << "computed dJ(1)(2, 2) value is " << (dJ(1)(2, 2))
          << " and the expected value is " << expected_dzdz(1) << ". ";
    }

    if (!test18) {
      msg << "computed dJ(2)(2, 2) value is " << (dJ(2)(2, 2))
          << " and the expected value is " << expected_dzdz(2) << ". ";
    }
  }

  return msg.str();
}

} // namespace CakeTests
} // namespace MultiPatch

/**
 * Runs all tests pertaining to the cake patch
 */
extern "C" void run_cake_tests() {
  DECLARE_CCTK_PARAMETERS

  using std::mt19937_64;
  using std::uniform_int_distribution;
  using std::uniform_real_distribution;

  using std::cos;
  using std::sin;

  using MultiPatch::Cake::get_owner_patch;
  using MultiPatch::Cake::patch_piece;
  using MultiPatch::Cake::piece_name;
  using MultiPatch::Cake::svec;

  using MultiPatch::dim;
  using MultiPatch::SetupCake;
  using MultiPatchTests::random_seed;

  using namespace MultiPatch::CakeTests;

  const auto ps = SetupCake();
  const auto &pt = ps.transformations;

  const auto r0 = pt.cake_inner_boundary_radius;
  const auto r1 = pt.cake_outer_boundary_radius;

  const auto cube_midpoint = r0 / 2.0;
  const auto patch_midpoint = r0 + (r1 - r0) / 2.0;

  CCTK_INFO("Running cake patch tests:");

  /*
   * Fixed point tests
   */
  const std::array<std::pair<svec, patch_piece>, 27> owner_test_data = {
      // Cartesian cube interior
      std::make_pair(svec{cube_midpoint, cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec{-cube_midpoint, cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec{cube_midpoint, -cube_midpoint, cube_midpoint},
                     patch_piece::cartesian),
      std::make_pair(svec{cube_midpoint, cube_midpoint, -cube_midpoint},
                     patch_piece::cartesian),

      // Cartesian cube corners
      std::make_pair(svec{r0, r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec{-r0, r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec{r0, -r0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec{r0, r0, -r0}, patch_piece::inner_boundary),

      // Middle of the +x patch
      std::make_pair(svec{r0, 0.0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec{patch_midpoint, 0.0, 0.0}, patch_piece::plus_x),
      std::make_pair(svec{r1, 0.0, 0.0}, patch_piece::outer_boundary),

      // Middle of the -x patch
      std::make_pair(svec{-r0, 0.0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec{-patch_midpoint, 0.0, 0.0}, patch_piece::minus_x),
      std::make_pair(svec{-r1, 0.0, 0.0}, patch_piece::outer_boundary),

      // Middle of the +y patch
      std::make_pair(svec{0.0, r0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec{0.0, patch_midpoint, 0.0}, patch_piece::plus_y),
      std::make_pair(svec{0.0, r1, 0.0}, patch_piece::outer_boundary),

      // Middle of the -y patch
      std::make_pair(svec{0.0, -r0, 0.0}, patch_piece::inner_boundary),
      std::make_pair(svec{0.0, -patch_midpoint, 0.0}, patch_piece::minus_y),
      std::make_pair(svec{0.0, -r1, 0.0}, patch_piece::outer_boundary),

      // Middle of the +z patch
      std::make_pair(svec{0.0, 0.0, r0}, patch_piece::inner_boundary),
      std::make_pair(svec{0.0, 0.0, patch_midpoint}, patch_piece::plus_z),
      std::make_pair(svec{0.0, 0.0, r1}, patch_piece::outer_boundary),

      // Middle of the -z patch
      std::make_pair(svec{0.0, 0.0, -r0}, patch_piece::inner_boundary),
      std::make_pair(svec{0.0, 0.0, -patch_midpoint}, patch_piece::minus_z),
      std::make_pair(svec{0.0, 0.0, -r1}, patch_piece::outer_boundary),

      std::make_pair(svec{r1 + 1.0, r1 + 1.0, r1 + 1.0},
                     patch_piece::exterior)};

  // Tests if the patch owner is correct.
  for (const auto &data : owner_test_data) {
    CCTK_VINFO("  Patch owner test at point (%f, %f, %f) %s", data.first(0),
               data.first(1), data.first(2),
               patch_owner_test(pt, data.first, data.second).c_str());
  }

  /*
   * Random point tests
   */
  mt19937_64 engine(random_seed);

  uniform_real_distribution<CCTK_REAL> r_distrib(0.0, r1);
  uniform_real_distribution<CCTK_REAL> theta_distrib(0.0, M_PI);
  uniform_real_distribution<CCTK_REAL> phi_distrib(0.0, 2.0 * M_PI);

  uniform_real_distribution<CCTK_REAL> local_distrib(-1.0, 1.0);

  uniform_int_distribution<int> patch_distrib(0, 6);

  auto global_point = svec{0.0, 0.0, 0.0};
  auto local_point = svec{0.0, 0.0, 0.0};
  double r = 0.0, theta = 0.0, phi = 0.0;
  int patch = 0;

  // Tests if local2global(global2local(global)) == global
  for (int i = 0; i < test_repetitions; i++) {
    r = r_distrib(engine);
    theta = theta_distrib(engine);
    phi = phi_distrib(engine);

    global_point = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
                    r * cos(theta)};

    CCTK_VINFO("  local2global(global2local(global)) transformation test at "
               "point (%f, %f, %f) %s",
               global_point(0), global_point(1), global_point(2),
               global_identity_test(pt, global_point).c_str());
  }

  // Tests if global2local(local2global(local, patch)) == (local, patch)
  for (int i = 0; i < test_repetitions; i++) {
    local_point = {local_distrib(engine), local_distrib(engine),
                   local_distrib(engine)};

    patch = patch_distrib(engine);

    CCTK_VINFO("  global2local(local2global(local, patch)) transformation test "
               "at point (%f, %f, %f) patch %s %s",
               local_point(0), local_point(1), local_point(2),
               piece_name(static_cast<patch_piece>(patch)).c_str(),
               local_identity_test(pt, patch, local_point).c_str());
  }

  // Tests if local -> global jacobians are correct.
  for (int i = 0; i < test_repetitions; i++) {
    r = r_distrib(engine);
    theta = theta_distrib(engine);
    phi = phi_distrib(engine);

    global_point = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
                    r * cos(theta)};

    CCTK_VINFO("  local -> global jacobian test at point (%f, %f, %f) patch %s",
               global_point(0), global_point(1), global_point(2),
               jacobian_test(pt, patch, global_point).c_str());
  }

  // Tests if local -> global jacobian derivatives are correct.
  for (int i = 0; i < test_repetitions; i++) {
    r = r_distrib(engine);
    theta = theta_distrib(engine);
    phi = phi_distrib(engine);

    global_point = {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi),
                    r * cos(theta)};

    CCTK_VINFO("  local -> global jacobian derivative test at point (%f, %f, "
               "%f) patch %s",
               global_point(0), global_point(1), global_point(2),
               djacobian_test(pt, patch, global_point).c_str());
  }
}
