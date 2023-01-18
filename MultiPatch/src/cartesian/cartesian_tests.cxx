#include "multipatch.hxx"
#include "tests.hxx"

#include <cctk_Parameters.h>

#include <string>
#include <sstream>

namespace MultiPatch {
namespace CartesianTests {

std::string local2global(const PatchTransformations &pt,
                         const vec<CCTK_REAL, dim> &x) {

  using namespace MultiPatchTests;

  std::ostringstream msg;
  msg << "has ";

  const auto l2g = pt.local2global(pt, 0, x);
  const bool eq_x = isapprox(l2g(0), x(0));
  const bool eq_y = isapprox(l2g(1), x(1));
  const bool eq_z = isapprox(l2g(2), x(2));

  if (eq_x && eq_y && eq_z) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason:";

    if (!eq_x) {
      msg << " The result in the x direction is " << l2g(0) << ".";
    }

    if (!eq_y) {
      msg << " The result in the y direction is " << l2g(1) << ".";
    }

    if (!eq_z) {
      msg << " The result in the z direction is " << l2g(02) << ".";
    }
  }

  return msg.str();
}

std::string global2local(const PatchTransformations &pt,
                         const vec<CCTK_REAL, dim> &x) {

  using namespace MultiPatchTests;

  std::ostringstream msg;
  msg << "has ";

  const auto g2l = pt.global2local(pt, x);
  const bool eq_x = isapprox(std::get<1>(g2l)(0), x(0));
  const bool eq_y = isapprox(std::get<1>(g2l)(1), x(1));
  const bool eq_z = isapprox(std::get<1>(g2l)(2), x(2));

  if (eq_x && eq_y && eq_z) {
    msg << MultiPatchTests::PASSED;
  } else {
    msg << MultiPatchTests::FAILED << ". Reason:";

    if (!eq_x) {
      msg << " The result in the x direction is " << std::get<1>(g2l)(0) << ".";
    }

    if (!eq_y) {
      msg << " The result in the y direction is " << std::get<1>(g2l)(1) << ".";
    }

    if (!eq_z) {
      msg << " The result in the z direction is " << std::get<1>(g2l)(2) << ".";
    }
  }

  return msg.str();
}

} // namespace CartesianTests

} // namespace MultiPatch

extern "C" void run_cartesian_tests() {
  DECLARE_CCTK_PARAMETERS

  using std::mt19937_64;
  using std::uniform_real_distribution;

  using namespace MultiPatch::CartesianTests;

  using MultiPatch::dim;
  using MultiPatch::SetupCartesian;
  using MultiPatchTests::random_seed;

  using Arith::vec;

  const auto ps = SetupCartesian();
  const auto &pt = ps.transformations;

  mt19937_64 engine(random_seed);
  uniform_real_distribution<CCTK_REAL> x_distrib(pt.cartesian_xmin,
                                                 pt.cartesian_xmax);
  uniform_real_distribution<CCTK_REAL> y_distrib(pt.cartesian_ymin,
                                                 pt.cartesian_ymax);
  uniform_real_distribution<CCTK_REAL> z_distrib(pt.cartesian_zmin,
                                                 pt.cartesian_zmax);
  auto point = vec<CCTK_REAL, dim>{0.0, 0.0, 0.0};

  CCTK_INFO("Running cartesian patch tests:");

  for (int i = 0; i < test_repetitions; i++) {
    point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    CCTK_VINFO("  Local to global transformation test at point (%f, %f, %f) %s",
               point(0), point(1), point(2), local2global(pt, point).c_str());
  }

  for (int i = 0; i < test_repetitions; i++) {
    point = {x_distrib(engine), y_distrib(engine), z_distrib(engine)};
    CCTK_VINFO("  Global to local transformation test at point (%f, %f, %f) %s",
               point(0), point(1), point(2), global2local(pt, point).c_str());
  }
}
