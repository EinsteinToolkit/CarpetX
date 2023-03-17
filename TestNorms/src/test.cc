#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <iostream>

#include "loop.hxx"

namespace TestNorms {
using namespace std;

CCTK_REAL fun1d(const CCTK_REAL x, const CCTK_REAL dx, const bool avgx,
                const int order) {
  if (!avgx)
    return pow(x, order);
  else
    return (pow(x + dx / 2, order + 1) - pow(x - dx / 2, order + 1)) /
           ((order + 1) * dx);
}

// The grid stores the average values of the underlying function (x*y*z)**n
// which amounts to storing differences of the anti-derivative 1/(n+1)**3 *
// (x*y*z)**(n+1)
CCTK_REAL fun(const Loop::PointDesc &p, const bool avgx, const bool avgy,
              const bool avgz, const int order) {
  return fun1d(p.x, p.dx, avgx, order) * fun1d(p.y, p.dy, avgy, order) *
         fun1d(p.z, p.dz, avgz, order);
}

extern "C" void TestNorms_SetError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestNorms_SetError;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
        fabs(p.z) <= refined_radius) {
      regrid_error(p.I) = 1;
    } else {
      regrid_error(p.I) = 0;
    }
  });
}

extern "C" void TestNorms_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestNorms_Set;
  DECLARE_CCTK_PARAMETERS;

  const int order = 3; // an order that is higher than linear so that the value
                       // at the cell center is not the same as the average
                       // over the cell

  Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    gf000(p.I) = fun(p, 0, 0, 0, order);
  });
  Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    gf001(p.I) = fun(p, 0, 0, 1, order);
  });
  Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    gf010(p.I) = fun(p, 0, 1, 0, order);
  });
  Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    gf011(p.I) = fun(p, 0, 1, 1, order);
  });
  Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    gf100(p.I) = fun(p, 1, 0, 0, order);
  });
  Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    gf101(p.I) = fun(p, 1, 0, 1, order);
  });
  Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    gf110(p.I) = fun(p, 1, 1, 0, order);
  });
  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    gf111(p.I) = fun(p, 1, 1, 1, order);
  });
}

} // namespace TestNorms
