#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "loop.hxx"

namespace {
// w(t,x,y,z) = (a + b x + c y + d z) * t  -- linear in time and space.
static inline CCTK_REAL analytic(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                                 CCTK_REAL z, CCTK_REAL a, CCTK_REAL b,
                                 CCTK_REAL c, CCTK_REAL d) {
  return (a + b * x + c * y + d * z) * t;
}
} // namespace

extern "C" void TestSubcyclingProlong_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingProlong_Init;
  grid.loop_int<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    snapshot(p.I) =
        analytic(cctk_time, p.x, p.y, p.z, poly_a, poly_b, poly_c, poly_d);
  });
}

extern "C" void TestSubcyclingProlong_Update(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingProlong_Update;
  grid.loop_int<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    snapshot(p.I) =
        analytic(cctk_time, p.x, p.y, p.z, poly_a, poly_b, poly_c, poly_d);
  });
}

extern "C" void TestSubcyclingProlong_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_TestSubcyclingProlong_Error;
  // Read snapshot everywhere (interior + ghosts) so the C/F-boundary ghost
  // error is captured.
  grid.loop_all<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    error(p.I) = snapshot(p.I) - analytic(cctk_time, p.x, p.y, p.z, poly_a,
                                          poly_b, poly_c, poly_d);
  });
}

extern "C" void TestSubcyclingProlong_Sync(CCTK_ARGUMENTS) {
  // do nothing
}
