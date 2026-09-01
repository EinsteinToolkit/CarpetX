#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cstdlib>

namespace TestTimelevels {

bool found_error = false;

// Compare one time level of test_gf against an expected constant value
// everywhere in the local box (interior, boundary and ghosts), warning and
// recording a failure otherwise.
void check(const Loop::GridDescBaseDevice &grid, const int cctk_iteration,
          const Loop::GF3D2<const CCTK_REAL> &gf, const CCTK_REAL expected,
          const char *restrict const which) {
  int error_count = 0;
  grid.loop_all<0, 0, 0>(grid.nghostzones,
                        [&] CCTK_HOST(const Loop::PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              if (gf(p.I) != expected) {
                                ++error_count;
                              }
                            });

  if (error_count > 0) {
    CCTK_VWARN(CCTK_WARN_ALERT,
              "TestTimelevels: %s holds the wrong value in %d points at "
              "iteration %d (expected %g)",
              which, error_count, cctk_iteration, double(expected));
    found_error = true;
  }
}

extern "C" void TestTimelevels_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_Init;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        test_gf(p.I) = 0;
        test_gf_p(p.I) = 0;
        test_gf_p_p(p.I) = 0;
      });
}

extern "C" void TestTimelevels_Write(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_Write;

  const CCTK_REAL value = cctk_iteration;
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        test_gf(p.I) = value;
      });
}

extern "C" void TestTimelevels_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_Check;

  // Time levels 1 and 2 were seeded to zero at t = 0, and only ever change by
  // being rotated forward from time level 0, so their expected value clamps
  // at zero rather than going negative.
  check(grid, cctk_iteration, test_gf, CCTK_REAL(cctk_iteration),
       "time level 0");
  check(grid, cctk_iteration, test_gf_p,
       CCTK_REAL(std::max(cctk_iteration - 1, 0)), "time level 1");
  check(grid, cctk_iteration, test_gf_p_p,
       CCTK_REAL(std::max(cctk_iteration - 2, 0)), "time level 2");
}

extern "C" void TestTimelevels_Terminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Terminate;

  if (found_error) {
    CCTK_Exit(cctkGH, EXIT_FAILURE);
  }
}

} // namespace TestTimelevels
