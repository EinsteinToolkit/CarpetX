// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>

namespace TestTimelevels {

bool found_error = false;

// resize_gf and swap_gf are addressed by group index, not through the
// variables the schedule clauses name, because the driver's time level
// storage API works on whole groups.
int resize_gi() {
  const int gi = CCTK_GroupIndex("TestTimelevels::resize_gf");
  assert(gi >= 0);
  return gi;
}
int resize_vi() {
  const int vi = CCTK_VarIndex("TestTimelevels::resize_gf");
  assert(vi >= 0);
  return vi;
}

int swap_gi() {
  const int gi = CCTK_GroupIndex("TestTimelevels::swap_gf");
  assert(gi >= 0);
  return gi;
}

// The mark SwapInit writes into swap_gf time level tl. Offset away from the
// iteration numbers test_gf and resize_gf carry, so that reading the wrong
// buffer cannot accidentally produce the expected value.
CCTK_REAL swap_mark(const int tl) { return 10 + tl; }

// How many time levels resize_gf should have at this iteration
int expected_resize_timelevels(const int cctk_iteration) {
  DECLARE_CCTK_PARAMETERS;
  return resize_iteration >= 0 && cctk_iteration >= resize_iteration
             ? final_resize_timelevels
             : initial_resize_timelevels;
}

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

extern "C" void TestTimelevels_SetupStorage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_SetupStorage;
  DECLARE_CCTK_PARAMETERS;

  // WRAGH runs before any level exists, so this only records the count that
  // the GroupData constructor will then use for every level it builds
  const int old_ntls =
      CarpetX::SetGroupTimelevels(resize_gi(), initial_resize_timelevels);
  CCTK_VINFO("TestTimelevels: resize_gf allocated with %d time levels "
             "(was %d, %d declared)",
             int(initial_resize_timelevels), old_ntls,
             CCTK_DeclaredTimeLevelsGI(resize_gi()));
}

extern "C" void TestTimelevels_ResizeWrite(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_ResizeWrite;

  const CCTK_REAL value = cctk_iteration;
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        resize_gf(p.I) = value;
      });
}

extern "C" void TestTimelevels_Resize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Resize;
  DECLARE_CCTK_PARAMETERS;

  if (resize_iteration < 0 || cctk_iteration != resize_iteration)
    return;

  const int old_ntls =
      CarpetX::SetGroupTimelevels(resize_gi(), final_resize_timelevels);
  CCTK_VINFO("TestTimelevels: resize_gf reallocated from %d to %d time levels "
             "at iteration %d",
             old_ntls, int(final_resize_timelevels), cctk_iteration);
}

extern "C" void TestTimelevels_ResizeCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_ResizeCheck;

  const int gi = resize_gi();
  const int vi = resize_vi();
  const int ntls = CarpetX::GetGroupTimelevels(gi);
  const int expected_ntls = expected_resize_timelevels(cctk_iteration);

  if (ntls != expected_ntls) {
    CCTK_VWARN(CCTK_WARN_ALERT,
              "TestTimelevels: resize_gf has %d time levels allocated at "
              "iteration %d, expected %d",
              ntls, cctk_iteration, expected_ntls);
    found_error = true;
  }

  // Exactly the allocated time levels must be reachable through
  // cctkGH->data, which is sized from the interface.ccl declaration and so
  // has slots for the unallocated ones too
  for (int tl = 0; tl < CCTK_DeclaredTimeLevelsGI(gi); ++tl) {
    const void *restrict const ptr = CCTK_VarDataPtrI(cctkGH, tl, vi);
    if ((ptr == nullptr) != (tl >= ntls)) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                "TestTimelevels: resize_gf time level %d has a %s data "
                "pointer at iteration %d, with %d time levels allocated",
                tl, ptr ? "non-null" : "null", cctk_iteration, ntls);
      found_error = true;
    }
  }

  // Reallocating must leave time level 0 alone
  check(grid, cctk_iteration, resize_gf, CCTK_REAL(cctk_iteration),
       "resize_gf time level 0");
}

extern "C" void TestTimelevels_SwapInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_SwapInit;

  const CCTK_REAL mark0 = swap_mark(0);
  const CCTK_REAL mark1 = swap_mark(1);
  const CCTK_REAL mark2 = swap_mark(2);
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        swap_gf(p.I) = mark0;
        swap_gf_p(p.I) = mark1;
        swap_gf_p_p(p.I) = mark2;
      });
}

extern "C" void TestTimelevels_Swap(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Swap;

  CarpetX::SwapGroupTimelevels(swap_gi(), 0, 2);
}

extern "C" void TestTimelevels_SwapCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_SwapCheck;

  // TestTimelevels_Swap exchanged slots 0 and 2, so their marks have traded
  // places and slot 1 is untouched. Reading these at all also exercises the
  // driver's own validity check: the exchange has to carry the validity flags
  // along with the buffers, or CallFunction rejects this routine's input.
  check(grid, cctk_iteration, swap_gf, swap_mark(2), "swap_gf time level 0");
  check(grid, cctk_iteration, swap_gf_p, swap_mark(1), "swap_gf time level 1");
  check(grid, cctk_iteration, swap_gf_p_p, swap_mark(0),
       "swap_gf time level 2");
}

extern "C" void TestTimelevels_Terminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Terminate;

  if (found_error) {
    CCTK_Exit(cctkGH, EXIT_FAILURE);
  }
}

} // namespace TestTimelevels
