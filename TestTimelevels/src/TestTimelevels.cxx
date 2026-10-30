// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdlib>

namespace TestTimelevels {

// Atomic because the local mode checks run once per tile, possibly on several
// threads at once
std::atomic<int> num_failed_checks{0};

// The groups below are addressed by group index or by name, not through the
// variables the schedule clauses name, because the group storage API works on
// whole groups.
int group_index(const char *const groupname) {
  const int gi = CCTK_GroupIndex(groupname);
  assert(gi >= 0);
  return gi;
}

const char *const resize_group = "TestTimelevels::resize_gf";
const char *const swap_group = "TestTimelevels::swap_gf";
const char *const storage_array_group = "TestTimelevels::storage_array";
const char *const storage_scalar_group = "TestTimelevels::storage_scalar";
const char *const cycle_array_group = "TestTimelevels::cycle_array";
const char *const cycle_scalar_group = "TestTimelevels::cycle_scalar";

// The grid array and grid scalar groups the storage cycle resizes, alongside
// resize_gf
const char *const resized_array_groups[] = {storage_array_group,
                                            storage_scalar_group};

// The grid array and grid scalar groups kept at their declared time levels,
// whose history has to survive every rotation
const char *const cycled_array_groups[] = {cycle_array_group,
                                           cycle_scalar_group};

// The mark SwapInit writes into swap_gf time level tl. Offset away from the
// iteration numbers test_gf and resize_gf carry, so that reading the wrong
// buffer cannot accidentally produce the expected value.
CCTK_REAL swap_mark(const int tl) { return 10 + tl; }

// How many time levels the resized groups should have at this iteration
int expected_resize_timelevels(const int cctk_iteration) {
  DECLARE_CCTK_PARAMETERS;
  if (regrow_iteration >= 0 && cctk_iteration >= regrow_iteration) {
    return regrow_resize_timelevels;
  }
  if (resize_iteration >= 0 && cctk_iteration >= resize_iteration) {
    return final_resize_timelevels;
  }
  return initial_resize_timelevels;
}

// Drive one group to exactly ntls time levels through the Cactus storage API.
// An increase never shrinks and a decrease never grows, so both calls are
// needed to move in either direction; the increase reports the count that was
// allocated before either of them acted.
int set_timelevels(const cGH *restrict const cctkGH, const int gi,
                   const int ntls) {
  int old_ntls = -1;
  CCTK_GroupStorageIncrease(cctkGH, 1, &gi, &ntls, &old_ntls);
  CCTK_GroupStorageDecrease(cctkGH, 1, &gi, &ntls, nullptr);
  return old_ntls;
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
    ++num_failed_checks;
  }
}

// Number of elements in one time level of one variable of a grid array or grid
// scalar group. Arrays and scalars have no ghost zones, so lsh is the whole
// extent, and a scalar has no dimensions at all.
int array_size(const cGH *restrict const cctkGH, const int gi) {
  cGroupDynamicData data;
  const int ierr = CCTK_GroupDynamicData(cctkGH, gi, &data);
  assert(not ierr);
  int size = 1;
  for (int d = 0; d < data.dim; ++d) {
    size *= data.lsh[d];
  }
  return size;
}

// Fill one time level of a grid array or grid scalar group with a constant
void write_array(const cGH *restrict const cctkGH, const char *const groupname,
                 const int tl, const CCTK_REAL value) {
  const int gi = group_index(groupname);
  CCTK_REAL *restrict const ptr = static_cast<CCTK_REAL *>(
      CCTK_VarDataPtrI(cctkGH, tl, CCTK_FirstVarIndexI(gi)));
  assert(ptr);
  const int size = array_size(cctkGH, gi);
  for (int i = 0; i < size; ++i) {
    ptr[i] = value;
  }
}

// The value one time level of a grid array or grid scalar group holds. Every
// element holds the same value, which check_array verifies.
CCTK_REAL read_array(const cGH *restrict const cctkGH,
                     const char *const groupname, const int tl) {
  const CCTK_REAL *restrict const ptr =
      static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
          cctkGH, tl, CCTK_FirstVarIndexI(group_index(groupname))));
  assert(ptr);
  return ptr[0];
}

// Compare one time level of a grid array or grid scalar group against an
// expected constant value, warning and recording a failure otherwise
void check_array(const cGH *restrict const cctkGH, const char *const groupname,
                 const int tl, const CCTK_REAL expected,
                 const char *restrict const where) {
  const int gi = group_index(groupname);
  const CCTK_REAL *restrict const ptr = static_cast<const CCTK_REAL *>(
      CCTK_VarDataPtrI(cctkGH, tl, CCTK_FirstVarIndexI(gi)));
  assert(ptr);
  const int size = array_size(cctkGH, gi);

  int error_count = 0;
  for (int i = 0; i < size; ++i) {
    if (ptr[i] != expected) {
      ++error_count;
    }
  }

  if (error_count > 0) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestTimelevels: %s time level %d holds the wrong value in %d "
               "points at iteration %d in %s (expected %g)",
               groupname, tl, error_count, cctkGH->cctk_iteration, where,
               double(expected));
    ++num_failed_checks;
  }
}

// Exactly the allocated time levels must be reachable through cctkGH->data,
// which is sized from the interface.ccl declaration and so has slots for the
// unallocated ones too
void check_data_pointers(const cGH *restrict const cctkGH,
                         const char *const groupname, const int ntls,
                         const char *restrict const where) {
  const int gi = group_index(groupname);
  const int vi = CCTK_FirstVarIndexI(gi);
  for (int tl = 0; tl < CCTK_DeclaredTimeLevelsGI(gi); ++tl) {
    const void *restrict const ptr = CCTK_VarDataPtrI(cctkGH, tl, vi);
    if ((ptr == nullptr) != (tl >= ntls)) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "TestTimelevels: %s time level %d has a %s data pointer at "
                 "iteration %d in %s, with %d time levels allocated",
                 groupname, tl, ptr ? "non-null" : "null",
                 cctkGH->cctk_iteration, where, ntls);
      ++num_failed_checks;
    }
  }
}

// Every way of asking how many time levels a group has must agree, and a valid
// group must always report that it has storage: CarpetX never releases one
// entirely.
void check_storage(const cGH *restrict const cctkGH,
                   const char *const groupname, const int ntls) {
  const int gi = group_index(groupname);
  const int declared_ntls = CCTK_DeclaredTimeLevelsGI(gi);

  const auto wrong = [&](const char *const what, const int got,
                         const int expected) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestTimelevels: %s reports %s = %d at iteration %d, "
               "expected %d",
               groupname, what, got, cctkGH->cctk_iteration, expected);
    ++num_failed_checks;
  };

  const int active = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
  if (active != ntls) {
    wrong("CCTK_ActiveTimeLevelsGI", active, ntls);
  }

  cGroupDynamicData data;
  const int ierr = CCTK_GroupDynamicData(cctkGH, gi, &data);
  assert(not ierr);
  if (data.activetimelevels != ntls) {
    wrong("cGroupDynamicData::activetimelevels", data.activetimelevels, ntls);
  }

  const int has_storage_by_index = CCTK_QueryGroupStorageI(cctkGH, gi);
  if (has_storage_by_index != 1) {
    wrong("CCTK_QueryGroupStorageI", has_storage_by_index, 1);
  }
  const int has_storage_by_name = CCTK_QueryGroupStorage(cctkGH, groupname);
  if (has_storage_by_name != 1) {
    wrong("CCTK_QueryGroupStorage", has_storage_by_name, 1);
  }

  const int max_ntls = CCTK_MaxActiveTimeLevelsGI(cctkGH, gi);
  if (max_ntls != declared_ntls) {
    wrong("CCTK_MaxActiveTimeLevelsGI", max_ntls, declared_ntls);
  }

  // An increase of zero is a pure query: it reports the active count without
  // changing it
  const int no_increase = 0;
  int status = -1;
  const int queried =
      CCTK_GroupStorageIncrease(cctkGH, 1, &gi, &no_increase, &status);
  if (queried != ntls) {
    wrong("CCTK_GroupStorageIncrease by 0", queried, ntls);
  }
  if (status != ntls) {
    wrong("CCTK_GroupStorageIncrease by 0 status", status, ntls);
  }

  // A decrease never grows a group, so asking for the declared count back is a
  // no-op whenever fewer time levels are allocated
  status = -1;
  const int unchanged =
      CCTK_GroupStorageDecrease(cctkGH, 1, &gi, &declared_ntls, &status);
  if (unchanged != ntls) {
    wrong("CCTK_GroupStorageDecrease to the declared count", unchanged, ntls);
  }
  if (status != ntls) {
    wrong("CCTK_GroupStorageDecrease to the declared count status", status,
          ntls);
  }
  const int still_active = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
  if (still_active != ntls) {
    wrong("CCTK_ActiveTimeLevelsGI after a decrease above the active count",
          still_active, ntls);
  }

  // Disabling storage floors a group at one time level instead of releasing
  // it, because time level 0 holds live data. Guarded on the count the group
  // actually has, not the expected one: a decrease from more than one time
  // level would shrink the group and hide a wrong count from every later
  // check.
  if (active == 1) {
    const int previous = CCTK_DisableGroupStorage(cctkGH, groupname);
    if (previous != 1) {
      wrong("CCTK_DisableGroupStorage", previous, 1);
    }
    const int floored = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
    if (floored != 1) {
      wrong("CCTK_ActiveTimeLevelsGI after CCTK_DisableGroupStorage", floored,
            1);
    }
  }
}

// The grid array and grid scalar data and data pointers, as seen through one
// cGH. The cycled groups carry their whole history; the resized groups lose
// theirs to the reallocation and only keep time level 0.
void check_arrays(const cGH *restrict const cctkGH,
                  const char *restrict const where) {
  const int cctk_iteration = cctkGH->cctk_iteration;

  for (const char *const groupname : cycled_array_groups) {
    const int gi = group_index(groupname);
    const int ntls = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
    check_data_pointers(cctkGH, groupname, ntls, where);
    for (int tl = 0; tl < ntls; ++tl) {
      check_array(cctkGH, groupname, tl, CCTK_REAL(cctk_iteration - tl), where);
    }
  }

  const int ntls = expected_resize_timelevels(cctk_iteration);
  for (const char *const groupname : resized_array_groups) {
    check_data_pointers(cctkGH, groupname, ntls, where);
    check_array(cctkGH, groupname, 0, CCTK_REAL(cctk_iteration), where);
  }
}

extern "C" void TestTimelevels_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_Init;

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
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
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { test_gf(p.I) = value; });
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

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      observed_test_tl1(p.I) = test_gf_p(p.I);
                                      observed_test_tl2(p.I) = test_gf_p_p(p.I);
                                    });
}

extern "C" void TestTimelevels_SetupStorage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_SetupStorage;
  DECLARE_CCTK_PARAMETERS;

  // WRAGH runs before any level exists, so this only records the count that
  // the GroupData constructor will then use for every level it builds
  const int gi = group_index(resize_group);
  const int old_ntls = set_timelevels(cctkGH, gi, initial_resize_timelevels);
  CCTK_VINFO("TestTimelevels: resize_gf allocated with %d time levels "
             "(was %d, %d declared)",
             int(initial_resize_timelevels), old_ntls,
             CCTK_DeclaredTimeLevelsGI(gi));
}

extern "C" void TestTimelevels_ResizeWrite(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_ResizeWrite;

  const CCTK_REAL value = cctk_iteration;
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { resize_gf(p.I) = value; });
}

extern "C" void TestTimelevels_Resize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Resize;

  const int gi = group_index(resize_group);
  const int ntls = expected_resize_timelevels(cctk_iteration);
  if (CCTK_ActiveTimeLevelsGI(cctkGH, gi) == ntls) {
    return;
  }

  const int old_ntls = set_timelevels(cctkGH, gi, ntls);
  CCTK_VINFO("TestTimelevels: resize_gf reallocated from %d to %d time levels "
             "at iteration %d",
             old_ntls, ntls, cctk_iteration);
}

extern "C" void TestTimelevels_ResizeCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_ResizeCheck;

  const int gi = group_index(resize_group);
  const int ntls = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
  const int expected_ntls = expected_resize_timelevels(cctk_iteration);

  if (ntls != expected_ntls) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestTimelevels: resize_gf has %d time levels allocated at "
               "iteration %d, expected %d",
               ntls, cctk_iteration, expected_ntls);
    ++num_failed_checks;
  }

  check_data_pointers(cctkGH, resize_group, ntls, "local mode");

  // Reallocating must leave time level 0 alone
  check(grid, cctk_iteration, resize_gf, CCTK_REAL(cctk_iteration),
        "resize_gf time level 0");
}

extern "C" void TestTimelevels_SwapInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestTimelevels_SwapInit;

  const CCTK_REAL mark0 = swap_mark(0);
  const CCTK_REAL mark1 = swap_mark(1);
  const CCTK_REAL mark2 = swap_mark(2);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      swap_gf(p.I) = mark0;
                                      swap_gf_p(p.I) = mark1;
                                      swap_gf_p_p(p.I) = mark2;
                                    });
}

extern "C" void TestTimelevels_Swap(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Swap;

  CarpetX::SwapGroupTimelevels(group_index(swap_group), 0, 2);
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

  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      observed_swap_tl1(p.I) = swap_gf_p(p.I);
                                      observed_swap_tl2(p.I) = swap_gf_p_p(p.I);
                                    });
}

extern "C" void TestTimelevels_ArraySetupStorage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArraySetupStorage;
  DECLARE_CCTK_PARAMETERS;

  for (const char *const groupname : resized_array_groups) {
    const int gi = group_index(groupname);
    const int old_ntls = set_timelevels(cctkGH, gi, initial_resize_timelevels);
    CCTK_VINFO("TestTimelevels: %s allocated with %d time levels "
               "(was %d, %d declared)",
               groupname, int(initial_resize_timelevels), old_ntls,
               CCTK_DeclaredTimeLevelsGI(gi));
  }
}

extern "C" void TestTimelevels_ArrayInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArrayInit;

  // Time level tl holds cctk_iteration - tl. Seeding every time level here
  // makes that invariant hold from the first iteration on, since each rotation
  // moves every value one level back.
  for (const char *const groupname : cycled_array_groups) {
    const int gi = group_index(groupname);
    for (int tl = 0; tl < CCTK_ActiveTimeLevelsGI(cctkGH, gi); ++tl) {
      write_array(cctkGH, groupname, tl, CCTK_REAL(cctk_iteration - tl));
    }
  }

  for (const char *const groupname : resized_array_groups) {
    write_array(cctkGH, groupname, 0, CCTK_REAL(cctk_iteration));
  }
}

extern "C" void TestTimelevels_ArrayWrite(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArrayWrite;

  for (const char *const groupname : cycled_array_groups) {
    write_array(cctkGH, groupname, 0, CCTK_REAL(cctk_iteration));
  }
  for (const char *const groupname : resized_array_groups) {
    write_array(cctkGH, groupname, 0, CCTK_REAL(cctk_iteration));
  }
}

extern "C" void TestTimelevels_ArrayResize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArrayResize;

  const int ntls = expected_resize_timelevels(cctk_iteration);
  for (const char *const groupname : resized_array_groups) {
    const int gi = group_index(groupname);
    if (CCTK_ActiveTimeLevelsGI(cctkGH, gi) == ntls) {
      continue;
    }

    const int old_ntls = set_timelevels(cctkGH, gi, ntls);
    CCTK_VINFO("TestTimelevels: %s reallocated from %d to %d time levels at "
               "iteration %d",
               groupname, old_ntls, ntls, cctk_iteration);
  }
}

extern "C" void TestTimelevels_StorageCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_StorageCheck;

  const int ntls = expected_resize_timelevels(cctk_iteration);
  check_storage(cctkGH, resize_group, ntls);
  for (const char *const groupname : resized_array_groups) {
    check_storage(cctkGH, groupname, ntls);
  }
  for (const char *const groupname : cycled_array_groups) {
    check_storage(cctkGH, groupname,
                  CCTK_DeclaredTimeLevelsGI(group_index(groupname)));
  }
}

extern "C" void TestTimelevels_ArrayCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArrayCheck;

  check_arrays(cctkGH, "global mode");

  *cycle_array_tl0 = read_array(cctkGH, cycle_array_group, 0);
  *cycle_array_tl1 = read_array(cctkGH, cycle_array_group, 1);
  *cycle_array_tl2 = read_array(cctkGH, cycle_array_group, 2);
  *cycle_scalar_tl0 = read_array(cctkGH, cycle_scalar_group, 0);
  *cycle_scalar_tl1 = read_array(cctkGH, cycle_scalar_group, 1);
  *cycle_scalar_tl2 = read_array(cctkGH, cycle_scalar_group, 2);
  *storage_array_tl0 = read_array(cctkGH, storage_array_group, 0);
  *storage_scalar_tl0 = read_array(cctkGH, storage_scalar_group, 0);
}

extern "C" void TestTimelevels_ArrayCheckLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_ArrayCheckLocal;

  // The local mode cGHs are deep copies, so they only hold the right pointers
  // if every cached cGH was re-pointed when the storage changed
  check_arrays(cctkGH, "local mode");
}

extern "C" void TestTimelevels_Record(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Record;

  *resize_timelevels =
      CCTK_ActiveTimeLevelsGI(cctkGH, group_index(resize_group));
  *cycle_array_timelevels =
      CCTK_ActiveTimeLevelsGI(cctkGH, group_index(cycle_array_group));
  *cycle_scalar_timelevels =
      CCTK_ActiveTimeLevelsGI(cctkGH, group_index(cycle_scalar_group));
  *storage_array_timelevels =
      CCTK_ActiveTimeLevelsGI(cctkGH, group_index(storage_array_group));
  *storage_scalar_timelevels =
      CCTK_ActiveTimeLevelsGI(cctkGH, group_index(storage_scalar_group));
  *num_errors = num_failed_checks;
}

extern "C" void TestTimelevels_Terminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestTimelevels_Terminate;

  if (num_failed_checks > 0) {
    CCTK_Exit(cctkGH, EXIT_FAILURE);
  }
}

} // namespace TestTimelevels
