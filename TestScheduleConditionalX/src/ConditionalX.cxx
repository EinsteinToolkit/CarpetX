// Conditional READS:/WRITES: clauses and DECLARE_CCTK_ARGUMENTSX_<func>.
//
// Two things only this macro can check, which is why this thorn lives in
// CarpetX rather than in cactustest:
//
//   - the GF3D2 views exist at all for a variable named by a single branch.
//     Which branch is taken is decided at startup, long after this file is
//     compiled, so the macro has to declare the union of every branch.
//
//   - const propagates correctly into the view. rdwr.pl marks a variable const
//     only when no branch writes it, and passes that through
//     CCTK_CENTERING_GF, so a variable written by any branch must come out as
//     a mutable view. The static_asserts below are the whole check; they are
//     evaluated at compile time and do not depend on the parameters.

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdlib>
#include <type_traits>

#include <util_String.h>

extern "C" void TestScheduleConditionalX_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_TestScheduleConditionalX_Init;

  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        condVarA(p.I) = 1.0;
        condVarB(p.I) = 2.0;
        condVarC(p.I) = 3.0;
      });
}

extern "C" void TestScheduleConditionalX_Check(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_TestScheduleConditionalX_Check;
  DECLARE_CCTK_PARAMETERS;

  // condVarA is READ by every branch and written by none, so it must be a
  // const view. condVarB and condVarC are each written by one branch, so they
  // must be mutable even when that branch is not the one taken.
  static_assert(std::is_same_v<std::remove_const_t<decltype(condVarA)>,
                               Loop::GF3D2<const CCTK_REAL> >,
                "condVarA is only READ, so its view must be const");
  static_assert(std::is_same_v<std::remove_const_t<decltype(condVarB)>,
                               Loop::GF3D2<CCTK_REAL> >,
                "condVarB is WRITTEN by a branch, so its view must be mutable");
  static_assert(std::is_same_v<std::remove_const_t<decltype(condVarC)>,
                               Loop::GF3D2<CCTK_REAL> >,
                "condVarC is WRITTEN by a branch, so its view must be mutable");

  // Run-time half: a variable belonging to a branch that was not taken is a
  // NULL pointer under presync_mode = "presync-only".
  *VarAAccessX = cctk_ptr_condVarA != nullptr;
  *VarBAccessX = cctk_ptr_condVarB != nullptr;
  *VarCAccessX = cctk_ptr_condVarC != nullptr;

  // Honour the WRITES clauses of whichever branch is active. Each says
  // (everywhere), so every point has to be written.
  if (cond_a) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          condVarB(p.I) = 4.0;
        });
  }
  if (CCTK_EQUALS(cond_mode, "one")) {
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          condVarC(p.I) = 5.0;
        });
  }
}

extern "C" void TestScheduleConditionalX_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_TestScheduleConditionalX_Output;
  DECLARE_CCTK_PARAMETERS;

  if (!*out_dir || CCTK_MyProc(cctkGH) != 0)
    return;

  char *filename;
  FILE *fh;
  if (Util_asprintf(&filename, "%s/BranchAccess.asc", out_dir) < 0 ||
      CCTK_CreateDirectory(0777, out_dir) < 0 ||
      (fh = fopen(filename, "w")) == nullptr)
    CCTK_VERROR("Could not open file '%s' for output", filename);

  fprintf(fh, "# Branch accessibility: 1 = variable accessible, 0 = NULL\n");
  fprintf(fh, "# Columns: condVarA condVarB condVarC\n");
  fprintf(fh, "%d %d %d\n", (int)*VarAAccessX, (int)*VarBAccessX,
          (int)*VarCAccessX);

  fclose(fh);
  free(filename);
}
