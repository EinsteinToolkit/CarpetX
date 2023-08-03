#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <sstream>

extern "C" void TestArrayGroup_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_Initialize;

  const int imax = 5;
  const int jmax = 6;
  const int nmax = 4;

  // Initialize test array
  for (int n = 0; n < nmax; n++) {
    for (int j = 0; j < jmax; j++) {
      for (int i = 0; i < imax; i++) {
        const int index = i + j * imax + n * imax * jmax;
        test1[index] = 1 + i * j * n;
        test2[index] = 1 + 7 * i * j * n;
        test3[index] = 1 + 13 * i * j * n;
      }
    }
  }

  // Initialize test grid function
  int gi = CCTK_GroupIndex("TestArrayGroup::test_gf");
  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  std::array<int, 3> lsh;
  ierr = CCTK_GrouplshGI(cctkGH, group.dim, lsh.data(), gi);
  assert(!ierr);

  const int i0 = cctk_lbnd[0];
  const int j0 = cctk_lbnd[1];
  const int k0 = cctk_lbnd[2];
  for (int k = 0; k < lsh[2]; k++) {
    for (int j = 0; j < lsh[1]; j++) {
      for (int i = 0; i < lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        test_gf[index] = (i0 + i) * (j0 + j) * (k0 + k);
      }
    }
  }

  *test_scalar = 1;
}

extern "C" void TestArrayGroup_Compare(cGH *cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_Compare;

  const int imax = 5;
  const int jmax = 6;
  const int nmax = 4;

  int error_count[3] = {0, 0, 0};

  for (int n = 0; n < nmax; n++) {
    std::ostringstream vname1;
    vname1 << "TestArrayGroup::test1[" << n << "]";
    CCTK_REAL *const var1 =
        (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, vname1.str().c_str());

    std::ostringstream vname2;
    vname2 << "TestArrayGroup::test2[" << n << "]";
    CCTK_REAL *const var2 =
        (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, vname2.str().c_str());

    std::ostringstream vname3;
    vname3 << "TestArrayGroup::test3[" << n << "]";
    CCTK_REAL *const var3 =
        (CCTK_REAL *)CCTK_VarDataPtr(cctkGH, 0, vname3.str().c_str());

    for (int j = 0; j < jmax; j++) {
      for (int i = 0; i < imax; i++) {
        const int index = i + j * imax;
        if (var1[index] != 1 + i * j * n)
          error_count[0] += 1;
        if (var2[index] != 1 + 7 * i * j * n)
          error_count[1] += 1;
        if (var3[index] != 1 + 13 * i * j * n)
          error_count[2] += 1;
      }
    }
  }

  if (error_count[0] > 0) {
    const int size = nmax * jmax * imax;
    CCTK_VERROR("TestArrayGroup: grid array test1 failed in %d of %d elements",
                error_count[0], size);
  }

  if (error_count[1] > 0) {
    const int size = nmax * jmax * imax;
    CCTK_VERROR("TestArrayGroup: grid array test2 failed in %d of %d elements",
                error_count[1], size);
  }

  if (error_count[2] > 0) {
    const int size = nmax * jmax * imax;
    CCTK_VERROR("TestArrayGroup: grid array test3 failed in %d of %d elements",
                error_count[2], size);
  }
}
