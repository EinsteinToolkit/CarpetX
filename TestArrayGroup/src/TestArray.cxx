#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <sstream>

using namespace std::literals::complex_literals;

static bool found_error = false;

extern "C" void CCTK_FNAME(TestArrayGroup_FoundError)(void) {
  found_error = true;
}

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
        test1[index] = CCTK_REAL(1 + (i + 1) * (j + 1) * (n + 1));
        test2[index] = CCTK_REAL(1 + 7 * (i + 1) * (j + 1) * (n + 1));
        test3[index] = CCTK_REAL(1 + 13 * (i + 1) * (j + 1) * (n + 1));
        test_int1[index] = CCTK_INT(2 + (i + 1) * (j + 1) * (n + 1));
        test_int2[index] = CCTK_INT(2 + 7 * (i + 1) * (j + 1) * (n + 1));
        test_int3[index] = CCTK_INT(2 + 13 * (i + 1) * (j + 1) * (n + 1));
        test_complex1[index] =
            CCTK_COMPLEX(3 + (i + 1) * (j + 1) * (n + 1)) + CCTK_COMPLEX(1i);
        test_complex2[index] =
            CCTK_COMPLEX(3 + 7 * (i + 1) * (j + 1) * (n + 1)) +
            CCTK_COMPLEX(1i);
        test_complex3[index] =
            CCTK_COMPLEX(3 + 13 * (i + 1) * (j + 1) * (n + 1)) +
            CCTK_COMPLEX(1i);
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
        test_gf[index] = (i0 + i + 1) * (j0 + j + 1) * (k0 + k + 1);
      }
    }
  }

  *test_scalar = CCTK_REAL(1);
  *test_scalar_int = CCTK_INT(2);
  *test_scalar_complex = CCTK_COMPLEX(3);
}

extern "C" void TestArrayGroup_Compare(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS_TestArrayGroup_Compare;

  // test array group
  const int imax = 5;
  const int jmax = 6;
  const int nmax = 4;

  const char *arrays[] = {"test", "test_int", "test_complex"};
  const int iscale[] = {1, 7, 13};
  for (char const *const array : arrays) {
    for (unsigned int suffix = 0; suffix < 3; ++suffix) {
      int array_error_count = 0;
      bool is_contiguous = true;

      std::ostringstream vname0;
      vname0 << "TestArrayGroup::" << array << (suffix + 1) << "[" << 0 << "]";
      void *const varptr0 = CCTK_VarDataPtr(cctkGH, 0, vname0.str().c_str());

      for (int n = 0; n < nmax; n++) {
        std::ostringstream vname;
        vname << "TestArrayGroup::" << array << (suffix + 1) << "[" << n << "]";
        void *const varptr = CCTK_VarDataPtr(cctkGH, 0, vname.str().c_str());
        const int vartype = CCTK_VarTypeI(CCTK_VarIndex(vname.str().c_str()));

        const int size = jmax * imax;
        is_contiguous = ((char *)varptr - (char *)varptr0) /
                            (CCTK_VarTypeSize(vartype) * size) ==
                        n;
        for (int j = 0; j < jmax; j++) {
          for (int i = 0; i < imax; i++) {
            const int index = i + j * imax;
            int ival = iscale[suffix] * (i + 1) * (j + 1) * (n + 1);
            bool is_equal = [&]() {
              switch (vartype) {
              case CCTK_VARIABLE_REAL:
                return static_cast<CCTK_REAL *>(varptr)[index] ==
                       CCTK_REAL(1 + ival);
              case CCTK_VARIABLE_INT:
                return static_cast<CCTK_INT *>(varptr)[index] ==
                       CCTK_INT(2 + ival);
              case CCTK_VARIABLE_COMPLEX:
                return static_cast<CCTK_COMPLEX *>(varptr)[index] ==
                       CCTK_COMPLEX(3 + ival) + CCTK_COMPLEX(1i);
              default:
                assert(0 && "Unexpected variable type");
                return false; // notreached
              }
            }();
            if (!is_equal)
              array_error_count += 1;
          }
        }

        if (array_error_count > 0) {
          CCTK_VWARN(
              CCTK_WARN_ALERT,
              "TestArrayGroup: grid array %s failed in %d of %d elements",
              vname.str().c_str(), array_error_count, size);
          found_error = true;
        }
        if (!is_contiguous) {
          CCTK_VWARN(CCTK_WARN_ALERT,
                     "TestArrayGroup: grid array %s is not contiguous in the "
                     "grid array vector",
                     vname.str().c_str());
          found_error = true;
        }
      }
    }
  }

  // test grid function
  int gf_gi = CCTK_GroupIndex("TestArrayGroup::test_gf");
  cGroup gf_group;
  int ierr = CCTK_GroupData(gf_gi, &gf_group);
  assert(!ierr);
  std::array<int, 3> lsh;
  ierr = CCTK_GrouplshGI(cctkGH, gf_group.dim, lsh.data(), gf_gi);
  assert(!ierr);

  const int i0 = cctk_lbnd[0];
  const int j0 = cctk_lbnd[1];
  const int k0 = cctk_lbnd[2];

  int gf_error_count = 0;

  for (int k = 0; k < lsh[2]; k++) {
    for (int j = 0; j < lsh[1]; j++) {
      for (int i = 0; i < lsh[0]; i++) {
        const int index = CCTK_GFINDEX3D(cctkGH, i, j, k);
        if (test_gf[index] != (i0 + i + 1) * (j0 + j + 1) * (k0 + k + 1))
          gf_error_count += 1;
      }
    }
  }

  if (gf_error_count > 0) {
    const int size = lsh[2] * lsh[1] * lsh[0];
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestArrayGroup: grid function test failed in %d of %d elements",
               gf_error_count, size);
    found_error = true;
  }

  // test grid scalar
  if (*test_scalar != CCTK_REAL(1)) {
    CCTK_VWARN(CCTK_WARN_ALERT, "TestArrayGroup: CCTK_REAL grid scalar failed");
    found_error = true;
  }

  if (*test_scalar_int != CCTK_INT(2)) {
    CCTK_VWARN(CCTK_WARN_ALERT, "TestArrayGroup: CCTK_INT grid scalar failed");
    found_error = true;
  }

  if (*test_scalar_complex != CCTK_COMPLEX(3)) {
    CCTK_VWARN(CCTK_WARN_ALERT,
               "TestArrayGroup: CCTK_COMPLEX grid scalar failed");
    found_error = true;
  }
}

extern "C" void TestArrayGroup_Terminate(CCTK_ARGUMENTS) {
  // this should be caugt by output files, but since we test the array code
  // itself, report an error explicitly
  if (found_error)
    CCTK_Exit(cctkGH, EXIT_FAILURE);
}
