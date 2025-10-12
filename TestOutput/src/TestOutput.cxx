#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include <complex>

namespace TestOutput {

using namespace std::literals;

extern "C" void TestOutput_SetVarsLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_SetVarsLocal;

  grid.loop_int<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    gf(p.I) = p.z * 10000 + p.y * 100 + p.x;
  });
}

extern "C" void TestOutput_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_Sync;
}

extern "C" void TestOutput_SetVarsGlobal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_SetVarsGlobal;

  *sc = CCTK_REAL(42);
  *sc_int = CCTK_INT(17);
  *sc_complex = CCTK_COMPLEX(1 + M_PI) * CCTK_COMPLEX(1i);

  for (int i = 0; i < 10; ++i) {
    a1[i] = CCTK_REAL(i);
    a1_int[i] = i + CCTK_INT(1);
    a1_complex[i] = CCTK_COMPLEX(i + 2) + CCTK_COMPLEX(1i);
  }

  for (int j = 0; j < 8; ++j)
    for (int i = 0; i < 9; ++i) {
      a2[j * 9 + i] = CCTK_REAL(100 * j + i);
      a2_int[j * 9 + i] = CCTK_INT(100 * j + i + 1);
      a2_complex[j * 9 + i] = CCTK_COMPLEX(100 * j + i + 2) + CCTK_COMPLEX(1i);
    }

  for (int k = 0; k < 5; ++k)
    for (int j = 0; j < 6; ++j)
      for (int i = 0; i < 7; ++i) {
        a3[(k * 6 + j) * 7 + i] = CCTK_REAL(10000 * k + 100 * j + i);
        a3_int[(k * 6 + j) * 7 + i] = CCTK_INT(10000 * k + 100 * j + i + 1);
        a3_complex[(k * 6 + j) * 7 + i] =
            CCTK_COMPLEX(10000 * k + 100 * j + i + 2) + CCTK_COMPLEX(1i);
      }
}

extern "C" void TestOutput_CopyVarsLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_CopyVarsLocal;

  grid.loop_all<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    gf_tmp(p.I) = gf(p.I);
  });
}

extern "C" void TestOutput_UpdateVarsLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_UpdateVarsLocal;

  grid.loop_int<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    CCTK_REAL lap = 0;
    for (int d = 0; d < 3; ++d)
      lap += gf_tmp(p.I - p.DI[d]) - 2 * gf_tmp(p.I) + gf_tmp(p.I + p.DI[d]);
    gf(p.I) = gf_tmp(p.I) + lap / 4;
  });
}

extern "C" void TestOutput_UpdateVarsGlobal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_UpdateVarsGlobal;

  *sc += CCTK_REAL(1);
  *sc_int += CCTK_INT(1);
  *sc_complex += CCTK_COMPLEX(1) + CCTK_COMPLEX(1i);

  for (int i = 0; i < 10; ++i) {
    a1[i] += CCTK_REAL(1);
    a1_int[i] += CCTK_INT(1);
    a1_complex[i] += CCTK_COMPLEX(1) + CCTK_COMPLEX(1i);
  }

  for (int j = 0; j < 8; ++j)
    for (int i = 0; i < 9; ++i) {
      a2[j * 9 + i] += CCTK_REAL(1);
      a2_int[j * 9 + i] += CCTK_INT(1);
      a2_complex[j * 9 + i] += CCTK_COMPLEX(1) + CCTK_COMPLEX(1i);
    }

  for (int k = 0; k < 5; ++k)
    for (int j = 0; j < 6; ++j)
      for (int i = 0; i < 7; ++i) {
        a3[(k * 6 + j) * 7 + i] += CCTK_REAL(1);
        a3_int[(k * 6 + j) * 7 + i] += CCTK_INT(1);
        a3_complex[(k * 6 + j) * 7 + i] += CCTK_COMPLEX(1) + CCTK_COMPLEX(1i);
      }
}

} // namespace TestOutput
