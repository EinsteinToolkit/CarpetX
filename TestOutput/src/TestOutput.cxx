#include <loop_device.hxx>

namespace TestOutput {

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

  *sc = 42;

  for (int i = 0; i < 10; ++i)
    a1[i] = i;

  for (int j = 0; j < 8; ++j)
    for (int i = 0; i < 9; ++i)
      a2[j * 9 + i] = 100 * j + i;

  for (int k = 0; k < 5; ++k)
    for (int j = 0; j < 6; ++j)
      for (int i = 0; i < 7; ++i)
        a3[(k * 5 + j) * 7 + i] = 10000 * k + 100 * j + i;
}

extern "C" void TestOutput_UpdateVarsLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_UpdateVarsLocal;

  grid.loop_int<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    CCTK_REAL lap = 0;
    for (int d = 0; d < 3; ++d)
      lap += gf(p.I - p.DI[d]) - 2 * gf(p.I) + gf(p.I + p.DI[d]);
    gf(p.I) += lap / 4;
  });
}

extern "C" void TestOutput_UpdateVarsGlobal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_UpdateVarsGlobal;

  *sc += 1;

  for (int i = 0; i < 10; ++i)
    a1[i] += 1;

  for (int j = 0; j < 8; ++j)
    for (int i = 0; i < 9; ++i)
      a2[j * 9 + i] += 1;

  for (int k = 0; k < 5; ++k)
    for (int j = 0; j < 6; ++j)
      for (int i = 0; i < 7; ++i)
        a3[(k * 5 + j) * 7 + i] += 1;
}

} // namespace TestOutput
