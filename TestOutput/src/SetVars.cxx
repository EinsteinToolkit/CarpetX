#include <loop_device.hxx>

namespace TestOutput {

extern "C" void TestOutput_SetVarsLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestOutput_SetVarsLocal;

  grid.loop_all<0, 0, 0>(grid.nghostzones, [=](const Loop::PointDesc &p) {
    gf(p.I) = p.k * 10000 + p.j * 100 + p.i;
  });
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

} // namespace TestOutput
