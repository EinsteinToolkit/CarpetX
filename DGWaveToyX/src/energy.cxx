#include <defs.hxx>
#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace DGWaveToyX {
using namespace Arith;

extern "C" void DGWaveToyX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGWaveToyX_Energy;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    eps(p.I) = (pown(ft(p.I), 2) + pown(fx(p.I), 2) + pown(fy(p.I), 2) +
                pown(fz(p.I), 2));
  });
}

} // namespace DGWaveToyX
