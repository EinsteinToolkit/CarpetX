#include <defs.hxx>
#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace DGWaveToyX {
using namespace Arith;

extern "C" void DGWaveToyX_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGWaveToyX_Init;
  DECLARE_CCTK_PARAMETERS;

  using std::sqrt;
  const CCTK_REAL kx = 2 * M_PI;
  const CCTK_REAL ky = 2 * M_PI;
  const CCTK_REAL kz = 2 * M_PI;
  const CCTK_REAL omega = sqrt(pown(kx, 2) + pown(ky, 2) + pown(kz, 2));

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL t = cctk_time;
    const CCTK_REAL x = coordx(p.I);
    const CCTK_REAL y = coordy(p.I);
    const CCTK_REAL z = coordz(p.I);

    using std::cos, std::sin;
    u(p.I) = cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);

    ft(p.I) = -omega * sin(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
    fx(p.I) = -kx * cos(omega * t) * sin(kx * x) * cos(ky * y) * cos(kz * z);
    fy(p.I) = -ky * cos(omega * t) * cos(kx * x) * sin(ky * y) * cos(kz * z);
    fz(p.I) = -kz * cos(omega * t) * cos(kx * x) * cos(ky * y) * sin(kz * z);
  });
}

} // namespace DGWaveToyX
