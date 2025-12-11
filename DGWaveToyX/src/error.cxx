#include <defs.hxx>
#include <dg.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace DGWaveToyX {
using namespace DG;

extern "C" void DGWaveToyX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGWaveToyX_Error;
  DECLARE_CCTK_PARAMETERS;

  using std::sqrt;
  const CCTK_REAL kx = 2 * M_PI;
  const CCTK_REAL ky = 2 * M_PI;
  const CCTK_REAL kz = 2 * M_PI;
  const CCTK_REAL omega = sqrt(pown(kx, 2) + pown(ky, 2) + pown(kz, 2));

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL t = cctk_time;
    const CCTK_REAL x = coordx(p.I);
    const CCTK_REAL y = coordy(p.I);
    const CCTK_REAL z = coordz(p.I);

    using std::cos, std::sin;
    u_error(p.I) =
        u(p.I) - cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);

    ft_error(p.I) = ft(p.I) - -omega * sin(omega * t) * cos(kx * x) *
                                  cos(ky * y) * cos(kz * z);
    fx_error(p.I) = fx(p.I) - -kx * cos(omega * t) * sin(kx * x) * cos(ky * y) *
                                  cos(kz * z);
    fy_error(p.I) = fy(p.I) - -ky * cos(omega * t) * cos(kx * x) * sin(ky * y) *
                                  cos(kz * z);
    fz_error(p.I) = fz(p.I) - -kz * cos(omega * t) * cos(kx * x) * cos(ky * y) *
                                  sin(kz * z);
  });
}

} // namespace DGWaveToyX
