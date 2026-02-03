#include <dg.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace DGWaveToyX {
using namespace DG;

using Loop::dim;

extern "C" void DGWaveToyX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGWaveToyX_RHS;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
    u_rhs(p.I) = ft(p.I);
  });

  switch (dg_npoints) {

  case 2: {
    constexpr int N = 2;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz, p);

      ft_rhs(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs(p.I) = dft[0];
      fy_rhs(p.I) = dft[1];
      fz_rhs(p.I) = dft[2];
    });
    break;
  }

  case 4: {
    constexpr int N = 4;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz, p);

      ft_rhs(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs(p.I) = dft[0];
      fy_rhs(p.I) = dft[1];
      fz_rhs(p.I) = dft[2];
    });
    break;
  }

  case 8: {
    constexpr int N = 8;
    grid.loop_int<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const vect<CCTK_REAL, dim> dft = get_derivs<N>(grid, ft, p);
      const vect<CCTK_REAL, dim> dfx = get_derivs<N>(grid, fx, p);
      const vect<CCTK_REAL, dim> dfy = get_derivs<N>(grid, fy, p);
      const vect<CCTK_REAL, dim> dfz = get_derivs<N>(grid, fz, p);

      ft_rhs(p.I) = dfx[0] + dfy[1] + dfz[2];
      fx_rhs(p.I) = dft[0];
      fy_rhs(p.I) = dft[1];
      fz_rhs(p.I) = dft[2];
    });
    break;
  }

  default:
    CCTK_ERROR("Unsupported DG element size");
  }
}

} // namespace DGWaveToyX
