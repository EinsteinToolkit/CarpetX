#include "physics.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Initial1;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g1{
      GF3D2<const CCTK_REAL>(layout1, gxx),
      GF3D2<const CCTK_REAL>(layout1, gxy),
      GF3D2<const CCTK_REAL>(layout1, gxz),
      GF3D2<const CCTK_REAL>(layout1, gyy),
      GF3D2<const CCTK_REAL>(layout1, gyz),
      GF3D2<const CCTK_REAL>(layout1, gzz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_K1{
      GF3D2<const CCTK_REAL>(layout1, kxx),
      GF3D2<const CCTK_REAL>(layout1, kxy),
      GF3D2<const CCTK_REAL>(layout1, kxz),
      GF3D2<const CCTK_REAL>(layout1, kyy),
      GF3D2<const CCTK_REAL>(layout1, kyz),
      GF3D2<const CCTK_REAL>(layout1, kzz)};

  const GF3D2<const CCTK_REAL> gf_alp1(layout1, alp);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta1{
      GF3D2<const CCTK_REAL>(layout1, betax),
      GF3D2<const CCTK_REAL>(layout1, betay),
      GF3D2<const CCTK_REAL>(layout1, betaz)};

  const GF3D2<CCTK_REAL> gf_chi1(layout1, chi);

  const smat<GF3D2<CCTK_REAL>, 3> gf_gammat1{
      GF3D2<CCTK_REAL>(layout1, gammatxx), GF3D2<CCTK_REAL>(layout1, gammatxy),
      GF3D2<CCTK_REAL>(layout1, gammatxz), GF3D2<CCTK_REAL>(layout1, gammatyy),
      GF3D2<CCTK_REAL>(layout1, gammatyz), GF3D2<CCTK_REAL>(layout1, gammatzz)};

  const GF3D2<CCTK_REAL> gf_Kh1(layout1, Kh);

  const smat<GF3D2<CCTK_REAL>, 3> gf_At1{
      GF3D2<CCTK_REAL>(layout1, Atxx), GF3D2<CCTK_REAL>(layout1, Atxy),
      GF3D2<CCTK_REAL>(layout1, Atxz), GF3D2<CCTK_REAL>(layout1, Atyy),
      GF3D2<CCTK_REAL>(layout1, Atyz), GF3D2<CCTK_REAL>(layout1, Atzz)};

  const GF3D2<CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec<GF3D2<CCTK_REAL>, 3> gf_betaG1{GF3D2<CCTK_REAL>(layout1, betaGx),
                                           GF3D2<CCTK_REAL>(layout1, betaGy),
                                           GF3D2<CCTK_REAL>(layout1, betaGz)};

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_Initial1::initial1");
#endif
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);

        // Load
        const smat<vreal, 3> g = gf_g1(mask, index1, one<smat<int, 3> >()());
        const smat<vreal, 3> K = gf_K1(mask, index1);
        const vreal alp = gf_alp1(mask, index1, 1);
        const vec<vreal, 3> beta = gf_beta1(mask, index1);

        // Calculate Z4c variables (all except Gammat)
        const vreal detg = calc_det(g);
        const smat<vreal, 3> gu = calc_inv(g, detg);

        const vreal chi = 1 / cbrt(detg);

        const smat<vreal, 3> gammat([&](int a, int b)
                                        ARITH_INLINE { return chi * g(a, b); });

        const vreal trK = sum_symm<3>(
            [&](int x, int y) ARITH_INLINE { return gu(x, y) * K(x, y); });

        const vreal Theta = 0;

        const vreal Kh = trK - 2 * Theta;

        const smat<vreal, 3> At([&](int a, int b) ARITH_INLINE {
          return chi * (K(a, b) - trK / 3 * g(a, b));
        });

        const vreal alphaG = alp - 1;

        const vec<vreal, 3> betaG([&](int a) ARITH_INLINE { return beta(a); });

        // Store
        gf_chi1.store(mask, index1, chi);
        gf_gammat1.store(mask, index1, gammat);
        gf_Kh1.store(mask, index1, Kh);
        gf_At1.store(mask, index1, At);
        gf_Theta1.store(mask, index1, Theta);
        gf_alphaG1.store(mask, index1, alphaG);
        gf_betaG1.store(mask, index1, betaG);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4c
