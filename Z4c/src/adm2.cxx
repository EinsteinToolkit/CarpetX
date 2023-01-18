#include <cctk.h>

#ifdef __CUDACC__
// Disable CCTK_DEBUG since the debug information takes too much
// parameter space to launch the kernels
#ifdef CCTK_DEBUG
#undef CCTK_DEBUG
#endif
#endif

#include "derivs.hxx"
#include "physics.hxx"
#include "z4c_vars.hxx"

#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef __CUDACC__
#include <nvToolsExt.h>
#endif

#include <cmath>

namespace Z4c {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Z4c_ADM2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < deriv_order / 2 + 1)
      CCTK_VERROR("Need at least %d ghost zones", deriv_order / 2 + 1);

  //

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout1(cctkGH, indextype);
  const GF3D5layout layout0(imin, imax);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gammat1{
      GF3D2<const CCTK_REAL>(layout1, gammatxx),
      GF3D2<const CCTK_REAL>(layout1, gammatxy),
      GF3D2<const CCTK_REAL>(layout1, gammatxz),
      GF3D2<const CCTK_REAL>(layout1, gammatyy),
      GF3D2<const CCTK_REAL>(layout1, gammatyz),
      GF3D2<const CCTK_REAL>(layout1, gammatzz)};

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_At1{
      GF3D2<const CCTK_REAL>(layout1, Atxx),
      GF3D2<const CCTK_REAL>(layout1, Atxy),
      GF3D2<const CCTK_REAL>(layout1, Atxz),
      GF3D2<const CCTK_REAL>(layout1, Atyy),
      GF3D2<const CCTK_REAL>(layout1, Atyz),
      GF3D2<const CCTK_REAL>(layout1, Atzz)};

  const vec<GF3D2<const CCTK_REAL>, 3> gf_Gamt1{
      GF3D2<const CCTK_REAL>(layout1, Gamtx),
      GF3D2<const CCTK_REAL>(layout1, Gamty),
      GF3D2<const CCTK_REAL>(layout1, Gamtz)};

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_betaG1{
      GF3D2<const CCTK_REAL>(layout1, betaGx),
      GF3D2<const CCTK_REAL>(layout1, betaGy),
      GF3D2<const CCTK_REAL>(layout1, betaGz)};

  //

  const int ntmps = 154;
  GF3D5vector<CCTK_REAL> tmps(layout0, ntmps);
  int itmp = 0;

  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(tmps(itmp++)); };
  const auto make_vec = [&](const auto &f) {
    return vec<result_of_t<decltype(f)()>, 3>([&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return smat<result_of_t<decltype(f)()>, 3>([&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };
  const auto make_vec_vec_gf = [&]() { return make_vec(make_vec_gf); };
  const auto make_vec_mat_gf = [&]() { return make_vec(make_mat_gf); };
  const auto make_mat_vec_gf = [&]() { return make_mat(make_vec_gf); };
  const auto make_mat_mat_gf = [&]() { return make_mat(make_mat_gf); };

  const GF3D5<CCTK_REAL> gf_chi0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dchi0(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddchi0(make_mat_gf());
  calc_derivs2(cctkGH, gf_chi1, gf_chi0, gf_dchi0, gf_ddchi0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_gammat0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dgammat0(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddgammat0(make_mat_mat_gf());
  calc_derivs2(cctkGH, gf_gammat1, gf_gammat0, gf_dgammat0, gf_ddgammat0,
               layout0);

  const GF3D5<CCTK_REAL> gf_Kh0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dKh0(make_vec_gf());
  calc_derivs(cctkGH, gf_Kh1, gf_Kh0, gf_dKh0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_At0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dAt0(make_mat_vec_gf());
  calc_derivs(cctkGH, gf_At1, gf_At0, gf_dAt0, layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_Gamt0(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dGamt0(make_vec_vec_gf());
  calc_derivs(cctkGH, gf_Gamt1, gf_Gamt0, gf_dGamt0, layout0);

  const GF3D5<CCTK_REAL> gf_Theta0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dTheta0(make_vec_gf());
  calc_derivs(cctkGH, gf_Theta1, gf_Theta0, gf_dTheta0, layout0);

  const GF3D5<CCTK_REAL> gf_alphaG0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dalphaG0(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddalphaG0(make_mat_gf());
  calc_derivs2(cctkGH, gf_alphaG1, gf_alphaG0, gf_dalphaG0, gf_ddalphaG0,
               layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_betaG0(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dbetaG0(make_vec_vec_gf());
  const vec<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddbetaG0(make_vec_mat_gf());
  calc_derivs2(cctkGH, gf_betaG1, gf_betaG0, gf_dbetaG0, gf_ddbetaG0, layout0);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  //

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_eTti1{
      GF3D2<const CCTK_REAL>(layout1, eTtx),
      GF3D2<const CCTK_REAL>(layout1, eTty),
      GF3D2<const CCTK_REAL>(layout1, eTtz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_eTij1{
      GF3D2<const CCTK_REAL>(layout1, eTxx),
      GF3D2<const CCTK_REAL>(layout1, eTxy),
      GF3D2<const CCTK_REAL>(layout1, eTxz),
      GF3D2<const CCTK_REAL>(layout1, eTyy),
      GF3D2<const CCTK_REAL>(layout1, eTyz),
      GF3D2<const CCTK_REAL>(layout1, eTzz)};

  //

  const smat<GF3D2<CCTK_REAL>, 3> gf_dtk1{
      GF3D2<CCTK_REAL>(layout1, dtkxx), GF3D2<CCTK_REAL>(layout1, dtkxy),
      GF3D2<CCTK_REAL>(layout1, dtkxz), GF3D2<CCTK_REAL>(layout1, dtkyy),
      GF3D2<CCTK_REAL>(layout1, dtkyz), GF3D2<CCTK_REAL>(layout1, dtkzz)};

  const GF3D2<CCTK_REAL> gf_dt2alp1(layout1, dt2alp);

  const vec<GF3D2<CCTK_REAL>, 3> gf_dt2beta1{
      GF3D2<CCTK_REAL>(layout1, dt2betax), GF3D2<CCTK_REAL>(layout1, dt2betay),
      GF3D2<CCTK_REAL>(layout1, dt2betaz)};

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);
#ifdef __CUDACC__
  const nvtxRangeId_t range = nvtxRangeStartA("Z4c_ADM2::adm2");
#endif
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index1(layout1, p.I);
        const GF3D5index index0(layout0, p.I);

        // load and calculate
        const z4c_vars<vreal> vars(
            set_Theta_zero, kappa1, kappa2, f_mu_L, f_mu_S, eta, //
            gf_chi0(mask, index0, 1), gf_dchi0(mask, index0),
            gf_ddchi0(mask, index0), //
            gf_gammat0(mask, index0, one<smat<int, 3> >()()),
            gf_dgammat0(mask, index0), gf_ddgammat0(mask, index0), //
            gf_Kh0(mask, index0), gf_dKh0(mask, index0),           //
            gf_At0(mask, index0), gf_dAt0(mask, index0),           //
            gf_Gamt0(mask, index0), gf_dGamt0(mask, index0),       //
            gf_Theta0(mask, index0), gf_dTheta0(mask, index0),     //
            gf_alphaG0(mask, index0), gf_dalphaG0(mask, index0),
            gf_ddalphaG0(mask, index0), //
            gf_betaG0(mask, index0), gf_dbetaG0(mask, index0),
            gf_ddbetaG0(mask, index0), //
            gf_eTtt1(mask, index1), gf_eTti1(mask, index1),
            gf_eTij1(mask, index1));

        // Store
        gf_dtk1.store(mask, index1, vars.K_rhs);
        gf_dt2alp1.store(mask, index1, vars.dtalpha_rhs);
        gf_dt2beta1.store(mask, index1, vars.dtbeta_rhs);
      });
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace Z4c
