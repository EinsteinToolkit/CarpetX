#include <fixmath.hxx>
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
#include "weyl_vars.hxx"

#include <defs.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <rten.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Weyl {
using namespace Arith;
using namespace Loop;
using namespace std;

extern "C" void Weyl_Weyl(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Weyl_Weyl;
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
  const GF3D5layout layout5(cctkGH, indextype);

  const smat<GF3D2<const CCTK_REAL>, 3> gf_gamma1{
      GF3D2<const CCTK_REAL>(layout1, gxx),
      GF3D2<const CCTK_REAL>(layout1, gxy),
      GF3D2<const CCTK_REAL>(layout1, gxz),
      GF3D2<const CCTK_REAL>(layout1, gyy),
      GF3D2<const CCTK_REAL>(layout1, gyz),
      GF3D2<const CCTK_REAL>(layout1, gzz)};

  const GF3D2<const CCTK_REAL> gf_alpha1(layout1, alp);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta1{
      GF3D2<const CCTK_REAL>(layout1, betax),
      GF3D2<const CCTK_REAL>(layout1, betay),
      GF3D2<const CCTK_REAL>(layout1, betaz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_k1{
      GF3D2<const CCTK_REAL>(layout1, kxx),
      GF3D2<const CCTK_REAL>(layout1, kxy),
      GF3D2<const CCTK_REAL>(layout1, kxz),
      GF3D2<const CCTK_REAL>(layout1, kyy),
      GF3D2<const CCTK_REAL>(layout1, kyz),
      GF3D2<const CCTK_REAL>(layout1, kzz)};

  const GF3D2<const CCTK_REAL> gf_dtalpha1(layout1, dtalp);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_dtbeta1{
      GF3D2<const CCTK_REAL>(layout1, dtbetax),
      GF3D2<const CCTK_REAL>(layout1, dtbetay),
      GF3D2<const CCTK_REAL>(layout1, dtbetaz)};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_dtk1{
      GF3D2<const CCTK_REAL>(layout1, dtkxx),
      GF3D2<const CCTK_REAL>(layout1, dtkxy),
      GF3D2<const CCTK_REAL>(layout1, dtkxz),
      GF3D2<const CCTK_REAL>(layout1, dtkyy),
      GF3D2<const CCTK_REAL>(layout1, dtkyz),
      GF3D2<const CCTK_REAL>(layout1, dtkzz)};

  const GF3D2<const CCTK_REAL> gf_dt2alpha1(layout1, dt2alp);

  const vec<GF3D2<const CCTK_REAL>, 3> gf_dt2beta1{
      GF3D2<const CCTK_REAL>(layout1, dt2betax),
      GF3D2<const CCTK_REAL>(layout1, dt2betay),
      GF3D2<const CCTK_REAL>(layout1, dt2betaz)};

  //

  const int ntmps = 146;
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

  const smat<GF3D5<CCTK_REAL>, 3> gf_gamma0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dgamma0(make_mat_vec_gf());
  const smat<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddgamma0(make_mat_mat_gf());
  calc_derivs2(cctkGH, gf_gamma1, gf_gamma0, gf_dgamma0, gf_ddgamma0, layout0);

  const GF3D5<CCTK_REAL> gf_alpha0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_dalpha0(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, 3> gf_ddalpha0(make_mat_gf());
  calc_derivs2(cctkGH, gf_alpha1, gf_alpha0, gf_dalpha0, gf_ddalpha0, layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_beta0(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dbeta0(make_vec_vec_gf());
  const vec<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddbeta0(make_vec_mat_gf());
  calc_derivs2(cctkGH, gf_beta1, gf_beta0, gf_dbeta0, gf_ddbeta0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_k0(make_mat_gf());
  const smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dk0(make_mat_vec_gf());
  calc_derivs(cctkGH, gf_k1, gf_k0, gf_dk0, layout0);

  const GF3D5<CCTK_REAL> gf_dtalpha0(make_gf());
  const vec<GF3D5<CCTK_REAL>, 3> gf_ddtalpha0(make_vec_gf());
  calc_derivs(cctkGH, gf_dtalpha1, gf_dtalpha0, gf_ddtalpha0, layout0);

  const vec<GF3D5<CCTK_REAL>, 3> gf_dtbeta0(make_vec_gf());
  const vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_ddtbeta0(make_vec_vec_gf());
  calc_derivs(cctkGH, gf_dtbeta1, gf_dtbeta0, gf_ddtbeta0, layout0);

  const smat<GF3D5<CCTK_REAL>, 3> gf_dtk0(make_mat_gf());
  calc_copy(cctkGH, gf_dtk1, gf_dtk0, layout0);

  if (itmp != ntmps)
    CCTK_VERROR("Wrong number of temporary variables: ntmps=%d itmp=%d", ntmps,
                itmp);
  itmp = -1;

  //

  const smat<GF3D2<CCTK_REAL>, 4> gf_g41{
      GF3D2<CCTK_REAL>(layout1, g4tt), GF3D2<CCTK_REAL>(layout1, g4tx),
      GF3D2<CCTK_REAL>(layout1, g4ty), GF3D2<CCTK_REAL>(layout1, g4tz),
      GF3D2<CCTK_REAL>(layout1, g4xx), GF3D2<CCTK_REAL>(layout1, g4xy),
      GF3D2<CCTK_REAL>(layout1, g4xz), GF3D2<CCTK_REAL>(layout1, g4yy),
      GF3D2<CCTK_REAL>(layout1, g4yz), GF3D2<CCTK_REAL>(layout1, g4zz)};

  // const GF3D2<CCTK_REAL> gf_Gamma4ttt1(layout1, Gamma4ttt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ttx1(layout1, Gamma4ttx);
  // const GF3D2<CCTK_REAL> gf_Gamma4tty1(layout1, Gamma4tty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ttz1(layout1, Gamma4ttz);
  // const GF3D2<CCTK_REAL> gf_Gamma4txx1(layout1, Gamma4txx);
  // const GF3D2<CCTK_REAL> gf_Gamma4txy1(layout1, Gamma4txy);
  // const GF3D2<CCTK_REAL> gf_Gamma4txz1(layout1, Gamma4txz);
  // const GF3D2<CCTK_REAL> gf_Gamma4tyy1(layout1, Gamma4tyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4tyz1(layout1, Gamma4tyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4tzz1(layout1, Gamma4tzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4xtt1(layout1, Gamma4xtt);
  // const GF3D2<CCTK_REAL> gf_Gamma4xtx1(layout1, Gamma4xtx);
  // const GF3D2<CCTK_REAL> gf_Gamma4xty1(layout1, Gamma4xty);
  // const GF3D2<CCTK_REAL> gf_Gamma4xtz1(layout1, Gamma4xtz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxx1(layout1, Gamma4xxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxy1(layout1, Gamma4xxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4xxz1(layout1, Gamma4xxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xyy1(layout1, Gamma4xyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4xyz1(layout1, Gamma4xyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4xzz1(layout1, Gamma4xzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4ytt1(layout1, Gamma4ytt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ytx1(layout1, Gamma4ytx);
  // const GF3D2<CCTK_REAL> gf_Gamma4yty1(layout1, Gamma4yty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ytz1(layout1, Gamma4ytz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxx1(layout1, Gamma4yxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxy1(layout1, Gamma4yxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4yxz1(layout1, Gamma4yxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yyy1(layout1, Gamma4yyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4yyz1(layout1, Gamma4yyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4yzz1(layout1, Gamma4yzz);

  // const GF3D2<CCTK_REAL> gf_Gamma4ztt1(layout1, Gamma4ztt);
  // const GF3D2<CCTK_REAL> gf_Gamma4ztx1(layout1, Gamma4ztx);
  // const GF3D2<CCTK_REAL> gf_Gamma4zty1(layout1, Gamma4zty);
  // const GF3D2<CCTK_REAL> gf_Gamma4ztz1(layout1, Gamma4ztz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxx1(layout1, Gamma4zxx);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxy1(layout1, Gamma4zxy);
  // const GF3D2<CCTK_REAL> gf_Gamma4zxz1(layout1, Gamma4zxz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zyy1(layout1, Gamma4zyy);
  // const GF3D2<CCTK_REAL> gf_Gamma4zyz1(layout1, Gamma4zyz);
  // const GF3D2<CCTK_REAL> gf_Gamma4zzz1(layout1, Gamma4zzz);

  // const GF3D2<CCTK_REAL> gf_rm4txtx1(layout1, rm4txtx);
  // const GF3D2<CCTK_REAL> gf_rm4txty1(layout1, rm4txty);
  // const GF3D2<CCTK_REAL> gf_rm4txtz1(layout1, rm4txtz);
  // const GF3D2<CCTK_REAL> gf_rm4txxy1(layout1, rm4txxy);
  // const GF3D2<CCTK_REAL> gf_rm4txxz1(layout1, rm4txxz);
  // const GF3D2<CCTK_REAL> gf_rm4txyz1(layout1, rm4txyz);

  // const GF3D2<CCTK_REAL> gf_rm4tyty1(layout1, rm4tyty);
  // const GF3D2<CCTK_REAL> gf_rm4tytz1(layout1, rm4tytz);
  // const GF3D2<CCTK_REAL> gf_rm4tyxy1(layout1, rm4tyxy);
  // const GF3D2<CCTK_REAL> gf_rm4tyxz1(layout1, rm4tyxz);
  // const GF3D2<CCTK_REAL> gf_rm4tyyz1(layout1, rm4tyyz);

  // const GF3D2<CCTK_REAL> gf_rm4tztz1(layout1, rm4tztz);
  // const GF3D2<CCTK_REAL> gf_rm4tzxy1(layout1, rm4tzxy);
  // const GF3D2<CCTK_REAL> gf_rm4tzxz1(layout1, rm4tzxz);
  // const GF3D2<CCTK_REAL> gf_rm4tzyz1(layout1, rm4tzyz);

  // const GF3D2<CCTK_REAL> gf_rm4xyxy1(layout1, rm4xyxy);
  // const GF3D2<CCTK_REAL> gf_rm4xyxz1(layout1, rm4xyxz);
  // const GF3D2<CCTK_REAL> gf_rm4xyyz1(layout1, rm4xyyz);

  // const GF3D2<CCTK_REAL> gf_rm4xzxz1(layout1, rm4xzxz);
  // const GF3D2<CCTK_REAL> gf_rm4xzyz1(layout1, rm4xzyz);

  // const GF3D2<CCTK_REAL> gf_rm4yzyz1(layout1, rm4yzyz);

  // const GF3D2<CCTK_REAL> gf_r4tt1(layout1, r4tt);
  // const GF3D2<CCTK_REAL> gf_r4tx1(layout1, r4tx);
  // const GF3D2<CCTK_REAL> gf_r4ty1(layout1, r4ty);
  // const GF3D2<CCTK_REAL> gf_r4tz1(layout1, r4tz);
  // const GF3D2<CCTK_REAL> gf_r4xx1(layout1, r4xx);
  // const GF3D2<CCTK_REAL> gf_r4xy1(layout1, r4xy);
  // const GF3D2<CCTK_REAL> gf_r4xz1(layout1, r4xz);
  // const GF3D2<CCTK_REAL> gf_r4yy1(layout1, r4yy);
  // const GF3D2<CCTK_REAL> gf_r4yz1(layout1, r4yz);
  // const GF3D2<CCTK_REAL> gf_r4zz1(layout1, r4zz);

  // const GF3D2<CCTK_REAL> gf_rsc41(layout1, rsc4);

  // const GF3D2<CCTK_REAL> gf_c4txtx1(layout1, c4txtx);
  // const GF3D2<CCTK_REAL> gf_c4txty1(layout1, c4txty);
  // const GF3D2<CCTK_REAL> gf_c4txtz1(layout1, c4txtz);
  // const GF3D2<CCTK_REAL> gf_c4txxy1(layout1, c4txxy);
  // const GF3D2<CCTK_REAL> gf_c4txxz1(layout1, c4txxz);
  // const GF3D2<CCTK_REAL> gf_c4txyz1(layout1, c4txyz);

  // const GF3D2<CCTK_REAL> gf_c4tyty1(layout1, c4tyty);
  // const GF3D2<CCTK_REAL> gf_c4tytz1(layout1, c4tytz);
  // const GF3D2<CCTK_REAL> gf_c4tyxy1(layout1, c4tyxy);
  // const GF3D2<CCTK_REAL> gf_c4tyxz1(layout1, c4tyxz);
  // const GF3D2<CCTK_REAL> gf_c4tyyz1(layout1, c4tyyz);

  // const GF3D2<CCTK_REAL> gf_c4tztz1(layout1, c4tztz);
  // const GF3D2<CCTK_REAL> gf_c4tzxy1(layout1, c4tzxy);
  // const GF3D2<CCTK_REAL> gf_c4tzxz1(layout1, c4tzxz);
  // const GF3D2<CCTK_REAL> gf_c4tzyz1(layout1, c4tzyz);

  // const GF3D2<CCTK_REAL> gf_c4xyxy1(layout1, c4xyxy);
  // const GF3D2<CCTK_REAL> gf_c4xyxz1(layout1, c4xyxz);
  // const GF3D2<CCTK_REAL> gf_c4xyyz1(layout1, c4xyyz);

  // const GF3D2<CCTK_REAL> gf_c4xzxz1(layout1, c4xzxz);
  // const GF3D2<CCTK_REAL> gf_c4xzyz1(layout1, c4xzyz);

  // const GF3D2<CCTK_REAL> gf_c4yzyz1(layout1, c4yzyz);

  //

  // const vec<GF3D2<CCTK_REAL>, 4> gf_l1{
  //     GF3D2<CCTK_REAL>(layout1, lt), GF3D2<CCTK_REAL>(layout1, lx),
  //     GF3D2<CCTK_REAL>(layout1, ly), GF3D2<CCTK_REAL>(layout1, lz)};

  // const vec<GF3D2<CCTK_REAL>, 4> gf_n1{
  //     GF3D2<CCTK_REAL>(layout1, nt), GF3D2<CCTK_REAL>(layout1, nx),
  //     GF3D2<CCTK_REAL>(layout1, ny), GF3D2<CCTK_REAL>(layout1, nz)};

  // const vec<GF3D2<CCTK_REAL>, 4> gf_mre1{
  //     GF3D2<CCTK_REAL>(layout1, mret), GF3D2<CCTK_REAL>(layout1, mrex),
  //     GF3D2<CCTK_REAL>(layout1, mrey), GF3D2<CCTK_REAL>(layout1, mrez)};

  // const vec<GF3D2<CCTK_REAL>, 4> gf_mim1{
  //     GF3D2<CCTK_REAL>(layout1, mimt), GF3D2<CCTK_REAL>(layout1, mimx),
  //     GF3D2<CCTK_REAL>(layout1, mimy), GF3D2<CCTK_REAL>(layout1, mimz)};

  //

  // const GF3D2<CCTK_REAL> gf_Lambda1(layout1, Lambda);
  // const GF3D2<CCTK_REAL> gf_Phi00(layout1, Phi00);
  // const GF3D2<CCTK_REAL> gf_Phi111(layout1, Phi11);
  // const GF3D2<CCTK_REAL> gf_Phi221(layout1, Phi22);
  // const GF3D2<CCTK_REAL> gf_Phi10re1(layout1, Phi10re);
  // const GF3D2<CCTK_REAL> gf_Phi10im1(layout1, Phi10im);
  // const GF3D2<CCTK_REAL> gf_Phi20re1(layout1, Phi20re);
  // const GF3D2<CCTK_REAL> gf_Phi20im1(layout1, Phi20im);
  // const GF3D2<CCTK_REAL> gf_Phi21re1(layout1, Phi21re);
  // const GF3D2<CCTK_REAL> gf_Phi21im1(layout1, Phi21im);

  // //

  const GF3D2<CCTK_REAL> gf_Psi0re1(layout1, Psi0re);
  const GF3D2<CCTK_REAL> gf_Psi0im1(layout1, Psi0im);
  const GF3D2<CCTK_REAL> gf_Psi1re1(layout1, Psi1re);
  const GF3D2<CCTK_REAL> gf_Psi1im1(layout1, Psi1im);
  const GF3D2<CCTK_REAL> gf_Psi2re1(layout1, Psi2re);
  const GF3D2<CCTK_REAL> gf_Psi2im1(layout1, Psi2im);
  const GF3D2<CCTK_REAL> gf_Psi3re1(layout1, Psi3re);
  const GF3D2<CCTK_REAL> gf_Psi3im1(layout1, Psi3im);
  const GF3D2<CCTK_REAL> gf_Psi4re1(layout1, Psi4re);
  const GF3D2<CCTK_REAL> gf_Psi4im1(layout1, Psi4im);

  // //

  // const GF3D2<CCTK_REAL> gf_npkappare1(layout1, npkappare);
  // const GF3D2<CCTK_REAL> gf_npkappaim1(layout1, npkappaim);
  // const GF3D2<CCTK_REAL> gf_npsigmare1(layout1, npsigmare);
  // const GF3D2<CCTK_REAL> gf_npsigmaim1(layout1, npsigmaim);
  // const GF3D2<CCTK_REAL> gf_nprhore1(layout1, nprhore);
  // const GF3D2<CCTK_REAL> gf_nprhoim1(layout1, nprhoim);
  // const GF3D2<CCTK_REAL> gf_nptaure1(layout1, nptaure);
  // const GF3D2<CCTK_REAL> gf_nptauim1(layout1, nptauim);
  // const GF3D2<CCTK_REAL> gf_npepsilonre1(layout1, npepsilonre);
  // const GF3D2<CCTK_REAL> gf_npepsilonim1(layout1, npepsilonim);
  // const GF3D2<CCTK_REAL> gf_npbetare1(layout1, npbetare);
  // const GF3D2<CCTK_REAL> gf_npbetaim1(layout1, npbetaim);
  // const GF3D2<CCTK_REAL> gf_npalphare1(layout1, npalphare);
  // const GF3D2<CCTK_REAL> gf_npalphaim1(layout1, npalphaim);
  // const GF3D2<CCTK_REAL> gf_npgammare1(layout1, npgammare);
  // const GF3D2<CCTK_REAL> gf_npgammaim1(layout1, npgammaim);
  // const GF3D2<CCTK_REAL> gf_nppire1(layout1, nppire);
  // const GF3D2<CCTK_REAL> gf_nppiim1(layout1, nppiim);
  // const GF3D2<CCTK_REAL> gf_npmure1(layout1, npmure);
  // const GF3D2<CCTK_REAL> gf_npmuim1(layout1, npmuim);
  // const GF3D2<CCTK_REAL> gf_nplambdare1(layout1, nplambdare);
  // const GF3D2<CCTK_REAL> gf_nplambdaim1(layout1, nplambdaim);
  // const GF3D2<CCTK_REAL> gf_npnure1(layout1, npnure);
  // const GF3D2<CCTK_REAL> gf_npnuim1(layout1, npnuim);

  //

  const smat<GF3D5<CCTK_REAL>, 4> gf_g45(std::array<GF3D5<CCTK_REAL>, 10>{
      GF3D5<CCTK_REAL>(layout5, g4tt), GF3D5<CCTK_REAL>(layout5, g4tx),
      GF3D5<CCTK_REAL>(layout5, g4ty), GF3D5<CCTK_REAL>(layout5, g4tz),
      GF3D5<CCTK_REAL>(layout5, g4xx), GF3D5<CCTK_REAL>(layout5, g4xy),
      GF3D5<CCTK_REAL>(layout5, g4xz), GF3D5<CCTK_REAL>(layout5, g4yy),
      GF3D5<CCTK_REAL>(layout5, g4yz), GF3D5<CCTK_REAL>(layout5, g4zz)});

  const GF3D5<CCTK_REAL> gf_Psi0re5(layout5, Psi0re);
  const GF3D5<CCTK_REAL> gf_Psi0im5(layout5, Psi0im);
  const GF3D5<CCTK_REAL> gf_Psi1re5(layout5, Psi1re);
  const GF3D5<CCTK_REAL> gf_Psi1im5(layout5, Psi1im);
  const GF3D5<CCTK_REAL> gf_Psi2re5(layout5, Psi2re);
  const GF3D5<CCTK_REAL> gf_Psi2im5(layout5, Psi2im);
  const GF3D5<CCTK_REAL> gf_Psi3re5(layout5, Psi3re);
  const GF3D5<CCTK_REAL> gf_Psi3im5(layout5, Psi3im);
  const GF3D5<CCTK_REAL> gf_Psi4re5(layout5, Psi4re);
  const GF3D5<CCTK_REAL> gf_Psi4im5(layout5, Psi4im);

  //

  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;

  const Loop::GridDescBaseDevice grid(cctkGH);

#if 1

  noinline([&]() __attribute__((__flatten__, __hot__)) {
    grid.loop_int_device<0, 0, 0, vsize>(
        grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
          const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
          const GF3D2index index1(layout1, p.I);
          const GF3D5index index0(layout0, p.I);
          const GF3D5index index5(layout5, p.I);

          // Load and calculate

          const auto id3 = one<smat<int, 3> >()();

          const vec<vreal, 3> coord3(
              [&](int d) { return p.X[d] + iota<vreal>() * p.DX[d]; });

          const weyl_vars<vreal> vars(
              cctk_time, coord3, //
              gf_gamma0(mask, index0, id3), gf_alpha0(mask, index0, 1),
              gf_beta0(mask, index0), //
              gf_k0(mask, index0), gf_dtalpha0(mask, index0),
              gf_dtbeta0(mask, index0), //
              gf_dgamma0(mask, index0), gf_dalpha0(mask, index0),
              gf_dbeta0(mask, index0), //
              gf_dtk1(mask, index1), gf_dt2alpha1(mask, index1),
              gf_dt2beta1(mask, index1), //
              gf_dk0(mask, index0), gf_ddtalpha0(mask, index0),
              gf_ddtbeta0(mask, index0), //
              gf_ddgamma0(mask, index0), gf_ddalpha0(mask, index0),
              gf_ddbeta0(mask, index0));

          // Store
          gf_g45.store(mask, index5, vars.g);

          // gf_Gamma4ttt1(p.I) = vars.Gamma(0)(0, 0);
          // gf_Gamma4ttx1(p.I) = vars.Gamma(0)(0, 1);
          // gf_Gamma4tty1(p.I) = vars.Gamma(0)(0, 2);
          // gf_Gamma4ttz1(p.I) = vars.Gamma(0)(0, 3);
          // gf_Gamma4txx1(p.I) = vars.Gamma(0)(1, 1);
          // gf_Gamma4txy1(p.I) = vars.Gamma(0)(1, 2);
          // gf_Gamma4txz1(p.I) = vars.Gamma(0)(1, 3);
          // gf_Gamma4tyy1(p.I) = vars.Gamma(0)(2, 2);
          // gf_Gamma4tyz1(p.I) = vars.Gamma(0)(2, 3);
          // gf_Gamma4tzz1(p.I) = vars.Gamma(0)(3, 3);

          // gf_Gamma4xtt1(p.I) = vars.Gamma(1)(0, 0);
          // gf_Gamma4xtx1(p.I) = vars.Gamma(1)(0, 1);
          // gf_Gamma4xty1(p.I) = vars.Gamma(1)(0, 2);
          // gf_Gamma4xtz1(p.I) = vars.Gamma(1)(0, 3);
          // gf_Gamma4xxx1(p.I) = vars.Gamma(1)(1, 1);
          // gf_Gamma4xxy1(p.I) = vars.Gamma(1)(1, 2);
          // gf_Gamma4xxz1(p.I) = vars.Gamma(1)(1, 3);
          // gf_Gamma4xyy1(p.I) = vars.Gamma(1)(2, 2);
          // gf_Gamma4xyz1(p.I) = vars.Gamma(1)(2, 3);
          // gf_Gamma4xzz1(p.I) = vars.Gamma(1)(3, 3);

          // gf_Gamma4ytt1(p.I) = vars.Gamma(2)(0, 0);
          // gf_Gamma4ytx1(p.I) = vars.Gamma(2)(0, 1);
          // gf_Gamma4yty1(p.I) = vars.Gamma(2)(0, 2);
          // gf_Gamma4ytz1(p.I) = vars.Gamma(2)(0, 3);
          // gf_Gamma4yxx1(p.I) = vars.Gamma(2)(1, 1);
          // gf_Gamma4yxy1(p.I) = vars.Gamma(2)(1, 2);
          // gf_Gamma4yxz1(p.I) = vars.Gamma(2)(1, 3);
          // gf_Gamma4yyy1(p.I) = vars.Gamma(2)(2, 2);
          // gf_Gamma4yyz1(p.I) = vars.Gamma(2)(2, 3);
          // gf_Gamma4yzz1(p.I) = vars.Gamma(2)(3, 3);

          // gf_Gamma4ztt1(p.I) = vars.Gamma(3)(0, 0);
          // gf_Gamma4ztx1(p.I) = vars.Gamma(3)(0, 1);
          // gf_Gamma4zty1(p.I) = vars.Gamma(3)(0, 2);
          // gf_Gamma4ztz1(p.I) = vars.Gamma(3)(0, 3);
          // gf_Gamma4zxx1(p.I) = vars.Gamma(3)(1, 1);
          // gf_Gamma4zxy1(p.I) = vars.Gamma(3)(1, 2);
          // gf_Gamma4zxz1(p.I) = vars.Gamma(3)(1, 3);
          // gf_Gamma4zyy1(p.I) = vars.Gamma(3)(2, 2);
          // gf_Gamma4zyz1(p.I) = vars.Gamma(3)(2, 3);
          // gf_Gamma4zzz1(p.I) = vars.Gamma(3)(3, 3);

          // gf_rm4txtx1(p.I) = vars.Rm(0, 1)(0, 1);
          // gf_rm4txty1(p.I) = vars.Rm(0, 1)(0, 2);
          // gf_rm4txtz1(p.I) = vars.Rm(0, 1)(0, 3);
          // gf_rm4txxy1(p.I) = vars.Rm(0, 1)(1, 2);
          // gf_rm4txxz1(p.I) = vars.Rm(0, 1)(1, 3);
          // gf_rm4txyz1(p.I) = vars.Rm(0, 1)(2, 3);

          // gf_rm4tyty1(p.I) = vars.Rm(0, 2)(0, 2);
          // gf_rm4tytz1(p.I) = vars.Rm(0, 2)(0, 3);
          // gf_rm4tyxy1(p.I) = vars.Rm(0, 2)(1, 2);
          // gf_rm4tyxz1(p.I) = vars.Rm(0, 2)(1, 3);
          // gf_rm4tyyz1(p.I) = vars.Rm(0, 2)(2, 3);

          // gf_rm4tztz1(p.I) = vars.Rm(0, 3)(0, 3);
          // gf_rm4tzxy1(p.I) = vars.Rm(0, 3)(1, 2);
          // gf_rm4tzxz1(p.I) = vars.Rm(0, 3)(1, 3);
          // gf_rm4tzyz1(p.I) = vars.Rm(0, 3)(2, 3);

          // gf_rm4xyxy1(p.I) = vars.Rm(1, 2)(1, 2);
          // gf_rm4xyxz1(p.I) = vars.Rm(1, 2)(1, 3);
          // gf_rm4xyyz1(p.I) = vars.Rm(1, 2)(2, 3);

          // gf_rm4xzxz1(p.I) = vars.Rm(1, 3)(1, 3);
          // gf_rm4xzyz1(p.I) = vars.Rm(1, 3)(2, 3);

          // gf_rm4yzyz1(p.I) = vars.Rm(2, 3)(2, 3);

          // vars.R.store(gf_r4tt1, gf_r4tx1, gf_r4ty1, gf_r4tz1, gf_r4xx1,
          // gf_r4xy1,
          //              gf_r4xz1, gf_r4yy1, gf_r4yz1, gf_r4zz1, p.I);

          // gf_rsc41(p.I) = vars.Rsc;

          // gf_c4txtx1(p.I) = vars.C(0, 1)(0, 1);
          // gf_c4txty1(p.I) = vars.C(0, 1)(0, 2);
          // gf_c4txtz1(p.I) = vars.C(0, 1)(0, 3);
          // gf_c4txxy1(p.I) = vars.C(0, 1)(1, 2);
          // gf_c4txxz1(p.I) = vars.C(0, 1)(1, 3);
          // gf_c4txyz1(p.I) = vars.C(0, 1)(2, 3);

          // gf_c4tyty1(p.I) = vars.C(0, 2)(0, 2);
          // gf_c4tytz1(p.I) = vars.C(0, 2)(0, 3);
          // gf_c4tyxy1(p.I) = vars.C(0, 2)(1, 2);
          // gf_c4tyxz1(p.I) = vars.C(0, 2)(1, 3);
          // gf_c4tyyz1(p.I) = vars.C(0, 2)(2, 3);

          // gf_c4tztz1(p.I) = vars.C(0, 3)(0, 3);
          // gf_c4tzxy1(p.I) = vars.C(0, 3)(1, 2);
          // gf_c4tzxz1(p.I) = vars.C(0, 3)(1, 3);
          // gf_c4tzyz1(p.I) = vars.C(0, 3)(2, 3);

          // gf_c4xyxy1(p.I) = vars.C(1, 2)(1, 2);
          // gf_c4xyxz1(p.I) = vars.C(1, 2)(1, 3);
          // gf_c4xyyz1(p.I) = vars.C(1, 2)(2, 3);

          // gf_c4xzxz1(p.I) = vars.C(1, 3)(1, 3);
          // gf_c4xzyz1(p.I) = vars.C(1, 3)(2, 3);

          // gf_c4yzyz1(p.I) = vars.C(2, 3)(2, 3);

          // gf_l1.store(mask, index1, vars.l);
          // gf_n1.store(mask, index1, vars.n);
          // gf_mre1.store(mask, index1,
          //               fmap([](const auto &x) { return real(x); }, vars.m));
          // gf_mim1.store(mask, index1,
          //               fmap([](const auto &x) { return imag(x); }, vars.m));

          // gf_Lambda1(p.I) = vars.Lambda;
          // gf_Phi00(p.I) = vars.Phi00;
          // gf_Phi111(p.I) = vars.Phi11;
          // gf_Phi221(p.I) = vars.Phi22;
          // gf_Phi10re1(p.I) = real(vars.Phi10);
          // gf_Phi10im1(p.I) = imag(vars.Phi10);
          // gf_Phi20re1(p.I) = real(vars.Phi20);
          // gf_Phi20im1(p.I) = imag(vars.Phi20);
          // gf_Phi21re1(p.I) = real(vars.Phi21);
          // gf_Phi21im1(p.I) = imag(vars.Phi21);

          gf_Psi0re5.store(mask, index5, real(vars.Psi0));
          gf_Psi0im5.store(mask, index5, imag(vars.Psi0));
          gf_Psi1re5.store(mask, index5, real(vars.Psi1));
          gf_Psi1im5.store(mask, index5, imag(vars.Psi1));
          gf_Psi2re5.store(mask, index5, real(vars.Psi2));
          gf_Psi2im5.store(mask, index5, imag(vars.Psi2));
          gf_Psi3re5.store(mask, index5, real(vars.Psi3));
          gf_Psi3im5.store(mask, index5, imag(vars.Psi3));
          gf_Psi4re5.store(mask, index5, real(vars.Psi4));
          gf_Psi4im5.store(mask, index5, imag(vars.Psi4));

          // gf_npkappare1(p.I) = real(vars.npkappa);
          // gf_npkappaim1(p.I) = imag(vars.npkappa);
          // gf_npsigmare1(p.I) = real(vars.npsigma);
          // gf_npsigmaim1(p.I) = imag(vars.npsigma);
          // gf_nprhore1(p.I) = real(vars.nprho);
          // gf_nprhoim1(p.I) = imag(vars.nprho);
          // gf_nptaure1(p.I) = real(vars.nptau);
          // gf_nptauim1(p.I) = imag(vars.nptau);
          // gf_npepsilonre1(p.I) = real(vars.npepsilon);
          // gf_npepsilonim1(p.I) = imag(vars.npepsilon);
          // gf_npbetare1(p.I) = real(vars.npbeta);
          // gf_npbetaim1(p.I) = imag(vars.npbeta);
          // gf_npalphare1(p.I) = real(vars.npalpha);
          // gf_npalphaim1(p.I) = imag(vars.npalpha);
          // gf_npgammare1(p.I) = real(vars.npgamma);
          // gf_npgammaim1(p.I) = imag(vars.npgamma);
          // gf_nppire1(p.I) = real(vars.nppi);
          // gf_nppiim1(p.I) = imag(vars.nppi);
          // gf_npmure1(p.I) = real(vars.npmu);
          // gf_npmuim1(p.I) = imag(vars.npmu);
          // gf_nplambdare1(p.I) = real(vars.nplambda);
          // gf_nplambdaim1(p.I) = imag(vars.nplambda);
          // gf_npnure1(p.I) = real(vars.npnu);
          // gf_npnuim1(p.I) = imag(vars.npnu);
        });
  });

#elif 0

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Load and calculate

        const vec3<CCTK_REAL> coord3{p.x, p.y, p.z};

        const weyl_vars<CCTK_REAL> vars(
            cctk_time, coord3, //
            gf_gamma0(layout0, p.I), gf_alpha0(layout0, p.I),
            gf_beta0(layout0, p.I), //
            gf_k0(layout0, p.I), gf_dtalpha0(layout0, p.I),
            gf_dtbeta0(layout0, p.I), //
            gf_dgamma0(layout0, p.I), gf_dalpha0(layout0, p.I),
            gf_dbeta0(layout0, p.I), //
            gf_dtk1(p.I), gf_dt2alpha1(p.I),
            gf_dt2beta1(p.I), //
            gf_dk0(layout0, p.I), gf_ddtalpha0(layout0, p.I),
            gf_ddtbeta0(layout0, p.I), //
            gf_ddgamma0(layout0, p.I), gf_ddalpha0(layout0, p.I),
            gf_ddbeta0(layout0, p.I));

        // Store
        vars.g.store(gf_g4tt1, gf_g4tx1, gf_g4ty1, gf_g4tz1, gf_g4xx1, gf_g4xy1,
                     gf_g4xz1, gf_g4yy1, gf_g4yz1, gf_g4zz1, p.I);

        // gf_Gamma4ttt1(p.I) = vars.Gamma(0)(0, 0);
        // gf_Gamma4ttx1(p.I) = vars.Gamma(0)(0, 1);
        // gf_Gamma4tty1(p.I) = vars.Gamma(0)(0, 2);
        // gf_Gamma4ttz1(p.I) = vars.Gamma(0)(0, 3);
        // gf_Gamma4txx1(p.I) = vars.Gamma(0)(1, 1);
        // gf_Gamma4txy1(p.I) = vars.Gamma(0)(1, 2);
        // gf_Gamma4txz1(p.I) = vars.Gamma(0)(1, 3);
        // gf_Gamma4tyy1(p.I) = vars.Gamma(0)(2, 2);
        // gf_Gamma4tyz1(p.I) = vars.Gamma(0)(2, 3);
        // gf_Gamma4tzz1(p.I) = vars.Gamma(0)(3, 3);

        // gf_Gamma4xtt1(p.I) = vars.Gamma(1)(0, 0);
        // gf_Gamma4xtx1(p.I) = vars.Gamma(1)(0, 1);
        // gf_Gamma4xty1(p.I) = vars.Gamma(1)(0, 2);
        // gf_Gamma4xtz1(p.I) = vars.Gamma(1)(0, 3);
        // gf_Gamma4xxx1(p.I) = vars.Gamma(1)(1, 1);
        // gf_Gamma4xxy1(p.I) = vars.Gamma(1)(1, 2);
        // gf_Gamma4xxz1(p.I) = vars.Gamma(1)(1, 3);
        // gf_Gamma4xyy1(p.I) = vars.Gamma(1)(2, 2);
        // gf_Gamma4xyz1(p.I) = vars.Gamma(1)(2, 3);
        // gf_Gamma4xzz1(p.I) = vars.Gamma(1)(3, 3);

        // gf_Gamma4ytt1(p.I) = vars.Gamma(2)(0, 0);
        // gf_Gamma4ytx1(p.I) = vars.Gamma(2)(0, 1);
        // gf_Gamma4yty1(p.I) = vars.Gamma(2)(0, 2);
        // gf_Gamma4ytz1(p.I) = vars.Gamma(2)(0, 3);
        // gf_Gamma4yxx1(p.I) = vars.Gamma(2)(1, 1);
        // gf_Gamma4yxy1(p.I) = vars.Gamma(2)(1, 2);
        // gf_Gamma4yxz1(p.I) = vars.Gamma(2)(1, 3);
        // gf_Gamma4yyy1(p.I) = vars.Gamma(2)(2, 2);
        // gf_Gamma4yyz1(p.I) = vars.Gamma(2)(2, 3);
        // gf_Gamma4yzz1(p.I) = vars.Gamma(2)(3, 3);

        // gf_Gamma4ztt1(p.I) = vars.Gamma(3)(0, 0);
        // gf_Gamma4ztx1(p.I) = vars.Gamma(3)(0, 1);
        // gf_Gamma4zty1(p.I) = vars.Gamma(3)(0, 2);
        // gf_Gamma4ztz1(p.I) = vars.Gamma(3)(0, 3);
        // gf_Gamma4zxx1(p.I) = vars.Gamma(3)(1, 1);
        // gf_Gamma4zxy1(p.I) = vars.Gamma(3)(1, 2);
        // gf_Gamma4zxz1(p.I) = vars.Gamma(3)(1, 3);
        // gf_Gamma4zyy1(p.I) = vars.Gamma(3)(2, 2);
        // gf_Gamma4zyz1(p.I) = vars.Gamma(3)(2, 3);
        // gf_Gamma4zzz1(p.I) = vars.Gamma(3)(3, 3);

        // gf_rm4txtx1(p.I) = vars.Rm(0, 1)(0, 1);
        // gf_rm4txty1(p.I) = vars.Rm(0, 1)(0, 2);
        // gf_rm4txtz1(p.I) = vars.Rm(0, 1)(0, 3);
        // gf_rm4txxy1(p.I) = vars.Rm(0, 1)(1, 2);
        // gf_rm4txxz1(p.I) = vars.Rm(0, 1)(1, 3);
        // gf_rm4txyz1(p.I) = vars.Rm(0, 1)(2, 3);

        // gf_rm4tyty1(p.I) = vars.Rm(0, 2)(0, 2);
        // gf_rm4tytz1(p.I) = vars.Rm(0, 2)(0, 3);
        // gf_rm4tyxy1(p.I) = vars.Rm(0, 2)(1, 2);
        // gf_rm4tyxz1(p.I) = vars.Rm(0, 2)(1, 3);
        // gf_rm4tyyz1(p.I) = vars.Rm(0, 2)(2, 3);

        // gf_rm4tztz1(p.I) = vars.Rm(0, 3)(0, 3);
        // gf_rm4tzxy1(p.I) = vars.Rm(0, 3)(1, 2);
        // gf_rm4tzxz1(p.I) = vars.Rm(0, 3)(1, 3);
        // gf_rm4tzyz1(p.I) = vars.Rm(0, 3)(2, 3);

        // gf_rm4xyxy1(p.I) = vars.Rm(1, 2)(1, 2);
        // gf_rm4xyxz1(p.I) = vars.Rm(1, 2)(1, 3);
        // gf_rm4xyyz1(p.I) = vars.Rm(1, 2)(2, 3);

        // gf_rm4xzxz1(p.I) = vars.Rm(1, 3)(1, 3);
        // gf_rm4xzyz1(p.I) = vars.Rm(1, 3)(2, 3);

        // gf_rm4yzyz1(p.I) = vars.Rm(2, 3)(2, 3);

        // vars.R.store(gf_r4tt1, gf_r4tx1, gf_r4ty1, gf_r4tz1, gf_r4xx1,
        // gf_r4xy1,
        //              gf_r4xz1, gf_r4yy1, gf_r4yz1, gf_r4zz1, p.I);

        // gf_rsc41(p.I) = vars.Rsc;

        // gf_c4txtx1(p.I) = vars.C(0, 1)(0, 1);
        // gf_c4txty1(p.I) = vars.C(0, 1)(0, 2);
        // gf_c4txtz1(p.I) = vars.C(0, 1)(0, 3);
        // gf_c4txxy1(p.I) = vars.C(0, 1)(1, 2);
        // gf_c4txxz1(p.I) = vars.C(0, 1)(1, 3);
        // gf_c4txyz1(p.I) = vars.C(0, 1)(2, 3);

        // gf_c4tyty1(p.I) = vars.C(0, 2)(0, 2);
        // gf_c4tytz1(p.I) = vars.C(0, 2)(0, 3);
        // gf_c4tyxy1(p.I) = vars.C(0, 2)(1, 2);
        // gf_c4tyxz1(p.I) = vars.C(0, 2)(1, 3);
        // gf_c4tyyz1(p.I) = vars.C(0, 2)(2, 3);

        // gf_c4tztz1(p.I) = vars.C(0, 3)(0, 3);
        // gf_c4tzxy1(p.I) = vars.C(0, 3)(1, 2);
        // gf_c4tzxz1(p.I) = vars.C(0, 3)(1, 3);
        // gf_c4tzyz1(p.I) = vars.C(0, 3)(2, 3);

        // gf_c4xyxy1(p.I) = vars.C(1, 2)(1, 2);
        // gf_c4xyxz1(p.I) = vars.C(1, 2)(1, 3);
        // gf_c4xyyz1(p.I) = vars.C(1, 2)(2, 3);

        // gf_c4xzxz1(p.I) = vars.C(1, 3)(1, 3);
        // gf_c4xzyz1(p.I) = vars.C(1, 3)(2, 3);

        // gf_c4yzyz1(p.I) = vars.C(2, 3)(2, 3);

        // vars.l.store(gf_lt1, gf_lx1, gf_ly1, gf_lz1, p.I);
        // vars.n.store(gf_nt1, gf_nx1, gf_ny1, gf_nz1, p.I);
        // gf_mret1(p.I) = real(vars.m(0));
        // gf_mrex1(p.I) = real(vars.m(1));
        // gf_mrey1(p.I) = real(vars.m(2));
        // gf_mrez1(p.I) = real(vars.m(3));
        // gf_mimt1(p.I) = imag(vars.m(0));
        // gf_mimx1(p.I) = imag(vars.m(1));
        // gf_mimy1(p.I) = imag(vars.m(2));
        // gf_mimz1(p.I) = imag(vars.m(3));

        // gf_Lambda1(p.I) = vars.Lambda;
        // gf_Phi00(p.I) = vars.Phi00;
        // gf_Phi111(p.I) = vars.Phi11;
        // gf_Phi221(p.I) = vars.Phi22;
        // gf_Phi10re1(p.I) = real(vars.Phi10);
        // gf_Phi10im1(p.I) = imag(vars.Phi10);
        // gf_Phi20re1(p.I) = real(vars.Phi20);
        // gf_Phi20im1(p.I) = imag(vars.Phi20);
        // gf_Phi21re1(p.I) = real(vars.Phi21);
        // gf_Phi21im1(p.I) = imag(vars.Phi21);

        gf_Psi0re1(p.I) = real(vars.Psi0);
        gf_Psi0im1(p.I) = imag(vars.Psi0);
        gf_Psi1re1(p.I) = real(vars.Psi1);
        gf_Psi1im1(p.I) = imag(vars.Psi1);
        gf_Psi2re1(p.I) = real(vars.Psi2);
        gf_Psi2im1(p.I) = imag(vars.Psi2);
        gf_Psi3re1(p.I) = real(vars.Psi3);
        gf_Psi3im1(p.I) = imag(vars.Psi3);
        gf_Psi4re1(p.I) = real(vars.Psi4);
        gf_Psi4im1(p.I) = imag(vars.Psi4);

        // gf_npkappare1(p.I) = real(vars.npkappa);
        // gf_npkappaim1(p.I) = imag(vars.npkappa);
        // gf_npsigmare1(p.I) = real(vars.npsigma);
        // gf_npsigmaim1(p.I) = imag(vars.npsigma);
        // gf_nprhore1(p.I) = real(vars.nprho);
        // gf_nprhoim1(p.I) = imag(vars.nprho);
        // gf_nptaure1(p.I) = real(vars.nptau);
        // gf_nptauim1(p.I) = imag(vars.nptau);
        // gf_npepsilonre1(p.I) = real(vars.npepsilon);
        // gf_npepsilonim1(p.I) = imag(vars.npepsilon);
        // gf_npbetare1(p.I) = real(vars.npbeta);
        // gf_npbetaim1(p.I) = imag(vars.npbeta);
        // gf_npalphare1(p.I) = real(vars.npalpha);
        // gf_npalphaim1(p.I) = imag(vars.npalpha);
        // gf_npgammare1(p.I) = real(vars.npgamma);
        // gf_npgammaim1(p.I) = imag(vars.npgamma);
        // gf_nppire1(p.I) = real(vars.nppi);
        // gf_nppiim1(p.I) = imag(vars.nppi);
        // gf_npmure1(p.I) = real(vars.npmu);
        // gf_npmuim1(p.I) = imag(vars.npmu);
        // gf_nplambdare1(p.I) = real(vars.nplambda);
        // gf_nplambdaim1(p.I) = imag(vars.nplambda);
        // gf_npnure1(p.I) = real(vars.npnu);
        // gf_npnuim1(p.I) = imag(vars.npnu);
      });

#else

  const mat4<GF3D5<CCTK_REAL> > gf_g0(make_mat_gf());
  const mat4<vec4<GF3D5<CCTK_REAL> > > gf_dg0(make_mat_vec_gf());
  const mat4<smat4<GF3D5<CCTK_REAL> > > gf_ddg0(make_mat_mat_gf());

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const GF3D5index layout0, p.I(layout0, p.I);

        // Load and calculate

        const vec3<CCTK_REAL> coord3{p.x, p.y, p.z};

        const weyl_vars<CCTK_REAL> vars(
            cctk_time, coord3, //
            gf_gamma0(layout0, p.I), gf_alpha0(layout0, p.I),
            gf_beta0(layout0, p.I), //
            gf_k0(layout0, p.I), gf_dtalpha0(layout0, p.I),
            gf_dtbeta0(layout0, p.I), //
            gf_dgamma0(layout0, p.I), gf_dalpha0(layout0, p.I),
            gf_dbeta0(layout0, p.I), //
            gf_dtk1(p.I), gf_dt2alpha1(p.I),
            gf_dt2beta1(p.I), //
            gf_dk0(layout0, p.I), gf_ddtalpha0(layout0, p.I),
            gf_ddtbeta0(layout0, p.I), //
            gf_ddgamma0(layout0, p.I), gf_ddalpha0(layout0, p.I),
            gf_ddbeta0(layout0, p.I));

        // Store
        for (int a = 0; a < 4; ++a)
          for (int b = a; b < 4; ++b)
            gf_g0(a, b)(layout0, p.I) = vars.g(a, b);
        for (int a = 0; a < 4; ++a)
          for (int b = a; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
              gf_dg0(a, b)(c)(layout0, p.I) = vars.dg(a, b)(c);
        for (int a = 0; a < 4; ++a)
          for (int b = a; b < 4; ++b)
            for (int c = 0; c < 4; ++c)
              for (int d = c; d < 4; ++d)
                gf_ddg0(a, b)(c, d)(layout0, p.I) = vars.ddg(a, b)(c, d);
      });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const GF3D5index layout0, p.I(layout0, p.I);

        // Load and calculate

        const vec4<CCTK_REAL> coord{cctk_time, p.x, p.y, p.z};

        const weyl_vars<CCTK_REAL> vars(
            coord, //
            mat4<CCTK_REAL>(
                [&](int a, int b) { return gf_g0(a, b)(layout0, p.I); }),
            mat4<vec4<CCTK_REAL> >([&](int a, int b) {
              return vec4<CCTK_REAL>(
                  [&](int c) { return gf_dg0(a, b)(c)(layout0, p.I); });
            }),
            mat4<smat4<CCTK_REAL> >([&](int a, int b) {
              return mat4<CCTK_REAL>([&](int c, int d) {
                return gf_ddg0(a, b)(c, d)(layout0, p.I);
              });
            }));

        // Store
        vars.g.store(gf_g4tt1, gf_g4tx1, gf_g4ty1, gf_g4tz1, gf_g4xx1, gf_g4xy1,
                     gf_g4xz1, gf_g4yy1, gf_g4yz1, gf_g4zz1, p.I);

        // gf_Gamma4ttt1(p.I) = vars.Gamma(0)(0, 0);
        // gf_Gamma4ttx1(p.I) = vars.Gamma(0)(0, 1);
        // gf_Gamma4tty1(p.I) = vars.Gamma(0)(0, 2);
        // gf_Gamma4ttz1(p.I) = vars.Gamma(0)(0, 3);
        // gf_Gamma4txx1(p.I) = vars.Gamma(0)(1, 1);
        // gf_Gamma4txy1(p.I) = vars.Gamma(0)(1, 2);
        // gf_Gamma4txz1(p.I) = vars.Gamma(0)(1, 3);
        // gf_Gamma4tyy1(p.I) = vars.Gamma(0)(2, 2);
        // gf_Gamma4tyz1(p.I) = vars.Gamma(0)(2, 3);
        // gf_Gamma4tzz1(p.I) = vars.Gamma(0)(3, 3);

        // gf_Gamma4xtt1(p.I) = vars.Gamma(1)(0, 0);
        // gf_Gamma4xtx1(p.I) = vars.Gamma(1)(0, 1);
        // gf_Gamma4xty1(p.I) = vars.Gamma(1)(0, 2);
        // gf_Gamma4xtz1(p.I) = vars.Gamma(1)(0, 3);
        // gf_Gamma4xxx1(p.I) = vars.Gamma(1)(1, 1);
        // gf_Gamma4xxy1(p.I) = vars.Gamma(1)(1, 2);
        // gf_Gamma4xxz1(p.I) = vars.Gamma(1)(1, 3);
        // gf_Gamma4xyy1(p.I) = vars.Gamma(1)(2, 2);
        // gf_Gamma4xyz1(p.I) = vars.Gamma(1)(2, 3);
        // gf_Gamma4xzz1(p.I) = vars.Gamma(1)(3, 3);

        // gf_Gamma4ytt1(p.I) = vars.Gamma(2)(0, 0);
        // gf_Gamma4ytx1(p.I) = vars.Gamma(2)(0, 1);
        // gf_Gamma4yty1(p.I) = vars.Gamma(2)(0, 2);
        // gf_Gamma4ytz1(p.I) = vars.Gamma(2)(0, 3);
        // gf_Gamma4yxx1(p.I) = vars.Gamma(2)(1, 1);
        // gf_Gamma4yxy1(p.I) = vars.Gamma(2)(1, 2);
        // gf_Gamma4yxz1(p.I) = vars.Gamma(2)(1, 3);
        // gf_Gamma4yyy1(p.I) = vars.Gamma(2)(2, 2);
        // gf_Gamma4yyz1(p.I) = vars.Gamma(2)(2, 3);
        // gf_Gamma4yzz1(p.I) = vars.Gamma(2)(3, 3);

        // gf_Gamma4ztt1(p.I) = vars.Gamma(3)(0, 0);
        // gf_Gamma4ztx1(p.I) = vars.Gamma(3)(0, 1);
        // gf_Gamma4zty1(p.I) = vars.Gamma(3)(0, 2);
        // gf_Gamma4ztz1(p.I) = vars.Gamma(3)(0, 3);
        // gf_Gamma4zxx1(p.I) = vars.Gamma(3)(1, 1);
        // gf_Gamma4zxy1(p.I) = vars.Gamma(3)(1, 2);
        // gf_Gamma4zxz1(p.I) = vars.Gamma(3)(1, 3);
        // gf_Gamma4zyy1(p.I) = vars.Gamma(3)(2, 2);
        // gf_Gamma4zyz1(p.I) = vars.Gamma(3)(2, 3);
        // gf_Gamma4zzz1(p.I) = vars.Gamma(3)(3, 3);

        // gf_rm4txtx1(p.I) = vars.Rm(0, 1)(0, 1);
        // gf_rm4txty1(p.I) = vars.Rm(0, 1)(0, 2);
        // gf_rm4txtz1(p.I) = vars.Rm(0, 1)(0, 3);
        // gf_rm4txxy1(p.I) = vars.Rm(0, 1)(1, 2);
        // gf_rm4txxz1(p.I) = vars.Rm(0, 1)(1, 3);
        // gf_rm4txyz1(p.I) = vars.Rm(0, 1)(2, 3);

        // gf_rm4tyty1(p.I) = vars.Rm(0, 2)(0, 2);
        // gf_rm4tytz1(p.I) = vars.Rm(0, 2)(0, 3);
        // gf_rm4tyxy1(p.I) = vars.Rm(0, 2)(1, 2);
        // gf_rm4tyxz1(p.I) = vars.Rm(0, 2)(1, 3);
        // gf_rm4tyyz1(p.I) = vars.Rm(0, 2)(2, 3);

        // gf_rm4tztz1(p.I) = vars.Rm(0, 3)(0, 3);
        // gf_rm4tzxy1(p.I) = vars.Rm(0, 3)(1, 2);
        // gf_rm4tzxz1(p.I) = vars.Rm(0, 3)(1, 3);
        // gf_rm4tzyz1(p.I) = vars.Rm(0, 3)(2, 3);

        // gf_rm4xyxy1(p.I) = vars.Rm(1, 2)(1, 2);
        // gf_rm4xyxz1(p.I) = vars.Rm(1, 2)(1, 3);
        // gf_rm4xyyz1(p.I) = vars.Rm(1, 2)(2, 3);

        // gf_rm4xzxz1(p.I) = vars.Rm(1, 3)(1, 3);
        // gf_rm4xzyz1(p.I) = vars.Rm(1, 3)(2, 3);

        // gf_rm4yzyz1(p.I) = vars.Rm(2, 3)(2, 3);

        // vars.R.store(gf_r4tt1, gf_r4tx1, gf_r4ty1, gf_r4tz1, gf_r4xx1,
        // gf_r4xy1,
        //              gf_r4xz1, gf_r4yy1, gf_r4yz1, gf_r4zz1, p.I);

        // gf_rsc41(p.I) = vars.Rsc;

        // gf_c4txtx1(p.I) = vars.C(0, 1)(0, 1);
        // gf_c4txty1(p.I) = vars.C(0, 1)(0, 2);
        // gf_c4txtz1(p.I) = vars.C(0, 1)(0, 3);
        // gf_c4txxy1(p.I) = vars.C(0, 1)(1, 2);
        // gf_c4txxz1(p.I) = vars.C(0, 1)(1, 3);
        // gf_c4txyz1(p.I) = vars.C(0, 1)(2, 3);

        // gf_c4tyty1(p.I) = vars.C(0, 2)(0, 2);
        // gf_c4tytz1(p.I) = vars.C(0, 2)(0, 3);
        // gf_c4tyxy1(p.I) = vars.C(0, 2)(1, 2);
        // gf_c4tyxz1(p.I) = vars.C(0, 2)(1, 3);
        // gf_c4tyyz1(p.I) = vars.C(0, 2)(2, 3);

        // gf_c4tztz1(p.I) = vars.C(0, 3)(0, 3);
        // gf_c4tzxy1(p.I) = vars.C(0, 3)(1, 2);
        // gf_c4tzxz1(p.I) = vars.C(0, 3)(1, 3);
        // gf_c4tzyz1(p.I) = vars.C(0, 3)(2, 3);

        // gf_c4xyxy1(p.I) = vars.C(1, 2)(1, 2);
        // gf_c4xyxz1(p.I) = vars.C(1, 2)(1, 3);
        // gf_c4xyyz1(p.I) = vars.C(1, 2)(2, 3);

        // gf_c4xzxz1(p.I) = vars.C(1, 3)(1, 3);
        // gf_c4xzyz1(p.I) = vars.C(1, 3)(2, 3);

        // gf_c4yzyz1(p.I) = vars.C(2, 3)(2, 3);

        // vars.l.store(gf_lt1, gf_lx1, gf_ly1, gf_lz1, p.I);
        // vars.n.store(gf_nt1, gf_nx1, gf_ny1, gf_nz1, p.I);
        // gf_mret1(p.I) = real(vars.m(0));
        // gf_mrex1(p.I) = real(vars.m(1));
        // gf_mrey1(p.I) = real(vars.m(2));
        // gf_mrez1(p.I) = real(vars.m(3));
        // gf_mimt1(p.I) = imag(vars.m(0));
        // gf_mimx1(p.I) = imag(vars.m(1));
        // gf_mimy1(p.I) = imag(vars.m(2));
        // gf_mimz1(p.I) = imag(vars.m(3));

        // gf_Lambda1(p.I) = vars.Lambda;
        // gf_Phi00(p.I) = vars.Phi00;
        // gf_Phi111(p.I) = vars.Phi11;
        // gf_Phi221(p.I) = vars.Phi22;
        // gf_Phi10re1(p.I) = real(vars.Phi10);
        // gf_Phi10im1(p.I) = imag(vars.Phi10);
        // gf_Phi20re1(p.I) = real(vars.Phi20);
        // gf_Phi20im1(p.I) = imag(vars.Phi20);
        // gf_Phi21re1(p.I) = real(vars.Phi21);
        // gf_Phi21im1(p.I) = imag(vars.Phi21);

        gf_Psi0re1(p.I) = real(vars.Psi0);
        gf_Psi0im1(p.I) = imag(vars.Psi0);
        gf_Psi1re1(p.I) = real(vars.Psi1);
        gf_Psi1im1(p.I) = imag(vars.Psi1);
        gf_Psi2re1(p.I) = real(vars.Psi2);
        gf_Psi2im1(p.I) = imag(vars.Psi2);
        gf_Psi3re1(p.I) = real(vars.Psi3);
        gf_Psi3im1(p.I) = imag(vars.Psi3);
        gf_Psi4re1(p.I) = real(vars.Psi4);
        gf_Psi4im1(p.I) = imag(vars.Psi4);

        // gf_npkappare1(p.I) = real(vars.npkappa);
        // gf_npkappaim1(p.I) = imag(vars.npkappa);
        // gf_npsigmare1(p.I) = real(vars.npsigma);
        // gf_npsigmaim1(p.I) = imag(vars.npsigma);
        // gf_nprhore1(p.I) = real(vars.nprho);
        // gf_nprhoim1(p.I) = imag(vars.nprho);
        // gf_nptaure1(p.I) = real(vars.nptau);
        // gf_nptauim1(p.I) = imag(vars.nptau);
        // gf_npepsilonre1(p.I) = real(vars.npepsilon);
        // gf_npepsilonim1(p.I) = imag(vars.npepsilon);
        // gf_npbetare1(p.I) = real(vars.npbeta);
        // gf_npbetaim1(p.I) = imag(vars.npbeta);
        // gf_npalphare1(p.I) = real(vars.npalpha);
        // gf_npalphaim1(p.I) = imag(vars.npalpha);
        // gf_npgammare1(p.I) = real(vars.npgamma);
        // gf_npgammaim1(p.I) = imag(vars.npgamma);
        // gf_nppire1(p.I) = real(vars.nppi);
        // gf_nppiim1(p.I) = imag(vars.nppi);
        // gf_npmure1(p.I) = real(vars.npmu);
        // gf_npmuim1(p.I) = imag(vars.npmu);
        // gf_nplambdare1(p.I) = real(vars.nplambda);
        // gf_nplambdaim1(p.I) = imag(vars.nplambda);
        // gf_npnure1(p.I) = real(vars.npnu);
        // gf_npnuim1(p.I) = imag(vars.npnu);
      });

#endif
}

} // namespace Weyl
