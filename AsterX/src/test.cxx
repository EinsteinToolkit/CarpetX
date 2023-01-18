#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "utils.hxx"

namespace AsterX {

using namespace Arith;

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  const smat<CCTK_REAL, 3> g{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const CCTK_REAL detg = -1.0;
  const smat<CCTK_REAL, 3> invg{1.0, -3.0, 2.0, 3.0, -1.0, 0.0};
  const vec<CCTK_REAL, 3> v_up{0.07, 0.08, 0.09};
  const CCTK_REAL tiny = 1e-14;

  {
    assert(calc_det(g) == detg);
    CCTK_VINFO("Test calc_det succeeded");
  }

  {
    const smat<CCTK_REAL, 3> invg_test = calc_inv(g, detg);
    assert(invg_test(0, 0) == 1.0);
    assert(invg_test(0, 1) == -3.0);
    assert(invg_test(0, 2) == 2.0);
    assert(invg_test(1, 1) == 3.0);
    assert(invg_test(1, 2) == -1.0);
    assert(invg_test(2, 2) == 0.0);
    assert(invg_test.elts[0] == 1.0);
    assert(invg_test.elts[1] == -3.0);
    assert(invg_test.elts[2] == 2.0);
    assert(invg_test.elts[3] == 3.0);
    assert(invg_test.elts[4] == -1.0);
    assert(invg_test.elts[5] == 0.0);
    assert(invg_test == invg);
    CCTK_VINFO("Test calc_inv succeeded");
  }

  {
    // const vec<CCTK_REAL, 3> v_dn([&](int i) ARITH_INLINE {
    //   return sum<3>([&](int j) ARITH_INLINE { return g(i, j) * v_up(j); });
    // });
    const vec<CCTK_REAL, 3> v_dn = calc_contraction(g, v_up);
    assert(v_dn(0) - 0.5 < tiny);
    assert(v_dn(1) - 0.91 < tiny);
    assert(v_dn(2) - 1.15 < tiny);
    CCTK_VINFO("Test calc_contraction of smat and vec succeeded");

    // const CCTK_REAL v2(
    //     sum<3>([&](int i) ARITH_INLINE { return v_up(i) * v_dn(i); }));
    const CCTK_REAL v2 = calc_contraction(v_up, v_dn);
    assert(v2 - 0.2113 < tiny);
    CCTK_VINFO("Test calc_contraction of vec and vec succeeded");

    const vec<CCTK_REAL, 3> vcv = calc_cross_product(v_dn, v_up);
    assert(vcv(0) - (-0.0101) < tiny);
    assert(vcv(1) - 0.0355 < tiny);
    assert(vcv(2) - (-0.0237) < tiny);
    CCTK_VINFO("Test calc_cross_product of vec and vec succeeded");

    const CCTK_REAL wlorentz = calc_wlorentz(v_up, v_dn);
    assert(wlorentz - 1.0 / sqrt(1.0 - v2) < tiny);
    CCTK_VINFO("Test calc_wlorentz of vec and vec succeeded");

    const CCTK_REAL two = 2.0;
    const vec<CCTK_REAL, 3> v_sum = v_dn * two + v_up / two;
    assert(v_sum(0) - 1.035 < tiny);
    assert(v_sum(1) - 1.86 < tiny);
    assert(v_sum(2) - 2.345 < tiny);
    CCTK_VINFO("Test sum/divide by scalar of vecs succeeded");

    const smat<CCTK_REAL, 3> vs3(
        [&](int i, int j) ARITH_INLINE { return v_dn(i) * v_dn(j); });
    assert(vs3(0, 0) - 0.5 * 0.5 < tiny);
    assert(vs3(0, 1) - 0.5 * 0.91 < tiny);
    assert(vs3(0, 2) - 0.5 * 1.15 < tiny);
    assert(vs3(1, 1) - 0.91 * 0.91 < tiny);
    assert(vs3(1, 2) - 0.91 * 1.15 < tiny);
    assert(vs3(2, 2) - 1.15 * 1.15 < tiny);
    CCTK_VINFO("Test tensor product of vec -> smat succeeded");
  }

  {
    const vec<CCTK_REAL, 3> v0{1.0, 2.0, 3.0};
    const array<CCTK_REAL, 3> v1_array = {4.0, 5.0, 6.0};
    const vec<CCTK_REAL, 3> v1{v1_array};
    const vec<vec<CCTK_REAL, 3>, 3> vv{v0, v1, {7.0, 8.0, 9.0}};
    assert(vv(0)(0) == 1.0);
    assert(vv(0)(1) == 2.0);
    assert(vv(0)(2) == 3.0);
    assert(vv(1)(0) == 4.0);
    assert(vv(1)(1) == 5.0);
    assert(vv(1)(2) == 6.0);
    assert(vv(2)(0) == 7.0);
    assert(vv(2)(1) == 8.0);
    assert(vv(2)(2) == 9.0);

    const vec<vec<CCTK_REAL, 3>, 3> vv2([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 3>([&](int j)
                                   ARITH_INLINE { return vv(i)(j) + 1.0; });
    });
    assert(vv2(0)(0) == 2.0);
    assert(vv2(0)(1) == 3.0);
    assert(vv2(0)(2) == 4.0);
    assert(vv2(1)(0) == 5.0);
    assert(vv2(1)(1) == 6.0);
    assert(vv2(1)(2) == 7.0);
    assert(vv2(2)(0) == 8.0);
    assert(vv2(2)(1) == 9.0);
    assert(vv2(2)(2) == 10.0);

    const vec<vec<CCTK_REAL, 3>, 3> vv3 = vv + vv2;
    assert(vv3(0)(0) == 3.0);
    assert(vv3(0)(1) == 5.0);
    assert(vv3(0)(2) == 7.0);
    assert(vv3(1)(0) == 9.0);
    assert(vv3(1)(1) == 11.0);
    assert(vv3(1)(2) == 13.0);
    assert(vv3(2)(0) == 15.0);
    assert(vv3(2)(1) == 17.0);
    assert(vv3(2)(2) == 19.0);

    const vec<vec<CCTK_REAL, 3>, 3> vv3_trans = calc_transpose(vv3);
    assert(vv3_trans(0)(0) == 3.0);
    assert(vv3_trans(1)(0) == 5.0);
    assert(vv3_trans(2)(0) == 7.0);
    assert(vv3_trans(0)(1) == 9.0);
    assert(vv3_trans(1)(1) == 11.0);
    assert(vv3_trans(2)(1) == 13.0);
    assert(vv3_trans(0)(2) == 15.0);
    assert(vv3_trans(1)(2) == 17.0);
    assert(vv3_trans(2)(2) == 19.0);

    CCTK_VINFO("Test vec of vecs succeeded");
  }

  {
    const vec<smat<CCTK_REAL, 3>, 2> vs{g, invg};
    assert(vs(0)(0, 0) == 1.0);
    assert(vs(0)(0, 1) == 2.0);
    assert(vs(0)(0, 2) == 3.0);
    assert(vs(0)(1, 1) == 4.0);
    assert(vs(0)(1, 2) == 5.0);
    assert(vs(0)(2, 2) == 6.0);
    assert(vs(1)(0, 0) == 1.0);
    assert(vs(1)(0, 1) == -3.0);
    assert(vs(1)(0, 2) == 2.0);
    assert(vs(1)(1, 1) == 3.0);
    assert(vs(1)(1, 2) == -1.0);
    assert(vs(1)(2, 2) == 0.0);

    const vec<smat<CCTK_REAL, 3>, 2> vs2([&](int k) ARITH_INLINE {
      return smat<CCTK_REAL, 3>([&](int i, int j)
                                    ARITH_INLINE { return vs(k)(i, j); });
    });
    assert(vs2(0)(0, 0) == 1.0);
    assert(vs2(0)(0, 1) == 2.0);
    assert(vs2(0)(0, 2) == 3.0);
    assert(vs2(0)(1, 1) == 4.0);
    assert(vs2(0)(1, 2) == 5.0);
    assert(vs2(0)(2, 2) == 6.0);
    assert(vs2(1)(0, 0) == 1.0);
    assert(vs2(1)(0, 1) == -3.0);
    assert(vs2(1)(0, 2) == 2.0);
    assert(vs2(1)(1, 1) == 3.0);
    assert(vs2(1)(1, 2) == -1.0);
    assert(vs2(1)(2, 2) == 0.0);

    CCTK_VINFO("Test vec of smats succeeded");
  }

  {
    const vect<CCTK_REAL, 3> v0{1.0, 2.0, 3.0};
    const vect<CCTK_REAL, 3> v1{4.0, 5.0, 6.0};
    const vect<CCTK_REAL, 3> v2 = v0 + v1;
    assert(v2[0] == 5.0);
    assert(v2[1] == 7.0);
    assert(v2[2] == 9.0);
    CCTK_VINFO("Test vect succeeded");
  }
}

} // namespace AsterX
