#include <loop_device.hxx>

#include <defs.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>
#include <derivs.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace TestDerivs {
using namespace Loop;
using namespace std;

template <typename T, typename VT>
constexpr void poly(const T kxx, const T kxy, const T kyz, const VT x,
                    const T y, const T z, VT &u) {
  u = kxx * x * x + kxy * x * y + kyz * y * z;
}

extern "C" void TestDerivs_SetError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDerivs_SetError;
  DECLARE_CCTK_PARAMETERS;

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
        fabs(p.z) <= refined_radius) {
      regrid_error(p.I) = 1;
    } else {
      regrid_error(p.I) = 0;
    }
  });
}

extern "C" void TestDerivs_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDerivs_Set;
  DECLARE_CCTK_PARAMETERS;

  using vreal = Arith::simd<CCTK_REAL>;
  using vbool = Arith::simdl<CCTK_REAL>;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const vbool mask = Arith::mask_for_loop_tail<vbool>(p.i, p.imax);
        const vreal x0 = p.x + Arith::iota<vreal>() * p.dx;
        const CCTK_REAL y0 = p.y;
        const CCTK_REAL z0 = p.z;
        vreal u0;
        poly(kxx, kxy, kyz, Arith::cos(x0), std::sin(y0), std::sin(z0), u0);

        chi.store(mask, p.I, u0);
      });
}

extern "C" void TestDerivs_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDerivs_Sync;

  // do nothing
}

extern "C" void TestDerivs_CalcDerivs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestDerivs_CalcDerivs;
  DECLARE_CCTK_PARAMETERS;

  using Arith::vect, Arith::vec, Arith::smat;

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};
  vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 1: with ghost zones, suffix 0: without ghost zones
  const GF3D2layout layout2(cctkGH, indextype);
  const GF3D5layout layout5(imin, imax);

  const GF3D2<const CCTK_REAL> gf2_chi(layout2, chi);

  /* Define temp tail vars */
  constexpr int nvars = 10;
  GF3D5vector<CCTK_REAL> vars(layout5, nvars);

  int ivar = 0;
  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(vars(ivar++)); };
  const auto make_vec = [&](const auto &f) {
    return vec<result_of_t<decltype(f)()>, dim>([&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return smat<result_of_t<decltype(f)()>, dim>([&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };

  const GF3D5<CCTK_REAL> t5_chi(make_gf());
  const vec<GF3D5<CCTK_REAL>, dim> t5_dchi(make_vec_gf());
  const smat<GF3D5<CCTK_REAL>, dim> t5_ddchi(make_mat_gf());

  const GridDescBaseDevice grid(cctkGH);
  const vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  Derivs::calc_derivs2<0, 0, 0>(t5_chi, t5_dchi, t5_ddchi, layout5, grid,
                                gf2_chi, dx, deriv_order);
  ivar = -1;

  const vec<GF3D2<CCTK_REAL>, dim> gf_dchi{GF3D2<CCTK_REAL>(layout2, dxchi),
                                           GF3D2<CCTK_REAL>(layout2, dychi),
                                           GF3D2<CCTK_REAL>(layout2, dzchi)};
  const smat<GF3D2<CCTK_REAL>, dim> gf_ddchi{
      GF3D2<CCTK_REAL>(layout2, dxxchi), GF3D2<CCTK_REAL>(layout2, dxychi),
      GF3D2<CCTK_REAL>(layout2, dxzchi), GF3D2<CCTK_REAL>(layout2, dyychi),
      GF3D2<CCTK_REAL>(layout2, dyzchi), GF3D2<CCTK_REAL>(layout2, dzzchi)};

  typedef Arith::simd<CCTK_REAL> vreal;
  typedef Arith::simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = tuple_size_v<vreal>;
  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = Arith::mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index2(layout2, p.I);
        const GF3D5index index5(layout5, p.I);

        gf_dchi.store(mask, index2, t5_dchi(mask, index5));
        gf_ddchi.store(mask, index2, t5_ddchi(mask, index5));
      });

#if CCTK_DEBUG
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const auto sinx = std::sin(p.x);
        const auto siny = std::sin(p.y);
        const auto sinz = std::sin(p.z);
        const auto cosx = std::cos(p.x);
        const auto cosy = std::cos(p.y);
        const auto cosz = std::cos(p.z);
        const auto dxxdchi =
            -2 * kxx * cosx * cosx + 2 * kxx * sinx * sinx - kxy * cosx * siny;
        const auto dxydchi = -kxy * cosy * sinx;
        const auto dxzdchi = 0.0;
        const auto dyydchi = -kxy * cosx * siny - kyz * siny * sinz;
        const auto dyzdchi = kyz * cosy * cosz;
        const auto dzzdchi = -kyz * siny * sinz;

        const auto tiny = 5e-3;
        if (dxxdchi - gf_ddchi(0, 0)(p.I) > tiny ||
            dxydchi - gf_ddchi(0, 1)(p.I) > tiny ||
            dxzdchi - gf_ddchi(0, 2)(p.I) > tiny ||
            dyydchi - gf_ddchi(1, 1)(p.I) > tiny ||
            dyzdchi - gf_ddchi(1, 2)(p.I) > tiny ||
            dzzdchi - gf_ddchi(2, 2)(p.I) > tiny) {
#ifndef SYCL_LANGUAGE_VERSION
          printf("ddxx = %f\n", dxxdchi - gf_ddchi(0, 0)(p.I));
          printf("ddxy = %f\n", dxydchi - gf_ddchi(0, 1)(p.I));
          printf("ddxz = %f\n", dxzdchi - gf_ddchi(0, 2)(p.I));
          printf("ddyy = %f\n", dyydchi - gf_ddchi(1, 1)(p.I));
          printf("ddyz = %f\n", dyzdchi - gf_ddchi(1, 2)(p.I));
          printf("ddzz = %f\n", dzzdchi - gf_ddchi(2, 2)(p.I));
#endif
          assert(0);
        }
      });
#endif
}

} // namespace TestDerivs
