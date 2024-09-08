#include <loop_device.hxx>

#include <defs.hxx>
#include <derivs.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace TestDerivs {
using namespace Loop;

template <typename T, typename VT>
constexpr VT poly(const T kxx, const T kxy, const T kyz, const VT x, const T y,
                  const T z) {
  return kxx * x * x + kxy * x * y + kyz * y * z;
}

template <typename T>
constexpr void poly_derivs(const T kxx, const T kxy, const T kyz, const T x,
                           const T y, const T z, Arith::vec<T, dim> &du,
                           Arith::smat<T, dim> &ddu) {
  const T sinx = std::sin(x);
  const T siny = std::sin(y);
  const T sinz = std::sin(z);
  const T cosx = std::cos(x);
  const T cosy = std::cos(y);
  const T cosz = std::cos(z);
  du(0) = -2 * kxx * cosx * sinx - kxy * sinx * siny;
  du(1) = kxy * cosx * cosy + kyz * cosy * sinz;
  du(2) = kyz * cosz * siny;
  ddu(0, 0) =
      -2 * kxx * cosx * cosx + 2 * kxx * sinx * sinx - kxy * cosx * siny;
  ddu(0, 1) = -kxy * cosy * sinx;
  ddu(0, 2) = 0.0;
  ddu(1, 1) = -kxy * cosx * siny - kyz * siny * sinz;
  ddu(1, 2) = kyz * cosy * cosz;
  ddu(2, 2) = -kyz * siny * sinz;
}

template <typename T>
constexpr Arith::vec<T, dim> beta_profile(const T x, const T y, const T z) {
  const T vx = 1 * std::sin(x);
  const T vy = 2 * std::sin(y);
  const T vz = 3 * std::sin(z);
  const T norm = std::sqrt(vx * vx + vy * vy + vz * vz);
  return Arith::vec<T, dim>({vx / norm, vy / norm, vy / norm});
}

template <typename T>
constexpr T poly_diss(const T kxx, const T kxy, const T kyz, const T x,
                      const T y, const T z, const int diss_order,
                      const Arith::vect<T, dim> &dx) {
  const T sinx = std::sin(x);
  const T siny = std::sin(y);
  const T sinz = std::sin(z);
  const T cosx = std::cos(x);
  int coeff = 0;

  switch (diss_order) {
  case 4:
    coeff = 8;
    break;
  case 6:
    coeff = -32;
    break;
  case 8:
    coeff = 128;
    break;
  case 10:
    coeff = -512;
    break;
  default:
    assert(0);
  }

  const int abc = std::abs(coeff);
  const int sig = (abc == coeff) ? 1 : -1;
  return sig * ((abc * kxx * cosx * cosx - abc * kxx * sinx * sinx +
                 kxy * cosx * siny) /
                    dx[0] +
                (kxy * cosx * siny + kyz * siny * sinz) / dx[1] +
                (kyz * siny * sinz) / dx[2]);
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
        chi.store(
            mask, p.I,
            poly(kxx, kxy, kyz, cos(x0), std::sin(y0), std::sin(z0)));
      });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const Arith::vec<CCTK_REAL, dim> betas = beta_profile(p.x, p.y, p.z);
        betax(p.I) = betas(0);
        betay(p.I) = betas(1);
        betaz(p.I) = betas(2);
      });
}

extern "C" void TestDerivs_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDerivs_Sync;
  // do nothing
}

extern "C" void TestDerivs_CalcDerivs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestDerivs_CalcDerivs;
  DECLARE_CCTK_PARAMETERS;

  const std::array<int, dim> indextype = {0, 0, 0};
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  Arith::vect<int, dim> imin, imax;
  GridDescBase(cctkGH).box_int<0, 0, 0>(nghostzones, imin, imax);
  // Suffix 2: with ghost zones, suffix 5: without ghost zones
  const GF3D2layout layout2(cctkGH, indextype);
  const GF3D5layout layout5(imin, imax);
  const GF3D2<const CCTK_REAL> gf2_chi(layout2, chi);

  /* Define temp tail vars */
  constexpr int nvars = 10;
  GF3D5vector<CCTK_REAL> vars(layout5, nvars);

  int ivar = 0;
  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(vars(ivar++)); };
  const auto make_vec = [&](const auto &f) {
    return Arith::vec<std::result_of_t<decltype(f)()>, dim>(
        [&](int) { return f(); });
  };
  const auto make_mat = [&](const auto &f) {
    return Arith::smat<std::result_of_t<decltype(f)()>, dim>(
        [&](int, int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };
  const auto make_mat_gf = [&]() { return make_mat(make_gf); };

  const GF3D5<CCTK_REAL> t5_chi(make_gf());
  const Arith::vec<GF3D5<CCTK_REAL>, dim> t5_dchi(make_vec_gf());
  const Arith::smat<GF3D5<CCTK_REAL>, dim> t5_ddchi(make_mat_gf());

  const GridDescBaseDevice grid(cctkGH);
  const Arith::vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  Derivs::calc_derivs2<0, 0, 0>(t5_chi, t5_dchi, t5_ddchi, layout5, grid,
                                gf2_chi, dx, deriv_order);
  ivar = -1;

  const Arith::vec<GF3D2<CCTK_REAL>, dim> gf_dchi{
      GF3D2<CCTK_REAL>(layout2, dxchi), GF3D2<CCTK_REAL>(layout2, dychi),
      GF3D2<CCTK_REAL>(layout2, dzchi)};
  const Arith::smat<GF3D2<CCTK_REAL>, dim> gf_ddchi{
      GF3D2<CCTK_REAL>(layout2, dxxchi), GF3D2<CCTK_REAL>(layout2, dxychi),
      GF3D2<CCTK_REAL>(layout2, dxzchi), GF3D2<CCTK_REAL>(layout2, dyychi),
      GF3D2<CCTK_REAL>(layout2, dyzchi), GF3D2<CCTK_REAL>(layout2, dzzchi)};
  const GF3D2<CCTK_REAL> gf_chi_diss(layout2, chi_diss);

  const Arith::vec<GF3D2<const CCTK_REAL>, 3> gf_beta{
      GF3D2<const CCTK_REAL>(layout2, betax),
      GF3D2<const CCTK_REAL>(layout2, betay),
      GF3D2<const CCTK_REAL>(layout2, betaz)};
  const GF3D2<CCTK_REAL> gf_chi_upwind(layout2, chi_upwind);

  typedef Arith::simd<CCTK_REAL> vreal;
  typedef Arith::simdl<CCTK_REAL> vbool;
  constexpr size_t vsize = std::tuple_size_v<vreal>;

  vreal (*calc_deriv_upwind)(const GF3D2<const CCTK_REAL> &, const vbool &,
                             const Arith::vect<int, dim> &,
                             const Arith::vect<CCTK_REAL, dim> &,
                             const Arith::vec<vreal, dim> &);
  vreal (*calc_diss)(const GF3D2<const CCTK_REAL> &, const vbool &,
                     const Arith::vect<int, dim> &,
                     const Arith::vect<CCTK_REAL, dim> &);

  switch (deriv_order) {
  case 2: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<2>;
    calc_diss = &Derivs::calc_diss<2>;
    break;
  }
  case 4: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<4>;
    calc_diss = &Derivs::calc_diss<4>;
    break;
  }
  case 6: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<6>;
    calc_diss = &Derivs::calc_diss<6>;
    break;
  }
  case 8: {
    calc_deriv_upwind = &Derivs::calc_deriv_upwind<8>;
    calc_diss = &Derivs::calc_diss<8>;
    break;
  }
  default:
    assert(0);
  }

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = Arith::mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D2index index2(layout2, p.I);
        const GF3D5index index5(layout5, p.I);
        gf_dchi.store(mask, index2, t5_dchi(mask, index5));
        gf_ddchi.store(mask, index2, t5_ddchi(mask, index5));
        gf_chi_diss.store(mask, index2, calc_diss(gf2_chi, mask, p.I, dx));
        gf_chi_upwind.store(
            mask, index2,
            calc_deriv_upwind(gf2_chi, mask, p.I, dx, gf_beta(mask, index2)));
      });

#if CCTK_DEBUG
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        Arith::vec<CCTK_REAL, dim> dchi;
        Arith::smat<CCTK_REAL, dim> ddchi;
        poly_derivs(kxx, kxy, kyz, p.x, p.y, p.z, dchi, ddchi);

        const auto tiny = 5e-3;
        if (dchi(0) - gf_dchi(0)(p.I) > tiny ||
            dchi(1) - gf_dchi(1)(p.I) > tiny ||
            dchi(2) - gf_dchi(2)(p.I) > tiny ||
            ddchi(0, 0) - gf_ddchi(0, 0)(p.I) > tiny ||
            ddchi(0, 1) - gf_ddchi(0, 1)(p.I) > tiny ||
            ddchi(0, 2) - gf_ddchi(0, 2)(p.I) > tiny ||
            ddchi(1, 1) - gf_ddchi(1, 1)(p.I) > tiny ||
            ddchi(1, 2) - gf_ddchi(1, 2)(p.I) > tiny ||
            ddchi(2, 2) - gf_ddchi(2, 2)(p.I) > tiny) {
#ifndef SYCL_LANGUAGE_VERSION
          printf("dx = %f\n", dchi(0) - gf_dchi(0)(p.I));
          printf("dy = %f\n", dchi(1) - gf_dchi(1)(p.I));
          printf("dz = %f\n", dchi(2) - gf_dchi(2)(p.I));
          printf("ddxx = %f\n", ddchi(0, 0) - gf_ddchi(0, 0)(p.I));
          printf("ddxy = %f\n", ddchi(0, 1) - gf_ddchi(0, 1)(p.I));
          printf("ddxz = %f\n", ddchi(0, 2) - gf_ddchi(0, 2)(p.I));
          printf("ddyy = %f\n", ddchi(1, 1) - gf_ddchi(1, 1)(p.I));
          printf("ddyz = %f\n", ddchi(1, 2) - gf_ddchi(1, 2)(p.I));
          printf("ddzz = %f\n", ddchi(2, 2) - gf_ddchi(2, 2)(p.I));
#endif
          assert(0);
        }
      });
#endif
}

extern "C" void TestDerivs_CalcError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestDerivs_CalcError;
  DECLARE_CCTK_PARAMETERS;

  const Arith::vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones, [=] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        Arith::vec<CCTK_REAL, dim> dchi;
        Arith::smat<CCTK_REAL, dim> ddchi;
        poly_derivs(kxx, kxy, kyz, p.x, p.y, p.z, dchi, ddchi);

        dxchi_error(p.I) = dxchi(p.I) - dchi(0);
        dychi_error(p.I) = dychi(p.I) - dchi(1);
        dzchi_error(p.I) = dzchi(p.I) - dchi(2);
        dxxchi_error(p.I) = dxxchi(p.I) - ddchi(0, 0);
        dxychi_error(p.I) = dxychi(p.I) - ddchi(0, 1);
        dxzchi_error(p.I) = dxzchi(p.I) - ddchi(0, 2);
        dyychi_error(p.I) = dyychi(p.I) - ddchi(1, 1);
        dyzchi_error(p.I) = dyzchi(p.I) - ddchi(1, 2);
        dzzchi_error(p.I) = dzzchi(p.I) - ddchi(2, 2);

        const int diss_order = deriv_order + 2;
        const CCTK_REAL diss =
            Arith::pown(-1, diss_order / 2.0 - 1.0) /
            Arith::pown(2, diss_order) *
            poly_diss(kxx, kxy, kyz, p.x, p.y, p.z, diss_order, dx);
        chi_diss_error(p.I) = diss - chi_diss(p.I);

        const Arith::vec<CCTK_REAL, dim> betas = beta_profile(p.x, p.y, p.z);
        const CCTK_REAL upwind = Arith::sum<dim>(
            [&](int i) ARITH_INLINE { return betas(i) * dchi(i); });
        chi_upwind_error(p.I) = upwind - chi_upwind(p.I);
      });
}

} // namespace TestDerivs
