#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>


// these are all wrong except that the memory access patterns should be about
// correct

namespace RegCheck {

// using wrappers
template <typename T, int d>
CCTK_DEVICE inline T calc_1st_deriv(const Loop::GF3D2<T> &gf,
                                    const Loop::PointDesc &p) {
  auto varm = gf(p.I - p.DI[d]);
  auto varp = gf(p.I + p.DI[d]);
  return (varp - varm)/(T(2)*p.dx);
}

template <typename T, int d>
CCTK_DEVICE inline T calc_2nd_deriv(const Loop::GF3D2<T> &gf,
                                    const Loop::PointDesc &p) {
  auto varm = gf(p.I - p.DI[d]);
  auto var0 = gf(p.I);
  auto varp = gf(p.I + p.DI[d]);
  return (varp - T(2)*var0 +  varm)/(p.dx*p.dx);
}

template <typename T, int d1, int d2>
CCTK_DEVICE inline T calc_2nd_deriv(const Loop::GF3D2<T> &gf,
                                    const Loop::PointDesc &p) {
  auto varm = (gf(p.I - p.DI[d1] + p.DI[d2]) - gf(p.I - p.DI[d1] - p.DI[d2]))/(T(2)*p.dx);
  auto varp = (gf(p.I + p.DI[d1] + p.DI[d2]) - gf(p.I + p.DI[d1] - p.DI[d2]))/(T(2)*p.dx);
  return ((varp - varm)/(T(2)*p.dx));
}

// for vectors
template <typename T, int d>
CCTK_DEVICE inline T calc_1st_deriv(const Loop::GF3D2<T> &gf,
                                    const Loop::PointDesc &p,
				    const int n) {
  auto varm = gf.ptr[n*gf.layout.np + gf.layout.linear(p.I - p.DI[d])];
  auto varp = gf.ptr[n*gf.layout.np + gf.layout.linear(p.I + p.DI[d])];
  return (varp - varm)/(T(2)*p.dx);
}

extern "C" void RegCheck_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RegCheck_Init;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        var1(p.I) = CCTK_REAL(0);
        var2(p.I) = CCTK_REAL(0);
        var3(p.I) = CCTK_REAL(0);
        var4(p.I) = CCTK_REAL(0);
        var5(p.I) = CCTK_REAL(0);
        var6(p.I) = CCTK_REAL(0);
        var7(p.I) = CCTK_REAL(0);
        var8(p.I) = CCTK_REAL(0);
        var9(p.I) = CCTK_REAL(0);
        var10(p.I) = CCTK_REAL(0);
        var11(p.I) = CCTK_REAL(0);
        var12(p.I) = CCTK_REAL(0);
        var13(p.I) = CCTK_REAL(0);
        var14(p.I) = CCTK_REAL(0);
        var15(p.I) = CCTK_REAL(0);
        var16(p.I) = CCTK_REAL(0);
        var17(p.I) = CCTK_REAL(0);
        var18(p.I) = CCTK_REAL(0);
        var19(p.I) = CCTK_REAL(0);
        var20(p.I) = CCTK_REAL(0);
        var21(p.I) = CCTK_REAL(0);
        var22(p.I) = CCTK_REAL(0);
        var23(p.I) = CCTK_REAL(0);
        var24(p.I) = CCTK_REAL(0);
        var25(p.I) = CCTK_REAL(0);
        var26(p.I) = CCTK_REAL(0);
        var27(p.I) = CCTK_REAL(0);
        var28(p.I) = CCTK_REAL(0);
	for(int n = 0 ; n < 28 ; n++) {
	  var.ptr[n*var.layout.np + var.layout.linear(p.I)] = CCTK_REAL(0);
	}
  });
}

extern "C" void RegCheck_DummyCalc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RegCheck_DummyCalc;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        rhsvar1(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var1, p); 
        rhsvar2(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var2, p); 
        rhsvar3(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var3, p); 
        rhsvar4(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var4, p); 
        rhsvar5(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var5, p); 
        rhsvar6(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var6, p); 
        rhsvar7(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var7, p); 
        rhsvar8(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var8, p); 
        rhsvar9(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var9, p); 
        rhsvar10(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var10, p); 
        rhsvar11(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var11, p); 
        rhsvar12(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var12, p); 
        rhsvar13(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var13, p); 
        rhsvar14(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var14, p); 
        rhsvar15(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var15, p); 
        rhsvar16(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var16, p); 
        rhsvar17(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var17, p); 
        rhsvar18(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var18, p); 
        rhsvar19(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var19, p); 
        rhsvar20(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var20, p); 
        rhsvar21(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var21, p); 
        rhsvar22(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var22, p); 
        rhsvar23(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var23, p); 
        rhsvar24(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var24, p); 
        rhsvar25(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var25, p); 
        rhsvar26(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var26, p); 
        rhsvar27(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var27, p); 
        rhsvar28(p.I) = calc_1st_deriv<CCTK_REAL, 0>(var28, p); 
  });
}

extern "C" void RegCheck_DummyCalcVect(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_RegCheck_DummyCalcVect;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
	for(int n = 0 ; n < 28 ; n++) {
	  rhsvar.ptr[n*var.layout.np + var.layout.linear(p.I)] = 
            calc_1st_deriv<CCTK_REAL, 0>(var, p, n); 
	}
  });
}

} // namespace RegCheck
