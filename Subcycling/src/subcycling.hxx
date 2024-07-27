#ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
#define CARPETX_SUBCYCLING_SUBCYCLING_HXX

#include <loop_device.hxx>

#include <sum.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

namespace Subcycling {
using namespace Arith;

/* copyed from CarpetX/src/driver.cxx */
inline array<int, Loop::dim> get_group_indextype(const int gi) {
  DECLARE_CCTK_PARAMETERS;

  assert(gi >= 0);

  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  array<CCTK_INT, Loop::dim> index;

  // The CST stage doesn't look for the `index` tag, and
  // `CCTK_ARGUMENTSX_...` would thus ignore it
  int iret = Util_TableGetIntArray(tags, Loop::dim, index.data(), "index");
  if (iret != UTIL_ERROR_TABLE_NO_SUCH_KEY)
    CCTK_VERROR(
        "The grid function group %s has a tag `index=...`. This is not "
        "supported any more; use a `CENTERING{...}` declaration instead.",
        CCTK_FullGroupName(gi));

  // Use the centering table
  const int centering = CCTK_GroupCenteringTableI(gi);
  assert(centering >= 0);
  iret = Util_TableGetIntArray(centering, Loop::dim, index.data(), "centering");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // Default: vertex-centred
    index = {0, 0, 0};
  } else if (iret >= 0) {
    assert(iret == Loop::dim);
  } else {
    assert(0);
  }

  // Convert to index type
  array<int, Loop::dim> indextype;
  for (int d = 0; d < Loop::dim; ++d)
    indextype[d] = index[d];

  return indextype;
}

/**
 * \brief Calculate Ys ghost points for fine grid using Ks on coarse grid
 *        using eq(39) of Peter McCorquodale and Phillip Colella,
 *        Applied Mathematics and Computational Science, 6(1):1-25, 2011.
 *
 * \param Yf        RK substage Ys on the fine side to be interperated into
 *                  the ghost zones
 * \param kcs       RK ks on the coarset side
 * \param u0        u at t0
 * \param dtc       Time step size on coarse side
 * \param xsi       which substep on fine level during a coarse time
 *                  step.  For an AMR simulation with subcycling and a
 *                  refinement ratio of 2, the number is either 0 or 0.5,
 *                  denoting the first and second substep, respectively.
 * \param stage     RK stage number starting from 1
 */
template <int RKSTAGES>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcYfFromKcs(const Loop::GridDescBaseDevice &grid,
              const array<int, Loop::dim> &indextype,
              const Loop::GF3D2<CCTK_REAL> &Yf,
              const Loop::GF3D2<const CCTK_REAL> &u0,
              array<const Loop::GF3D2<const CCTK_REAL>, RKSTAGES> &kcs,
              const Loop::GF3D2<const CCTK_REAL> &isrmbndry,
              const CCTK_REAL dtc, const CCTK_REAL xsi, const CCTK_INT stage) {
  assert(stage > 0 && stage <= 4);

  CCTK_REAL r = 0.5; // ratio between coarse and fine cell size (2 to 1 MR case)
  CCTK_REAL xsi2 = xsi * xsi;
  CCTK_REAL xsi3 = xsi2 * xsi;
  // coefficients for U
  CCTK_REAL b1 = xsi - CCTK_REAL(1.5) * xsi2 + CCTK_REAL(2. / 3.) * xsi3;
  CCTK_REAL b2 = xsi2 - CCTK_REAL(2. / 3.) * xsi3;
  CCTK_REAL b3 = b2;
  CCTK_REAL b4 = CCTK_REAL(-0.5) * xsi2 + CCTK_REAL(2. / 3.) * xsi3;
  // coefficients for Ut
  CCTK_REAL c1 = CCTK_REAL(1.) - CCTK_REAL(3.) * xsi + CCTK_REAL(2.) * xsi2;
  CCTK_REAL c2 = CCTK_REAL(2.) * xsi - CCTK_REAL(2.) * xsi2;
  CCTK_REAL c3 = c2;
  CCTK_REAL c4 = -xsi + CCTK_REAL(2.) * xsi2;
  // coefficients for Utt
  CCTK_REAL d1 = CCTK_REAL(-3.) + CCTK_REAL(4.) * xsi;
  CCTK_REAL d2 = CCTK_REAL(2.) - CCTK_REAL(4.) * xsi;
  CCTK_REAL d3 = d2;
  CCTK_REAL d4 = CCTK_REAL(-1.) + CCTK_REAL(4.) * xsi;
  // coefficients for Uttt
  constexpr CCTK_REAL e1 = CCTK_REAL(4.);
  constexpr CCTK_REAL e2 = CCTK_REAL(-4.);
  constexpr CCTK_REAL e3 = CCTK_REAL(-4.);
  constexpr CCTK_REAL e4 = CCTK_REAL(4.);

  if (stage == 1) {
    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            CCTK_REAL k1 = kcs[0](p.I);
            CCTK_REAL k2 = kcs[1](p.I);
            CCTK_REAL k3 = kcs[2](p.I);
            CCTK_REAL k4 = kcs[3](p.I);
            CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
            Yf(p.I) = u0(p.I) + dtc * uu;
          }
        });
  } else if (stage == 2) {
    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            CCTK_REAL k1 = kcs[0](p.I);
            CCTK_REAL k2 = kcs[1](p.I);
            CCTK_REAL k3 = kcs[2](p.I);
            CCTK_REAL k4 = kcs[3](p.I);
            CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
            CCTK_REAL ut = c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4;
            Yf(p.I) = u0(p.I) + dtc * (uu + CCTK_REAL(0.5) * r * ut);
          }
        });
  } else if (stage == 3 || stage == 4) {
    CCTK_REAL r2 = r * r;
    CCTK_REAL r3 = r2 * r;
    CCTK_REAL at = (stage == 3) ? CCTK_REAL(0.5) * r : r;
    CCTK_REAL att = (stage == 3) ? CCTK_REAL(0.25) * r2 : CCTK_REAL(0.5) * r2;
    CCTK_REAL attt =
        (stage == 3) ? CCTK_REAL(0.0625) * r3 : CCTK_REAL(0.125) * r3;
    CCTK_REAL ak = (stage == 3) ? CCTK_REAL(-4.) : CCTK_REAL(4.);

    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            CCTK_REAL k1 = kcs[0](p.I);
            CCTK_REAL k2 = kcs[1](p.I);
            CCTK_REAL k3 = kcs[2](p.I);
            CCTK_REAL k4 = kcs[3](p.I);
            CCTK_REAL uu = b1 * k1 + b2 * k2 + b3 * k3 + b4 * k4;
            CCTK_REAL ut = c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4;
            CCTK_REAL utt = d1 * k1 + d2 * k2 + d3 * k3 + d4 * k4;
            CCTK_REAL uttt = e1 * k1 + e2 * k2 + e3 * k3 + e4 * k4;
            Yf(p.I) = u0(p.I) + dtc * (uu + at * ut + att * utt +
                                       attt * (uttt + ak * (k3 - k2)));
          }
        });
  }
}

/* Varlist version */
template <int RKSTAGES>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcYfFromKcs(CCTK_ARGUMENTS, vector<int> &Yfs, vector<int> &u0s,
              const array<vector<int>, RKSTAGES> &kcss, const CCTK_REAL dtc,
              const CCTK_REAL xsi, const CCTK_INT stage) {

  const Loop::GridDescBaseDevice grid(cctkGH);
  const int tl = 0;
  // TODO: we need different centering types of flag for refinement boundary,
  // maybe make it a group
  const int isrmbndry_0 =
      CCTK_FirstVarIndexI(CCTK_GroupIndex("Subcycling::isrmbndry"));
  const Loop::GF3D2<const CCTK_REAL> isrmbndry(
      Loop::GF3D2layout(cctkGH, array<int, Loop::dim>{0, 0, 0}),
      static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, isrmbndry_0 + 0)));
  for (size_t i = 0; i < Yfs.size(); ++i) {
    const int nvars = CCTK_NumVarsInGroupI(Yfs[i]);
    const array<int, Loop::dim> indextype = get_group_indextype(Yfs[i]);
    const Loop::GF3D2layout layout(cctkGH, indextype);

    const int Yf_0 = CCTK_FirstVarIndexI(Yfs[i]);
    const int u0_0 = CCTK_FirstVarIndexI(u0s[i]);
    for (int vi = 0; vi < nvars; vi++) {
      const Loop::GF3D2<CCTK_REAL> Yf(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, Yf_0 + vi)));
      const Loop::GF3D2<const CCTK_REAL> u0(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, u0_0 + vi)));
      switch (RKSTAGES) {
      case 4: {
        array<const Loop::GF3D2<const CCTK_REAL>, 4> kcs{
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                            cctkGH, tl, CCTK_FirstVarIndexI(kcss[0][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                            cctkGH, tl, CCTK_FirstVarIndexI(kcss[1][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                            cctkGH, tl, CCTK_FirstVarIndexI(kcss[2][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout,
                static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, tl, CCTK_FirstVarIndexI(kcss[3][i]) + vi)))};
        CalcYfFromKcs<4>(grid, indextype, Yf, u0, kcs, isrmbndry, dtc, xsi,
                         stage);
        break;
      }
      default: {
        CCTK_ERROR("Unsupported RK stages with subcycling");
        break;
      }
      }
    }
  }
}

CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
SetK(const Loop::GridDescBaseDevice &grid,
     const array<int, Loop::dim> &indextype, const Loop::GF3D2<CCTK_REAL> &K,
     const Loop::GF3D2<const CCTK_REAL> &rhs) {
  grid.loop_device_idx<Loop::where_t::interior>(
      indextype, grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { K(p.I) = rhs(p.I); });
}

/* Varlist version */
template <int RKSTAGES>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
SetK(CCTK_ARGUMENTS, const array<vector<int>, RKSTAGES> &kss, vector<int> &rhss,
     const CCTK_INT stage) {
  assert(stage > 0 && stage <= 4);
  const Loop::GridDescBaseDevice grid(cctkGH);
  const int tl = 0;
  for (size_t i = 0; i < rhss.size(); ++i) {
    const int nvars = CCTK_NumVarsInGroupI(rhss[i]);
    const array<int, Loop::dim> indextype = get_group_indextype(rhss[i]);
    const Loop::GF3D2layout layout(cctkGH, indextype);
    const int rhs_0 = CCTK_FirstVarIndexI(rhss[i]);
    const int K_0 = CCTK_FirstVarIndexI(kss[stage - 1][i]);
    for (int vi = 0; vi < nvars; vi++) {
      const Loop::GF3D2<CCTK_REAL> K(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, K_0 + vi)));
      const Loop::GF3D2<const CCTK_REAL> rhs(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, rhs_0 + vi)));
      SetK(grid, indextype, K, rhs);
    }
  }
}

} // namespace Subcycling

#endif // #ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
