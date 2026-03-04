#ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
#define CARPETX_SUBCYCLING_SUBCYCLING_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>

#include "utils.hxx"

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
 * \brief Compute fine-grid ghost points for Ys using the prolongated Ks from
 *        the coarse-grid
 *
 * \param Yf        RK substage values (Ys) at the fine-grid ghost points.
 * \param kcs       RK substage derivatives (Ks) prolongated from the
 *                  coarse-grid to the fine-grid ghost points.
 * \param u0        Field values at the initial time t0.
 * \param dtc       Time step size on the coarse-grid.
 * \param xsi       Substep position within a coarse time step on the fine grid.
 *                  In an AMR simulation with subcycling and a refinement ratio
 *                  of 2, this value is either 0 (first substep) or 0.5 (second
 *                  substep).
 * \param stage     RK stage number (starting from 1).
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

  constexpr CCTK_REAL r =
      0.5; // ratio between coarse and fine cell size (2 to 1 MR case)

  const CCTK_REAL xsi2 = xsi * xsi;
  const CCTK_REAL xsi3 = xsi2 * xsi;

  // Coefficients for the dense output formulas (U, Ut, Utt, Uttt)
  const array<CCTK_REAL, RKSTAGES> b = {
      xsi - 1.5 * xsi2 + (2. / 3.) * xsi3, // b1
      xsi2 - (2. / 3.) * xsi3,             // b2
      xsi2 - (2. / 3.) * xsi3,             // b3
      -0.5 * xsi2 + (2. / 3.) * xsi3       // b4
  };

  const array<CCTK_REAL, RKSTAGES> bt = {
      1.0 - 3.0 * xsi + 2.0 * xsi2, // bt1
      2.0 * xsi - 2.0 * xsi2,       // bt2
      2.0 * xsi - 2.0 * xsi2,       // bt3
      -xsi + 2.0 * xsi2             // bt4
  };

  const array<CCTK_REAL, RKSTAGES> btt = {
      -3.0 + 4.0 * xsi, // btt1
      2.0 - 4.0 * xsi,  // btt2
      2.0 - 4.0 * xsi,  // btt3
      -1.0 + 4.0 * xsi  // btt4
  };

  constexpr array<CCTK_REAL, RKSTAGES> bttt = {
      4.0,  // bttt1
      -4.0, // bttt2
      -4.0, // bttt3
      4.0   // bttt4
  };

  // Create and launch the appropriate lambda based on the stage.
  if (stage == 1) {
    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            const array<CCTK_REAL, RKSTAGES> k = {kcs[0](p.I), kcs[1](p.I),
                                                  kcs[2](p.I), kcs[3](p.I)};
            const CCTK_REAL uu =
                b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];
            Yf(p.I) = u0(p.I) + dtc * uu;
          }
        });
  } else if (stage == 2) {
    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            const array<CCTK_REAL, RKSTAGES> k = {kcs[0](p.I), kcs[1](p.I),
                                                  kcs[2](p.I), kcs[3](p.I)};
            const CCTK_REAL uu =
                b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];
            const CCTK_REAL ut =
                bt[0] * k[0] + bt[1] * k[1] + bt[2] * k[2] + bt[3] * k[3];
            Yf(p.I) = u0(p.I) + dtc * (uu + 0.5 * r * ut);
          }
        });
  } else { // stage 3 or stage 4
    const CCTK_REAL r2 = r * r;
    const CCTK_REAL r3 = r2 * r;
    const CCTK_REAL at = (stage == 3) ? 0.5 * r : r;
    const CCTK_REAL att = (stage == 3) ? 0.25 * r2 : 0.5 * r2;
    const CCTK_REAL attt = (stage == 3) ? 0.0625 * r3 : 0.125 * r3;
    const CCTK_REAL ak = (stage == 3) ? -4.0 : 4.0;

    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (isrmbndry(p.I)) {
            const array<CCTK_REAL, RKSTAGES> k = {kcs[0](p.I), kcs[1](p.I),
                                                  kcs[2](p.I), kcs[3](p.I)};
            const CCTK_REAL uu =
                b[0] * k[0] + b[1] * k[1] + b[2] * k[2] + b[3] * k[3];
            const CCTK_REAL ut =
                bt[0] * k[0] + bt[1] * k[1] + bt[2] * k[2] + bt[3] * k[3];
            const CCTK_REAL utt =
                btt[0] * k[0] + btt[1] * k[1] + btt[2] * k[2] + btt[3] * k[3];
            const CCTK_REAL uttt = bttt[0] * k[0] + bttt[1] * k[1] +
                                   bttt[2] * k[2] + bttt[3] * k[3];
            Yf(p.I) = u0(p.I) + dtc * (uu + at * ut + att * utt +
                                       attt * (uttt + ak * (k[2] - k[1])));
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
  const auto isrmbndry_map = construct_index_map();
  const int tl = 0;

  // Helper to find the isrmbndry gf for a given indextype
  auto get_isrmbndry =
      [&](const array<int, Loop::dim> &indextype,
          const Loop::GF3D2layout &layout) -> Loop::GF3D2<const CCTK_REAL> {
    auto it = isrmbndry_map.find(indextype);
    if (it != isrmbndry_map.end()) {
      return {layout, static_cast<CCTK_REAL *>(
                          CCTK_VarDataPtrI(cctkGH, tl, it->second))};
    }
    CCTK_ERROR("get_isrmbndry: No matching indextype found.");
  };

  for (size_t i = 0; i < Yfs.size(); ++i) {
    const int nvars = CCTK_NumVarsInGroupI(Yfs[i]);
    const array<int, Loop::dim> indextype = get_group_indextype(Yfs[i]);
    const Loop::GF3D2layout layout(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> isrmbndry =
        get_isrmbndry(indextype, layout);

    const int Yf_0 = CCTK_FirstVarIndexI(Yfs[i]);
    const int u0_0 = CCTK_FirstVarIndexI(u0s[i]);
    for (int vi = 0; vi < nvars; ++vi) {
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
SetStateInterior(const Loop::GridDescBaseDevice &grid,
                 const array<int, Loop::dim> &indextype,
                 const Loop::GF3D2<CCTK_REAL> &u,
                 const Loop::GF3D2<const CCTK_REAL> &var) {
  grid.loop_device_idx<Loop::where_t::interior>(
      indextype, grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { u(p.I) = var(p.I); });
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
    const int k_0 = CCTK_FirstVarIndexI(kss[stage - 1][i]);
    for (int vi = 0; vi < nvars; ++vi) {
      const Loop::GF3D2<CCTK_REAL> K(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, k_0 + vi)));
      const Loop::GF3D2<const CCTK_REAL> rhs(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, rhs_0 + vi)));
      SetStateInterior(grid, indextype, K, rhs);
    }
  }
}

/* Varlist version */
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
SetOld(CCTK_ARGUMENTS, vector<int> &olds, vector<int> &vars) {
  const Loop::GridDescBaseDevice grid(cctkGH);
  const int tl = 0;

  // Loop over groups
  for (size_t i = 0; i < vars.size(); ++i) {
    const int nvars = CCTK_NumVarsInGroupI(vars[i]);
    const array<int, Loop::dim> indextype = get_group_indextype(vars[i]);
    const Loop::GF3D2layout layout(cctkGH, indextype);

    const int var_0 = CCTK_FirstVarIndexI(vars[i]);
    const int old_0 = CCTK_FirstVarIndexI(olds[i]);
    for (int vi = 0; vi < nvars; ++vi) {
      const Loop::GF3D2<CCTK_REAL> old(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, old_0 + vi)));
      const Loop::GF3D2<const CCTK_REAL> var(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, var_0 + vi)));
      SetStateInterior(grid, indextype, old, var);
    }
  }
}

} // namespace Subcycling

#endif // #ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
