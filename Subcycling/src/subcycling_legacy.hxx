#ifndef CARPETX_SUBCYCLING_SUBCYCLING_LEGACY_HXX
#define CARPETX_SUBCYCLING_SUBCYCLING_LEGACY_HXX

#include <loop_device.hxx>

#include "CarpetX/CarpetX/src/driver.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <AMReX_Array4.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace Subcycling {
using namespace Arith;
using std::array;
using std::vector;

/**
 * \brief Compute fine-grid ghost points for Ys using the prolongated Ks from
 *        the coarse-grid
 *
 * \param Yf        RK substage values (Ys) at the fine-grid ghost points.
 * \param kcs       RK substage derivatives (Ks) prolongated from the
 *                  coarse-grid to the fine-grid ghost points.
 * \param u0        Field values at the initial time t0.
 * \param cf_mask   Coarse-fine ghost mask (truthy at CF ghosts).
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
              const amrex::Array4<const int> &cf_mask, const CCTK_REAL dtc,
              const CCTK_REAL xsi, const CCTK_INT stage) {
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

  // Cactus PointDesc indices are local to the FAB (start at 0 for the first
  // allocated cell). The Array4 cf_mask uses global indices;
  // amrex::lbound(cf_mask) holds the global index of that same first
  // allocated cell.
  const amrex::Dim3 cf_lo = amrex::lbound(cf_mask);
  const int cf_ox = cf_lo.x;
  const int cf_oy = cf_lo.y;
  const int cf_oz = cf_lo.z;

  // Create and launch the appropriate lambda based on the stage.
  if (stage == 1) {
    grid.loop_device_idx<Loop::where_t::ghosts>(
        indextype, grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (cf_mask(p.i + cf_ox, p.j + cf_oy, p.k + cf_oz)) {
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
          if (cf_mask(p.i + cf_ox, p.j + cf_oy, p.k + cf_oz)) {
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
          if (cf_mask(p.i + cf_ox, p.j + cf_oy, p.k + cf_oz)) {
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

namespace detail {

/* Shared varlist implementation: Yf written at Yf_tl, u0 read at u0_tl.
   kcs are scratch groups read at tl=0. The coarse-fine ghost mask is fetched
   per-group from CarpetX::LevelData::get_cf_mask. */
template <int RKSTAGES>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void CalcYfFromKcs_impl(
    CCTK_ARGUMENTS, vector<int> &Yfs, const int Yf_tl, vector<int> &u0s,
    const int u0_tl, const array<vector<int>, RKSTAGES> &kcss,
    const CCTK_REAL dtc, const CCTK_REAL xsi, const CCTK_INT stage) {

  const Loop::GridDescBaseDevice grid(cctkGH);

  const int patch = cctkGH->cctk_patch;
  const int level = cctkGH->cctk_level;
  const int component = cctkGH->cctk_component;
  const auto &leveldata =
      CarpetX::ghext->patchdata.at(patch).leveldata.at(level);

  for (size_t i = 0; i < Yfs.size(); ++i) {
    const int nvars = CCTK_NumVarsInGroupI(Yfs[i]);
    const auto &groupdata = *leveldata.groupdata.at(Yfs[i]);
    const array<int, Loop::dim> &indextype = groupdata.indextype;
    const Loop::GF3D2layout layout(cctkGH, indextype);

    amrex::iMultiFab *const cf_mfab =
        leveldata.get_cf_mask(groupdata.indextype, groupdata.nghostzones);
    // No coarse-fine ghosts to fill at level 0 or when subcycling is disabled.
    if (cf_mfab == nullptr)
      continue;
    const amrex::Array4<const int> cf_mask = cf_mfab->const_array(component);

    const int Yf_0 = CCTK_FirstVarIndexI(Yfs[i]);
    const int u0_0 = CCTK_FirstVarIndexI(u0s[i]);
    for (int vi = 0; vi < nvars; ++vi) {
      const Loop::GF3D2<CCTK_REAL> Yf(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, Yf_tl, Yf_0 + vi)));
      const Loop::GF3D2<const CCTK_REAL> u0(
          layout,
          static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, u0_tl, u0_0 + vi)));
      switch (RKSTAGES) {
      case 4: {
        array<const Loop::GF3D2<const CCTK_REAL>, 4> kcs{
            Loop::GF3D2<const CCTK_REAL>(
                layout,
                static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, /*tl=*/0, CCTK_FirstVarIndexI(kcss[0][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout,
                static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, /*tl=*/0, CCTK_FirstVarIndexI(kcss[1][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout,
                static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, /*tl=*/0, CCTK_FirstVarIndexI(kcss[2][i]) + vi))),
            Loop::GF3D2<const CCTK_REAL>(
                layout,
                static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                    cctkGH, /*tl=*/0, CCTK_FirstVarIndexI(kcss[3][i]) + vi)))};
        CalcYfFromKcs<4>(grid, indextype, Yf, u0, kcs, cf_mask, dtc, xsi,
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

} // namespace detail

/* Varlist version (legacy two-group form): Yf and u0 are separate groups,
   both at tl=0. Used by TestSubcyclingMC. */
template <int RKSTAGES>
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
CalcYfFromKcs(CCTK_ARGUMENTS, vector<int> &Yfs, vector<int> &u0s,
              const array<vector<int>, RKSTAGES> &kcss, const CCTK_REAL dtc,
              const CCTK_REAL xsi, const CCTK_INT stage) {
  detail::CalcYfFromKcs_impl<RKSTAGES>(CCTK_PASS_CTOC, Yfs, /*Yf_tl=*/0, u0s,
                                       /*u0_tl=*/0, kcss, dtc, xsi, stage);
}

} // namespace Subcycling

#endif // #ifndef CARPETX_SUBCYCLING_SUBCYCLING_LEGACY_HXX
