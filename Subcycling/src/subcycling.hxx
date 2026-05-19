#ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
#define CARPETX_SUBCYCLING_SUBCYCLING_HXX

#include "CarpetX/CarpetX/src/driver.hxx"

#include <cctk.h>

#include <AMReX_FabArray.H>      // MultiArray4, MultiFab::arrays()
#include <AMReX_GpuContainers.H> // amrex::GpuArray
#include <AMReX_IntVect.H>       // amrex::IntVect
#include <AMReX_MFParallelFor.H> // amrex::ParallelFor(MF, IntVect, ncomp, F)
#include <AMReX_MultiFab.H>      // amrex::MultiFab

#include <array>
#include <cassert>
#include <vector>

namespace Subcycling {

/**
 * \brief MF-level fused entry point. Dispatches one amrex::ParallelFor per
 *        evolved group per level, fused across all local boxes and all
 *        components of the group.
 *
 * \param leveldata Per-level data (provides groupdata, get_cf_mask).
 * \param Yfs       Cactus group indices receiving the dense-output value.
 * \param Yf_tl     Time level of Yfs to write into.
 * \param u0s       Cactus group indices providing u(t0).
 * \param u0_tl     Time level of u0s to read from.
 * \param kcss      Per-stage Cactus group indices for RK substage derivatives.
 * \param dtc       Coarse-grid time step size.
 * \param xsi       Substep position within the coarse step (0 or 0.5).
 * \param stage     RK stage number (1..RKSTAGES).
 */
template <int RKSTAGES>
CCTK_HOST inline void CalcYfFromKcs_MFlevel(
    CarpetX::GHExt::PatchData::LevelData &leveldata,
    const std::vector<int> &Yfs, const int Yf_tl, const std::vector<int> &u0s,
    const int u0_tl, const std::array<std::vector<int>, RKSTAGES> &kcss,
    const CCTK_REAL dtc, const CCTK_REAL xsi, const CCTK_INT stage) {
  static_assert(RKSTAGES == 4,
                "CalcYfFromKcs_MFlevel only supports RKSTAGES == 4");
  assert(stage > 0 && stage <= 4);

  constexpr CCTK_REAL r =
      0.5; // ratio between coarse and fine cell size (2 to 1 MR case)

  const CCTK_REAL xsi2 = xsi * xsi;
  const CCTK_REAL xsi3 = xsi2 * xsi;

  // Coefficients for the dense output formulas (U, Ut, Utt, Uttt)
  const std::array<CCTK_REAL, 4> b = {
      xsi - 1.5 * xsi2 + (2. / 3.) * xsi3, // b1
      xsi2 - (2. / 3.) * xsi3,             // b2
      xsi2 - (2. / 3.) * xsi3,             // b3
      -0.5 * xsi2 + (2. / 3.) * xsi3       // b4
  };
  const std::array<CCTK_REAL, 4> bt = {
      1.0 - 3.0 * xsi + 2.0 * xsi2, // bt1
      2.0 * xsi - 2.0 * xsi2,       // bt2
      2.0 * xsi - 2.0 * xsi2,       // bt3
      -xsi + 2.0 * xsi2             // bt4
  };
  const std::array<CCTK_REAL, 4> btt = {
      -3.0 + 4.0 * xsi, // btt1
      2.0 - 4.0 * xsi,  // btt2
      2.0 - 4.0 * xsi,  // btt3
      -1.0 + 4.0 * xsi  // btt4
  };
  constexpr std::array<CCTK_REAL, 4> bttt = {
      4.0,  // bttt1
      -4.0, // bttt2
      -4.0, // bttt3
      4.0   // bttt4
  };

  for (size_t i = 0; i < Yfs.size(); ++i) {
    const auto &groupdata = *leveldata.groupdata.at(Yfs[i]);
    const auto &u0_groupdata = *leveldata.groupdata.at(u0s[i]);
    const int nvars = groupdata.numvars;
    assert(u0_groupdata.numvars == nvars);

    amrex::iMultiFab *const cf_mfab =
        leveldata.get_cf_mask(groupdata.indextype, groupdata.nghostzones);
    // No coarse-fine ghosts to fill at level 0 or when subcycling is disabled.
    if (cf_mfab == nullptr)
      continue;

    amrex::MultiFab &Yf_mf = *groupdata.mfab.at(Yf_tl);
    const amrex::MultiFab &u0_mf = *u0_groupdata.mfab.at(u0_tl);

    auto Yf_arrs = Yf_mf.arrays();
    auto u0_arrs = u0_mf.const_arrays();
    auto cf_arrs = cf_mfab->const_arrays();

    amrex::GpuArray<amrex::MultiArray4<const CCTK_REAL>, 4> kcs_arrs;
    for (int s = 0; s < 4; ++s) {
      const auto &kcs_groupdata = *leveldata.groupdata.at(kcss[s][i]);
      assert(kcs_groupdata.numvars == nvars);
      kcs_arrs[s] = kcs_groupdata.mfab.at(0)->const_arrays();
    }

    const amrex::IntVect ng(groupdata.nghostzones[0], groupdata.nghostzones[1],
                            groupdata.nghostzones[2]);

    if (stage == 1) {
      amrex::ParallelFor(
          Yf_mf, ng, nvars,
          [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
            if (cf_arrs[b_](i, j, k)) {
              const std::array<CCTK_REAL, 4> kk = {
                  kcs_arrs[0][b_](i, j, k, n), kcs_arrs[1][b_](i, j, k, n),
                  kcs_arrs[2][b_](i, j, k, n), kcs_arrs[3][b_](i, j, k, n)};
              const CCTK_REAL uu =
                  b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2] + b[3] * kk[3];
              Yf_arrs[b_](i, j, k, n) = u0_arrs[b_](i, j, k, n) + dtc * uu;
            }
          });
    } else if (stage == 2) {
      amrex::ParallelFor(
          Yf_mf, ng, nvars,
          [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
            if (cf_arrs[b_](i, j, k)) {
              const std::array<CCTK_REAL, 4> kk = {
                  kcs_arrs[0][b_](i, j, k, n), kcs_arrs[1][b_](i, j, k, n),
                  kcs_arrs[2][b_](i, j, k, n), kcs_arrs[3][b_](i, j, k, n)};
              const CCTK_REAL uu =
                  b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2] + b[3] * kk[3];
              const CCTK_REAL ut =
                  bt[0] * kk[0] + bt[1] * kk[1] + bt[2] * kk[2] + bt[3] * kk[3];
              Yf_arrs[b_](i, j, k, n) =
                  u0_arrs[b_](i, j, k, n) + dtc * (uu + 0.5 * r * ut);
            }
          });
    } else { // stage 3 or stage 4
      const CCTK_REAL r2 = r * r;
      const CCTK_REAL r3 = r2 * r;
      const CCTK_REAL at = (stage == 3) ? 0.5 * r : r;
      const CCTK_REAL att = (stage == 3) ? 0.25 * r2 : 0.5 * r2;
      const CCTK_REAL attt = (stage == 3) ? 0.0625 * r3 : 0.125 * r3;
      const CCTK_REAL ak = (stage == 3) ? -4.0 : 4.0;
      amrex::ParallelFor(
          Yf_mf, ng, nvars,
          [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
            if (cf_arrs[b_](i, j, k)) {
              const std::array<CCTK_REAL, 4> kk = {
                  kcs_arrs[0][b_](i, j, k, n), kcs_arrs[1][b_](i, j, k, n),
                  kcs_arrs[2][b_](i, j, k, n), kcs_arrs[3][b_](i, j, k, n)};
              const CCTK_REAL uu =
                  b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2] + b[3] * kk[3];
              const CCTK_REAL ut =
                  bt[0] * kk[0] + bt[1] * kk[1] + bt[2] * kk[2] + bt[3] * kk[3];
              const CCTK_REAL utt = btt[0] * kk[0] + btt[1] * kk[1] +
                                    btt[2] * kk[2] + btt[3] * kk[3];
              const CCTK_REAL uttt = bttt[0] * kk[0] + bttt[1] * kk[1] +
                                     bttt[2] * kk[2] + bttt[3] * kk[3];
              Yf_arrs[b_](i, j, k, n) =
                  u0_arrs[b_](i, j, k, n) +
                  dtc * (uu + at * ut + att * utt +
                         attt * (uttt + ak * (kk[2] - kk[1])));
            }
          });
    }
  }
}

} // namespace Subcycling

#endif // #ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
