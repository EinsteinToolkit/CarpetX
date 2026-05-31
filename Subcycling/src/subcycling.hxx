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
 * Band-native: the RK substage derivatives are read from each evolved group's
 * coarse-fine consumer bands (zero-ghost MultiFabs over the level's cf-ghost
 * region, filled by band->band prolongation from the parent's source bands).
 * The fine old state u(t0) is staged onto a band of the same geometry, the
 * dense-output value Yf is computed in place on that band, and Yf is scattered
 * into the fine state's ghost halo. Iterating the band's own BoxArray makes the
 * old cf-mask gate redundant: the band covers exactly the cf-ghost cells.
 *
 * \param leveldata Per-level data (provides groupdata / bands).
 * \param Yfs       Cactus group indices receiving the dense-output value.
 * \param Yf_tl     Time level of Yfs to write into.
 * \param dtc       Coarse-grid time step size.
 * \param xsi       Substep position within the coarse step (0 or 0.5).
 * \param stage     RK stage number (1..RKSTAGES).
 */
template <int RKSTAGES>
CCTK_HOST inline void
CalcYfFromKcs_MFlevel(CarpetX::GHExt::PatchData::LevelData &leveldata,
                      const std::vector<int> &Yfs, const int Yf_tl,
                      const CCTK_REAL dtc, const CCTK_REAL xsi,
                      const CCTK_INT stage) {
  static_assert(RKSTAGES == 3 || RKSTAGES == 4,
                "CalcYfFromKcs_MFlevel only supports RKSTAGES == 3 or 4");
  assert(stage > 0 && stage <= RKSTAGES);

  constexpr CCTK_REAL r =
      0.5; // ratio between coarse and fine cell size (2 to 1 MR case)

  const CCTK_REAL xsi2 = xsi * xsi;

  const auto &geom = CarpetX::ghext->patchdata.at(leveldata.patch)
                         .amrcore->Geom(leveldata.level);

  for (size_t i = 0; i < Yfs.size(); ++i) {
    const auto &groupdata = *leveldata.groupdata.at(Yfs[i]);
    const int nvars = groupdata.numvars;

    // No coarse-fine bands at level 0 or when subcycling is disabled.
    if (!groupdata.ks_consumer_band[0])
      continue;

    amrex::MultiFab &Yf_mf = *groupdata.mfab.at(Yf_tl);

    // All bands of this group share one BoxArray/DistributionMapping (the
    // level's cf-ghost region). The Yf scratch band uses the same.
    const amrex::MultiFab &cb0 = *groupdata.ks_consumer_band[0];
    const amrex::BoxArray &ba = cb0.boxArray();
    const amrex::DistributionMapping &dm = cb0.DistributionMap();

    const amrex::IntVect ng(groupdata.nghostzones[0], groupdata.nghostzones[1],
                            groupdata.nghostzones[2]);

    // Stage the fine old state u(t_n) from old_consumer_band (same BA/DM, so a
    // direct Copy suffices); Yf is then computed in place on yf_band.
    amrex::MultiFab yf_band(ba, dm, nvars, 0);
    assert(groupdata.old_consumer_band);
    amrex::MultiFab::Copy(yf_band, *groupdata.old_consumer_band, 0, 0, nvars,
                          0);

    auto yf_arrs = yf_band.arrays();
    amrex::GpuArray<amrex::MultiArray4<const CCTK_REAL>, RKSTAGES> kcs_arrs;
    for (int s = 0; s < RKSTAGES; ++s) {
      assert(groupdata.ks_consumer_band[s]);
      kcs_arrs[s] = groupdata.ks_consumer_band[s]->const_arrays();
    }

    if constexpr (RKSTAGES == 3) {
      // AMReX RK3 (SSPRK3) dense-output coefficients
      // (AMReX_FillPatcher.H:474-527), a degree-2 Taylor expansion over
      // S_old + k1, k2, k3 (no uttt term).
      const std::array<CCTK_REAL, 3> b = {
          xsi - (5. / 6.) * xsi2, // b1
          (1. / 6.) * xsi2,       // b2
          (2. / 3.) * xsi2        // b3
      };
      const std::array<CCTK_REAL, 3> bt = {
          1.0 - (5. / 3.) * xsi, // bt1
          (1. / 3.) * xsi,       // bt2
          (4. / 3.) * xsi        // bt3
      };
      constexpr std::array<CCTK_REAL, 3> btt = {
          -5. / 3., // btt1
          1. / 3.,  // btt2
          4. / 3.   // btt3
      };

      if (stage == 1) {
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
              const std::array<CCTK_REAL, 3> kk = {kcs_arrs[0][b_](i, j, k, n),
                                                   kcs_arrs[1][b_](i, j, k, n),
                                                   kcs_arrs[2][b_](i, j, k, n)};
              const CCTK_REAL uu = b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2];
              yf_arrs[b_](i, j, k, n) += dtc * uu;
            });
      } else if (stage == 2) {
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
              const std::array<CCTK_REAL, 3> kk = {kcs_arrs[0][b_](i, j, k, n),
                                                   kcs_arrs[1][b_](i, j, k, n),
                                                   kcs_arrs[2][b_](i, j, k, n)};
              const CCTK_REAL uu = b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2];
              const CCTK_REAL ut =
                  bt[0] * kk[0] + bt[1] * kk[1] + bt[2] * kk[2];
              // note r*ut (not 0.5*r*ut as in the RK4 stage 2)
              yf_arrs[b_](i, j, k, n) += dtc * (uu + r * ut);
            });
      } else { // stage 3
        const CCTK_REAL r2 = r * r;
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
              const std::array<CCTK_REAL, 3> kk = {kcs_arrs[0][b_](i, j, k, n),
                                                   kcs_arrs[1][b_](i, j, k, n),
                                                   kcs_arrs[2][b_](i, j, k, n)};
              const CCTK_REAL uu = b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2];
              const CCTK_REAL ut =
                  bt[0] * kk[0] + bt[1] * kk[1] + bt[2] * kk[2];
              const CCTK_REAL utt =
                  btt[0] * kk[0] + btt[1] * kk[1] + btt[2] * kk[2];
              yf_arrs[b_](i, j, k, n) +=
                  dtc * (uu + 0.5 * r * ut + 0.25 * r2 * utt);
            });
      }
    } else { // RKSTAGES == 4
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

      if (stage == 1) {
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
              const std::array<CCTK_REAL, 4> kk = {
                  kcs_arrs[0][b_](i, j, k, n), kcs_arrs[1][b_](i, j, k, n),
                  kcs_arrs[2][b_](i, j, k, n), kcs_arrs[3][b_](i, j, k, n)};
              const CCTK_REAL uu =
                  b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2] + b[3] * kk[3];
              yf_arrs[b_](i, j, k, n) += dtc * uu;
            });
      } else if (stage == 2) {
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
              const std::array<CCTK_REAL, 4> kk = {
                  kcs_arrs[0][b_](i, j, k, n), kcs_arrs[1][b_](i, j, k, n),
                  kcs_arrs[2][b_](i, j, k, n), kcs_arrs[3][b_](i, j, k, n)};
              const CCTK_REAL uu =
                  b[0] * kk[0] + b[1] * kk[1] + b[2] * kk[2] + b[3] * kk[3];
              const CCTK_REAL ut =
                  bt[0] * kk[0] + bt[1] * kk[1] + bt[2] * kk[2] + bt[3] * kk[3];
              yf_arrs[b_](i, j, k, n) += dtc * (uu + 0.5 * r * ut);
            });
      } else { // stage 3 or stage 4
        const CCTK_REAL r2 = r * r;
        const CCTK_REAL r3 = r2 * r;
        const CCTK_REAL at = (stage == 3) ? 0.5 * r : r;
        const CCTK_REAL att = (stage == 3) ? 0.25 * r2 : 0.5 * r2;
        const CCTK_REAL attt = (stage == 3) ? 0.0625 * r3 : 0.125 * r3;
        const CCTK_REAL ak = (stage == 3) ? -4.0 : 4.0;
        amrex::ParallelFor(
            yf_band, amrex::IntVect{0}, nvars,
            [=] AMREX_GPU_DEVICE(int b_, int i, int j, int k, int n) noexcept {
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
              yf_arrs[b_](i, j, k, n) +=
                  dtc * (uu + at * ut + att * utt +
                         attt * (uttt + ak * (kk[2] - kk[1])));
            });
      }
    }

    // Wait for the device kernel before the scatter reads yf_band.
    amrex::Gpu::synchronize();

    // Scatter the band (cf-ghost values of Yf) into the fine state's ghost
    // halo, mirroring FillPatch_Prolongate's final ParallelCopy.
    Yf_mf.ParallelCopy(yf_band, 0, 0, nvars, amrex::IntVect{0}, ng,
                       amrex::Periodicity::NonPeriodic());
  }
}

} // namespace Subcycling

#endif // #ifndef CARPETX_SUBCYCLING_SUBCYCLING_HXX
