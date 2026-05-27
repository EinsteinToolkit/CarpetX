#ifndef CARPETX_CARPETX_FILLPATCH_HXX
#define CARPETX_CARPETX_FILLPATCH_HXX

#include "driver.hxx"
#include "task_manager.hxx"

#include <functional>

namespace CarpetX {

// Sync
void FillPatch_Sync(task_manager &tasks2,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    amrex::MultiFab &mfab, const amrex::Geometry &geom);

// Prolongate ghosts from coarse level, optionally with same-level sync.
// When do_sync=true, also performs FillBoundary (same-level ghost exchange).
// When do_sync=false, only performs coarse-to-fine interpolation.
//
// Optional time-blend: when cmfab_old != nullptr and w_new != 1, the coarse
// patch is filled as w_new*cmfab + (1-w_new)*cmfab_old before interpolation.
// This is used to time-interpolate the coarse source when the fine level is
// mid-subcycle (misaligned in time with the coarse level).
void FillPatch_Prolongate(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    bool do_sync, const amrex::MultiFab *cmfab_old = nullptr,
    CCTK_REAL w_new = 1);

// Prolongate and sync ghosts (same-level exchange + coarse-to-fine
// interpolation)
inline void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    const amrex::MultiFab *cmfab_old = nullptr, CCTK_REAL w_new = 1) {
  FillPatch_Prolongate(tasks2, tasks3, groupdata, coarsegroupdata, mfab, cmfab,
                       fgeom, cgeom, mapper, bcrecs, /*do_sync=*/true,
                       cmfab_old, w_new);
}

// Prolongate only (coarse-to-fine interpolation, no same-level exchange)
inline void FillPatch_ProlongateOnly(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs) {
  FillPatch_Prolongate(tasks2, tasks3, groupdata, coarsegroupdata, mfab, cmfab,
                       fgeom, cgeom, mapper, bcrecs, /*do_sync=*/false);
}

// Band-to-band prolongation for the subcycling RK k-stages: interpolate the
// coarse level's source band into the fine level's consumer band, reusing the
// same coarse-to-fine machinery as FillPatch_Prolongate. The source band must
// cover the coarse patch (fpc.ba_crse_patch) and the consumer band must be
// allocated over fpc.ba_fine_patch / fpc.dm_patch (see build_bands).
void FillPatch_ProlongateToBand(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &consumer_band, const amrex::MultiFab &source_band,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and
// synced. ("InterpFromCoarseLevel")
void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

// ("FillPatchTwoLevels")
void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::MultiFab &fmfab, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom, amrex::Interpolater *mapper,
    const amrex::Vector<amrex::BCRec> &bcrecs);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_FILLPATCH_HXX
