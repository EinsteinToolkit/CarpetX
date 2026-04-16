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
void FillPatch_Prolongate(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    bool do_sync);

// Prolongate and sync ghosts (same-level exchange + coarse-to-fine
// interpolation)
inline void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs) {
  FillPatch_Prolongate(tasks2, tasks3, groupdata, coarsegroupdata, mfab, cmfab,
                       fgeom, cgeom, mapper, bcrecs, /*do_sync=*/true);
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
