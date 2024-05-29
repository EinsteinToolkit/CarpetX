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

// Prolongate (but do not sync) ghosts. Expects coarse mfab synced (?)
// (but not necessarily ghost-prolongated).
void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
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
