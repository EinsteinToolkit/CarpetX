#ifndef CARPETX_CARPETX_FILLPATCH_HXX
#define CARPETX_CARPETX_FILLPATCH_HXX

#include "driver.hxx"
#include "task_manager.hxx"

#include <functional>

namespace CarpetX {

// Sync
//
// `bc_pass` selects which faces the boundary-condition call this function
// submits is responsible for. It is captured BY VALUE into the deferred task
// (as `tl` already is at schedule.cxx), because the closure runs long after
// this function has returned.
void FillPatch_Sync(task_manager &tasks2,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    amrex::MultiFab &mfab, const amrex::Geometry &geom,
                    bc_pass_t bc_pass = bc_pass_t::all);

// Prolongate (but do not sync) ghosts. Expects coarse mfab synced (?)
// (but not necessarily ghost-prolongated).
//
// `bc_pass` applies to the two calls that fill the REAL `mfab` (the
// no-coarser-data early return, and the post-prolongation one). The third call
// in this function's body boundary-fills the coarse TEMPORARY that
// `FillPatchInterp` then reads; that temporary is never followed by a second
// pass, so it always gets `bc_pass_t::all` and the parameter must not be
// blanket-forwarded to it.
void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &fgeom, const amrex::Geometry &cgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs,
    bc_pass_t bc_pass = bc_pass_t::all);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and
// synced. ("InterpFromCoarseLevel")
//
// No `bc_pass` parameter, deliberately: both boundary-condition calls in this
// function and in `FillPatch_RemakeLevel` take the default `all`. The coarse
// temporaries need a complete fill for the same reason as above, and the two
// real-`mfab` calls sit on the regrid path, which is out of scope here.
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
