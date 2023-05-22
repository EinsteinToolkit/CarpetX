#ifndef CARPETX_CARPETX_FILLPATCH_HXX
#define CARPETX_CARPETX_FILLPATCH_HXX

#include "driver.hxx"

#include <functional>

namespace CarpetX {

// Sync
std::function<std::function<void()>()>
FillPatch_Sync(const GHExt::PatchData::LevelData::GroupData &groupdata,
               amrex::MultiFab &mfab, const amrex::Geometry &geom);

// Prolongate (but do not sync) ghosts. Expects coarse mfab synced (but not
// necessarily ghost-prolongated).
std::function<std::function<void()>()> FillPatch_ProlongateGhosts(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and synced.
// ("InterpFromCoarseLevel")
void FillPatch_NewLevel(const GHExt::PatchData::LevelData::GroupData &groupdata,
                        amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
                        const amrex::Geometry &cgeom,
                        const amrex::Geometry &fgeom,
                        amrex::Interpolater *mapper,
                        const amrex::Vector<amrex::BCRec> &bcrecs);

// ("FillPatchTwoLevels")
void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::MultiFab &fmfab, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom, amrex::Interpolater *mapper,
    const amrex::Vector<amrex::BCRec> &bcrecs);

} // namespace CarpetX

#endif // #ifndef CARPETX_CARPETX_FILLPATCH_HXX
