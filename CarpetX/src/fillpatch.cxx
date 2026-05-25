#include "fillpatch.hxx"
#include "schedule.hxx"

#include <utility>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Version.H>

namespace CarpetX {

using namespace amrex;
#if AMREX_RELEASE_NUMBER >= 240500
using namespace amrex::detail;
#endif

// The code in this file is written in a "coroutine style"
// <https://en.wikipedia.org/wiki/Coroutine>. That is, each function
// returns another function that describes what to do next. This
// allows the caller to interleave many function calls, for example to
// schedule many calls to `MPI_Irecv` and `MPI_Isend` simultaneously.
//
// This programming style is obviously quite tedious. C++20 will have
// special support for this via `co_yield` etc.
// <https://en.cppreference.com/w/cpp/coroutine>, and the functions in
// this file will then look like normal functions.
//
// Coroutines were popularized in the "Modula" language in the 1980s.
// Welcome to the future, C++, you're only 40 years behind.

void FillPatch_Sync(task_manager &tasks2,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    MultiFab &mfab, const Geometry &geom) {
  assert(!groupdata.mfab.empty());
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  tasks2.submit_serially([&groupdata, &mfab]() {
    mfab.FillBoundary_finish();
    groupdata.apply_boundary_conditions(mfab);
  });
}

void FillPatch_Prolongate(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &fgeom,
    const Geometry &cgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs, const bool do_sync,
    const MultiFab *const cmfab_old, const CCTK_REAL w_new) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return;

  // Time-blend: when an old coarse snapshot is provided and the weight is not
  // unity, the coarse patch is filled as w_new*cmfab + (1-w_new)*cmfab_old
  // before the spatial interpolation.
  const bool do_blend = cmfab_old && w_new != CCTK_REAL(1);

  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      mfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  // Same-level ghost exchange (only when do_sync)
  if (do_sync)
    mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                             fgeom.periodicity());

  if (fpc.ba_crse_patch.empty()) {
    // There is no coarser level for our boundaries, i.e. there is no
    // prolongation.

    // Early-exit cleanup (only when do_sync)
    if (do_sync) {
      tasks2.submit_serially([&groupdata, &mfab]() {
        mfab.FillBoundary_finish();
        groupdata.apply_boundary_conditions(mfab);
      });
    }
    return;
  }

  // Prolongate from the next coarser level. Apply the boundary
  // conditions after the prolongation is done (because symmetry
  // boundary conditions might require prolongated points).

  // Copy parts of coarse grid into temporary buffer
  MultiFab *const mfab_crse_patch_ptr =
      new MultiFab(make_mf_crse_patch<MultiFab>(fpc, ncomps));
  MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;
  mf_set_domain_bndry(mfab_crse_patch, cgeom);

  // This is not local
  mfab_crse_patch.ParallelCopy_nowait(
      cmfab, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
      mfab_crse_patch.nGrowVect(), cgeom.periodicity());

  // Optionally copy the old coarse snapshot into a sibling buffer, to be
  // blended with the new one before interpolation. Same lifetime/delete
  // pattern as mfab_crse_patch; its copy is overlapped with the one above.
  MultiFab *mfab_crse_patch_old_ptr = nullptr;
  if (do_blend) {
    mfab_crse_patch_old_ptr =
        new MultiFab(make_mf_crse_patch<MultiFab>(fpc, ncomps));
    mf_set_domain_bndry(*mfab_crse_patch_old_ptr, cgeom);
    mfab_crse_patch_old_ptr->ParallelCopy_nowait(
        *cmfab_old, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
        mfab_crse_patch_old_ptr->nGrowVect(), cgeom.periodicity());
  }

  tasks2.submit_serially([&tasks3, &groupdata, &coarsegroupdata, &mfab, &cgeom,
                          &fgeom, mapper, &bcrecs, &fpc, mfab_crse_patch_ptr,
                          mfab_crse_patch_old_ptr, w_new, do_sync]() {
    const IntVect &nghosts = mfab.nGrowVect();
    const int ncomps = mfab.nComp();
    const IntVect ratio{2, 2, 2};
    MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

    // Finish same-level sync (only when do_sync)
    if (do_sync)
      mfab.FillBoundary_finish();

    // Finish copying parts of coarse grid into temporary buffer
    mfab_crse_patch.ParallelCopy_finish();

    // Blend the new and old coarse snapshots in place, in physical time, before
    // applying boundary conditions: applying the (linear) coarse BC after the
    // blend is equivalent to blending the BCs and keeps a single code path.
    if (mfab_crse_patch_old_ptr) {
      MultiFab &mfab_crse_patch_old = *mfab_crse_patch_old_ptr;
      mfab_crse_patch_old.ParallelCopy_finish();
      const CCTK_REAL w_old = CCTK_REAL(1) - w_new;
      // dst = w_new*new + w_old*old  (elementwise, in place)
      MultiFab::LinComb(mfab_crse_patch, w_new, mfab_crse_patch, 0, w_old,
                        mfab_crse_patch_old, 0, 0, ncomps, IntVect{0});
      delete mfab_crse_patch_old_ptr;
    }

    coarsegroupdata.apply_boundary_conditions(mfab_crse_patch);

    MultiFab *const mfab_fine_patch_ptr =
        new MultiFab(make_mf_fine_patch<MultiFab>(fpc, ncomps));
    MultiFab &mfab_fine_patch = *mfab_fine_patch_ptr;

    // Interpolate coarse buffer into fine buffer (in space, local)
    FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps,
                    IntVect{0} /* don't add any new ghosts */, cgeom, fgeom,
                    grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
                    ratio, mapper, bcrecs, 0);

    // Copy fine buffer into destination
    mfab.ParallelCopy_nowait(
        mfab_fine_patch, 0, 0, ncomps,
        IntVect{0} /* don't use any ghosts from the buffer */, nghosts);

    delete mfab_crse_patch_ptr;

    tasks3.submit_serially([&groupdata, &mfab, mfab_fine_patch_ptr]() {
      // Finish copying fine buffer into destination
      mfab.ParallelCopy_finish();

      // Apply symmetry and boundary conditions
      groupdata.apply_boundary_conditions(mfab);

      delete mfab_fine_patch_ptr;
    });
  });
}

void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &cgeom,
    const Geometry &fgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  // const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const BoxArray &ba = mfab.boxArray();
  const DistributionMapping &dm = mfab.DistributionMap();

  const IndexType &ixtype = ba.ixType();
  assert(ixtype == cmfab.boxArray().ixType());

  // Suffix `_g` is for "with ghosts added"
  Box fdomain_g(amrex::convert(fgeom.Domain(), mfab.ixType()));
  for (int d = 0; d < dim; ++d)
    if (fgeom.isPeriodic(d))
      fdomain_g.grow(d, nghosts[d]);

  const int nboxes = ba.size();
  BoxArray cba_g(nboxes);
  for (int i = 0; i < nboxes; ++i) {
    Box box = amrex::convert(amrex::grow(ba[i], nghosts), ixtype);
    box &= fdomain_g;
    cba_g.set(i, coarsener.doit(box));
  }
  MultiFab cmfab_g(cba_g, dm, ncomps, 0);
  mf_set_domain_bndry(cmfab_g, cgeom);

  cmfab_g.ParallelCopy(cmfab, 0, 0, ncomps, cgeom.periodicity());

  coarsegroupdata.apply_boundary_conditions(cmfab_g);

  FillPatchInterp(mfab, 0, cmfab_g, 0, ncomps, nghosts, cgeom, fgeom, fdomain_g,
                  ratio, mapper, bcrecs, 0);

  groupdata.apply_boundary_conditions(mfab);
}

void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const MultiFab &fmfab,
    const Geometry &cgeom, const Geometry &fgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      fmfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  if (!fpc.ba_crse_patch.empty()) {
    MultiFab mfab_crse_patch = make_mf_crse_patch<MultiFab>(fpc, ncomps);
    mf_set_domain_bndry(mfab_crse_patch, cgeom);

    mfab_crse_patch.ParallelCopy(
        cmfab, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
        mfab_crse_patch.nGrowVect(), cgeom.periodicity());
    coarsegroupdata.apply_boundary_conditions(mfab_crse_patch);

    MultiFab mfab_fine_patch = make_mf_fine_patch<MultiFab>(fpc, ncomps);

    // In space, local
    FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps,
                    IntVect{0} /* don't add any new ghosts */, cgeom, fgeom,
                    grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
                    ratio, mapper, bcrecs, 0);

    mfab.ParallelCopy_nowait(
        mfab_fine_patch, 0, 0, ncomps,
        IntVect{0} /* don't use any ghosts from the buffer */, nghosts);
    mfab.ParallelCopy_finish();
  }

  mfab.ParallelCopy(fmfab, 0, 0, ncomps, IntVect{0} /* don't use old ghosts */,
                    nghosts, fgeom.periodicity());
  groupdata.apply_boundary_conditions(mfab);
}

// Band-to-band prolongation for the subcycling RK k-stages. This reuses the
// same TheFPinfo + FillPatchInterp + apply_boundary_conditions primitives as
// FillPatch_Prolongate / FillPatch_RemakeLevel, but reads its coarse source
// from the coarse level's source band and writes the interpolated result into
// the fine level's consumer band, instead of reading the full coarse mfab
// interior and scattering into the fine mfab's ghost halo. It is synchronous
// (modelled on FillPatch_RemakeLevel) because the destination is a small
// pre-allocated band, so the coroutine/task overlap of FillPatch_Prolongate is
// unnecessary. The time-blend (cmfab_old) path is intentionally absent.
//
// Preconditions (guaranteed by build_bands using the same fpc):
//   - source_band covers fpc.ba_crse_patch (its cells are valid interior cells
//     written by setks, so the ParallelCopy into the coarse patch is complete);
//   - consumer_band is allocated over fpc.ba_fine_patch / fpc.dm_patch (==
//     make_mf_fine_patch geometry), so FillPatchInterp stays local.
void FillPatch_ProlongateToBand(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &consumer_band, const MultiFab &source_band, const Geometry &fgeom,
    const Geometry &cgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
  assert(!groupdata.mfab.empty());
  assert(!coarsegroupdata.mfab.empty());

  const MultiFab &mfab = *groupdata.mfab.at(0);
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return;

  const int ncomps = consumer_band.nComp();
  const IntVect ratio{2, 2, 2};
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      mfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  if (fpc.ba_crse_patch.empty())
    // There is no coarser level for our boundaries, i.e. there is no
    // prolongation.
    return;

  // Copy the coarse source band into the coarse patch buffer. ParallelCopy
  // redistributes the source band's cells onto fpc.dm_patch.
  MultiFab mfab_crse_patch = make_mf_crse_patch<MultiFab>(fpc, ncomps);
  mf_set_domain_bndry(mfab_crse_patch, cgeom);
  mfab_crse_patch.ParallelCopy(
      source_band, 0, 0, ncomps, IntVect{0} /* don't use source-band ghosts */,
      mfab_crse_patch.nGrowVect(), cgeom.periodicity());
  coarsegroupdata.apply_boundary_conditions(mfab_crse_patch);

  // Interpolate the coarse patch directly into the caller's consumer band.
  FillPatchInterp(consumer_band, 0, mfab_crse_patch, 0, ncomps,
                  IntVect{0} /* don't add any new ghosts */, cgeom, fgeom,
                  grow(convert(fgeom.Domain(), mfab.ixType()), nghosts), ratio,
                  mapper, bcrecs, 0);
}

} // namespace CarpetX
