#include "fillpatch.hxx"

#include <utility>

#include <AMReX_FillPatchUtil.H>

namespace CarpetX {

using namespace amrex;

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

namespace {
using task0 = void;
using task1 = std::function<task0()>;
using task2 = std::function<task1()>;
task0 do_nothing1() {}
task1 do_nothing2() { return do_nothing1; }
} // namespace

task2 FillPatch_Sync(const GHExt::PatchData::LevelData::GroupData &groupdata,
                     MultiFab &mfab, const Geometry &geom) {
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  return [&groupdata, &mfab]() -> task1 {
    mfab.FillBoundary_finish();
    groupdata.apply_boundary_conditions(mfab);
    return do_nothing1;
  };
}

task2 FillPatch_ProlongateGhosts(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const MultiFab &cmfab, const Geometry &cgeom, const Geometry &fgeom,
    Interpolater *const mapper, const Vector<BCRec> &bcrecs) {
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return do_nothing2;

  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      mfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  // Synchronize
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           fgeom.periodicity());

  if (fpc.ba_crse_patch.empty()) {
    // There is no coarser level for our boundaries, i.e. there is no
    // prolongation. Apply the boundary conditions right away.

    return [&groupdata, &mfab]() -> task1 {
      // Finish synchronizing
      mfab.FillBoundary_finish();

      // Apply symmetry and boundary conditions
      groupdata.apply_boundary_conditions(mfab);

      return do_nothing1;
    };
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

  return [&groupdata, &mfab, &cgeom, &fgeom, mapper, &bcrecs, &fpc,
          mfab_crse_patch_ptr]() -> task1 {
    const IntVect &nghosts = mfab.nGrowVect();
    const int ncomps = mfab.nComp();
    const IntVect ratio{2, 2, 2};
    MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

    // Finish synchronizing
    mfab.FillBoundary_finish();

    // Finish copying parts of coarse grid into temporary buffer
    mfab_crse_patch.ParallelCopy_finish();
    groupdata.apply_boundary_conditions(mfab_crse_patch);

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

    return [&groupdata, &mfab, mfab_fine_patch_ptr]() {
      // Finish copying fine buffer into destination
      mfab.ParallelCopy_finish();

      // Apply symmetry and boundary conditions
      groupdata.apply_boundary_conditions(mfab);

      delete mfab_fine_patch_ptr;
    };
  };
}

void FillPatch_NewLevel(const GHExt::PatchData::LevelData::GroupData &groupdata,
                        MultiFab &mfab, const MultiFab &cmfab,
                        const Geometry &cgeom, const Geometry &fgeom,
                        Interpolater *const mapper,
                        const Vector<BCRec> &bcrecs) {

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

  groupdata.apply_boundary_conditions(cmfab_g);

  FillPatchInterp(mfab, 0, cmfab_g, 0, ncomps, nghosts, cgeom, fgeom, fdomain_g,
                  ratio, mapper, bcrecs, 0);

  groupdata.apply_boundary_conditions(mfab);
}

void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const MultiFab &cmfab, const MultiFab &fmfab, const Geometry &cgeom,
    const Geometry &fgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs) {
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
    groupdata.apply_boundary_conditions(mfab_crse_patch);

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

} // namespace CarpetX
