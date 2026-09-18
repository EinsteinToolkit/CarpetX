#include "fillpatch.hxx"
#include "schedule.hxx"

#include <utility>

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Version.H>

namespace CarpetX {

using namespace amrex;
#if AMREX_RELEASE_NUMBER >= 240500
using namespace amrex::detail;
#endif

////////////////////////////////////////////////////////////////////////////////
//
// C-AMR -- THE CONTRACT THAT REPLACES THE BLANKET REFUSAL OF MESH REFINEMENT
// ON A MULTI-PATCH GRID.
//
//   C-AMR.  On a patch system with more than one patch, a refined level may
//   exist only where no prolongation reads coarse data from outside the patch
//   that owns it.  Concretely: for every `level > 0` fill, every box of the
//   coarse temporary must be contained in the coarse level's domain in every
//   direction whose face carries `symmetry_t::interpatch`.
//
// The refusal this replaces was a PARAMCHECK in `driver.cxx` that declined
// `num_patches() > 1 && max_num_levels > 1` outright.  It is deleted in the
// same commit: a driver that can state the condition it needs should enforce
// that condition, not the whole configuration class that might violate it.
//
// WHY THE PREDICATE IS AT THE COARSE TEMPORARY AND NOT AT THE FINE BOX.  The
// hole is not where a reader expects it, so it is written out here once:
//
//   1. `FabArrayBase::TheFPinfo` builds `ba_crse_patch` by growing each fine
//      box by `nghost`, CLIPPING that to the fine domain, and then applying
//      `Interpolater::CoarseBox` -- which grows AGAIN, by the prolongation
//      stencil.  The second grow is not clipped by anything.
//   2. `mf_set_domain_bndry(mf, cgeom)` is `FabArray::setDomainBndry`: it
//      NaN-fills, in every box, the complement of the coarse domain.  All
//      three of these temporaries carry ZERO ghost zones, so that complement
//      is the whole story.
//   3. `apply_boundary_conditions(..., bc_pass_t::all)` then writes the
//      temporary's PHYSICAL outer faces.  An interpatch face carries
//      `boundary_t::none`, so `all` writes NOTHING there, in every pass --
//      see the dispatch note in `boundaries_impl.hxx`.
//   4. `FillPatchInterp` reads what is left, which is the NaN from step 2.
//
// So a refinement box that is "well inside the patch" by eye can still violate
// the contract, if it comes within `ceil(nghost/2) + stencil` coarse cells of
// an interpatch face.  That margin is a NUMBER -- it depends on
// `blocking_factor`, `max_grid_size`, AMReX's error-buffer width and the tag
// set -- and it is QUANTISED: on the geometry this driver is used with it
// steps from +5 coarse cells of clearance to -1 with no intermediate value,
// over about 0.19 M of refinement-box travel, because the refined region snaps
// to multiples of `blocking_factor / 2`.  A guard on a position, or a warning
// when a trend looks bad, would therefore be useless.  This one evaluates the
// predicate itself, on the actual `BoxArray`, at the fill.
//
// WHY IT IS HERE AND NOT AT PARAMCHECK.  At PARAMCHECK neither the box array
// nor the patch symmetry table exists.  The predicate needs both.
//
// WHY A GUARD IS NEEDED AT ALL, RATHER THAN THE DETECTOR THAT IS ALREADY
// THERE.  The only thing on this path that notices is
// `assert(isfinite(crse(i,j,k)))` in `prolongate_3d_rf2_impl.hxx`, and that
// assert sits inside `#ifdef CCTK_DEBUG`.  One parameter file, one knob -- the
// build -- and the same 10088 NaN cells in the coarse temporary: the debug
// build aborts inside the prolongation, and the optimized build prints `Done.`
// and exits 0.  A C-AMR violation is silent in the build people run.
//
// WHAT THIS DOES NOT CLAIM.  It does not fix the hole.  The interpatch faces
// of a coarse temporary are still written by nobody; what changes is that the
// configuration which reads them is refused by name, before the first
// prolongation, under a contract that is stated rather than hoped for.
// Filling those faces from the coarse level's already-interpolated interpatch
// ghosts is a different and larger change.
//
// THE SITE IS SERIAL, AND THAT IS LOAD-BEARING: `CCTK_VERROR` raised from
// inside an `omp parallel` region corrupts the heap.  All three call sites are
// reached from serial code -- `FillPatch_ProlongateGhosts` from a
// `tasks1.submit_serially` closure run by `run_tasks_serially`, the other two
// from the per-group loops in `CactusAmrCore::MakeNewLevelFromCoarse` and
// `::RemakeLevel`.  Rather than rest on that reading, the one-shot report
// below prints `omp_in_parallel`, so the claim is measured in every log that
// reaches a refined level rather than asserted in a comment.
//
// COST.  Nothing per cell and nothing per point: `nboxes` integer box
// operations per `mf_set_domain_bndry` call, and such a call happens only when
// a `level > 0` fill happens at all.  At `max_num_levels = 1` no site here is
// ever entered, so this is inert -- not cheap, inert -- on a single-level run.
//
////////////////////////////////////////////////////////////////////////////////

namespace {

enum class camr_site_t { prolongate_ghosts = 0, new_level = 1, remake_level = 2 };

const char *camr_site_name(const camr_site_t site) {
  switch (site) {
  case camr_site_t::prolongate_ghosts:
    return "FillPatch_ProlongateGhosts";
  case camr_site_t::new_level:
    return "FillPatch_NewLevel";
  case camr_site_t::remake_level:
    return "FillPatch_RemakeLevel";
  }
  return "unknown";
}

int camr_in_parallel() {
#ifdef _OPENMP
  return omp_in_parallel();
#else
  return 0;
#endif
}

// Evaluate C-AMR on one coarse temporary, immediately after
// `mf_set_domain_bndry` has NaN-filled it and before anything prolongates from
// it.  `fine_level` is the level being filled; `coarsegroupdata` is the group
// data of the level being read.
void check_camr_contract(
    const camr_site_t site,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    const int fine_level, const MultiFab &mfab_crse, const Geometry &cgeom) {
  DECLARE_CCTK_PARAMETERS;

  const char *const sitename = camr_site_name(site);
  const int patch = coarsegroupdata.patch;
  const int clevel = coarsegroupdata.level;
  const char *const groupname = coarsegroupdata.groupname.c_str();

  // Two INDEPENDENT one-shot flags per site.  The first makes "the guard was
  // dispatched and declined" a measurement rather than an absence of evidence,
  // and is the only output this guard produces on a run that satisfies the
  // contract.  The second bounds the volume in `warn` mode.  They must not be
  // the same flag: a site whose first fill holds and whose tenth violates
  // would then report the holding one and swallow the violation.  The three
  // call sites are serial (see the note above), so plain statics are enough.
  static bool announced_ok[3] = {false, false, false};
  static bool announced_violation[3] = {false, false, false};
  const int isite = int(site);

  const auto &symmetries = ghext->patchdata.at(patch).symmetries;
  bool have_interpatch = false;
  for (int f = 0; f < 2; ++f)
    for (int d = 0; d < dim; ++d)
      have_interpatch |= symmetries[f][d] == symmetry_t::interpatch;

  if (!have_interpatch) {
    // A single-patch grid, `patch_system = "Cartesian"` included: six physical
    // outer faces and no interpatch face anywhere.  Mesh refinement on such a
    // grid works today and this guard must not touch it.  The predicate is the
    // patch's symmetry table, never `CCTK_IsFunctionAliased`, because a
    // Cartesian patch system aliases the multipatch functions and is still one
    // patch.
    if (!announced_ok[isite]) {
      announced_ok[isite] = true;
      CCTK_VINFO("C-AMR is not applicable at %s: patch %d has no interpatch "
                 "face, so no prolongation here can read across one "
                 "(level %d -> %d, group %s, omp_in_parallel=%d)",
                 sitename, patch, clevel, fine_level, groupname,
                 camr_in_parallel());
    }
    return;
  }

  // `setDomainBndry`'s own domain box, verbatim: the coarse geometry's domain
  // converted to THIS multifab's index type, grown by a full domain length in
  // each periodic direction.  The index type matters -- CarpetX groups are
  // vertex centred by default, so a 32-cell coarse domain is [0,32] here.
  Box domain_box = amrex::convert(cgeom.Domain(), mfab_crse.boxArray().ixType());
  for (int d = 0; d < dim; ++d)
    if (cgeom.isPeriodic(d))
      domain_box.grow(d, domain_box.length(d));

  const BoxArray &ba = mfab_crse.boxArray();
  const int nboxes = ba.size();

  int min_clear = INT_MAX;
  int bind_box = -1, bind_d = -1, bind_f = -1;
  long long nan_ip = 0;
  int nesc_ip = 0;

  for (int b = 0; b < nboxes; ++b) {
    const Box box = ba[b];
    // Clip only the INTERPATCH faces: cells outside the domain across a
    // physical outer face are written by `apply_boundary_conditions`, cells
    // outside across an interpatch face are written by nobody.
    Box ipclip = box;
    for (int f = 0; f < 2; ++f) {
      for (int d = 0; d < dim; ++d) {
        if (symmetries[f][d] != symmetry_t::interpatch)
          continue;
        const int clear = f == 0 ? box.smallEnd(d) - domain_box.smallEnd(d)
                                 : domain_box.bigEnd(d) - box.bigEnd(d);
        if (clear < min_clear) {
          min_clear = clear;
          bind_box = b;
          bind_d = d;
          bind_f = f;
        }
        if (f == 0)
          ipclip.setSmall(d, std::max(ipclip.smallEnd(d), domain_box.smallEnd(d)));
        else
          ipclip.setBig(d, std::min(ipclip.bigEnd(d), domain_box.bigEnd(d)));
      }
    }
    const long long escaped =
        (long long)box.numPts() - (ipclip.ok() ? (long long)ipclip.numPts() : 0);
    nan_ip += escaped;
    if (escaped > 0)
      ++nesc_ip;
  }

  if (nboxes == 0 || min_clear == INT_MAX) {
    // Nothing to prolongate from at this site on this process.
    if (!announced_ok[isite]) {
      announced_ok[isite] = true;
      CCTK_VINFO("C-AMR at %s: no coarse temporary box on this process "
                 "(patch %d, level %d -> %d, group %s, omp_in_parallel=%d)",
                 sitename, patch, clevel, fine_level, groupname,
                 camr_in_parallel());
    }
    return;
  }

  if (min_clear >= 0) {
    if (!announced_ok[isite]) {
      announced_ok[isite] = true;
      // This reports THIS call, which is the first fill at this site on this
      // process -- it is a witness that the predicate was evaluated, not a
      // claim about the whole run.  Every later call is evaluated too; only
      // this line is once.
      CCTK_VINFO("C-AMR holds at %s (first fill at this site on this "
                 "process): patch %d, level %d -> %d, group %s, %d box(es), "
                 "smallest interpatch clearance in this call %d coarse "
                 "cell(s) at direction %d face %d; omp_in_parallel=%d",
                 sitename, patch, clevel, fine_level, groupname, nboxes,
                 min_clear, bind_d, bind_f, camr_in_parallel());
    }
    return;
  }

  // Violated.  Buffered stdout is discarded by `CCTK_VERROR` and by a failing
  // assert, and everything above this point that a reader would want is on
  // stdout.
  const Box &box = ba[bind_box];
  const bool warn_only = CCTK_EQUALS(multipatch_amr_contract, "warn");
  std::fflush(nullptr);

  // In `warn` mode this is a diagnosis aid, not a mode of operation: report
  // once per site and then stay quiet, so that a run kept alive on purpose is
  // still readable.  The counts are this process's.
  if (warn_only && announced_violation[isite])
    return;
  announced_violation[isite] = true;

  char msg[2400];
  std::snprintf(
      msg, sizeof msg,
      "C-AMR is VIOLATED at %s. The coarse temporary for group %s on patch %d, "
      "level %d -> %d, box %d of %d, [%d:%d,%d:%d,%d:%d], escapes the coarse "
      "patch domain [%d:%d,%d:%d,%d:%d] by %d cell(s) across direction %d "
      "face %d, and that face carries symmetry_t::interpatch. "
      "mf_set_domain_bndry has NaN-filled %lld cell(s) in %d box(es) there; "
      "nothing writes an interpatch face of a coarse temporary, because such a "
      "face carries boundary_t::none and bc_pass_t::all skips it; and "
      "FillPatchInterp is about to prolongate from those cells. Mesh "
      "refinement on a multi-patch grid is supported only where this contract "
      "holds. Move the refined region at least %d coarse cell(s) further from "
      "the patch face (BoxInBox::radius_*, BoxInBox::position_*), or lower "
      "CarpetX::blocking_factor_{x,y,z}, which is what quantises the "
      "clearance, or raise the patch resolution -- keeping every patch's cell "
      "count a multiple of blocking_factor in every direction, which AMReX "
      "checks separately. (AMReX's error-buffer width is amr.n_error_buf and "
      "is NOT a CarpetX parameter; it is reachable only through "
      "CarpetX::amrex_parameters, whose first element is currently discarded "
      "by amrex::Initialize.) CarpetX::multipatch_amr_contract = \"warn\" "
      "downgrades this to a warning, for diagnosis only: the prolongation then "
      "reads NaN. The counts above are this process's.",
      sitename, groupname, patch, clevel, fine_level, bind_box, nboxes,
      box.smallEnd(0), box.bigEnd(0), box.smallEnd(1), box.bigEnd(1),
      box.smallEnd(2), box.bigEnd(2), domain_box.smallEnd(0),
      domain_box.bigEnd(0), domain_box.smallEnd(1), domain_box.bigEnd(1),
      domain_box.smallEnd(2), domain_box.bigEnd(2), -min_clear, bind_d, bind_f,
      nan_ip, nesc_ip, -min_clear);

  if (warn_only)
    CCTK_VWARN(CCTK_WARN_ALERT, "%s", msg);
  else
    CCTK_VERROR("%s", msg);
}

////////////////////////////////////////////////////////////////////////////////
//
// AMR-B6 / AMR-D8 INSTRUMENT -- WHAT THE FINAL BOUNDARY PASS OF EACH FILL SITE
// ACTUALLY REACHES.
//
// AMR-D8's claim is that `FillPatch_NewLevel` and `FillPatch_RemakeLevel` run
// `bc_pass_t::all` over the real `mfab` before `MultiPatch_Interpolate` has run
// on the new level, so an interpatch x outer corner there is written twice --
// once from an unpopulated interpatch ghost zone, then again by the
// corners-only pass in `regrid_interpatch_repair_apply`.  Whether that happens
// AT ALL is a number nobody had measured: `[P211]`, `[P225]` and `[P270]` say
// no rig in this project has ever entered `FillPatch_RemakeLevel`, and
// `check_camr_contract` above cannot answer it either, because at the remake
// site it is called only when the coarse temporary is non-empty while the
// boundary pass runs unconditionally.
//
// This counts, per call and without touching a single grid value:
//
//   * `nboxes`     -- local FABs of this MultiFab;
//   * `bnd_boxes`  -- how many of them fail `gdomain.contains(fab.box())`,
//                     which is EXACTLY the test `apply_boundary_conditions`
//                     makes before it constructs a `BoundaryCondition` at all
//                     (`driver.cxx`).  `bnd_boxes = 0` means the pass writes
//                     nothing whatsoever, not merely no corner;
//   * `corner_*`   -- regions and cells that are interpatch in one direction
//                     and carry a real outer condition in another.  This is
//                     AMR-D8's double-written set, and the predicate is
//                     transcribed from `boundaries_impl.hxx`, `!ip &&`
//                     included;
//   * `write_*`    -- regions and cells in which at least one direction
//                     survives the dispatch in `apply_on_face_symbc{x,y,z}`,
//                     i.e. does not map to `symmetry_t::none` /
//                     `boundary_t::none`.  A region that maps entirely to none
//                     is skipped by the `if constexpr` early-out in
//                     `apply_on_face_symbcxyz` and writes nothing, so counting
//                     regions alone would overstate the work;
//   * `crse_boxes` -- boxes in the coarse temporary this site built, so that
//                     "the site was not entered" and "the site was entered with
//                     nothing to prolongate from" are different readings.
//
// The region arithmetic is `BoundaryCondition::apply_on_face`'s, transcribed:
// same `imin`/`imax` (including the `+ !indextype[d]` that makes a
// vertex-centred domain one wider than the cell-centred one), same
// `bmin`/`bmax`, same emptiness test, and the same 26 normals -- with `INT,
// INT, INT` skipped, exactly as `BoundaryCondition::apply` skips it.  It is a
// transcription and not a call because the real thing runs inside
// `#pragma omp parallel` and writes data; this must be readable from a serial
// point and must write nothing.
//
// WHY NOT `bc_pass_census`.  That instrument exists (step B2 / B10) and counts
// the same classes, but it is `#ifdef CCTK_DEBUG` and it deliberately does not
// census `bc_pass_t::all` -- "`all` has nothing to census, it skips nothing".
// AMR-D8 is a question ABOUT the `all` pass, and a two-level moving-box rig
// costs about two minutes optimized against several hours in the debug build,
// so the measurement has to be readable in the build people run.  This
// therefore follows AMR-A4's `CAPYRX_LOG_CAMR` precedent (`[P220]`): always
// compiled, off unless asked.
//
// COST (R2).  Nothing at all unless `CARPETX_LOG_BCSITE` is set: one
// function-local static `int` compare per call.  When set, `nboxes` box
// containment tests plus 26 integer region computations per box that fails
// one, and one `std::cerr` line per call.  It reads no grid data, allocates
// nothing per cell, and is never entered from inside an `omp parallel` region
// (trap 11).
//
// `CARPETX_LOG_BCSITE=1` reports the two REGRID sites, which is AMR-D8's
// subject.  `=2` adds the two sync-path sites, which are what the
// `skip_interpatch_corners` / `interpatch_corners_only` partition already
// covers and which therefore serve as this instrument's own positive control:
// on a multipatch run they must print a NON-zero `corner_cells`, so a zero at
// the regrid sites is a measurement rather than a silent instrument
// (`[P135]`).
//
// OUTPUT GOES TO `std::cerr`, NOT `CCTK_VINFO` (`[P223]`): these counts are
// rank-local and the flesh reopens non-root stdout to the null device, so a
// `CCTK_VINFO` census would silently measure rank 0 alone.
//
////////////////////////////////////////////////////////////////////////////////

enum class bcsite_t {
  sync = 0,          // FillPatch_Sync, the level-0 branch
  prolongate = 1,    // FillPatch_ProlongateGhosts, the level>0 sync branch
  new_level = 2,     // FillPatch_NewLevel
  remake_level = 3,  // FillPatch_RemakeLevel
};

const char *bcsite_name(const bcsite_t site) {
  switch (site) {
  case bcsite_t::sync:
    return "FillPatch_Sync";
  case bcsite_t::prolongate:
    return "FillPatch_ProlongateGhosts";
  case bcsite_t::new_level:
    return "FillPatch_NewLevel";
  case bcsite_t::remake_level:
    return "FillPatch_RemakeLevel";
  }
  return "unknown";
}

// 0 = off; 1 = the two regrid sites; 2 = all four.
int bcsite_log_level() {
  static const int level = []() {
    const char *const env = std::getenv("CARPETX_LOG_BCSITE");
    if (!env)
      return 0;
    const int l = std::atoi(env);
    return l < 0 ? 0 : l;
  }();
  return level;
}

bool bcsite_wanted(const bcsite_t site) {
  const int level = bcsite_log_level();
  if (level <= 0)
    return false;
  if (site == bcsite_t::new_level || site == bcsite_t::remake_level)
    return true;
  return level >= 2;
}

void report_bc_site(const bcsite_t site,
                    const GHExt::PatchData::LevelData::GroupData &groupdata,
                    const MultiFab &mfab, const bc_pass_t bc_pass,
                    const int crse_boxes) {
  if (!bcsite_wanted(site))
    return;

  const int patch = groupdata.patch;
  const int level = groupdata.level;
  const Geometry &geom = ghext->patchdata.at(patch).amrcore->Geom(level);
  const auto &symmetries = ghext->patchdata.at(patch).symmetries;

  // `apply_boundary_conditions`' own two facts, verbatim.
  const bool all_periodic = geom.isAllPeriodic();
  Box gdomain = amrex::convert(geom.Domain(), mfab.boxArray().ixType());
  for (int d = 0; d < dim; ++d)
    if (geom.isPeriodic(d))
      gdomain.grow(d, mfab.nGrow(d));

  // `BoundaryCondition`'s own domain bounds, verbatim.
  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = geom.Domain().smallEnd(d);
    imax[d] = geom.Domain().bigEnd(d) + 1 + !groupdata.indextype[d];
  }

  bool has_interpatch = false, has_outerbc = false;
  for (int f = 0; f < 2; ++f)
    for (int d = 0; d < dim; ++d) {
      const bool ip = symmetries[f][d] == symmetry_t::interpatch;
      has_interpatch |= ip;
      has_outerbc |= !ip && groupdata.boundaries[f][d] != boundary_t::none;
    }

  int nboxes = 0, bnd_boxes = 0;
  long long corner_regions = 0, corner_cells = 0;
  long long write_regions = 0, write_cells = 0;

  const auto mfitinfo = amrex::MFItInfo().DisableDeviceSync();
  for (amrex::MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
    ++nboxes;
    const Box fbox = mfi.fabbox();
    if (all_periodic || gdomain.contains(fbox))
      continue;
    ++bnd_boxes;

    int dmin[dim], dmax[dim];
    for (int d = 0; d < dim; ++d) {
      dmin[d] = fbox.smallEnd(d);
      dmax[d] = fbox.bigEnd(d) + 1;
    }

    for (int ni = -1; ni <= 1; ++ni) {
      for (int nj = -1; nj <= 1; ++nj) {
        for (int nk = -1; nk <= 1; ++nk) {
          const int inormal[dim] = {ni, nj, nk};
          if (ni == 0 && nj == 0 && nk == 0)
            continue; // `BoundaryCondition::apply` skips <INT,INT,INT> too

          int bmin[dim], bmax[dim];
          bool empty = false;
          for (int d = 0; d < dim; ++d) {
            if (inormal[d] < 0) {
              bmin[d] = dmin[d];
              bmax[d] = std::min(dmax[d], imin[d]);
            } else if (inormal[d] > 0) {
              bmin[d] = std::max(dmin[d], imax[d]);
              bmax[d] = dmax[d];
            } else {
              bmin[d] = std::max(dmin[d], imin[d]);
              bmax[d] = std::min(dmax[d], imax[d]);
            }
            empty |= bmax[d] <= bmin[d];
          }
          if (empty)
            continue;

          long long ncells = 1;
          for (int d = 0; d < dim; ++d)
            ncells *= bmax[d] - bmin[d];

          // The corner predicate, transcribed from `boundaries_impl.hxx`.
          bool ip_here = false, obc_here = false, writes = false;
          for (int d = 0; d < dim; ++d) {
            if (inormal[d] == 0)
              continue;
            const int f = inormal[d] > 0;
            const symmetry_t sym = symmetries[f][d];
            const boundary_t bnd = groupdata.boundaries[f][d];
            const bool ip = sym == symmetry_t::interpatch;
            ip_here |= ip;
            obc_here |= !ip && bnd != boundary_t::none;
            // The dispatch's own "maps to nothing" condition, transcribed.
            const bool maps_to_none =
                (sym == symmetry_t::none && bnd == boundary_t::none) ||
                (sym == symmetry_t::interpatch &&
                 (bnd == boundary_t::none || bc_pass != bc_pass_t::all)) ||
                sym == symmetry_t::periodic;
            writes |= !maps_to_none;
          }
          if (ip_here && obc_here) {
            ++corner_regions;
            corner_cells += ncells;
          }
          if (writes) {
            ++write_regions;
            write_cells += ncells;
          }
        }
      }
    }
  }

  std::ostringstream line;
  line << "BCSITE site=" << bcsite_name(site) << " patch=" << patch
       << " level=" << level << " group=" << groupdata.groupname
       << " pass=" << bc_pass << " rank=" << amrex::ParallelDescriptor::MyProc()
       << " nboxes=" << nboxes << " bnd_boxes=" << bnd_boxes
       << " crse_boxes=" << crse_boxes
       << " corner_regions=" << corner_regions
       << " corner_cells=" << corner_cells
       << " write_regions=" << write_regions << " write_cells=" << write_cells
       << " has_interpatch=" << int(has_interpatch)
       << " has_outerbc=" << int(has_outerbc)
       << " omp_in_parallel=" << camr_in_parallel() << "\n";
  std::cerr << line.str();
}

} // namespace

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
                    MultiFab &mfab, const Geometry &geom,
                    const bc_pass_t bc_pass) {
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  // `bc_pass` is captured by value: this closure runs from `tasks2`, long after
  // this function has returned.
  tasks2.submit_serially([&groupdata, &mfab, bc_pass]() {
    mfab.FillBoundary_finish();
    // AMR-B6 instrument, `CARPETX_LOG_BCSITE >= 2`.  This site is the sync
    // path's level-0 branch and is ALREADY partitioned (step B2), so it is what
    // shows the counter can print a non-zero corner census on the same run in
    // which the regrid sites print zero.  There is no coarse temporary here.
    report_bc_site(bcsite_t::sync, groupdata, mfab, bc_pass, -1);
    groupdata.apply_boundary_conditions(mfab, bc_pass);
  });
}

void FillPatch_ProlongateGhosts(
    task_manager &tasks2, task_manager &tasks3,
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &fgeom,
    const Geometry &cgeom, Interpolater *const mapper,
    const Vector<BCRec> &bcrecs, const bc_pass_t bc_pass) {
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return;

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

    tasks2.submit_serially([&groupdata, &mfab, bc_pass]() {
      // Finish synchronizing
      mfab.FillBoundary_finish();

      // AMR-B6 instrument, `CARPETX_LOG_BCSITE >= 2`; `crse_boxes = 0` is the
      // reason this branch was taken.
      report_bc_site(bcsite_t::prolongate, groupdata, mfab, bc_pass, 0);

      // Apply symmetry and boundary conditions
      groupdata.apply_boundary_conditions(mfab, bc_pass);
    });
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
  check_camr_contract(camr_site_t::prolongate_ghosts, coarsegroupdata,
                      groupdata.level, mfab_crse_patch, cgeom);

  // This is not local
  mfab_crse_patch.ParallelCopy_nowait(
      cmfab, 0, 0, ncomps, IntVect{0} /* don't use coarse ghosts */,
      mfab_crse_patch.nGrowVect(), cgeom.periodicity());

  const int ncrse = fpc.ba_crse_patch.size();
  tasks2.submit_serially([&tasks3, &groupdata, &coarsegroupdata, &mfab, &cgeom,
                          &fgeom, mapper, &bcrecs, &fpc, mfab_crse_patch_ptr,
                          bc_pass, ncrse]() {
    const IntVect &nghosts = mfab.nGrowVect();
    const int ncomps = mfab.nComp();
    const IntVect ratio{2, 2, 2};
    MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

    // Finish synchronizing
    mfab.FillBoundary_finish();

    // Finish copying parts of coarse grid into temporary buffer
    mfab_crse_patch.ParallelCopy_finish();

    // The coarse TEMPORARY, not the real `mfab`: `FillPatchInterp` below reads
    // it and nothing runs a second BC pass over it, so it must be filled
    // completely. Deliberately `all`, never `bc_pass`.
    //
    // "Completely" stopped being true at step B7 on a multipatch grid: `all`
    // writes the temporary's physical outer faces and, since B7, not its
    // interpatch faces, which `mf_set_domain_bndry` above has NaN-filled and
    // which nothing else touches -- the interpolator only ever writes the real
    // `mfab`. `FillPatchInterp` then reads them. That hole is still here; what
    // closed is the route to it, at the `check_camr_contract` call above,
    // which refuses by name any configuration in which one of these boxes
    // reaches outside the coarse patch domain across an interpatch face. See
    // the C-AMR block at the top of this file.
    coarsegroupdata.apply_boundary_conditions(mfab_crse_patch, bc_pass_t::all);

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

    tasks3.submit_serially([&groupdata, &mfab, mfab_fine_patch_ptr, bc_pass,
                            ncrse]() {
      // Finish copying fine buffer into destination
      mfab.ParallelCopy_finish();

      // AMR-B6 instrument, `CARPETX_LOG_BCSITE >= 2`.
      report_bc_site(bcsite_t::prolongate, groupdata, mfab, bc_pass, ncrse);

      // Apply symmetry and boundary conditions
      groupdata.apply_boundary_conditions(mfab, bc_pass);

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
  check_camr_contract(camr_site_t::new_level, coarsegroupdata, groupdata.level,
                      cmfab_g, cgeom);

  cmfab_g.ParallelCopy(cmfab, 0, 0, ncomps, cgeom.periodicity());

  coarsegroupdata.apply_boundary_conditions(cmfab_g);

  FillPatchInterp(mfab, 0, cmfab_g, 0, ncomps, nghosts, cgeom, fgeom, fdomain_g,
                  ratio, mapper, bcrecs, 0);

  // AMR-B6 instrument, `CARPETX_LOG_BCSITE >= 1`.  AMR-D8's first site.
  report_bc_site(bcsite_t::new_level, groupdata, mfab, bc_pass_t::all, nboxes);
  // Deliberately `all` and not `skip_interpatch_corners`: see the paragraph
  // headed "WHAT IS NOT DISJOINT HERE" in `schedule.cxx`.  The line the
  // instrument above prints carries the three numbers that argument rests on
  // -- `bnd_boxes`, `has_outerbc` and `corner_cells`.
  groupdata.apply_boundary_conditions(mfab);
}

void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    const GHExt::PatchData::LevelData::GroupData &coarsegroupdata,
    MultiFab &mfab, const MultiFab &cmfab, const MultiFab &fmfab,
    const Geometry &cgeom, const Geometry &fgeom, Interpolater *const mapper,
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
    check_camr_contract(camr_site_t::remake_level, coarsegroupdata,
                        groupdata.level, mfab_crse_patch, cgeom);

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
  // AMR-B6 instrument, `CARPETX_LOG_BCSITE >= 1`.  AMR-D8's second site, and
  // the one that matters most: it is OUTSIDE the `!fpc.ba_crse_patch.empty()`
  // block above, so this line is printed whenever `FillPatch_RemakeLevel` runs
  // at all -- which is what separates `[P225]`'s "site 3 is silent" into "the
  // function was never called" and "it was called with an empty coarse
  // temporary".  `check_camr_contract` cannot make that distinction.
  report_bc_site(bcsite_t::remake_level, groupdata, mfab, bc_pass_t::all,
                 fpc.ba_crse_patch.size());
  // Deliberately `all`, for the same three reasons as the twin call in
  // `FillPatch_NewLevel` above; see `schedule.cxx`.
  groupdata.apply_boundary_conditions(mfab);
}

} // namespace CarpetX
