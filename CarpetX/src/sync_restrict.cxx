#include "schedule.hxx"
#include "driver.hxx"
#include "fillpatch.hxx"
#include "sync_restrict_internal.hxx"
#include "task_manager.hxx"
#include "timer.hxx"
#include "valid.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <AMReX_MultiFabUtil.H>

#include <cassert>
#include <sstream>
#include <vector>

namespace CarpetX {

// Forward declaration for the 3-arg Restrict wrapper, defined later in this
// TU. SyncGroupsByDirI's restrict_during_sync branch calls it, but to keep
// the top-down order (helpers → Sync* → Reflux → Restrict*) the definition
// appears below.
void Restrict(const cGH *cctkGH, int level, const std::vector<int> &groups);

// =======================================================================
// Sync helpers
// =======================================================================

bool sync_active = false; // Catch recursive calls

struct mark_sync_active {
  mark_sync_active() {
    if (sync_active)
      CCTK_ERROR(
          "Recursive call to SyncGroupsByDirI. Maybe you are syncing grid "
          "functions in the \"restrict\" bin while the parameter "
          "\"restrict_during_sync\" is true?");
    sync_active = true;
  }
  ~mark_sync_active() { sync_active = false; }
};

static void sync_log_groups(const char *label, int numgroups,
                            const int *groups0) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose) {
    std::ostringstream buf;
    for (int n = 0; n < numgroups; ++n) {
      if (n != 0)
        buf << ", ";
      buf << CCTK_FullGroupName(groups0[n]);
    }
#pragma omp critical
    CCTK_VINFO("%s %s", label, buf.str().c_str());
  }
}

static std::vector<int> sync_filter_groups(int numgroups, const int *groups0) {
  static const int gi_regrid_error =
      CCTK_GroupIndex("CarpetXRegrid::regrid_error");
  assert(gi_regrid_error >= 0);

  std::vector<int> groups;
  for (int n = 0; n < numgroups; ++n) {
    const int gi = groups0[n];
    if (ghext->active_timelevels.at(gi) == 0)
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;
    // Don't restrict the regridding error
    if (gi == gi_regrid_error)
      continue;
    groups.push_back(gi);
  }
  return groups;
}

static void sync_multipatch_postcheck(const cGH *cctkGH,
                                      const std::vector<int> &groups,
                                      const char *label) {
  static const bool have_multipatch_boundaries =
      CCTK_IsFunctionAliased("MultiPatch_Interpolate");

  if (have_multipatch_boundaries) {
    std::vector<CCTK_INT> cactusvarinds;
    for (int group : groups) {
      const auto &groupdata =
          *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(group);
      for (int var = 0; var < groupdata.numvars; ++var)
        cactusvarinds.push_back(groupdata.firstvarindex + var);
    }
    MultiPatch_Interpolate(cctkGH, cactusvarinds.size(), cactusvarinds.data());

    for (const int gi : groups) {
      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      assert(!groupdata0.mfab.empty());
      const nan_handling_t nan_handling = groupdata0.do_evolve
                                              ? nan_handling_t::forbid_nans
                                              : nan_handling_t::allow_nans;
      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls0 = groupdata0.mfab.size();
      const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

      for (int tl = 0; tl < sync_tl0; ++tl)
        for (int vi = 0; vi < groupdata0.numvars; ++vi)
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling, [label]() {
            return std::string(label) + " after syncing";
          });

    } // for gi

  } else {
    assert(ghext->num_patches() == 1);
  }
}

static std::vector<int> collect_restrictable_groups() {
  const int numgroups = CCTK_NumGroups();
  std::vector<int> groups;
  groups.reserve(numgroups);
  const auto &patchdata0 = ghext->patchdata.at(0);
  const auto &leveldata0 = patchdata0.leveldata.at(0);
  for (const auto &groupdataptr : leveldata0.groupdata) {
    // Only grid functions
    if (groupdataptr) {
      auto &restrict groupdata = *groupdataptr;
      // Only grid functions with storage
      if (groupdata.mfab.empty())
        continue;
      // Only grid functions with restriction enabled
      if (groupdata.do_restrict)
        groups.push_back(groupdata.groupindex);
    }
  }
  return groups;
}

// =======================================================================
// Sync entry points
// =======================================================================

static int
SyncGroupsByDirIProlongateOnly_impl(const cGH *restrict cctkGH, int numgroups,
                                    const int *groups0, const int *directions,
                                    const bool prolongate_on_same_iteration,
                                    const int tl_arg = -1) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH) || in_level_mode(cctkGH));

  mark_sync_active marked;

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  sync_log_groups("SyncGroupsProlongateOnly", numgroups, groups0);

  std::vector<int> groups = sync_filter_groups(numgroups, groups0);

  // We need to loop over groups, patches, and levels in a definite
  // order so that AMReX's communication pattern does not get
  // confused. Therefore all the loops here are serial. The only
  // parallelization happens within AMReX and within our boundary
  // conditions. This is not efficient.

  task_manager tasks1;
  task_manager tasks2;
  task_manager tasks3;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      // We always sync all directions.
      // tl_arg >= 0 syncs only that timelevel; tl_arg = -1 syncs every
      // timelevel except the oldest (or the only one, if ntls == 1).
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int tl_lo = (tl_arg < 0) ? 0 : tl_arg;
      const int tl_hi =
          (tl_arg < 0) ? (ntls > 1 ? ntls - 1 : ntls) : (tl_arg + 1);
      assert(tl_lo >= 0 && tl_hi <= ntls);

      if (leveldata.level == 0) {
        // Level 0 requires no interpolation, so return early
        return;
      }

      // For levels greater than 0, interpolate from the next coarser level

      const int level = leveldata.level;
      const auto &restrict coarseleveldata =
          ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);

      if (ghext->use_subcycling && prolongate_on_same_iteration) {
        if (leveldata.iteration != coarseleveldata.iteration)
          return;
      }

      auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
      assert(!coarsegroupdata.mfab.empty());
      assert(coarsegroupdata.numvars == groupdata.numvars);

      amrex::Interpolater *const interpolator = groupdata.interpolator;

      for (int tl = tl_lo; tl < tl_hi; ++tl) {

        tasks1.submit_serially([&tasks2, &tasks3, &leveldata, &groupdata,
                                &coarsegroupdata, interpolator, tl]() {
          FillPatch_ProlongateOnly(tasks2, tasks3, groupdata, coarsegroupdata,
                                   *groupdata.mfab.at(tl),
                                   *coarsegroupdata.mfab.at(tl),
                                   ghext->patchdata.at(leveldata.patch)
                                       .amrcore->Geom(leveldata.level),
                                   ghext->patchdata.at(leveldata.patch)
                                       .amrcore->Geom(leveldata.level - 1),
                                   interpolator, groupdata.bcrecs);
        });

      } // for tl
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();
  tasks3.run_tasks_serially();
  synchronize();

  sync_multipatch_postcheck(cctkGH, groups, "SyncGroupsByDirIProlongateOnly");

  assert(sync_active);

  return numgroups; // number of groups synchronized
}

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups0, const int *directions) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH) || in_level_mode(cctkGH));

  mark_sync_active marked;

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  sync_log_groups("SyncGroups", numgroups, groups0);

  std::vector<int> groups = sync_filter_groups(numgroups, groups0);

  // Skip groups that have valid ghosts and boundaries
  if (CCTK_EQUALS(presync_mode, "presync-only")) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      std::vector<int> new_groups;
      for (const int gi : groups) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        assert(!groupdata.mfab.empty());
        bool need_sync = false;
        for (int tl = 0; tl < int(groupdata.valid.size()); tl++) {
          if (need_sync)
            break;
          auto &timeleveldata = groupdata.valid.at(tl);
          for (int vi = 0; vi < int(timeleveldata.size()); vi++) {
            if (need_sync)
              break;
            valid_t have = groupdata.valid.at(tl).at(vi).get();
            if (!have.valid_ghosts || !have.valid_outer) {
              need_sync = true;
            }
          }
        }
        if (need_sync) {
          new_groups.push_back(gi);
        }
      }
      groups = new_groups;
    });
    if (groups.size() == 0) {
      return 0;
    }
  }

  if (restrict_during_sync) {
    active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
      if (leveldata.level < ghext->num_levels() - 1)
        Restrict(cctkGH, leveldata.level, groups);
    });
    // FIXME: cannot call POSTRESTRICT since this could contain a SYNC leading
    // to an infinite loop. This means that outer boundaries will be left
    // invalid after an implicit restrict
    // CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
  }

  static const bool have_multipatch_boundaries =
      CCTK_IsFunctionAliased("MultiPatch_Interpolate");

  // Check preconditions
  for (const int gi : groups) {
    const auto &patchdata0 = ghext->patchdata.at(0);
    const auto &leveldata0 = patchdata0.leveldata.at(0);
    const auto &groupdata0 = *leveldata0.groupdata.at(gi);
    assert(!groupdata0.mfab.empty());
    const nan_handling_t nan_handling = groupdata0.do_evolve
                                            ? nan_handling_t::forbid_nans
                                            : nan_handling_t::allow_nans;
    // We always sync all directions.
    // If there is more than one time level, then we don't sync the
    // oldest.
    // TODO: during evolution, sync only one time level
    const int ntls0 = groupdata0.mfab.size();
    const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      if (leveldata.level > 0) {

        const int level = leveldata.level;
        const auto &restrict coarseleveldata =
            ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(!coarsegroupdata.mfab.empty());
        assert(coarsegroupdata.numvars == groupdata.numvars);

        for (int tl = 0; tl < sync_tl0; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            error_if_invalid(coarsegroupdata, vi, tl, make_valid_int(), []() {
              return "SyncGroupsByDirI on coarse level before prolongation";
            });
          }
        } // for tl

      } // if leveldata.level > 0

      for (int tl = 0; tl < sync_tl0; ++tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          // Synchronization only uses the interior
          error_if_invalid(groupdata, vi, tl, make_valid_int(),
                           []() { return "SyncGroupsByDirI before syncing"; });
          groupdata.valid.at(tl).at(vi).set_invalid(make_valid_ghosts(), []() {
            return "SyncGroupsByDirI before syncing: "
                   "Mark ghost zones as invalid";
          });
        }
      } // for tl
    });

    active_levels_t active_fine_levels = *active_levels;
    using std::max;
    active_fine_levels.min_level = max(active_fine_levels.min_level, 1);
    for (int tl = 0; tl < sync_tl0; ++tl) {
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        check_valid_gf(active_fine_levels, gi, vi, tl, nan_handling, []() {
          return "SyncGroupsByDirI on coarse level before prolongation";
        });
        poison_invalid_gf(*active_levels, gi, vi, tl);
        check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                       []() { return "SyncGroupsByDirI before syncing"; });
      }
    } // for tl
  } // for gi

  // We need to loop over groups, patches, and levels in a definite
  // order so that AMReX's communication pattern does not get
  // confused. Therefore all the loops here are serial. The only
  // parallelization happens within AMReX and within our boundary
  // conditions. This is not efficient.

  task_manager tasks1;
  task_manager tasks2;
  task_manager tasks3;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      // const int level = leveldata.level;
      // const auto &restrict coarseleveldata =
      //     level == 0
      //         ? ghext->patchdata.at(leveldata.patch).leveldata.at(level)
      //         : ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
      // const int rhs_key_exists = Util_TableQueryValueInfo(
      //     CCTK_GroupTagsTableI(gi), nullptr, nullptr, "rhs");

      // const bool exchange_ghost_only =
      //     level == 0 ||
      //     (rhs_key_exists && leveldata.iteration !=
      //     coarseleveldata.iteration);

      // if (exchange_ghost_only) {
      if (leveldata.level == 0) {
        // Copy from adjacent boxes on same level

        for (int tl = 0; tl < sync_tl; ++tl) {
          tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
            FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                           ghext->patchdata.at(leveldata.patch)
                               .amrcore->Geom(leveldata.level));
          });
        } // for tl

      } else { // if leveldata.level > 0
        // Copy from adjacent boxes on same level, and interpolate
        // from next coarser level

        const int level = leveldata.level;
        const auto &restrict coarseleveldata =
            ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(!coarsegroupdata.mfab.empty());
        assert(coarsegroupdata.numvars == groupdata.numvars);

        amrex::Interpolater *const interpolator = groupdata.interpolator;

        for (int tl = 0; tl < sync_tl; ++tl) {

          tasks1.submit_serially([&tasks2, &tasks3, &leveldata, &groupdata,
                                  &coarsegroupdata, interpolator, tl]() {
            FillPatch_ProlongateGhosts(tasks2, tasks3, groupdata,
                                       coarsegroupdata, *groupdata.mfab.at(tl),
                                       *coarsegroupdata.mfab.at(tl),
                                       ghext->patchdata.at(leveldata.patch)
                                           .amrcore->Geom(leveldata.level),
                                       ghext->patchdata.at(leveldata.patch)
                                           .amrcore->Geom(leveldata.level - 1),
                                       interpolator, groupdata.bcrecs);
          });

        } // for tl

      } // if leveldata.level > 0
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();
  tasks3.run_tasks_serially();
  synchronize();

  // Check postconditions
  for (const int gi : groups) {
    const auto &patchdata0 = ghext->patchdata.at(0);
    const auto &leveldata0 = patchdata0.leveldata.at(0);
    const auto &groupdata0 = *leveldata0.groupdata.at(gi);
    assert(!groupdata0.mfab.empty());
    const nan_handling_t nan_handling = groupdata0.do_evolve
                                            ? nan_handling_t::forbid_nans
                                            : nan_handling_t::allow_nans;
    // We always sync all directions.
    // If there is more than one time level, then we don't sync the
    // oldest.
    // TODO: during evolution, sync only one time level
    const int ntls0 = groupdata0.mfab.size();
    const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      for (int tl = 0; tl < sync_tl0; ++tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          groupdata.valid.at(tl).at(vi).set_ghosts(true, []() {
            return "SyncGroupsByDirI after syncing: "
                   "Mark ghost zones as valid";
          });
          if (groupdata.all_faces_have_symmetries_or_boundaries())
            groupdata.valid.at(tl).at(vi).set_outer(true, []() {
              return "SyncGroupsByDirI after syncing: "
                     "Mark outer boundaries as valid";
            });
        }
      } // for tl
    });

    for (int tl = 0; tl < sync_tl0; ++tl) {
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        poison_invalid_gf(*active_levels, gi, vi, tl);
        // TODO: Check after applying multi-patch boundaries
        if (!have_multipatch_boundaries)
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling,
                         []() { return "SyncGroupsByDirI after syncing"; });
      }
    } // for tl
  } // for gi

  sync_multipatch_postcheck(cctkGH, groups, "SyncGroupsByDirI");

  assert(sync_active);

  return numgroups; // number of groups synchronized
}

int SyncGroupsByDirISubcycling(const cGH *restrict cctkGH, int numgroups,
                               const int *groups0, const int *directions) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH) || in_level_mode(cctkGH));

  mark_sync_active marked;

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  sync_log_groups("SyncGroupsBySubcycling", numgroups, groups0);

  std::vector<int> groups = sync_filter_groups(numgroups, groups0);

  static const bool have_multipatch_boundaries =
      CCTK_IsFunctionAliased("MultiPatch_Interpolate");

  // We need to loop over groups, patches, and levels in a definite
  // order so that AMReX's communication pattern does not get
  // confused. Therefore all the loops here are serial. The only
  // parallelization happens within AMReX and within our boundary
  // conditions. This is not efficient.

  task_manager tasks1;
  task_manager tasks2;
  task_manager tasks3;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      if (leveldata.level == 0) {
        // Copy from adjacent boxes on same level

        for (int tl = 0; tl < sync_tl; ++tl) {
          tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
            FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                           ghext->patchdata.at(leveldata.patch)
                               .amrcore->Geom(leveldata.level));
          });
        } // for tl

      } else { // if leveldata.level > 0

        const int level = leveldata.level;
        const auto &restrict coarseleveldata =
            ghext->patchdata.at(leveldata.patch).leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(!coarsegroupdata.mfab.empty());
        assert(coarsegroupdata.numvars == groupdata.numvars);

        const bool evolving_subiter =
            groupdata.do_evolve && leveldata.iteration > 0;

        if (evolving_subiter) {
          // Copy from adjacent boxes on same level only
          for (int tl = 0; tl < sync_tl; ++tl) {
            tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
              FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                             ghext->patchdata.at(leveldata.patch)
                                 .amrcore->Geom(leveldata.level));
            });
          } // for tl
        } else {
          // Copy from adjacent boxes on same level, and interpolate
          // from next coarser level
          amrex::Interpolater *const interpolator = groupdata.interpolator;

          // Time-blend gate: only a non-evolved group can reach this branch
          // while misaligned in time with its parent. When the fine level is
          // ahead of the coarse, blend coarse tl=0 (new) with coarse tl=1
          // (old) into the coarse patch before interpolating. Mirrors AMReX
          // beta=(time-t0)/(t1-t0): with t0 = t_new - cdt (old coarse time),
          // t1 = t_new (new coarse time), time = t_fin (fine time).
          const bool aligned =
              (leveldata.iteration == coarseleveldata.iteration);
          const rat64 cdt = coarseleveldata.delta_iteration;
          const rat64 t_new = coarseleveldata.iteration;
          const rat64 t_fin = leveldata.iteration;
          const CCTK_REAL w_new = CCTK_REAL((t_fin - t_new + cdt) / cdt);

          // The blend needs a valid old coarse snapshot (tl=1).
          bool old_valid = coarsegroupdata.mfab.size() >= 2;
          if (old_valid)
            for (int vi = 0; vi < coarsegroupdata.numvars; ++vi)
              old_valid = old_valid &&
                          coarsegroupdata.valid.at(1).at(vi).get().valid_int;

          for (int tl = 0; tl < sync_tl; ++tl) {
            // Only tl=0 is the "new" coarse snapshot whose old partner is tl=1.
            const bool do_blend = !aligned && old_valid && tl == 0;
            const amrex::MultiFab *const cmfab_old =
                do_blend ? coarsegroupdata.mfab.at(1).get() : nullptr;
            const CCTK_REAL tl_w_new = do_blend ? w_new : CCTK_REAL(1);
            tasks1.submit_serially([&tasks2, &tasks3, &leveldata, &groupdata,
                                    &coarsegroupdata, interpolator, tl,
                                    cmfab_old, tl_w_new]() {
              FillPatch_ProlongateGhosts(
                  tasks2, tasks3, groupdata, coarsegroupdata,
                  *groupdata.mfab.at(tl), *coarsegroupdata.mfab.at(tl),
                  ghext->patchdata.at(leveldata.patch)
                      .amrcore->Geom(leveldata.level),
                  ghext->patchdata.at(leveldata.patch)
                      .amrcore->Geom(leveldata.level - 1),
                  interpolator, groupdata.bcrecs, cmfab_old, tl_w_new);
            });
          } // for tl
        }

      } // if leveldata.level > 0
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();
  tasks3.run_tasks_serially();
  synchronize();

  // Check postconditions
  for (const int gi : groups) {
    const auto &patchdata0 = ghext->patchdata.at(0);
    const auto &leveldata0 = patchdata0.leveldata.at(0);
    const auto &groupdata0 = *leveldata0.groupdata.at(gi);
    assert(!groupdata0.mfab.empty());
    const nan_handling_t nan_handling = groupdata0.do_evolve
                                            ? nan_handling_t::forbid_nans
                                            : nan_handling_t::allow_nans;
    // We always sync all directions.
    // If there is more than one time level, then we don't sync the
    // oldest.
    // TODO: during evolution, sync only one time level
    const int ntls0 = groupdata0.mfab.size();
    const int sync_tl0 = ntls0 > 1 ? ntls0 - 1 : ntls0;

    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      const bool evolving_subiter =
          groupdata.do_evolve && leveldata.iteration > 0;

      for (int tl = 0; tl < sync_tl0; ++tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          if (!evolving_subiter) {
            groupdata.valid.at(tl).at(vi).set_ghosts(true, []() {
              return "SyncGroupsByDirISubcycling after syncing: "
                     "Mark ghost zones as valid";
            });
          }
          if (groupdata.all_faces_have_symmetries_or_boundaries())
            groupdata.valid.at(tl).at(vi).set_outer(true, []() {
              return "SyncGroupsByDirISubcycling after syncing: "
                     "Mark outer boundaries as valid";
            });
        }
      } // for tl
    });

    for (int tl = 0; tl < sync_tl0; ++tl) {
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        poison_invalid_gf(*active_levels, gi, vi, tl);
        // TODO: Check after applying multi-patch boundaries
        if (!have_multipatch_boundaries)
          check_valid_gf(*active_levels, gi, vi, tl, nan_handling, []() {
            return "SyncGroupsByDirISubcycling after syncing";
          });
      }
    } // for tl
  } // for gi

  sync_multipatch_postcheck(cctkGH, groups, "SyncGroupsByDirISubcycling");

  assert(sync_active);

  return numgroups; // number of groups synchronized
}

int SyncGroupsByDirIProlongateOnly(const cGH *restrict cctkGH, int numgroups,
                                   const int *groups0, const int *directions,
                                   const int tl) {
  return SyncGroupsByDirIProlongateOnly_impl(cctkGH, numgroups, groups0,
                                             directions, false, tl);
}

int SyncGroupsByDirIGhostOnly(const cGH *restrict cctkGH, int numgroups,
                              const int *groups0, const int *directions,
                              const int tl_arg) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH) || in_level_mode(cctkGH));

  mark_sync_active marked;

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  sync_log_groups("SyncGroupsGhostOnly", numgroups, groups0);

  std::vector<int> groups = sync_filter_groups(numgroups, groups0);

  // We need to loop over groups, patches, and levels in a definite
  // order so that AMReX's communication pattern does not get
  // confused. Therefore all the loops here are serial. The only
  // parallelization happens within AMReX and within our boundary
  // conditions. This is not efficient.

  task_manager tasks1;
  task_manager tasks2;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      assert(!groupdata.mfab.empty());

      // We always sync all directions.
      // tl_arg >= 0 syncs only that timelevel; tl_arg = -1 syncs every
      // timelevel except the oldest (or the only one, if ntls == 1).
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int tl_lo = (tl_arg < 0) ? 0 : tl_arg;
      const int tl_hi =
          (tl_arg < 0) ? (ntls > 1 ? ntls - 1 : ntls) : (tl_arg + 1);
      assert(tl_lo >= 0 && tl_hi <= ntls);

      // Copy from adjacent boxes on same level
      for (int tl = tl_lo; tl < tl_hi; ++tl) {
        tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
          FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                         ghext->patchdata.at(leveldata.patch)
                             .amrcore->Geom(leveldata.level));
        });
      } // for tl
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();

  sync_multipatch_postcheck(cctkGH, groups, "SyncGroupsByDirIGhostOnly");

  assert(sync_active);

  return numgroups; // number of groups synchronized
}

void ProlongateRestrictedGFs(const cGH *cctkGH) {
  const std::vector<int> groups = collect_restrictable_groups();
  SyncGroupsByDirIProlongateOnly_impl(cctkGH, groups.size(), groups.data(),
                                      nullptr, true);
}

// =======================================================================
// Reflux
// =======================================================================

void Reflux(const cGH *cctkGH, int level) {
  DECLARE_CCTK_PARAMETERS;

  if (!do_reflux)
    return;

  static Timer timer("Reflux");
  Interval interval(timer);

  for (const auto &patchdata : ghext->patchdata) {
    if (level + 1 < int(patchdata.leveldata.size())) {
      auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      for (int gi = 0; gi < int(leveldata.groupdata.size()); ++gi) {
        const int tl = 0;
        cGroup group;
        int ierr = CCTK_GroupData(gi, &group);
        assert(!ierr);

        if (group.grouptype != CCTK_GF)
          continue;

        auto &groupdata = *leveldata.groupdata.at(gi);
        if (groupdata.mfab.empty())
          continue;
        const auto &finegroupdata = *fineleveldata.groupdata.at(gi);
        assert(!finegroupdata.mfab.empty());
        const nan_handling_t nan_handling = groupdata.do_evolve
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        // If the group has associated fluxes
        if (finegroupdata.freg) {

          // Check coarse and fine data and fluxes are valid
          for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
            error_if_invalid(finegroupdata, vi, tl, make_valid_int(), []() {
              return "Reflux before refluxing: Fine level data";
            });
            error_if_invalid(groupdata, vi, tl, make_valid_int(), []() {
              return "Reflux before refluxing: Coarse level data";
            });
          }
          for (int d = 0; d < dim; ++d) {
            const int flux_gi = finegroupdata.fluxes.at(d);
            const auto &flux_finegroupdata =
                *fineleveldata.groupdata.at(flux_gi);
            const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
            for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
              error_if_invalid(
                  flux_finegroupdata, vi, tl, make_valid_int(), [&]() {
                    std::ostringstream buf;
                    buf << "Reflux: Fine level flux in direction " << d;
                    return buf.str();
                  });
              error_if_invalid(flux_groupdata, vi, tl, make_valid_int(), [&]() {
                std::ostringstream buf;
                buf << "Reflux: Coarse level flux in direction " << d;
                return buf.str();
              });
            }
          }

          for (int d = 0; d < dim; ++d) {
            const int flux_gi = finegroupdata.fluxes.at(d);
            const auto &flux_finegroupdata =
                *fineleveldata.groupdata.at(flux_gi);
            const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
            finegroupdata.freg->CrseInit(*flux_groupdata.mfab.at(tl), d, 0, 0,
                                         flux_groupdata.numvars, -1);
            finegroupdata.freg->FineAdd(*flux_finegroupdata.mfab.at(tl), d, 0,
                                        0, flux_finegroupdata.numvars, 1);
          }
          const amrex::Geometry &geom = patchdata.amrcore->Geom(level);
          finegroupdata.freg->Reflux(*groupdata.mfab.at(tl), 1.0, 0, 0,
                                     groupdata.numvars, geom);

          const active_levels_t active_levels(level, level + 1);
          for (int vi = 0; vi < finegroupdata.numvars; ++vi)
            check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
              return "Reflux after refluxing: Fine level data";
            });
        }
      } // for gi
    } // if level exists
  } // for patchdata
}

// =======================================================================
// Restrict
// =======================================================================

static void Restrict_impl(const cGH *cctkGH, int level,
                          const std::vector<int> &groups,
                          const bool do_validity_tracking) {
  DECLARE_CCTK_PARAMETERS;

  assert(do_restrict);

  static Timer timer("Restrict");
  Interval interval(timer);

  int gi_regrid_error = -1;
  gi_regrid_error = CCTK_GroupIndex("CarpetXRegrid::regrid_error");
  assert(gi_regrid_error >= 0);

  for (const auto &patchdata : ghext->patchdata) {
    const int patch = patchdata.patch;
    if (level + 1 < int(patchdata.leveldata.size())) {
      auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      const active_levels_t active_levels(level, level + 1, patch, patch + 1);
      const active_levels_t active_fine_levels(level + 1, level + 2, patch,
                                               patch + 1);

      for (const int gi : groups) {
        cGroup group;
        int ierr = CCTK_GroupData(gi, &group);
        assert(!ierr);

        assert(group.grouptype == CCTK_GF);

        auto &groupdata = *leveldata.groupdata.at(gi);
        assert(!groupdata.mfab.empty());
        const auto &finegroupdata = *fineleveldata.groupdata.at(gi);
        assert(!finegroupdata.mfab.empty());
        const amrex::IntVect reffact{2, 2, 2};
        const nan_handling_t nan_handling = groupdata.do_evolve
                                                ? nan_handling_t::forbid_nans
                                                : nan_handling_t::allow_nans;

        // Don't restrict the regridding error
        if (gi == gi_regrid_error)
          continue;

        // If there is more than one time level, then we don't restrict the
        // oldest.
        // TODO: during evolution, restrict only one time level
        int ntls = groupdata.mfab.size();
        int restrict_tl = ntls > 1 ? ntls - 1 : ntls;
        for (int tl = 0; tl < restrict_tl; ++tl) {

          for (int vi = 0; vi < groupdata.numvars; ++vi) {

            // Restriction only uses the interior
            error_if_invalid(finegroupdata, vi, tl, make_valid_int(), []() {
              return "Restrict on fine level before restricting";
            });
            poison_invalid_gf(active_fine_levels, gi, vi, tl);
            check_valid_gf(active_fine_levels, gi, vi, tl, nan_handling, []() {
              return "Restrict on fine level before restricting";
            });
            error_if_invalid(groupdata, vi, tl, make_valid_int(), []() {
              return "Restrict on coarse level before restricting";
            });
            poison_invalid_gf(active_levels, gi, vi, tl);
            check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
              return "Restrict on coarse level before restricting";
            });
          }

#if 1
          {
            static Timer timer("Restrict::average_down");
            Interval interval(timer);
#warning                                                                       \
    "TODO: Allow different restriction operators, and ensure this is conservative"
            // rank: 0: vertex, 1: edge, 2: face, 3: volume
            int rank = 0;
            for (int d = 0; d < dim; ++d)
              rank += groupdata.indextype.at(d);
            switch (rank) {
            case 0:
              average_down_nodal(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 1:
              average_down_edges(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 2:
              average_down_faces(*finegroupdata.mfab.at(tl),
                                 *groupdata.mfab.at(tl), reffact);
              break;
            case 3:
              average_down(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl),
                           0, groupdata.numvars, reffact);
              break;
            default:
              assert(0);
            }
          }
#endif

          if (do_validity_tracking) {
            // TODO: Also remember old why_valid for interior?
            for (int vi = 0; vi < groupdata.numvars; ++vi) {
              // Should we mark ghosts and maybe outer boundaries as
              // valid as well?
              groupdata.valid.at(tl).at(vi).set_invalid(
                  make_valid_outer() | make_valid_ghosts(),
                  []() { return "Restrict"; });
              poison_invalid_gf(active_levels, gi, vi, tl);
              check_valid_gf(active_levels, gi, vi, tl, nan_handling, []() {
                return "Restrict on coarse level after restricting";
              });
            }
          }

        } // for tl
      } // for gi
    } // if level exists
  } // for patchdata
}

void Restrict(const cGH *cctkGH, int level, const std::vector<int> &groups) {
  Restrict_impl(cctkGH, level, groups, /*do_validity_tracking=*/true);
}

void RestrictNoPoison(const cGH *cctkGH, int level, const std::vector<int> &groups) {
  Restrict_impl(cctkGH, level, groups, /*do_validity_tracking=*/false);
}

void Restrict(const cGH *cctkGH, int level) {
  const std::vector<int> groups = collect_restrictable_groups();
  if (ghext->use_subcycling) {
    RestrictNoPoison(cctkGH, level, groups);
  } else {
    Restrict(cctkGH, level, groups);
  }
}

} // namespace CarpetX
