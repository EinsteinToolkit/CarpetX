#include "solve.hxx"
#include <subcycling.hxx>

namespace ODESolvers {

constexpr int rkstages = 4;

extern "C" void ODESolvers_Solve_Subcycling(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_Solve_Subcycling;
  DECLARE_CCTK_PARAMETERS;

  static bool did_output = false;
  if (verbose || !did_output)
    CCTK_VINFO("Integrator is %s", method);
  did_output = true;

  static Timer timer("ODESolvers::Solve");
  Interval interval(timer);

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const int tl = 0;

  static Timer timer_setup("ODESolvers::Solve::setup");
  std::optional<Interval> interval_setup(timer_setup);

  statecomp_t var, rhs, old;
  std::array<statecomp_t, rkstages> ks;
  std::vector<int> var_groups, rhs_groups, dep_groups, old_groups;
  std::array<std::vector<int>, rkstages> ks_groups;
  int nvars = 0;
  bool do_accumulate_nvars = true;
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      // TODO: add support for evolving grid scalars
      if (groupdataptr == nullptr)
        continue;

      auto &groupdata = *groupdataptr;
      const int rhs_gi = get_group_rhs(groupdata.groupindex);
      const int old_gi = get_group_old(groupdata.groupindex);
      const auto &ks_gi = get_group_ks<int, rkstages>(groupdata.groupindex);
      if (rhs_gi >= 0 && old_gi < 0)
        CCTK_VERROR("Variable group \"%s\" declares rhs but is missing "
                    "required subcycling tag old=\"...\". Add old=\"...\" "
                    "to the group's TAGS in interface.ccl.",
                    CCTK_FullGroupName(groupdata.groupindex));
      if (rhs_gi >= 0) {
        for (int i = 0; i < rkstages; i++) {
          if (ks_gi[i] < 0)
            CCTK_VERROR("Variable group \"%s\" declares rhs but is missing "
                        "required subcycling tag ks=\"...\". Add "
                        "ks=\"k1 k2 k3 k4\" to the group's TAGS in "
                        "interface.ccl.",
                        CCTK_FullGroupName(groupdata.groupindex));
        }
      }
      if (rhs_gi >= 0) {
        assert(rhs_gi != groupdata.groupindex);
        auto &rhs_groupdata = *leveldata.groupdata.at(rhs_gi);
        assert(rhs_groupdata.numvars == groupdata.numvars);
        var.groupdatas.push_back(&groupdata);
        var.mfabs.push_back(groupdata.mfab.at(tl).get());
        rhs.groupdatas.push_back(&rhs_groupdata);
        rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
        auto &old_groupdata = *leveldata.groupdata.at(old_gi);
        old.groupdatas.push_back(&old_groupdata);
        old.mfabs.push_back(old_groupdata.mfab.at(tl).get());
        for (int i = 0; i < rkstages; i++) {
          auto &ki_groupdata = *leveldata.groupdata.at(ks_gi[i]);
          ks[i].groupdatas.push_back(&ki_groupdata);
          ks[i].mfabs.push_back(ki_groupdata.mfab.at(tl).get());
        }

        if (do_accumulate_nvars) {
          nvars += groupdata.numvars;
          var_groups.push_back(groupdata.groupindex);
          rhs_groups.push_back(rhs_gi);
          old_groups.push_back(old_gi);
          for (int i = 0; i < rkstages; i++) {
            ks_groups[i].push_back(ks_gi[i]);
          }
          const auto &dependents = get_group_dependents(groupdata.groupindex);
          dep_groups.insert(dep_groups.end(), dependents.begin(),
                            dependents.end());
        }
      }
    }
    do_accumulate_nvars = false;
  });
  if (verbose)
    CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

  {
    std::sort(var_groups.begin(), var_groups.end());
    const auto last = std::unique(var_groups.begin(), var_groups.end());
    assert(last == var_groups.end());
  }

  {
    std::sort(rhs_groups.begin(), rhs_groups.end());
    const auto last = std::unique(rhs_groups.begin(), rhs_groups.end());
    assert(last == rhs_groups.end());
  }

  {
    std::sort(old_groups.begin(), old_groups.end());
    const auto last = std::unique(old_groups.begin(), old_groups.end());
    assert(last == old_groups.end());
  }

  for (int i = 0; i < rkstages; i++) {
    std::sort(ks_groups[i].begin(), ks_groups[i].end());
    const auto last = std::unique(ks_groups[i].begin(), ks_groups[i].end());
    assert(last == ks_groups[i].end());
  }

  // Add RHS variables to dependent variables
  dep_groups.insert(dep_groups.end(), rhs_groups.begin(), rhs_groups.end());

  {
    std::sort(dep_groups.begin(), dep_groups.end());
    const auto last = std::unique(dep_groups.begin(), dep_groups.end());
    dep_groups.erase(last, dep_groups.end());
  }

  for (const int gi : var_groups)
    assert(std::find(dep_groups.begin(), dep_groups.end(), gi) ==
           dep_groups.end());
  for (const int gi : rhs_groups)
    assert(std::find(var_groups.begin(), var_groups.end(), gi) ==
           var_groups.end());

  interval_setup.reset();

  {
    static Timer timer_alloc_temps("ODESolvers::Solve::alloc_temps");
    Interval interval_alloc_temps(timer_alloc_temps);
    statecomp_t::init_tmp_mfabs();
  }

  const CCTK_REAL saved_time = cctkGH->cctk_time;
  const CCTK_REAL old_time = cctkGH->cctk_time - dt;

  static Timer timer_lincomb("ODESolvers::Solve::lincomb");
  static Timer timer_rhs("ODESolvers::Solve::rhs");
  static Timer timer_poststep("ODESolvers::Solve::poststep");

  const auto calcrhs = [&](const int n) {
    Interval interval_rhs(timer_rhs);
    if (verbose)
      CCTK_VINFO("Calculating RHS #%d at t=%g", n, double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");
    rhs.check_valid(make_valid_int(),
                    "ODESolvers after calling ODESolvers_RHS");
  };
  // t = t_0 + c
  // var = a_0 * var + \Sum_i a_i * var_i
  const auto calcupdate = [&](const int n, const CCTK_REAL c,
                              const CCTK_REAL a0, const auto &as,
                              const auto &vars) {
    {
      Interval interval_lincomb(timer_lincomb);
      if (verbose)
        CCTK_VINFO("Calculated new state #%d at t=%g", n,
                   double(cctkGH->cctk_time));
      statecomp_t::lincomb(var, a0, as, vars, make_valid_int());
      var.check_valid(make_valid_int(),
                      "ODESolvers after defining new state vector");
      mark_invalid(dep_groups);
    }
    {
      Interval interval_poststep(timer_poststep);
      *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + c;
    }
  };
  // calling ODESolvers_PostStep Group
  const auto calcpoststep = [&]() {
    if (!use_odesolvers_poststep_during_rksubsteps) {
      CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
    }
  };
  // calculate Ys from ks and old on the mesh refinement boundary
  const auto calcys_rmbnd = [&](const int stage) {
    if (verbose)
      CCTK_VINFO(
          "Fill refinement boundary ghost zones using Ys for stage #%d at t=%g",
          stage, double(cctkGH->cctk_time));

    active_levels->loop_parallel([&](int patch, int level, int index,
                                     int component, const cGH *local_cctkGH) {
      if (level == 0)
        return;

      const auto &patchdata = ghext->patchdata.at(patch);
      const auto &leveldata = patchdata.leveldata.at(level);
      const auto &prev_leveldata = patchdata.leveldata.at(level - 1);
      CCTK_REAL xsi =
          (leveldata.iteration == prev_leveldata.iteration) ? 0.5 : 0.0;
      if (stage == 5) {
        xsi += 0.5;
      }
      const int stage0 = (stage == 5 ? 1 : stage);
      update_cctkGH(const_cast<cGH *>(local_cctkGH), cctkGH);
      Subcycling::CalcYfFromKcs<rkstages>(const_cast<cGH *>(local_cctkGH),
                                          var_groups, old_groups, ks_groups,
                                          dt * 2, xsi, stage0);
    });
    synchronize();
    var.set_valid(make_valid_all());
  };
  // set ks in the interior which will be used for prolongation later
  const auto setks = [&](const int stage) {
    if (verbose)
      CCTK_VINFO(
          "Set interior Ks for stage #%d at t=%g, to be prolongated later",
          stage, double(cctkGH->cctk_time));
    active_levels->loop_parallel([&](int patch, int level, int index,
                                     int component, const cGH *local_cctkGH) {
      update_cctkGH(const_cast<cGH *>(local_cctkGH), cctkGH);
      Subcycling::SetK<rkstages>(const_cast<cGH *>(local_cctkGH), ks_groups,
                                 rhs_groups, stage);
    });
    synchronize();

    // apply boundary condition to account for mesh refinement overlapping the
    // outer boundary
    const int s = stage - 1;
    SyncGroupsByDirIGhostOnly(cctkGH, ks_groups[s].size(), ks_groups[s].data(),
                              nullptr);
  };
  // set old in the interior which will be used for prolongation later
  const auto setold = [&]() {
    if (verbose)
      CCTK_VINFO("Set interior old state at t=%g, to be prolongated later",
                 double(cctkGH->cctk_time));
    active_levels->loop_parallel([&](int patch, int level, int index,
                                     int component, const cGH *local_cctkGH) {
      update_cctkGH(const_cast<cGH *>(local_cctkGH), cctkGH);
      Subcycling::SetOld(const_cast<cGH *>(local_cctkGH), old_groups,
                         var_groups);
    });
    synchronize();

    // apply boundary condition to account for mesh refinement overlapping the
    // outer boundary
    SyncGroupsByDirIGhostOnly(cctkGH, old_groups.size(), old_groups.data(),
                              nullptr);
  };

  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;

  if (CCTK_EQUALS(method, "constant")) {

    // y1 = y0

    // do nothing

  } else if (CCTK_EQUALS(method, "RK4")) {

    // k1 = f(y0)
    // k2 = f(y0 + h/2 k1)
    // k3 = f(y0 + h/2 k2)
    // k4 = f(y0 + h k3)
    // y1 = y0 + h/6 k1 + h/3 k2 + h/3 k3 + h/6 k4

    // Sync OldState and Ks: prolongate old and ks from parent level which are
    // set in previous steps.
    // CallScheduleGroup(cctkGH, "ODESolvers_SyncKsOld");
    if (old_groups.size() > 0) {
      old.set_valid(make_valid_int()); // mark interior valid to work around
                                       // poison mechanism
      SyncGroupsByDirIProlongateOnly(cctkGH, old_groups.size(),
                                     old_groups.data(), nullptr);
      for (int s = 0; s < rkstages; ++s) {
        ks[s].set_valid(make_valid_int()); // mark interior valid to work around
                                           // poison mechanism
        SyncGroupsByDirIProlongateOnly(cctkGH, ks_groups[s].size(),
                                       ks_groups[s].data(), nullptr);
      }
    }

    // Grid functions used to fill the refinement boundary substate.
    // Temporary variables cannot be used for old values here because
    // they need to be accessed in subsequent CallScheduleGroup functions,
    // which do not yet support access to temporary variables.
    setold();

    calcys_rmbnd(1); // refinement boundary only

    // k1 = f(Y1)
    calcrhs(1);
    setks(1); // interior only
    calcupdate(1, dt / 2, 1.0, reals<1>{dt / 2}, states<1>{&rhs});
    calcys_rmbnd(2); // refinement boundary only
    calcpoststep();

    // k2 = f(Y2)
    calcrhs(2);
    setks(2); // interior only
    calcupdate(2, dt / 2, 0.0, reals<2>{1.0, dt / 2}, states<2>{&old, &rhs});
    calcys_rmbnd(3); // refinement boundary only
    calcpoststep();

    // k3 = f(Y3)
    calcrhs(3);
    setks(3); // interior only
    calcupdate(3, dt, 0.0, reals<2>{1.0, dt}, states<2>{&old, &rhs});
    calcys_rmbnd(4); // refinement boundary only
    calcpoststep();

    // k4 = f(Y4)
    calcrhs(4);
    setks(4); // interior only
    calcupdate(4, dt, 0.0, reals<5>{1.0, dt / 6, dt / 3, dt / 3, dt / 6},
               states<5>{&old, &ks[0], &ks[1], &ks[2], &ks[3]});
    calcys_rmbnd(5); // refinement boundary only
    calcpoststep();

    // In the interprocess_ghost_sync_during_substep case, the refinement
    // boundary is not synchronized at this point. Instead, we rely on the SYNCs
    // in CCTK_POSTRESTRICT to fill in these ghost zones. If the fine level is
    // not properly aligned with the coarse level in time, the filled-in values
    // may be incorrect. However, they will be corrected during calcys_rmbnd(1).

  } else {
    assert(0);
  }

  {
    static Timer timer_free_temps("ODESolvers::Solve::free_temps");
    Interval interval_free_temps(timer_free_temps);
    statecomp_t::free_tmp_mfabs();
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolvers
