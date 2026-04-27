#include "solve.hxx"
#include <subcycling.hxx>

namespace ODESolvers {

constexpr int rkstages = 4;

namespace {

struct solve_setup_t {
  statecomp_t var, rhs, old;
  std::array<statecomp_t, rkstages> ks;
  std::vector<int> var_groups, rhs_groups, dep_groups, old_groups;
  std::array<std::vector<int>, rkstages> ks_groups;
  int nvars = 0;
};

// Collect evolved groups and their old=/ks= companions into statecomp_t
// bundles. Operates on CarpetX::active_levels; no cGH needed.
solve_setup_t collect_solve_setup() {
  const int tl = 0;
  solve_setup_t s;
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
        s.var.groupdatas.push_back(&groupdata);
        s.var.mfabs.push_back(groupdata.mfab.at(tl).get());
        s.rhs.groupdatas.push_back(&rhs_groupdata);
        s.rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
        auto &old_groupdata = *leveldata.groupdata.at(old_gi);
        s.old.groupdatas.push_back(&old_groupdata);
        s.old.mfabs.push_back(old_groupdata.mfab.at(tl).get());
        for (int i = 0; i < rkstages; i++) {
          auto &ki_groupdata = *leveldata.groupdata.at(ks_gi[i]);
          s.ks[i].groupdatas.push_back(&ki_groupdata);
          s.ks[i].mfabs.push_back(ki_groupdata.mfab.at(tl).get());
        }

        if (do_accumulate_nvars) {
          s.nvars += groupdata.numvars;
          s.var_groups.push_back(groupdata.groupindex);
          s.rhs_groups.push_back(rhs_gi);
          s.old_groups.push_back(old_gi);
          for (int i = 0; i < rkstages; i++) {
            s.ks_groups[i].push_back(ks_gi[i]);
          }
          const auto &dependents = get_group_dependents(groupdata.groupindex);
          s.dep_groups.insert(s.dep_groups.end(), dependents.begin(),
                              dependents.end());
        }
      }
    }
    do_accumulate_nvars = false;
  });

  {
    std::sort(s.var_groups.begin(), s.var_groups.end());
    const auto last = std::unique(s.var_groups.begin(), s.var_groups.end());
    assert(last == s.var_groups.end());
  }

  {
    std::sort(s.rhs_groups.begin(), s.rhs_groups.end());
    const auto last = std::unique(s.rhs_groups.begin(), s.rhs_groups.end());
    assert(last == s.rhs_groups.end());
  }

  {
    std::sort(s.old_groups.begin(), s.old_groups.end());
    const auto last = std::unique(s.old_groups.begin(), s.old_groups.end());
    assert(last == s.old_groups.end());
  }

  for (int i = 0; i < rkstages; i++) {
    std::sort(s.ks_groups[i].begin(), s.ks_groups[i].end());
    const auto last = std::unique(s.ks_groups[i].begin(), s.ks_groups[i].end());
    assert(last == s.ks_groups[i].end());
  }

  // Add RHS variables to dependent variables
  s.dep_groups.insert(s.dep_groups.end(), s.rhs_groups.begin(),
                      s.rhs_groups.end());

  {
    std::sort(s.dep_groups.begin(), s.dep_groups.end());
    const auto last = std::unique(s.dep_groups.begin(), s.dep_groups.end());
    s.dep_groups.erase(last, s.dep_groups.end());
  }

  for (const int gi : s.var_groups)
    assert(std::find(s.dep_groups.begin(), s.dep_groups.end(), gi) ==
           s.dep_groups.end());
  for (const int gi : s.rhs_groups)
    assert(std::find(s.var_groups.begin(), s.var_groups.end(), gi) ==
           s.var_groups.end());

  return s;
}

} // namespace

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

  static Timer timer_setup("ODESolvers::Solve::setup");
  std::optional<Interval> interval_setup(timer_setup);

  auto setup = collect_solve_setup();
  auto &var = setup.var;
  auto &rhs = setup.rhs;
  auto &old = setup.old;
  auto &ks = setup.ks;
  auto &var_groups = setup.var_groups;
  auto &rhs_groups = setup.rhs_groups;
  auto &dep_groups = setup.dep_groups;
  auto &old_groups = setup.old_groups;
  auto &ks_groups = setup.ks_groups;
  const int nvars = setup.nvars;
  if (verbose)
    CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

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
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
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
      // mark interior valid to work around poison mechanism
      old.set_valid(make_valid_int());
      SyncGroupsByDirIProlongateOnly(cctkGH, old_groups.size(),
                                     old_groups.data(), nullptr);
      for (int s = 0; s < rkstages; ++s) {
        // mark interior valid to work around poison mechanism
        ks[s].set_valid(make_valid_int());
        SyncGroupsByDirIProlongateOnly(cctkGH, ks_groups[s].size(),
                                       ks_groups[s].data(), nullptr);
      }
    }

    // Grid functions used to fill the refinement boundary substate.
    // Temporary variables cannot be used for old values here because
    // they need to be accessed in subsequent CallScheduleGroup functions,
    // which do not yet support access to temporary variables.
    setold();

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

    // No calcys_rmbnd(1) here: refinement-boundary ghosts are kept aligned by
    // subcycling-aware POSTRESTRICT SYNCs. The post-recovery case is handled
    // by ODESolvers_Solve_Subcycling_Recovery at CCTK_CPINITIAL.

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

extern "C" void ODESolvers_Solve_Subcycling_Recovery(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_Solve_Subcycling_Recovery;
  DECLARE_CCTK_PARAMETERS;

  // Skip on fresh initialization; cctk_iteration > 0 only on recovery.
  if (cctk_iteration <= 0)
    return;

  if (verbose)
    CCTK_VINFO("Subcycling recovery: replaying calcys_rmbnd(1) prerequisites "
               "for the first post-recovery RK4 step");

  static Timer timer("ODESolvers::Solve_Subcycling_Recovery");
  Interval interval(timer);

  auto setup = collect_solve_setup();
  auto &var = setup.var;
  auto &old = setup.old;
  auto &ks = setup.ks;
  auto &var_groups = setup.var_groups;
  auto &old_groups = setup.old_groups;
  auto &ks_groups = setup.ks_groups;
  if (setup.nvars == 0)
    return;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  // Prolongate old and ks from the coarse level: the checkpoint restored only
  // GF interiors, so refinement-boundary ghosts must be rebuilt here.
  if (old_groups.size() > 0) {
    // mark interior valid to work around poison mechanism
    old.set_valid(make_valid_int());
    SyncGroupsByDirIProlongateOnly(cctkGH, old_groups.size(), old_groups.data(),
                                   nullptr);
    for (int s = 0; s < rkstages; ++s) {
      // mark interior valid to work around poison mechanism
      ks[s].set_valid(make_valid_int());
      SyncGroupsByDirIProlongateOnly(cctkGH, ks_groups[s].size(),
                                     ks_groups[s].data(), nullptr);
    }
  }

  // calcys_rmbnd(1) body: fill refinement-boundary ghosts of var_groups.
  // xsi=0.5 and stage0=1 are the values calcys_rmbnd(1) would have computed
  // at the start of the first post-recovery RK step.
  active_levels->loop_parallel([&](int patch, int level, int /*index*/,
                                   int /*component*/, const cGH *local_cctkGH) {
    if (level == 0)
      return;

    // When this level is time-aligned with its parent, SyncRestrictGFs fills
    // the refinement-boundary ghosts; skip CalcYfFromKcs to avoid overwriting
    // them with a stale interpolation from old/k_i.
    const auto &patchdata = ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    const auto &prev_leveldata = patchdata.leveldata.at(level - 1);
    if (leveldata.iteration == prev_leveldata.iteration)
      return;

    constexpr CCTK_REAL xsi = 0.5;
    constexpr int stage0 = 1;
    update_cctkGH(const_cast<cGH *>(local_cctkGH), cctkGH);
    Subcycling::CalcYfFromKcs<rkstages>(const_cast<cGH *>(local_cctkGH),
                                        var_groups, old_groups, ks_groups,
                                        dt * 2, xsi, stage0);
  });
  synchronize();
  var.set_valid(make_valid_all());
}

} // namespace ODESolvers
