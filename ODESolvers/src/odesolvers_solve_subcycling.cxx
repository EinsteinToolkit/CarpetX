#include "solve.hxx"
#include <subcycling.hxx>

#include <AMReX_MultiFab.H>

namespace ODESolvers {

constexpr int rkstages = 4;

namespace {

struct solve_setup_t {
  statecomp_t var, rhs, old;
  std::array<statecomp_t, rkstages> ks;
  std::vector<int> var_groups, rhs_groups, dep_groups;
  std::array<std::vector<int>, rkstages> ks_groups;
  int nvars = 0;
};

// Collect evolved groups and their ks= companions into statecomp_t bundles.
// The previous-step "old" anchor lives at tl=1 of the evolved group itself.
// Operates on CarpetX::active_levels; no cGH needed.
solve_setup_t collect_solve_setup() {
  solve_setup_t s;
  s.var.timelevel = 0;
  s.rhs.timelevel = 0;
  s.old.timelevel = 1;
  for (int i = 0; i < rkstages; i++)
    s.ks[i].timelevel = 0;
  bool do_accumulate_nvars = true;
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      // TODO: add support for evolving grid scalars
      if (groupdataptr == nullptr)
        continue;

      auto &groupdata = *groupdataptr;
      const int rhs_gi = get_group_rhs(groupdata.groupindex);
      const auto &ks_gi = get_group_ks<int, rkstages>(groupdata.groupindex);
      if (rhs_gi >= 0) {
        if (int(groupdata.mfab.size()) < 2)
          CCTK_VERROR("Variable group \"%s\" declares rhs but has TIMELEVELS<2."
                      " Bump TIMELEVELS>=2 in interface.ccl for use with the "
                      "subcycling integrator.",
                      CCTK_FullGroupName(groupdata.groupindex));
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
        s.var.mfabs.push_back(groupdata.mfab.at(0).get());
        s.rhs.groupdatas.push_back(&rhs_groupdata);
        s.rhs.mfabs.push_back(rhs_groupdata.mfab.at(0).get());
        // old aliases the evolved groupdata at tl=1
        s.old.groupdatas.push_back(&groupdata);
        s.old.mfabs.push_back(groupdata.mfab.at(1).get());
        for (int i = 0; i < rkstages; i++) {
          auto &ki_groupdata = *leveldata.groupdata.at(ks_gi[i]);
          s.ks[i].groupdatas.push_back(&ki_groupdata);
          s.ks[i].mfabs.push_back(ki_groupdata.mfab.at(0).get());
        }

        if (do_accumulate_nvars) {
          s.nvars += groupdata.numvars;
          s.var_groups.push_back(groupdata.groupindex);
          s.rhs_groups.push_back(rhs_gi);
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
                                          var_groups, ks_groups, dt * 2, xsi,
                                          stage0);
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
    const int s = stage - 1;
    active_levels->loop_coarse_to_fine([&](const auto &restrict leveldata) {
      for (size_t i = 0; i < rhs_groups.size(); ++i) {
        const auto &rhs_groupdata = *leveldata.groupdata.at(rhs_groups[i]);
        const auto &k_groupdata = *leveldata.groupdata.at(ks_groups[s][i]);
        auto &rhs_mf = *rhs_groupdata.mfab.at(0);
        auto &k_mf = *k_groupdata.mfab.at(0);
        assert(k_mf.ixType() == rhs_mf.ixType());
        assert(k_mf.nComp() == rhs_mf.nComp());
        assert(k_mf.nGrowVect().allGE(amrex::IntVect(0)));
        assert(rhs_mf.nGrowVect().allGE(amrex::IntVect(0)));
        amrex::MultiFab::Copy(k_mf, rhs_mf, 0, 0, k_mf.nComp(),
                              amrex::IntVect(0));
      }
    });
    synchronize();

    // apply boundary condition to account for mesh refinement overlapping the
    // outer boundary
    SyncGroupsByDirIGhostOnly(cctkGH, ks_groups[s].size(), ks_groups[s].data(),
                              nullptr);
  };
  // populate var(tl=0) from var(tl=1) over the full FAB (interior + outer +
  // ghosts). This is the start-of-substep state that the RHS reads.
  const auto init_substep = [&]() {
    if (verbose)
      CCTK_VINFO("Init substep: copy var(tl=1) -> var(tl=0) at t=%g",
                 double(cctkGH->cctk_time));
    statecomp_t::lincomb(var, 0.0, reals<1>{1.0}, states<1>{&old},
                         make_valid_all());
    var.set_valid(make_valid_all());
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

    // Initialize the substep working state: var(tl=0) <- var(tl=1).
    init_substep();

    // Sync OldState and Ks: prolongate the previous-step anchor (var at tl=1)
    // and the ks scratch groups from the parent level, which were set in
    // previous steps.
    if (var_groups.size() > 0) {
      // mark interior valid to work around poison mechanism
      old.set_valid(make_valid_int());
      SyncGroupsByDirIProlongateOnly(cctkGH, var_groups.size(),
                                     var_groups.data(), nullptr, /*tl=*/1);
      for (int s = 0; s < rkstages; ++s) {
        // mark interior valid to work around poison mechanism
        ks[s].set_valid(make_valid_int());
        SyncGroupsByDirIProlongateOnly(cctkGH, ks_groups[s].size(),
                                       ks_groups[s].data(), nullptr);
      }
    }

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
  auto &ks_groups = setup.ks_groups;
  if (setup.nvars == 0)
    return;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  // Prolongate old (var at tl=1) and ks from the coarse level: the checkpoint
  // restored only GF interiors, so refinement-boundary ghosts must be rebuilt
  // here.
  if (var_groups.size() > 0) {
    // Verify the IO layer populated tl=1 interior on recovery; this assertion
    // catches regressions where a checkpoint read fails to mark every
    // timelevel valid.
    old.check_valid(make_valid_int(),
                    "ODESolvers_Solve_Subcycling_Recovery requires tl=1 "
                    "interior to be populated by checkpoint recovery");
    SyncGroupsByDirIProlongateOnly(cctkGH, var_groups.size(), var_groups.data(),
                                   nullptr, /*tl=*/1);
    for (int s = 0; s < rkstages; ++s) {
      ks[s].check_valid(make_valid_int(),
                        "ODESolvers_Solve_Subcycling_Recovery requires ks "
                        "interior to be populated by checkpoint recovery");
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
                                        var_groups, ks_groups, dt * 2, xsi,
                                        stage0);
  });
  synchronize();
  var.set_valid(make_valid_all());
}

} // namespace ODESolvers
