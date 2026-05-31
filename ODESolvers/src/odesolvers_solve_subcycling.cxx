#include "solve.hxx"
#include <subcycling.hxx>

// For FillPatch_ProlongateToBand (band->band prolongation of the RK k-stages).
// TODO: Don't include files from other thorns; create a proper interface.
#include "../../CarpetX/src/fillpatch.hxx"

#include <AMReX_MultiFab.H>

namespace ODESolvers {

constexpr int max_num_rk_stages = 4;

namespace {

struct solve_setup_t {
  statecomp_t var, rhs;
  std::vector<int> var_groups, rhs_groups, dep_groups;
  int nvars = 0;
};

// Collect evolved groups into statecomp_t bundles. The old-state anchor is now
// a scratch copy of var(tl=0) made by the solver (no extra timelevel); the RK
// k-stages live as coarse-fine bands on each group's GroupData. Operates on
// CarpetX::active_levels; no cGH needed.
solve_setup_t collect_solve_setup() {
  solve_setup_t s;
  s.var.timelevel = 0;
  s.rhs.timelevel = 0;
  bool do_accumulate_nvars = true;
  assert(CarpetX::active_levels);
  CarpetX::active_levels->loop_serially([&](const auto &leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      // TODO: add support for evolving grid scalars
      if (groupdataptr == nullptr)
        continue;

      auto &groupdata = *groupdataptr;
      const int rhs_gi = get_group_rhs(groupdata.groupindex);
      if (rhs_gi >= 0) {
        assert(rhs_gi != groupdata.groupindex);
        auto &rhs_groupdata = *leveldata.groupdata.at(rhs_gi);
        assert(rhs_groupdata.numvars == groupdata.numvars);
        s.var.groupdatas.push_back(&groupdata);
        s.var.mfabs.push_back(groupdata.mfab.at(0).get());
        s.rhs.groupdatas.push_back(&rhs_groupdata);
        s.rhs.mfabs.push_back(rhs_groupdata.mfab.at(0).get());

        if (do_accumulate_nvars) {
          s.nvars += groupdata.numvars;
          s.var_groups.push_back(groupdata.groupindex);
          s.rhs_groups.push_back(rhs_gi);
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
  auto &var_groups = setup.var_groups;
  auto &rhs_groups = setup.rhs_groups;
  auto &dep_groups = setup.dep_groups;
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

  // Allocate the coarse-fine RK k-stage bands single-threaded. The source-band
  // geometry reads the next-finer level, so this must run once all levels
  // exist; like build_cf_mask it opens its own MFIter/OpenMP region and so must
  // not run inside a parallel consume. This is the band allocation site that
  // setks, the band->band prolongation, and the dense-output kernel rely on.
  active_levels->loop_serially([&](const auto &leveldata) {
    for (const int gi : var_groups) {
      const auto &gd = *leveldata.groupdata.at(gi);
      leveldata.build_bands(gd);
    }
  });

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
  // calculate Ys from the k-stage bands and old on the mesh refinement boundary
  const auto calcys_rmbnd = [&](const int stage) {
    if (verbose)
      CCTK_VINFO(
          "Fill refinement boundary ghost zones using Ys for stage #%d at t=%g",
          stage, double(cctkGH->cctk_time));

    active_levels->loop_coarse_to_fine([&](auto &leveldata) {
      const int level = leveldata.level;
      if (level == 0)
        return;

      const auto &patchdata = ghext->patchdata.at(leveldata.patch);
      const auto &prev_leveldata = patchdata.leveldata.at(level - 1);
      // Virtual end-of-step stage = num_rk_stages + 1 (RK4 -> 5, SSPRK3 -> 4).
      const int virtual_end = ghext->num_rk_stages + 1;
      CCTK_REAL xsi =
          (leveldata.iteration == prev_leveldata.iteration) ? 0.5 : 0.0;
      if (stage == virtual_end) {
        xsi += 0.5;
      }
      const int stage0 = (stage == virtual_end ? 1 : stage);
      if (ghext->num_rk_stages == 3)
        Subcycling::CalcYfFromKcs_MFlevel<3>(leveldata, var_groups, /*Yf_tl=*/0,
                                             dt * 2, xsi, stage0);
      else
        Subcycling::CalcYfFromKcs_MFlevel<4>(leveldata, var_groups, /*Yf_tl=*/0,
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
    const int s = stage - 1;
    active_levels->loop_coarse_to_fine([&](const auto &restrict leveldata) {
      // rhs_groups[i] and var_groups[i] are paired by sort order.
      for (size_t i = 0; i < rhs_groups.size(); ++i) {
        const auto &rhs_groupdata = *leveldata.groupdata.at(rhs_groups[i]);
        // The k-stage bands live on the evolved group's GroupData.
        const auto &groupdata = *leveldata.groupdata.at(var_groups[i]);
        // The finest level has no source band (no children to prolongate to).
        if (!groupdata.ks_source_band[s])
          continue;
        auto &rhs_mf = *rhs_groupdata.mfab.at(0);
        auto &src_band = *groupdata.ks_source_band[s];
        assert(src_band.ixType() == rhs_mf.ixType());
        assert(src_band.nComp() == rhs_mf.nComp());
        // Fill the source band's interior from the RHS interior. The band is
        // zero-ghost and feeds band->band prolongation, which pulls from valid
        // interior cells, so the old same-level FillBoundary is unnecessary.
        src_band.ParallelCopy(rhs_mf, 0, 0, src_band.nComp(), amrex::IntVect{0},
                              amrex::IntVect{0},
                              amrex::Periodicity::NonPeriodic());
      }
    });
    synchronize();
  };
  // Capture u(t_n) = var(tl=0) into each level's old_source_band (interior
  // only), for prolongation into the children's old_consumer_band. Like setks
  // but from var(tl=0) once per step. The finest level has no source band.
  const auto fill_old_source_band = [&]() {
    active_levels->loop_coarse_to_fine([&](const auto &restrict leveldata) {
      for (const int gi : var_groups) {
        const auto &groupdata = *leveldata.groupdata.at(gi);
        if (!groupdata.old_source_band)
          continue;
        auto &var_mf = *groupdata.mfab.at(0);
        auto &src_band = *groupdata.old_source_band;
        assert(src_band.ixType() == var_mf.ixType());
        assert(src_band.nComp() == var_mf.nComp());
        src_band.ParallelCopy(var_mf, 0, 0, src_band.nComp(), amrex::IntVect{0},
                              amrex::IntVect{0},
                              amrex::Periodicity::NonPeriodic());
      }
    });
    synchronize();
  };
  // Prolongate the coarse-fine bands from the parent (band->band), filling each
  // fine level's consumer bands: the RK k-stage bands (written by setks) and
  // the single old-state band (written by fill_old_source_band). Both feed
  // calcys_rmbnd.
  const auto prolongate_bands = [&]() {
    active_levels->loop_coarse_to_fine([&](const auto &restrict leveldata) {
      const int level = leveldata.level;
      if (level == 0)
        return;
      const auto &patchdata = ghext->patchdata.at(leveldata.patch);
      const auto &coarseleveldata = patchdata.leveldata.at(level - 1);
      // Prolongate the bands are redundant for the second fine iteration.
      if (leveldata.iteration == coarseleveldata.iteration)
        return;
      const auto &fgeom = patchdata.amrcore->Geom(level);
      const auto &cgeom = patchdata.amrcore->Geom(level - 1);
      for (size_t i = 0; i < var_groups.size(); ++i) {
        const auto &groupdata = *leveldata.groupdata.at(var_groups[i]);
        const auto &coarsegroupdata =
            *coarseleveldata.groupdata.at(var_groups[i]);
        amrex::Interpolater *const interpolator = groupdata.interpolator;
        for (int s = 0; s < max_num_rk_stages; ++s) {
          // Skip where either band is absent (level 0 has no consumer band, the
          // finest level no source band).
          if (!groupdata.ks_consumer_band[s] ||
              !coarsegroupdata.ks_source_band[s])
            continue;
          CarpetX::FillPatch_ProlongateToBand(
              groupdata, coarsegroupdata, *groupdata.ks_consumer_band[s],
              *coarsegroupdata.ks_source_band[s], fgeom, cgeom, interpolator,
              groupdata.bcrecs);
        }
        // Old-state band (single), reusing the same band->band helper.
        if (groupdata.old_consumer_band && coarsegroupdata.old_source_band)
          CarpetX::FillPatch_ProlongateToBand(
              groupdata, coarsegroupdata, *groupdata.old_consumer_band,
              *coarsegroupdata.old_source_band, fgeom, cgeom, interpolator,
              groupdata.bcrecs);
      }
    });
    synchronize();
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

    // Scratch copy of u(t_n) = var(tl=0), the RK4 interior anchor y0. At one
    // timelevel var(tl=0) holds the previous step's result, so no init copy is
    // needed (mirrors the non-subcycling solver).
    const auto old = var.copy(make_valid_all());

    // Capture u(t_n) into the old bands before the RK stages overwrite var,
    // then prolongate the parent's bands (old + k-stage) into this level's
    // consumer bands.
    if (var_groups.size() > 0) {
      fill_old_source_band();
      prolongate_bands();
    }

    // k1 = f(Y1)
    calcrhs(1);
    setks(1); // interior only
    const auto kaccum = rhs.copy(make_valid_int());
    calcupdate(1, dt / 2, 1.0, reals<1>{dt / 2}, states<1>{&rhs});
    calcys_rmbnd(2); // refinement boundary only
    calcpoststep();

    // k2 = f(Y2)
    calcrhs(2);
    setks(2); // interior only
    statecomp_t::lincomb(kaccum, 1.0, reals<1>{2.0}, states<1>{&rhs},
                         make_valid_int());
    calcupdate(2, dt / 2, 0.0, reals<2>{1.0, dt / 2}, states<2>{&old, &rhs});
    calcys_rmbnd(3); // refinement boundary only
    calcpoststep();

    // k3 = f(Y3)
    calcrhs(3);
    setks(3); // interior only
    statecomp_t::lincomb(kaccum, 1.0, reals<1>{2.0}, states<1>{&rhs},
                         make_valid_int());
    calcupdate(3, dt, 0.0, reals<2>{1.0, dt}, states<2>{&old, &rhs});
    calcys_rmbnd(4); // refinement boundary only
    calcpoststep();

    // k4 = f(Y4)
    calcrhs(4);
    setks(4); // interior only
    calcupdate(4, dt, 0.0, reals<3>{1.0, dt / 6, dt / 6},
               states<3>{&old, &kaccum, &rhs});
    calcys_rmbnd(5); // refinement boundary only
    calcpoststep();

    // No calcys_rmbnd(1) here: refinement-boundary ghosts are kept aligned by
    // subcycling-aware POSTRESTRICT SYNCs. The post-recovery case is handled
    // by ODESolvers_Solve_Subcycling_Recovery at CCTK_CPINITIAL.

  } else if (CCTK_EQUALS(method, "SSPRK3")) {

    // k1 = f(y0)
    // k2 = f(y0 + h k1)
    // k3 = f(y0 + h/4 k1 + h/4 k2)
    // y1 = y0 + h/6 k1 + h/6 k2 + 2/3 h k3

    assert(ghext->num_rk_stages == 3);

    // Scratch copy of u(t_n) = var(tl=0), the SSPRK3 interior anchor y0.
    const auto old = var.copy(make_valid_all());

    // Capture u(t_n) into the old bands before the RK stages overwrite var,
    // then prolongate the parent's bands (old + k-stage) into this level's
    // consumer bands.
    if (var_groups.size() > 0) {
      fill_old_source_band();
      prolongate_bands();
    }

    // k1 = f(Y1)
    calcrhs(1);
    setks(1); // interior only
    const auto k1 = rhs.copy(make_valid_int());
    calcupdate(1, dt, 1.0, reals<1>{dt}, states<1>{&rhs}); // var = y0 + dt*k1
    calcys_rmbnd(2); // refinement boundary only
    calcpoststep();

    // k2 = f(Y2)
    calcrhs(2);
    setks(2); // interior only
    const auto k2 = rhs.copy(make_valid_int());
    calcupdate(2, dt / 2, 0.0, reals<3>{1.0, dt / 4, dt / 4},
               states<3>{&old, &k1, &k2});
    calcys_rmbnd(3); // refinement boundary only
    calcpoststep();

    // k3 = f(Y3)
    calcrhs(3);
    setks(3); // interior only
    calcupdate(3, dt, 0.0, reals<4>{1.0, dt / 6, dt / 6, 2 * dt / 3},
               states<4>{&old, &k1, &k2, &rhs});
    calcys_rmbnd(4); // virtual end-of-step (num_rk_stages + 1)
    calcpoststep();

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
    CCTK_VINFO("Subcycling recovery: refilling refinement-boundary ghosts "
               "(spatial prolongation on time-aligned levels, dense output "
               "from restored consumer bands otherwise)");

  static Timer timer("ODESolvers::Solve_Subcycling_Recovery");
  Interval interval(timer);

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  auto setup = collect_solve_setup();
  auto &var = setup.var;
  auto &var_groups = setup.var_groups;
  if (setup.nvars == 0)
    return;

  // Refill each recovered fine level's refinement-boundary (cf) ghosts.
  // Synchronized levels (and old/synchronized checkpoints with no restored band
  // data) use spatial tl=0 prolongation. Unsynchronized levels carry mid-cycle
  // dense-output state in their restored consumer bands and reconstruct the
  // cf-ghosts the uninterrupted run's previous fine substep last wrote.
  if (var_groups.size() > 0) {
    var.check_valid(make_valid_int(),
                    "ODESolvers_Solve_Subcycling_Recovery requires the tl=0 "
                    "interior to be populated by checkpoint recovery");

    // Spatial prolongation fills every fine level's cf-ghosts; the
    // unsynchronized levels below are then overwritten with dense output.
    SyncGroupsByDirIProlongateOnly(cctkGH, var_groups.size(), var_groups.data(),
                                   nullptr, /*tl=*/0);

    active_levels->loop_coarse_to_fine([&](auto &restrict leveldata) {
      const int level = leveldata.level;
      if (level == 0)
        return;
      const auto &patchdata = ghext->patchdata.at(leveldata.patch);
      const auto &prev_leveldata = patchdata.leveldata.at(level - 1);
      // Time-aligned with the parent: spatial prolongation above is correct.
      if (leveldata.iteration == prev_leveldata.iteration)
        return;
      // Reconstruct only where restored consumer-band data is present (a band
      // rebuilt but never read falls back to the spatial path).
      bool have_bands = false;
      for (const int gi : var_groups) {
        const auto &groupdata = *leveldata.groupdata.at(gi);
        if (groupdata.ks_consumer_band[0] &&
            !groupdata.ks_consumer_band[0]->empty()) {
          have_bands = true;
          break;
        }
      }
      if (!have_bands)
        return;
      // Mirror the previous fine substep's calcys_rmbnd at the virtual
      // end-of-step: base offset 0.0 plus the +0.5 give xsi = 0.5, stage0 = 1,
      // dtc = dt*2.
      const CCTK_REAL xsi = 0.5;
      if (ghext->num_rk_stages == 3)
        Subcycling::CalcYfFromKcs_MFlevel<3>(leveldata, var_groups, /*Yf_tl=*/0,
                                             dt * 2, xsi, /*stage0=*/1);
      else
        Subcycling::CalcYfFromKcs_MFlevel<4>(leveldata, var_groups, /*Yf_tl=*/0,
                                             dt * 2, xsi, /*stage0=*/1);
    });
    synchronize();
    var.set_valid(make_valid_all());
  }
}

extern "C" void ODESolvers_CheckTimelevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_CheckTimelevels;
  DECLARE_CCTK_PARAMETERS;

  for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
    if (get_group_rhs(gi) < 0)
      continue; // not an ODE-evolved group
    const int ntls = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
    if (ntls >= 2)
      CCTK_VERROR("ODESolvers subcycling requires evolved groups to have a "
                  "single timelevel, but group \"%s\" has %d active "
                  "timelevels. Subcycling does not support timelevels >= 2 "
                  "for evolution variables.",
                  CCTK_FullGroupName(gi), ntls);
  }
}

} // namespace ODESolvers
