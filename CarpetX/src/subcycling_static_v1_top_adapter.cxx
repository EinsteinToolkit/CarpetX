#include "subcycling_static_v1_top_adapter.hxx"

#include "driver.hxx"
#include "hierarchy_stepper.hxx"
#include "schedule.hxx"
#include "subcycling_native_gate.hxx"
#include "subcycling_runtime_clock.hxx"
#include "subcycling_schedule_certification.hxx"
#include "subcycling_schedule_state.hxx"
#include "subcycling_scratch_adapter.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"
#include "subcycling_stage_spatial_preparer.hxx"
#include "subcycling_step_context.hxx"
#include "transaction_level_step_session.hxx"

#include <cctk.h>
#include <cctk_ActiveThorns.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>

#include <AMReX_IntVect.H>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {
namespace {

using LevelData = GHExt::PatchData::LevelData;
using GroupData = GHExt::PatchData::LevelData::GroupData;

step_clock_t exact_clock(const rat64 value) {
  return step_clock_t(value.num, value.den);
}

rat64 level_clock(const step_clock_t value) {
  return rat64(value.num, value.den);
}

bool ratio_is_two(const amrex::IntVect &ratio) noexcept {
  for (int direction = 0; direction < AMREX_SPACEDIM; ++direction)
    if (ratio[direction] != 2)
      return false;
  return true;
}

bool same_step_context(const StepContext &left,
                       const StepContext &right) noexcept {
  return left.level == right.level &&
         left.begin_clock == right.begin_clock &&
         left.end_clock == right.end_clock &&
         left.begin_time == right.begin_time &&
         left.end_time == right.end_time && left.method == right.method &&
         left.require_dense_output == right.require_dense_output &&
         left.endpoint_accepted_step == right.endpoint_accepted_step;
}

std::vector<int>
contract_groups(const SubcyclingGroupSchema &group_schema) {
  std::vector<int> groups;
  groups.reserve(2 * group_schema.ordered_group_pairs.size() +
                 group_schema.dependent_groups.size());
  for (const auto &pair : group_schema.ordered_group_pairs) {
    groups.push_back(pair.evolved_group);
    groups.push_back(pair.rhs_group);
  }
  groups.insert(groups.end(), group_schema.dependent_groups.begin(),
                group_schema.dependent_groups.end());
  std::sort(groups.begin(), groups.end());
  groups.erase(std::unique(groups.begin(), groups.end()), groups.end());
  return groups;
}

bool bounded_test_odesolvers2_schema(
    const SubcyclingGroupSchema &group_schema) {
  const int state = CCTK_GroupIndex("TestODESolvers2::state");
  const int rhs = CCTK_GroupIndex("TestODESolvers2::rhs");
  const int dependent =
      CCTK_GroupIndex("TestODESolvers2::gate_dependent");
  std::vector<int> expected_dependents{rhs, dependent};
  std::sort(expected_dependents.begin(), expected_dependents.end());
  const auto *const native_gate = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("native_subcycling_gate", "TestODESolvers2",
                        nullptr));
  return state >= 0 && rhs >= 0 && dependent >= 0 && native_gate != nullptr &&
         *native_gate == 0 &&
         group_schema.ordered_group_pairs.size() == 1 &&
         group_schema.ordered_group_pairs.front().evolved_group == state &&
         group_schema.ordered_group_pairs.front().rhs_group == rhs &&
         group_schema.dependent_groups == expected_dependents;
}

bool bounded_test_odesolvers2_active_thorns() {
  const int compiled_count = CCTK_NumCompiledThorns();
  if (compiled_count < 0)
    return false;
  std::vector<std::string> active_thorns;
  active_thorns.reserve(static_cast<std::size_t>(compiled_count));
  for (int index = 0; index < compiled_count; ++index) {
    const char *const thorn = CCTK_CompiledThorn(index);
    if (thorn == nullptr)
      return false;
    const int active = CCTK_IsThornActive(thorn);
    if (active < 0)
      return false;
    if (active != 0)
      active_thorns.emplace_back(thorn);
  }
  return has_exact_static_v1_test_odesolvers2_active_thorns(
      std::move(active_thorns));
}

bool cactus_group_is_gf_real(const int group_index) {
  cGroup group{};
  return group_index >= 0 && group_index < CCTK_NumGroups() &&
         CCTK_GroupData(group_index, &group) == 0 &&
         group.grouptype == CCTK_GF &&
         group.vartype == CCTK_VARIABLE_REAL && group.numvars > 0;
}

bool schema_is_supported(const GHExt &driver,
                         const SubcyclingGroupSchema &group_schema) {
  if (group_schema.ordered_group_pairs.empty())
    return false;

  const auto groups = contract_groups(group_schema);
  for (const int group : groups)
    if (!cactus_group_is_gf_real(group))
      return false;

  for (const auto &pair : group_schema.ordered_group_pairs) {
    cGroup evolved{};
    cGroup rhs{};
    if (pair.evolved_group == pair.rhs_group ||
        CCTK_GroupData(pair.evolved_group, &evolved) != 0 ||
        CCTK_GroupData(pair.rhs_group, &rhs) != 0 ||
        evolved.numvars != rhs.numvars || evolved.numtimelevels != 1)
      return false;
  }

  if (driver.patchdata.size() != 1 ||
      driver.patchdata.front().leveldata.size() != 2)
    return false;
  for (const auto &level : driver.patchdata.front().leveldata) {
    for (const int group : groups) {
      if (group < 0 ||
          static_cast<std::size_t>(group) >= level.groupdata.size() ||
          level.groupdata[static_cast<std::size_t>(group)] == nullptr)
        return false;
      const auto &data =
          *level.groupdata[static_cast<std::size_t>(group)];
      if (data.groupindex != group || data.mfab.empty() ||
          data.valid.empty() || data.mfab.size() != data.valid.size() ||
          data.numvars <= 0)
        return false;
      for (std::size_t tl = 0; tl < data.mfab.size(); ++tl)
        if (data.mfab[tl] == nullptr || !data.mfab[tl]->isDefined() ||
            data.mfab[tl]->nComp() != data.numvars ||
            data.valid[tl].size() !=
                static_cast<std::size_t>(data.numvars) ||
            !data.mfab[tl]->boxArray().ok())
          return false;
    }
    for (const auto &pair : group_schema.ordered_group_pairs) {
      const auto &evolved = *level.groupdata.at(
          static_cast<std::size_t>(pair.evolved_group));
      if (evolved.mfab.size() != 1 || evolved.valid.size() != 1 ||
          evolved.valid.front().size() !=
              static_cast<std::size_t>(evolved.numvars))
        return false;
    }
  }
  return true;
}

struct GroupContractIdentity {
  int group_index;
  const GroupData *group_identity;
  int numvars;
  bool do_restrict;
  std::size_t timelevel_count;
  std::size_t validity_timelevel_count;
  std::vector<const void *> storage_identities;
  std::vector<const void *> validity_identities;
};

struct LevelContractIdentity {
  const LevelData *level_identity;
  const cGH *patch_gh_identity;
  const amrex::FabArrayBase *fab_identity;
  rat64 delta_iteration;
  bool is_subcycling_level;
  std::vector<GroupContractIdentity> groups;
};

LevelContractIdentity
capture_level_contract(const LevelData &level,
                       const std::vector<int> &groups) {
  LevelContractIdentity result{&level, level.get_patch_cctkGH(), level.fab.get(),
                               level.delta_iteration,
                               level.is_subcycling_level, {}};
  result.groups.reserve(groups.size());
  for (const int group : groups) {
    const auto &data =
        *level.groupdata.at(static_cast<std::size_t>(group));
    std::vector<const void *> storage_identities;
    std::vector<const void *> validity_identities;
    storage_identities.reserve(data.mfab.size());
    validity_identities.reserve(data.valid.size());
    for (const auto &storage : data.mfab)
      storage_identities.push_back(storage.get());
    for (const auto &validity : data.valid)
      validity_identities.push_back(validity.data());
    const auto canonicalize = [](auto &identities) {
      std::sort(identities.begin(), identities.end(),
                std::less<const void *>{});
    };
    canonicalize(storage_identities);
    canonicalize(validity_identities);
    result.groups.push_back({group,
                             &data,
                             data.numvars,
                             data.do_restrict,
                             data.mfab.size(),
                             data.valid.size(),
                             std::move(storage_identities),
                             std::move(validity_identities)});
  }
  return result;
}

struct StaticV1Preflight {
  GHExt *driver;
  SubcyclingMethodContractSnapshot method_snapshot;
  std::int64_t hierarchy_epoch;
  double initial_physical_time;
  double base_delta_time;
  std::vector<int> restricted_evolved_groups;
  std::array<LevelContractIdentity, 2> level_contracts;
};

StaticV1Preflight preflight_static_v1(cGH &root,
                                      const bool recovering) {
  DECLARE_CCTK_PARAMETERS;

  if (!ghext)
    throw std::invalid_argument("static-v1 requires an active CarpetX GHExt");
  if (!in_global_mode(&root))
    throw std::invalid_argument(
        "static-v1 requires the root GH to be in global mode");

  auto method_snapshot =
      require_complete_registered_subcycling_method_contract();
  validate_registered_subcycling_method_contract(method_snapshot);
  if (!method_snapshot.group_schema)
    throw std::logic_error(
        "static-v1 complete method snapshot lost its group schema");

  const int patch_count = static_cast<int>(ghext->patchdata.size());
  const int level_count =
      patch_count == 1
          ? static_cast<int>(ghext->patchdata.front().leveldata.size())
          : 0;
  int spatial_refinement_ratio = -1;
  if (patch_count == 1 && level_count == 2 &&
      ghext->patchdata.front().amrcore != nullptr &&
      ratio_is_two(ghext->patchdata.front().amrcore->refRatio(0)))
    spatial_refinement_ratio = 2;

  bool level_clocks_zero = level_count == 2;
  bool level_one_subcycles = false;
  if (level_count == 2) {
    const auto &levels = ghext->patchdata.front().leveldata;
    level_clocks_zero =
        exact_clock(levels[0].iteration) == step_clock_t(0) &&
        exact_clock(levels[1].iteration) == step_clock_t(0);
    level_one_subcycles = levels[1].is_subcycling_level;
  }
  const bool supported_schema =
      schema_is_supported(*ghext, *method_snapshot.group_schema);
  const bool bounded_configuration =
      method_snapshot.contract.provider_id == "ODESolvers::method" &&
      bounded_test_odesolvers2_schema(*method_snapshot.group_schema) &&
      bounded_test_odesolvers2_active_thorns();

  validate_static_v1_policy_envelope(StaticV1PolicyEnvelope{
      patch_count,
      level_count,
      spatial_refinement_ratio,
      regrid_every,
      recovering,
      bool(do_reflux),
      bool(do_restrict),
      bool(restrict_during_sync),
      bool(poison_undefined_values),
      level_one_subcycles,
      root.cctk_iteration == 0,
      level_clocks_zero,
      true,
      supported_schema,
      bounded_configuration});

  const auto &levels = ghext->patchdata.front().leveldata;
  for (std::size_t level = 0; level < levels.size(); ++level) {
    const auto &data = levels[level];
    if (data.patch != 0 || data.level != static_cast<int>(level) ||
        data.fab == nullptr || data.get_patch_cctkGH() == nullptr ||
        data.fab->boxArray().d_numPts() <= 0)
      throw std::invalid_argument(
          "static-v1 requires initialized, non-empty level and patch-GH "
          "storage for levels zero and one");
  }
  validate_static_v1_clock_envelope(StaticV1ClockEnvelope{
      2,
      2,
      recovering,
      regrid_every != 0,
      step_clock_t(0),
      {exact_clock(levels[0].iteration),
       exact_clock(levels[1].iteration)},
      {exact_clock(levels[0].delta_iteration),
       exact_clock(levels[1].delta_iteration)},
      0,
      {0, 0}});

  if (!std::isfinite(root.cctk_time))
    throw std::invalid_argument(
        "static-v1 initial physical time must be finite");
  const double base_delta_time =
      static_v1_base_delta_time(root.cctk_delta_time, level_count);
  const auto groups = contract_groups(*method_snapshot.group_schema);
  std::vector<int> restricted_evolved_groups;
  restricted_evolved_groups.reserve(
      method_snapshot.group_schema->ordered_group_pairs.size());
  for (const auto &pair :
       method_snapshot.group_schema->ordered_group_pairs) {
    const auto &coarse_evolved = *levels[0].groupdata.at(
        static_cast<std::size_t>(pair.evolved_group));
    const auto &fine_evolved = *levels[1].groupdata.at(
        static_cast<std::size_t>(pair.evolved_group));
    if (!coarse_evolved.do_restrict || !fine_evolved.do_restrict)
      throw std::invalid_argument(
          "static-v1 requires restriction enabled for every evolved group");
    if (std::find(restricted_evolved_groups.begin(),
                  restricted_evolved_groups.end(),
                  pair.evolved_group) != restricted_evolved_groups.end())
      throw std::invalid_argument(
          "static-v1 evolved group schema contains a duplicate group");
    restricted_evolved_groups.push_back(pair.evolved_group);
  }

  const auto hierarchy_epoch =
      static_cast<std::int64_t>(CarpetX_GetEpoch());
  if (hierarchy_epoch < 0)
    throw std::invalid_argument(
        "static-v1 hierarchy epoch must be non-negative");

  return StaticV1Preflight{
      ghext.get(),
      std::move(method_snapshot),
      hierarchy_epoch,
      root.cctk_time,
      base_delta_time,
      std::move(restricted_evolved_groups),
      {capture_level_contract(levels[0], groups),
       capture_level_contract(levels[1], groups)}};
}

struct BeginLevelSnapshot {
  std::array<rat64, 2> clocks;
  rat64 endpoint_clock;
};

class StaticV1EvolutionAdapter final : public HierarchyEvolutionAdapter {
public:
  StaticV1EvolutionAdapter(cGH &root, StaticV1Preflight preflight)
      : root_(root), driver_(preflight.driver),
        method_snapshot_(std::move(preflight.method_snapshot)),
        hierarchy_epoch_(preflight.hierarchy_epoch),
        initial_physical_time_(preflight.initial_physical_time),
        base_delta_time_(preflight.base_delta_time),
        restricted_evolved_groups_(
            std::move(preflight.restricted_evolved_groups)),
        level_contracts_(std::move(preflight.level_contracts)) {
    schedule_registry_ = load_certified_local_schedule_registry();
    dense_provider_ = std::make_shared<const DenseOutputProvider>(
        method_snapshot_.contract.dense);
    dense_registry_.register_provider(dense_provider_);
    stage_preparer_ = std::make_unique<TwoLevelStageSpatialPreparer>(
        *driver_,
        StageSpatialScheduleOwnership::certified_local_no_global_sync);

    const LevelAdvanceConfig method{
        method_snapshot_.contract.dense.method,
        method_snapshot_.contract.dense.tableau_fingerprint};
    HierarchyStepperConfig config{
        {method, method}, step_clock_t(0), initial_physical_time_,
        base_delta_time_, 2, 0, {0, 0}};
    stepper_ =
        std::make_unique<HierarchyStepper>(std::move(config), dense_registry_);

    // Everything above is preflight/read-only. This is the single transition
    // from the legacy finest physical delta to the static-v1 coarse base delta.
    const auto initial_metadata = full_sync_root_runtime_clock(
        0, initial_physical_time_, base_delta_time_);
    apply_persistent_metadata(initial_metadata);
  }

  HierarchyAdvanceResult advance_one_epoch() {
    return stepper_->advance_one_epoch(*this);
  }

  double weighted_cells_per_epoch() const {
    validate_immutable_contract();
    const auto &levels = driver_->patchdata.front().leveldata;
    return levels[0].fab->boxArray().d_numPts() +
           2.0 * levels[1].fab->boxArray().d_numPts();
  }

  std::unique_ptr<LevelStepSession>
  begin_level_step(const StepContext &context,
                   const bool require_dense) override {
    const auto begin = validate_begin_context(context, require_dense);
    return with_candidate_scope(context, [&]() {
      const auto &schema = *method_snapshot_.group_schema;
      auto frame = cycle_then_capture_static_v1_level_state(
          [&] {
            ScheduleInternal::cycle_active_level_grid_function_timelevels(
                &root_);
          },
          [&] {
            return copy_canonical_tl0_collective(*driver_, 0,
                                                 context.level);
          });
      ScratchStateTransactionFactoryMetadata metadata{
          hierarchy_epoch_,
          1,
          2,
          static_cast<int>(context.endpoint_accepted_step),
          step_clock_t(1),
          base_delta_time_,
          context.level == 0 ? 1 : 2,
          false,
          schema.ordered_group_pairs,
          schema.dependent_groups,
          {},
          {}};
      auto transaction = ScratchStateTransactionFactory::create_native(
          *driver_, *schedule_registry_, std::move(frame),
          std::move(metadata));

      TransactionLevelEvolution evolution =
          [this, context](ScratchStateTransaction &transaction) {
            with_candidate_scope(context, [&] {
              const auto *const active_context = current_step_context();
              if (active_context == nullptr ||
                  !same_step_context(*active_context, context) ||
                  current_scratch_state_transaction() != &transaction)
                throw std::logic_error(
                    "static-v1 evolution lost its active step transaction");
              CCTK_Traverse(&root_, "CCTK_PRESTEP");
              CCTK_Traverse(&root_, "CCTK_EVOL");
            });
          };
      ValidatedLevelCommit commit = [this, context, begin] {
        auto &target = validate_commit_target(context, begin);
        // All throwing work is complete. This is the callback's only mutation
        // and there is deliberately no operation after it.
        target.iteration = begin.endpoint_clock;
      };
      return std::make_unique<TransactionLevelStepSession>(
          std::move(transaction), require_dense, std::move(evolution),
          std::move(commit));
    });
  }

  void prepare_stage(const StepContext &context,
                     const StagePoint &stage_point,
                     const DenseInterval *parent_dense) override {
    const auto *const active_context = current_step_context();
    auto *const transaction = current_scratch_state_transaction();
    if (active_context == nullptr ||
        !same_step_context(*active_context, context) || transaction == nullptr)
      throw std::logic_error(
          "static-v1 stage preparation has no active transaction");
    stage_preparer_->prepare(*transaction, context, stage_point, parent_dense);
  }

  void begin_synchronization(const int coarse_level, const int fine_level,
                             const hierarchy_time_t time) override {
    if (synchronization_active_ || coarse_level != 0 || fine_level != 1)
      throw std::logic_error(
          "static-v1 received an invalid synchronization pair");
    validate_immutable_contract();
    const auto &levels = driver_->patchdata.front().leveldata;
    if (exact_clock(levels[0].iteration) != time ||
        exact_clock(levels[1].iteration) != time)
      throw std::logic_error(
          "static-v1 levels are not aligned before synchronization");
    synchronization_active_ = true;
  }

  void synchronize_levels(const int coarse_level, const int fine_level,
                          const hierarchy_time_t time) override {
    if (!synchronization_active_ || coarse_level != 0 || fine_level != 1)
      throw std::logic_error(
          "static-v1 synchronization scope is not active");
    if (time.den != 1 || time.num < 0)
      throw std::logic_error(
          "static-v1 synchronization clock is not a coarse endpoint");
    const double physical_time =
        initial_physical_time_ + static_cast<double>(time) * base_delta_time_;
    const auto runtime = full_sync_root_runtime_clock(
        static_cast<std::uint64_t>(time.num), physical_time,
        base_delta_time_);
    const ScheduleInternal::LevelTimeMetadata<CCTK_REAL> metadata{
        runtime.iteration, runtime.end_time, runtime.base_delta_time,
        runtime.timefac};
    ScheduleInternal::ScopedLevelTimeMetadata<cGH, PropagateRootMetadata>
        synchronized(root_, metadata, PropagateRootMetadata{});
    ScheduleInternal::restrict_and_sync_static_v1_evolved_groups(
        &root_, restricted_evolved_groups_);
  }

  void end_synchronization(const int, const int,
                           const hierarchy_time_t) noexcept override {
    synchronization_active_ = false;
  }

  void run_sync_observers(const hierarchy_time_t time,
                          const double physical_time,
                          const std::uint64_t completed_epoch,
                          const bool stop_requested) override {
    static_cast<void>(stop_requested);
    validate_immutable_contract();
    const auto &levels = driver_->patchdata.front().leveldata;
    const auto &stepper_clocks = stepper_->clocks();
    if (stepper_clocks.size() != 2)
      throw std::logic_error(
          "static-v1 stepper lost its two-level clock inventory");
    validate_static_v1_clock_envelope(StaticV1ClockEnvelope{
        2,
        2,
        false,
        false,
        time,
        {exact_clock(levels[0].iteration),
         exact_clock(levels[1].iteration)},
        {exact_clock(levels[0].delta_iteration),
         exact_clock(levels[1].delta_iteration)},
        completed_epoch,
        {stepper_clocks[0].accepted_steps,
         stepper_clocks[1].accepted_steps}});

    const auto metadata = full_sync_root_runtime_clock(
        completed_epoch, physical_time, base_delta_time_);
    apply_persistent_metadata(metadata);

    ScheduleInternal::ScopedActiveLevels<active_levels_t> all_levels(
        active_levels, active_levels_t(0, 2, 0, 1));
    for (const auto action : static_v1_sync_observer_order()) {
      switch (action) {
      case StaticV1SyncObserverAction::cycle_synchronized_globals:
        ScheduleInternal::cycle_synchronized_global_timelevels(&root_);
        break;
      case StaticV1SyncObserverAction::postrestrict:
        CCTK_Traverse(&root_, "CCTK_POSTRESTRICT");
        break;
      case StaticV1SyncObserverAction::poststep:
        CCTK_Traverse(&root_, "CCTK_POSTSTEP");
        break;
      case StaticV1SyncObserverAction::analysis:
        CCTK_Traverse(&root_, "CCTK_ANALYSIS");
        break;
      case StaticV1SyncObserverAction::output:
        CCTK_OutputGH(&root_);
        break;
      case StaticV1SyncObserverAction::checkpoint:
        CCTK_Traverse(&root_, "CCTK_CHECKPOINT");
        break;
      }
    }
  }

private:
  struct PropagateRootMetadata {
    void operator()(cGH &root) const noexcept {
      ScheduleInternal::propagate_root_runtime_metadata(root);
    }
  };

  template <class Function>
  auto with_candidate_scope(const StepContext &context, Function &&function)
      -> std::invoke_result_t<Function> {
    const auto runtime =
        candidate_runtime_clock(context, base_delta_time_);
    const ScheduleInternal::LevelTimeMetadata<CCTK_REAL> metadata{
        runtime.iteration, runtime.end_time, runtime.base_delta_time,
        runtime.timefac};
    ScheduleInternal::ScopedActiveLevels<active_levels_t> requested_level(
        active_levels,
        active_levels_t(context.level, context.level + 1, 0, 1));
    ScheduleInternal::ScopedLevelTimeMetadata<cGH, PropagateRootMetadata>
        candidate(root_, metadata, PropagateRootMetadata{});
    return std::forward<Function>(function)();
  }

  void apply_persistent_metadata(
      const RuntimeClockMetadata &runtime) noexcept {
    ScheduleInternal::apply_level_time_metadata(
        root_, ScheduleInternal::LevelTimeMetadata<CCTK_REAL>{
                   runtime.iteration, runtime.end_time,
                   runtime.base_delta_time, runtime.timefac});
    ScheduleInternal::propagate_root_runtime_metadata(root_);
  }

  void validate_immutable_contract() const {
    if (!ghext || ghext.get() != driver_ ||
        static_cast<std::int64_t>(CarpetX_GetEpoch()) != hierarchy_epoch_ ||
        driver_->patchdata.size() != 1 ||
        driver_->patchdata.front().leveldata.size() != 2)
      throw std::runtime_error(
          "static-v1 hierarchy identity changed after preflight");
    validate_registered_subcycling_method_contract(method_snapshot_);

    const auto &levels = driver_->patchdata.front().leveldata;
    for (std::size_t level_index = 0; level_index < levels.size();
         ++level_index) {
      const auto &level = levels[level_index];
      const auto &expected = level_contracts_[level_index];
      if (&level != expected.level_identity ||
          level.get_patch_cctkGH() != expected.patch_gh_identity ||
          level.fab.get() != expected.fab_identity ||
          level.delta_iteration != expected.delta_iteration ||
          level.is_subcycling_level != expected.is_subcycling_level ||
          expected.groups.empty())
        throw std::runtime_error(
            "static-v1 LevelData identity or clock contract changed");
      for (const auto &group : expected.groups) {
        if (group.group_index < 0 ||
            static_cast<std::size_t>(group.group_index) >=
                level.groupdata.size() ||
            level.groupdata[static_cast<std::size_t>(group.group_index)] ==
                nullptr)
          throw std::runtime_error(
              "static-v1 group identity disappeared");
        const auto &current = *level.groupdata[static_cast<std::size_t>(
            group.group_index)];
        std::vector<const void *> storage_identities;
        std::vector<const void *> validity_identities;
        storage_identities.reserve(current.mfab.size());
        validity_identities.reserve(current.valid.size());
        for (const auto &storage : current.mfab)
          storage_identities.push_back(storage.get());
        for (const auto &validity : current.valid)
          validity_identities.push_back(validity.data());
        const auto canonicalize = [](auto &identities) {
          std::sort(identities.begin(), identities.end(),
                    std::less<const void *>{});
        };
        canonicalize(storage_identities);
        canonicalize(validity_identities);
        if (&current != group.group_identity ||
            current.numvars != group.numvars ||
            current.do_restrict != group.do_restrict ||
            current.mfab.size() != group.timelevel_count ||
            current.valid.size() != group.validity_timelevel_count ||
            current.mfab.empty() || current.valid.empty() ||
            storage_identities != group.storage_identities ||
            validity_identities != group.validity_identities)
          throw std::runtime_error(
              "static-v1 group storage contract changed");
      }
    }
  }

  BeginLevelSnapshot
  validate_begin_context(const StepContext &context,
                         const bool require_dense) const {
    validate_immutable_contract();
    if (context.level < 0 || context.level > 1 ||
        context.method != method_snapshot_.contract.dense.method ||
        require_dense != (context.level == 0) ||
        context.require_dense_output != require_dense)
      throw std::invalid_argument(
          "static-v1 level-step context differs from the frozen method");
    static_cast<void>(candidate_runtime_clock(context, base_delta_time_));

    const auto &levels = driver_->patchdata.front().leveldata;
    BeginLevelSnapshot begin{{levels[0].iteration, levels[1].iteration},
                             level_clock(context.end_clock)};
    if (exact_clock(begin.clocks[static_cast<std::size_t>(context.level)]) !=
        context.begin_clock)
      throw std::logic_error(
          "static-v1 active LevelData clock differs from the step begin");
    return begin;
  }

  LevelData &validate_commit_target(const StepContext &context,
                                    const BeginLevelSnapshot &begin) {
    validate_immutable_contract();
    auto &levels = driver_->patchdata.front().leveldata;
    if (levels[0].iteration != begin.clocks[0] ||
        levels[1].iteration != begin.clocks[1] ||
        exact_clock(levels[static_cast<std::size_t>(context.level)]
                        .iteration) != context.begin_clock ||
        exact_clock(begin.endpoint_clock) != context.end_clock)
      throw std::logic_error(
          "static-v1 LevelData clocks changed before accepted commit");
    return levels[static_cast<std::size_t>(context.level)];
  }

  cGH &root_;
  GHExt *const driver_;
  const SubcyclingMethodContractSnapshot method_snapshot_;
  const std::int64_t hierarchy_epoch_;
  const double initial_physical_time_;
  const double base_delta_time_;
  const std::vector<int> restricted_evolved_groups_;
  const std::array<LevelContractIdentity, 2> level_contracts_;
  DenseOutputRegistry dense_registry_;
  std::shared_ptr<const DenseOutputProvider> dense_provider_;
  std::unique_ptr<CertifiedLocalScheduleRegistry> schedule_registry_;
  std::unique_ptr<TwoLevelStageSpatialPreparer> stage_preparer_;
  std::unique_ptr<HierarchyStepper> stepper_;
  bool synchronization_active_{false};
};

} // namespace

int EvolveStaticV1(tFleshConfig *const config) {
  try {
    if (config == nullptr || config->nGHs == 0 || config->GH == nullptr ||
        config->GH[0] == nullptr)
      throw std::invalid_argument(
          "static-v1 requires a configured root grid hierarchy");
    cGH &root = *config->GH[0];
    auto preflight =
        preflight_static_v1(root, config->recovered != 0);
    StaticV1EvolutionAdapter adapter(root, std::move(preflight));

#pragma omp critical
    CCTK_VINFO("Starting static-v1 factor-two subcycling evolution...");

    const double weighted_cells = adapter.weighted_cells_per_epoch();
    double total_seconds = 0.0;
    std::uint64_t coarse_epochs = 0;
    while (!ScheduleInternal::evolution_is_done_at_sync(&root)) {
      const auto started = std::chrono::steady_clock::now();
      const auto result = adapter.advance_one_epoch();
      const auto finished = std::chrono::steady_clock::now();
      const double elapsed =
          std::chrono::duration<double>(finished - started).count();
      total_seconds += elapsed;
      ++coarse_epochs;
      const double rate = elapsed > 0.0 ? weighted_cells / elapsed : 0.0;
#pragma omp critical
      CCTK_VINFO("Static-v1 epoch: %llu   simulation time: %.17g   "
                 "weighted cells: %.0f   weighted cell updates/s: %.17g",
                 static_cast<unsigned long long>(result.epoch),
                 result.synchronized_physical_time, weighted_cells, rate);
    }

    if (coarse_epochs > 0 && total_seconds > 0.0) {
#pragma omp critical
      CCTK_VINFO("Static-v1 performance: coarse epochs=%llu   "
                 "evolution seconds=%.17g   average weighted cell updates/s="
                 "%.17g",
                 static_cast<unsigned long long>(coarse_epochs),
                 total_seconds,
                 coarse_epochs * weighted_cells / total_seconds);
    }
    return 0;
  } catch (const std::exception &error) {
    CCTK_VERROR("Static-v1 subcycling evolution failed closed: %s",
                error.what());
  } catch (...) {
    CCTK_ERROR(
        "Static-v1 subcycling evolution failed closed with an unknown error");
  }
  return -1;
}

} // namespace CarpetX
