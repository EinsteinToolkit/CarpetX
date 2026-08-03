#include "subcycling_stage_spatial_preparer.hxx"
#include "subcycling_stage_spatial_preparer_internal.hxx"

#include "driver.hxx"
#include "fillpatch.hxx"
#include "subcycling_dense_mfab_state.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"
#include "task_manager.hxx"

#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CarpetX {

// The driver-wide synchronization primitive is defined in schedule.cxx.
void synchronize();

namespace {

bool same_layout(const amrex::MultiFab &left,
                 const amrex::MultiFab &right) {
  return left.isDefined() && right.isDefined() &&
         left.boxArray() == right.boxArray() &&
         left.DistributionMap() == right.DistributionMap() &&
         left.nComp() == right.nComp() &&
         left.nGrowVect() == right.nGrowVect();
}

bool ratio_is_two(const amrex::IntVect &ratio) {
  for (int direction = 0; direction < AMREX_SPACEDIM; ++direction)
    if (ratio[direction] != 2)
      return false;
  return true;
}

} // namespace

class TwoLevelStageSpatialNativeBackend final
    : public detail::StageSpatialPreparerBackend {
  using GroupData = GHExt::PatchData::LevelData::GroupData;
  using StageAccess = ScratchStateTransactionFactory::StageSpatialAccess;

  struct Binding {
    int group_index;
    std::size_t frame_entry;
    GroupData *groupdata;
    amrex::MultiFab *target;
    std::vector<why_valid_t> *live_validity;
    std::vector<ScratchValidity> *scratch_validity;
  };

  struct Pending {
    StageSpatialTarget target;
    int level;
    step_clock_t stage_clock;
    std::vector<int> evolved_groups;
  };

public:
  TwoLevelStageSpatialNativeBackend(
      GHExt &ghext, const StageSpatialScheduleOwnership ownership)
      : ghext_(ghext), ownership_(ownership) {}

  StageSpatialMetadata
  inspect(ScratchStateTransaction *transaction) override {
    auto &native_transaction = require_transaction(transaction);
    const auto access = ScratchStateTransactionFactory::stage_spatial_access(
        native_transaction);
    const int patch_count = static_cast<int>(ghext_.patchdata.size());
    const int level_count =
        patch_count == 1
            ? static_cast<int>(ghext_.patchdata.front().leveldata.size())
            : 0;
    int spatial_ratio = 2;
    if (patch_count == 1 && level_count == 2) {
      const auto &patchdata = ghext_.patchdata.front();
      spatial_ratio =
          patchdata.amrcore != nullptr &&
                  ratio_is_two(patchdata.amrcore->refRatio(0))
              ? 2
              : -1;
    }
    return StageSpatialMetadata{
        patch_count,
        level_count,
        spatial_ratio,
        access.level,
        access.hierarchy_epoch,
        static_cast<std::int64_t>(CarpetX_GetEpoch()),
        ownership_ == StageSpatialScheduleOwnership::conflicting_global_sync,
        native_transaction.group_pairs().size()};
  }

  void prepare_level_zero(ScratchStateTransaction *transaction,
                          const StageSpatialTarget target,
                          const step_clock_t stage_clock) override {
    auto &native_transaction = require_transaction(transaction);
    require_no_pending();
    auto access = ScratchStateTransactionFactory::stage_spatial_access(
        native_transaction);
    if (access.level != 0)
      throw std::invalid_argument(
          "level-zero spatial preparation received another level");
    auto bindings = collect_bindings(access, native_transaction, target);
    invalidate_spatial(bindings);
    pending_ = Pending{target, 0, stage_clock, group_indices(bindings)};

    task_manager tasks2;
    const auto &patchdata = ghext_.patchdata.at(0);
    if (patchdata.amrcore == nullptr)
      throw std::runtime_error("stage spatial AMR core disappeared");
    const auto &geometry = patchdata.amrcore->Geom(0);
    for (auto &binding : bindings)
      FillPatch_Sync(tasks2, *binding.groupdata, *binding.target, geometry);
    synchronize();
    tasks2.run_tasks_serially();
    synchronize();
  }

  void prepare_level_one(ScratchStateTransaction *transaction,
                         const StageSpatialTarget target,
                         const step_clock_t stage_clock,
                         const step_clock_t parent_theta,
                         const DenseInterval &parent_dense) override {
    auto &native_transaction = require_transaction(transaction);
    require_no_pending();
    auto access = ScratchStateTransactionFactory::stage_spatial_access(
        native_transaction);
    if (access.level != 1)
      throw std::invalid_argument(
          "level-one spatial preparation received another level");
    auto bindings = collect_bindings(access, native_transaction, target);

    // The parent sample is fully materialized and compatibility-checked before
    // any fine target validity is invalidated.
    auto parent_state = evaluate_parent(access, bindings, parent_theta,
                                        parent_dense);
    validate_level_one_fill(bindings, *parent_state);
    invalidate_spatial(bindings);
    pending_ = Pending{target, 1, stage_clock, group_indices(bindings)};

    auto &patchdata = ghext_.patchdata.at(0);
    const auto &fine_geometry = patchdata.amrcore->Geom(1);
    const auto &coarse_geometry = patchdata.amrcore->Geom(0);
    task_manager tasks2;
    task_manager tasks3;
    for (std::size_t index = 0; index < bindings.size(); ++index) {
      auto &binding = bindings[index];
      auto &coarse_group = require_group(0, binding.group_index);
      FillPatch_ProlongateGhosts(
          tasks2, tasks3, *binding.groupdata, coarse_group, *binding.target,
          parent_state->multifab(index), fine_geometry, coarse_geometry,
          binding.groupdata->interpolator, binding.groupdata->bcrecs);
    }
    synchronize();
    tasks2.run_tasks_serially();
    synchronize();
    tasks3.run_tasks_serially();
    synchronize();
  }

  void promote(ScratchStateTransaction *transaction,
               const StageSpatialTarget target) override {
    auto &native_transaction = require_transaction(transaction);
    if (!pending_.has_value() || pending_->target != target)
      throw std::logic_error(
          "stage spatial validity promotion lacks its preparation");
    auto access = ScratchStateTransactionFactory::stage_spatial_access(
        native_transaction);
    if (access.level != pending_->level)
      throw std::runtime_error(
          "stage spatial level changed before validity promotion");
    auto bindings = collect_bindings(access, native_transaction, target);
    if (group_indices(bindings) != pending_->evolved_groups)
      throw std::runtime_error(
          "stage spatial group schema changed before validity promotion");

    promote_spatial(bindings);
    pending_.reset();
  }

  void fault(ScratchStateTransaction *transaction) noexcept override {
    pending_.reset();
    if (transaction != nullptr)
      ScratchStateTransactionFactory::fault_stage_spatial_preparation(
          *transaction);
  }

private:
  static ScratchStateTransaction &
  require_transaction(ScratchStateTransaction *transaction) {
    if (transaction == nullptr)
      throw std::invalid_argument(
          "native stage spatial preparation requires a transaction");
    return *transaction;
  }

  void require_no_pending() const {
    if (pending_.has_value())
      throw std::logic_error(
          "prior stage spatial validity has not been promoted");
  }

  GroupData &require_group(const int level, const int group_index) {
    if (ghext_.patchdata.size() != 1 || level < 0 || level > 1 ||
        static_cast<std::size_t>(level) >=
            ghext_.patchdata.front().leveldata.size())
      throw std::runtime_error("stage spatial hierarchy layout changed");
    auto &leveldata = ghext_.patchdata.front().leveldata.at(
        static_cast<std::size_t>(level));
    if (group_index < 0 ||
        static_cast<std::size_t>(group_index) >=
            leveldata.groupdata.size() ||
        leveldata.groupdata[static_cast<std::size_t>(group_index)] == nullptr)
      throw std::runtime_error("stage spatial group disappeared");
    return *leveldata.groupdata[static_cast<std::size_t>(group_index)];
  }

  static std::size_t find_frame_entry(const ScratchLevelFrame &frame,
                                      const int group_index) {
    for (std::size_t entry = 0; entry < frame.entry_count(); ++entry)
      if (frame.key(entry).group_index == group_index)
        return entry;
    throw std::runtime_error("stage spatial scratch group disappeared");
  }

  static const ScratchLiveEntrySnapshot &
  find_live_entry(const StageAccess &access, const int group_index) {
    const auto found = std::find_if(
        access.live_entries.begin(), access.live_entries.end(),
        [group_index](const ScratchLiveEntrySnapshot &entry) {
          return entry.key.group_index == group_index;
        });
    if (found == access.live_entries.end())
      throw std::runtime_error("stage spatial live group disappeared");
    return *found;
  }

  std::vector<Binding>
  collect_bindings(StageAccess &access,
                   const ScratchStateTransaction &transaction,
                   const StageSpatialTarget target) {
    if (access.working_frame == nullptr ||
        access.working_frame->hierarchy_epoch() != access.hierarchy_epoch ||
        access.working_frame->level() != access.level)
      throw std::runtime_error("stage spatial scratch frame changed");
    std::vector<Binding> bindings;
    bindings.reserve(transaction.group_pairs().size());
    for (const auto &pair : transaction.group_pairs()) {
      const int group_index = pair.evolved_group;
      auto &groupdata = require_group(access.level, group_index);
      if (groupdata.groupindex != group_index || groupdata.mfab.empty() ||
          groupdata.valid.empty() || groupdata.mfab.at(0) == nullptr)
        throw std::runtime_error(
            "stage spatial live TL0 schema changed");
      auto *const live = groupdata.mfab.at(0).get();
      auto &live_validity = groupdata.valid.at(0);
      const auto &snapshot = find_live_entry(access, group_index);
      const auto frame_entry =
          find_frame_entry(*access.working_frame, group_index);
      const auto &frame_key = access.working_frame->key(frame_entry);
      if (snapshot.key.hierarchy_epoch != access.hierarchy_epoch ||
          snapshot.key.patch != 0 || snapshot.key.level != access.level ||
          frame_key.hierarchy_epoch != access.hierarchy_epoch ||
          frame_key.patch != 0 || frame_key.level != access.level ||
          snapshot.group_identity != &groupdata ||
          snapshot.tl0_identity != &groupdata.mfab.at(0) ||
          snapshot.storage_identity != live || snapshot.multifab != live ||
          snapshot.validity_identity != &live_validity ||
          !snapshot.grid_function_real || live == nullptr ||
          !live->isDefined())
        throw std::runtime_error(
            "stage spatial live identity or key changed");

      amrex::MultiFab *target_mfab = nullptr;
      std::vector<why_valid_t> *target_live_validity = nullptr;
      std::vector<ScratchValidity> *target_scratch_validity = nullptr;
      if (target == StageSpatialTarget::primary_live_tl0) {
        target_mfab = live;
        target_live_validity = &live_validity;
      } else {
        target_mfab =
            &access.working_frame->mutable_multifab(frame_entry);
        target_scratch_validity =
            &access.working_frame->mutable_validity(frame_entry);
      }
      if (target_mfab == nullptr || !same_layout(*target_mfab, *live) ||
          target_mfab->nComp() != groupdata.numvars ||
          live_validity.size() !=
              static_cast<std::size_t>(groupdata.numvars) ||
          snapshot.validity.size() != live_validity.size() ||
          (target_scratch_validity != nullptr &&
           target_scratch_validity->size() != live_validity.size()))
        throw std::runtime_error(
            "stage spatial target layout or schema changed");

      for (std::size_t component = 0; component < live_validity.size();
           ++component) {
        const bool interior =
            target_live_validity != nullptr
                ? target_live_validity->at(component).get().valid_int
                : target_scratch_validity->at(component).interior;
        if (!interior)
          throw std::runtime_error(
              "stage spatial target interior is invalid");
      }
      bindings.push_back(Binding{group_index, frame_entry, &groupdata,
                                 target_mfab, target_live_validity,
                                 target_scratch_validity});
    }
    return bindings;
  }

  std::unique_ptr<OwnedMultiFabDenseState>
  evaluate_parent(const StageAccess &access,
                  const std::vector<Binding> &bindings,
                  const step_clock_t parent_theta,
                  const DenseInterval &parent_dense) {
    std::vector<DenseMFabView> views;
    views.reserve(bindings.size());
    for (const auto &binding : bindings) {
      auto &coarse_group = require_group(0, binding.group_index);
      if (coarse_group.mfab.empty() || coarse_group.mfab.at(0) == nullptr ||
          !coarse_group.mfab.at(0)->isDefined())
        throw std::runtime_error(
            "stage spatial parent live group disappeared");
      views.push_back(
          DenseMFabView{{access.hierarchy_epoch, 0, 0,
                         binding.group_index},
                        coarse_group.mfab.at(0).get()});
    }
    auto state = OwnedMultiFabDenseState::allocate_like(views);
    parent_dense.evaluate(parent_theta, *state);
    return state;
  }

  void validate_level_one_fill(
      const std::vector<Binding> &bindings,
      const OwnedMultiFabDenseState &parent_state) {
    if (ghext_.patchdata.size() != 1 ||
        ghext_.patchdata.front().leveldata.size() != 2 ||
        ghext_.patchdata.front().amrcore == nullptr ||
        !ratio_is_two(ghext_.patchdata.front().amrcore->refRatio(0)) ||
        parent_state.entry_count() != bindings.size())
      throw std::runtime_error(
          "stage spatial level-one hierarchy changed");
    for (std::size_t index = 0; index < bindings.size(); ++index) {
      const auto &binding = bindings[index];
      auto &coarse_group = require_group(0, binding.group_index);
      const auto &key = parent_state.key(index);
      if (key.patch != 0 || key.level != 0 ||
          key.group_index != binding.group_index ||
          coarse_group.numvars != binding.groupdata->numvars ||
          binding.groupdata->interpolator == nullptr ||
          !same_layout(parent_state.multifab(index),
                       *coarse_group.mfab.at(0)))
        throw std::runtime_error(
            "stage spatial parent sample or prolongation schema changed");
    }
  }

  static std::vector<int>
  group_indices(const std::vector<Binding> &bindings) {
    std::vector<int> result;
    result.reserve(bindings.size());
    for (const auto &binding : bindings)
      result.push_back(binding.group_index);
    return result;
  }

  static void invalidate_spatial(std::vector<Binding> &bindings) {
    for (auto &binding : bindings) {
      if (binding.live_validity != nullptr) {
        for (auto &validity : *binding.live_validity) {
          validity.set_outer(false, [] {
            return "TwoLevelStageSpatialPreparer fill pending";
          });
          validity.set_ghosts(false, [] {
            return "TwoLevelStageSpatialPreparer fill pending";
          });
        }
      } else {
        for (auto &validity : *binding.scratch_validity) {
          validity.outer = false;
          validity.ghosts = false;
        }
      }
    }
  }

  static void promote_spatial(std::vector<Binding> &bindings) {
    for (auto &binding : bindings) {
      if (binding.live_validity != nullptr) {
        for (auto &validity : *binding.live_validity) {
          validity.set_outer(true, [] {
            return "TwoLevelStageSpatialPreparer fill completed";
          });
          validity.set_ghosts(true, [] {
            return "TwoLevelStageSpatialPreparer fill completed";
          });
        }
      } else {
        for (auto &validity : *binding.scratch_validity) {
          validity.outer = true;
          validity.ghosts = true;
        }
      }
    }
  }

  GHExt &ghext_;
  StageSpatialScheduleOwnership ownership_;
  std::optional<Pending> pending_;
};

TwoLevelStageSpatialPreparer::TwoLevelStageSpatialPreparer(
    GHExt &ghext, const StageSpatialScheduleOwnership schedule_ownership)
    : TwoLevelStageSpatialPreparer(
          std::make_unique<TwoLevelStageSpatialNativeBackend>(
              ghext, schedule_ownership)) {}

StageSpatialPreparationReceipt TwoLevelStageSpatialPreparer::prepare(
    ScratchStateTransaction &transaction, const StepContext &context,
    const StagePoint &stage_point, const DenseInterval *parent_dense) {
  return prepare_impl(&transaction, context, stage_point, parent_dense);
}

} // namespace CarpetX
