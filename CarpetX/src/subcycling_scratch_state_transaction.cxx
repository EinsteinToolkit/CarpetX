#include "subcycling_scratch_state_transaction.hxx"

#include "subcycling_dense_control_builder.hxx"
#include "subcycling_dense_mfab_state.hxx"
#include "subcycling_dense_output.hxx"
#include "subcycling_dense_stencil.hxx"
#ifdef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
#define CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
#endif
#include "subcycling_scratch_local_gh.hxx"
#include "subcycling_scratch_schedule_executor.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"

#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
#include "driver.hxx"
#include <cctk_Parameters.h>
extern "C" CCTK_INT CarpetX_GetEpoch(void);
#endif

#include <AMReX_MultiFab.H>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>

namespace CarpetX {
namespace {

std::atomic<std::uint64_t> next_transaction_owner{1};

bool same_key(const ScratchGFKey &left, const ScratchGFKey &right) noexcept {
  return std::tie(left.hierarchy_epoch, left.level, left.patch,
                  left.group_index) ==
         std::tie(right.hierarchy_epoch, right.level, right.patch,
                  right.group_index);
}

bool same_layout(const amrex::MultiFab &left,
                 const amrex::MultiFab &right) noexcept {
  return left.boxArray() == right.boxArray() &&
         left.DistributionMap() == right.DistributionMap() &&
         left.nComp() == right.nComp() &&
         left.nGrowVect() == right.nGrowVect();
}

bool finite_positive(const double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool close_time(const double left, const double right,
                const double begin, const double end) noexcept {
  const auto scale =
      std::max({1.0, std::abs(left), std::abs(right), std::abs(begin),
                std::abs(end)});
  const auto tolerance =
      16.0 * std::numeric_limits<double>::epsilon() * scale;
  return std::abs(left - right) <= tolerance;
}

std::size_t find_group(const ScratchLevelFrame &frame, const int group) {
  for (std::size_t entry = 0; entry < frame.entry_count(); ++entry)
    if (frame.key(entry).group_index == group)
      return entry;
  throw std::invalid_argument("scratch transaction group is absent");
}

DenseSampleConstraint sample_constraint(const ScratchDenseSampleRef &sample) {
  return {sample.theta,
          sample.kind == ScratchDenseSampleKind::value
              ? DenseSampleKind::value
              : DenseSampleKind::scaled_derivative};
}

bool same_constraint(const DenseSampleConstraint &left,
                     const DenseSampleConstraint &right) noexcept {
  return left.theta == right.theta && left.kind == right.kind;
}

} // namespace

ScratchStateToken::ScratchStateToken() noexcept = default;

ScratchStateToken::ScratchStateToken(const std::uint64_t owner,
                                     const std::uint64_t state,
                                     const std::int64_t epoch,
                                     const std::uint64_t schema,
                                     const ScratchStateKind kind) noexcept
    : owner_(owner), state_(state), epoch_(epoch), schema_(schema),
      kind_(kind) {}

ScratchStateToken::ScratchStateToken(ScratchStateToken &&other) noexcept
    : owner_(other.owner_), state_(other.state_), epoch_(other.epoch_),
      schema_(other.schema_), kind_(other.kind_) {
  other.reset();
}

ScratchStateToken &
ScratchStateToken::operator=(ScratchStateToken &&other) noexcept {
  if (this == &other)
    return *this;
  owner_ = other.owner_;
  state_ = other.state_;
  epoch_ = other.epoch_;
  schema_ = other.schema_;
  kind_ = other.kind_;
  other.reset();
  return *this;
}

bool ScratchStateToken::valid() const noexcept {
  return owner_ != 0 && state_ != 0 && epoch_ >= 0 && schema_ != 0;
}

void ScratchStateToken::reset() noexcept {
  owner_ = 0;
  state_ = 0;
  epoch_ = -1;
  schema_ = 0;
  kind_ = ScratchStateKind::evolved;
}

struct ScratchStateTransaction::Storage {
  struct PairIndex {
    ScratchGroupPair groups;
    std::size_t evolved_entry;
    std::size_t rhs_entry;
  };

  struct StateRecord {
    ScratchStateKind kind;
    ScratchLevelFrame frame;
  };

  struct ExpectedLiveEntry {
    ScratchLiveEntrySnapshot snapshot;
  };

  std::uint64_t owner{next_transaction_owner.fetch_add(1)};
  std::uint64_t schema{1};
  std::int64_t hierarchy_epoch{-1};
  std::vector<ScratchGroupPair> group_pairs;
  std::vector<int> dependent_groups;
  std::vector<PairIndex> pair_indices;
  std::vector<std::size_t> invalidated_entries;
  int level{-1};
  int level_iteration{-1};
  step_clock_t base_delta_clock{};
  double base_delta_time{0.0};
  int time_refinement_factor{0};
  std::function<std::int64_t()> epoch_reader;
  std::vector<ScratchLiveEntryReader> live_readers;
  std::vector<ExpectedLiveEntry> expected_live;
  std::vector<std::unique_ptr<StateRecord>> states;
  std::size_t rhs_evaluations{0};
  bool faulted{false};
  bool discarded{false};
  bool dense_committed_once{false};
  std::shared_ptr<const DenseInterval> committed_dense;

  // Declaration order is intentional: executor is destroyed before the
  // binding whose lease it owns.
#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
  std::unique_ptr<ScratchLocalLevelBinding> binding;
  std::unique_ptr<CertifiedScratchStageExecutor> executor;
#else
  std::optional<ScratchLevelFrame> test_frame;
  ScratchStateTransactionFactory::TestScheduleExecutor test_executor;
#endif
  ScratchLevelFrame *working_frame{nullptr};

  void release_execution_storage() noexcept {
#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
    executor.reset();
    binding.reset();
#else
    test_executor = {};
    test_frame.reset();
#endif
    working_frame = nullptr;
  }

  void fault_and_discard() noexcept {
    faulted = true;
    discarded = true;
    states.clear();
    committed_dense.reset();
    release_execution_storage();
  }
};

class ScratchStateTransactionCore final {
public:
  static void require_available(ScratchStateTransaction::Storage &);
  static ScratchStateTransaction::Storage::StateRecord &
  require_state(ScratchStateTransaction::Storage &,
                const ScratchStateToken &);
  static const ScratchStateTransaction::Storage::StateRecord &
  require_state(const ScratchStateTransaction::Storage &,
                const ScratchStateToken &);
  static ScratchStateToken publish_state(ScratchStateTransaction::Storage &,
                                         ScratchStateKind,
                                         ScratchLevelFrame);
  static std::vector<ScratchLiveEntrySnapshot>
  read_and_validate_live(ScratchStateTransaction::Storage &);
  static ScratchLevelFrame capture_live_frame(
      ScratchStateTransaction::Storage &, ScratchStateKind);
  static ScratchLevelFrame capture_working_frame(
      ScratchStateTransaction::Storage &, ScratchStateKind);
  static ScratchStageCoordinates stage_coordinates(
      ScratchStateTransaction::Storage &, const StepContext &,
      const StagePoint &);
  static void validate_factory_metadata(
      const ScratchLevelFrame &,
      const ScratchStateTransactionFactoryMetadata &,
      std::vector<ScratchStateTransaction::Storage::PairIndex> &,
      std::vector<std::size_t> &,
      std::vector<ScratchStateTransaction::Storage::ExpectedLiveEntry> &);
};

namespace {

} // namespace

void ScratchStateTransactionCore::require_available(
    ScratchStateTransaction::Storage &storage) {
  if (storage.discarded || storage.working_frame == nullptr)
    throw std::logic_error("scratch transaction is discarded");
  if (storage.faulted)
    throw std::logic_error("scratch transaction is faulted");
  try {
    if (!storage.epoch_reader ||
        storage.epoch_reader() != storage.hierarchy_epoch) {
      storage.fault_and_discard();
      throw std::runtime_error("scratch transaction hierarchy epoch changed");
    }
  } catch (...) {
    storage.fault_and_discard();
    throw;
  }
}

ScratchStateTransaction::Storage::StateRecord &
ScratchStateTransactionCore::require_state(
    ScratchStateTransaction::Storage &storage,
    const ScratchStateToken &token) {
  require_available(storage);
  if (!token.valid() || token.owner_ != storage.owner ||
      token.epoch_ != storage.hierarchy_epoch ||
      token.schema_ != storage.schema || token.state_ > storage.states.size())
    throw std::invalid_argument("scratch state token does not belong here");
  auto *const record = storage.states[token.state_ - 1].get();
  if (record == nullptr || record->kind != token.kind_)
    throw std::invalid_argument("scratch state token metadata is stale");
  return *record;
}

const ScratchStateTransaction::Storage::StateRecord &
ScratchStateTransactionCore::require_state(
    const ScratchStateTransaction::Storage &storage,
    const ScratchStateToken &token) {
  return ScratchStateTransactionCore::require_state(
      const_cast<ScratchStateTransaction::Storage &>(storage), token);
}

ScratchStateToken ScratchStateTransactionCore::publish_state(
    ScratchStateTransaction::Storage &storage, const ScratchStateKind kind,
    ScratchLevelFrame frame) {
  auto record = std::make_unique<ScratchStateTransaction::Storage::StateRecord>(
      ScratchStateTransaction::Storage::StateRecord{kind, std::move(frame)});
  storage.states.push_back(std::move(record));
  return ScratchStateToken(storage.owner, storage.states.size(),
                           storage.hierarchy_epoch, storage.schema, kind);
}

std::vector<ScratchLiveEntrySnapshot>
ScratchStateTransactionCore::read_and_validate_live(
    ScratchStateTransaction::Storage &storage) {
  require_available(storage);
  if (!storage.epoch_reader ||
      storage.epoch_reader() != storage.hierarchy_epoch)
    throw std::runtime_error("scratch live hierarchy epoch changed");
  if (storage.live_readers.size() != storage.expected_live.size())
    throw std::runtime_error("scratch live reader schema changed");

  std::vector<ScratchLiveEntrySnapshot> snapshots;
  snapshots.reserve(storage.live_readers.size());
  for (std::size_t entry = 0; entry < storage.live_readers.size(); ++entry) {
    if (!storage.live_readers[entry])
      throw std::runtime_error("scratch live reader is unavailable");
    const auto observed = storage.live_readers[entry]();
    const auto &expected = storage.expected_live[entry].snapshot;
    if (!same_key(observed.key, expected.key) ||
        observed.level_identity != expected.level_identity ||
        observed.group_identity != expected.group_identity ||
        observed.tl0_identity != expected.tl0_identity ||
        observed.storage_identity != expected.storage_identity ||
        observed.validity_identity != expected.validity_identity ||
        observed.multifab != expected.multifab ||
        observed.grid_function_real != expected.grid_function_real ||
        !observed.grid_function_real || !observed.restore ||
        observed.multifab == nullptr ||
        !observed.multifab->isDefined() ||
        !same_layout(*observed.multifab,
                     storage.working_frame->multifab(entry)) ||
        observed.validity.size() !=
            storage.working_frame->validity(entry).size())
      throw std::runtime_error("scratch live identity or schema changed");
    snapshots.push_back(observed);
  }
  return snapshots;
}

ScratchLevelFrame ScratchStateTransactionCore::capture_live_frame(
    ScratchStateTransaction::Storage &storage,
    const ScratchStateKind kind) {
  const auto snapshots = read_and_validate_live(storage);
  UncertifiedScratchLevelManifest manifest{storage.hierarchy_epoch,
                                           storage.level, {}};
  std::vector<ScratchGFView> views;
  manifest.entries.reserve(snapshots.size());
  views.reserve(snapshots.size());
  for (const auto &snapshot : snapshots) {
    manifest.entries.push_back({snapshot.key, snapshot.storage_identity});
    views.push_back({snapshot.key, 0, snapshot.storage_identity,
                     snapshot.multifab, &snapshot.validity});
  }
  auto frame = ScratchLevelFrame::copy_tl0(manifest, views);
  if (kind == ScratchStateKind::raw_rhs) {
    for (const auto &pair : storage.pair_indices) {
      frame.copy_interior_for_transaction(pair.evolved_entry, frame,
                                          pair.rhs_entry);
      frame.copy_interior_validity_for_transaction(
          pair.evolved_entry, frame, pair.rhs_entry, true);
    }
  }
  return frame;
}

ScratchLevelFrame ScratchStateTransactionCore::capture_working_frame(
    ScratchStateTransaction::Storage &storage,
    const ScratchStateKind kind) {
  require_available(storage);
  auto frame = storage.working_frame->clone_for_transaction();
  if (kind == ScratchStateKind::raw_rhs) {
    for (const auto &pair : storage.pair_indices) {
      frame.copy_interior_for_transaction(pair.evolved_entry, frame,
                                          pair.rhs_entry);
      frame.copy_interior_validity_for_transaction(
          pair.evolved_entry, frame, pair.rhs_entry, true);
    }
  }
  return frame;
}

ScratchStageCoordinates ScratchStateTransactionCore::stage_coordinates(
    ScratchStateTransaction::Storage &storage,
    const StepContext &context, const StagePoint &stage_point) {
  require_available(storage);
  if (current_step_context() != &context ||
      current_scratch_state_transaction() == nullptr)
    throw std::logic_error(
        "scratch schedule requires its active StepContext scope");
  if (context.level != storage.level)
    throw std::invalid_argument(
        "scratch transaction StepContext level differs");
  if (stage_point.stage_index <= 0 || stage_point.stage_count <= 0 ||
      stage_point.stage_index > stage_point.stage_count)
    throw std::invalid_argument("scratch stage index/count are invalid");
  switch (stage_point.kind) {
  case StageKind::primary:
  case StageKind::fractional:
    break;
  case StageKind::endpoint_probe:
    if (stage_point.stage_index != 1 || stage_point.stage_count != 1)
      throw std::invalid_argument(
          "scratch endpoint probe must use the singleton stage index");
    break;
  default:
    throw std::invalid_argument("scratch stage kind is invalid");
  }
  if (stage_point.stage_fraction < step_clock_t(0) ||
      stage_point.stage_fraction > step_clock_t(1))
    throw std::out_of_range(
        "scratch exact stage fraction is outside [0,1]");
  if (!std::isfinite(stage_point.stage_time))
    throw std::invalid_argument("scratch stage time must be finite");
  const auto level_delta_clock =
      storage.base_delta_clock / storage.time_refinement_factor;
  if (context.end_clock - context.begin_clock != level_delta_clock)
    throw std::invalid_argument("scratch stage exact clock differs");
  const auto level_delta_time =
      storage.base_delta_time /
      static_cast<double>(storage.time_refinement_factor);
  if (!close_time(context.end_time - context.begin_time, level_delta_time,
                  context.begin_time, context.end_time))
    throw std::invalid_argument("scratch stage floating clock differs");
  const auto scale =
      std::max({1.0, std::abs(context.begin_time),
                std::abs(context.end_time)});
  const auto tolerance =
      16.0 * std::numeric_limits<double>::epsilon() * scale;
  if (stage_point.stage_time < context.begin_time - tolerance ||
      stage_point.stage_time > context.end_time + tolerance)
    throw std::out_of_range("scratch stage time is outside the step");
  const auto exact_stage_clock =
      context.begin_clock +
      stage_point.stage_fraction * (context.end_clock - context.begin_clock);
  if (exact_stage_clock < context.begin_clock ||
      exact_stage_clock > context.end_clock)
    throw std::out_of_range("scratch exact stage clock is outside the step");
  const auto exact_fraction = static_cast<double>(stage_point.stage_fraction);
  const auto expected_stage_time =
      context.begin_time +
      exact_fraction * (context.end_time - context.begin_time);
  if (!std::isfinite(expected_stage_time) ||
      std::abs(stage_point.stage_time - expected_stage_time) > tolerance)
    throw std::invalid_argument(
        "scratch stage time differs from its exact stage fraction");
  return {storage.level_iteration, stage_point.stage_time,
          storage.base_delta_clock,
          storage.base_delta_time, storage.time_refinement_factor};
}

void ScratchStateTransactionCore::validate_factory_metadata(
    const ScratchLevelFrame &frame,
    const ScratchStateTransactionFactoryMetadata &metadata,
    std::vector<ScratchStateTransaction::Storage::PairIndex> &pair_indices,
    std::vector<std::size_t> &invalidated_entries,
    std::vector<ScratchStateTransaction::Storage::ExpectedLiveEntry>
        &expected_live) {
  if (metadata.hierarchy_epoch < 0 ||
      frame.hierarchy_epoch() != metadata.hierarchy_epoch)
    throw std::invalid_argument("scratch transaction epoch differs");
  if (metadata.patch_count != 1 || metadata.level_count < 1 ||
      metadata.level_count > 2 || frame.level() < 0 || frame.level() > 1 ||
      frame.level() >= metadata.level_count)
    throw std::invalid_argument(
        "scratch transaction requires patch zero and level zero or one");
  const int expected_time_refinement_factor = 1 << frame.level();
  if (metadata.level_iteration < 0 ||
      !(step_clock_t(0) < metadata.base_delta_clock) ||
      !finite_positive(metadata.base_delta_time) ||
      metadata.time_refinement_factor != expected_time_refinement_factor ||
      metadata.poison_undefined_values)
    throw std::invalid_argument(
        "scratch transaction time or configuration is unsupported");
  if (!metadata.epoch_reader ||
      metadata.epoch_reader() != metadata.hierarchy_epoch)
    throw std::invalid_argument("scratch transaction epoch reader differs");
  if (metadata.group_pairs.empty() ||
      metadata.live_entry_readers.size() != frame.entry_count())
    throw std::invalid_argument("scratch transaction schema is incomplete");

  std::unordered_set<int> evolved_groups;
  std::unordered_set<int> rhs_groups;
  pair_indices.reserve(metadata.group_pairs.size());
  for (std::size_t pair_ordinal = 0;
       pair_ordinal < metadata.group_pairs.size(); ++pair_ordinal) {
    const auto &pair = metadata.group_pairs[pair_ordinal];
    if (pair.evolved_group < 0 || pair.rhs_group < 0 ||
        pair.evolved_group == pair.rhs_group ||
        (pair_ordinal != 0 &&
         metadata.group_pairs[pair_ordinal - 1].evolved_group >=
             pair.evolved_group) ||
        !evolved_groups.insert(pair.evolved_group).second ||
        !rhs_groups.insert(pair.rhs_group).second)
      throw std::invalid_argument("scratch transaction group pairs overlap");
    const auto evolved = find_group(frame, pair.evolved_group);
    const auto rhs = find_group(frame, pair.rhs_group);
    if (!same_layout(frame.multifab(evolved), frame.multifab(rhs)) ||
        frame.validity(evolved).size() != frame.validity(rhs).size())
      throw std::invalid_argument("scratch transaction pair layouts differ");
    pair_indices.push_back({pair, evolved, rhs});
  }
  for (const auto evolved : evolved_groups)
    if (rhs_groups.count(evolved) != 0)
      throw std::invalid_argument(
          "scratch transaction evolved and RHS sets overlap");

  std::unordered_set<int> invalidation_groups;
  for (const auto group : metadata.dependent_groups) {
    if (group < 0 || !invalidation_groups.insert(group).second)
      throw std::invalid_argument(
          "scratch transaction invalidation groups overlap");
    if (evolved_groups.count(group) != 0)
      throw std::invalid_argument(
          "scratch transaction cannot invalidate an evolved group");
    invalidated_entries.push_back(find_group(frame, group));
  }
  for (const auto rhs : rhs_groups)
    if (invalidation_groups.count(rhs) == 0)
      throw std::invalid_argument(
          "scratch transaction RHS lacks invalidation coverage");

  std::unordered_set<const void *> level_identities;
  std::unordered_set<const void *> group_identities;
  std::unordered_set<const void *> tl0_identities;
  std::unordered_set<const void *> storage_identities;
  std::unordered_set<const void *> validity_identities;
  std::unordered_set<const amrex::MultiFab *> multifabs;
  expected_live.reserve(frame.entry_count());
  for (std::size_t entry = 0; entry < frame.entry_count(); ++entry) {
    if (frame.key(entry).patch != 0 ||
        frame.key(entry).level != frame.level())
      throw std::invalid_argument(
          "scratch transaction frame key differs from requested level");
    if (!metadata.live_entry_readers[entry])
      throw std::invalid_argument("scratch transaction live reader is empty");
    const auto snapshot = metadata.live_entry_readers[entry]();
    if (!same_key(snapshot.key, frame.key(entry)) ||
        snapshot.key.patch != 0 || snapshot.key.level != frame.level() ||
        snapshot.level_identity == nullptr ||
        snapshot.group_identity == nullptr || snapshot.tl0_identity == nullptr ||
        snapshot.storage_identity == nullptr ||
        snapshot.validity_identity == nullptr ||
        snapshot.multifab == nullptr ||
        !snapshot.grid_function_real || !snapshot.restore ||
        !snapshot.multifab->isDefined() ||
        !same_layout(*snapshot.multifab, frame.multifab(entry)) ||
        snapshot.validity.size() != frame.validity(entry).size())
      throw std::invalid_argument("scratch transaction live schema differs");
    level_identities.insert(snapshot.level_identity);
    if (!group_identities.insert(snapshot.group_identity).second ||
        !tl0_identities.insert(snapshot.tl0_identity).second ||
        !storage_identities.insert(snapshot.storage_identity).second ||
        !validity_identities.insert(snapshot.validity_identity).second ||
        !multifabs.insert(snapshot.multifab).second)
      throw std::invalid_argument("scratch transaction live identities overlap");
    expected_live.push_back({snapshot});
  }
  if (level_identities.size() != 1)
    throw std::invalid_argument(
        "scratch transaction live entries span multiple levels");
}

ScratchStateTransaction::ScratchStateTransaction(
    std::unique_ptr<Storage> storage)
    : storage_(std::move(storage)) {
  if (storage_ == nullptr)
    throw std::invalid_argument("scratch transaction storage is missing");
}

ScratchStateTransaction::~ScratchStateTransaction() = default;

const std::vector<ScratchGroupPair> &
ScratchStateTransaction::group_pairs() const noexcept {
  return storage_->group_pairs;
}

const std::vector<int> &
ScratchStateTransaction::dependent_groups() const noexcept {
  return storage_->dependent_groups;
}

ScratchStateToken ScratchStateTransaction::capture_live_evolved() {
  try {
    return ScratchStateTransactionCore::publish_state(
        *storage_, ScratchStateKind::evolved,
        ScratchStateTransactionCore::capture_live_frame(
            *storage_, ScratchStateKind::evolved));
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

ScratchStateToken ScratchStateTransaction::capture_live_rhs() {
  try {
    return ScratchStateTransactionCore::publish_state(
        *storage_, ScratchStateKind::raw_rhs,
        ScratchStateTransactionCore::capture_live_frame(
            *storage_, ScratchStateKind::raw_rhs));
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

ScratchStateToken ScratchStateTransaction::capture_scratch_evolved() {
  return ScratchStateTransactionCore::publish_state(
      *storage_, ScratchStateKind::evolved,
      ScratchStateTransactionCore::capture_working_frame(
          *storage_, ScratchStateKind::evolved));
}

ScratchStateToken ScratchStateTransaction::capture_scratch_rhs() {
  return ScratchStateTransactionCore::publish_state(
      *storage_, ScratchStateKind::raw_rhs,
      ScratchStateTransactionCore::capture_working_frame(
          *storage_, ScratchStateKind::raw_rhs));
}

ScratchStateToken
ScratchStateTransaction::clone_state(const ScratchStateToken &state) {
  const auto &source = ScratchStateTransactionCore::require_state(*storage_, state);
  return ScratchStateTransactionCore::publish_state(
      *storage_, source.kind, source.frame.clone_for_transaction());
}

ScratchStateKind
ScratchStateTransaction::state_kind(const ScratchStateToken &state) const {
  return ScratchStateTransactionCore::require_state(std::as_const(*storage_),
                                                     state)
      .kind;
}

bool ScratchStateTransaction::state_valid(
    const ScratchStateToken &state, const ScratchStateRegion region) const {
  const auto &record = ScratchStateTransactionCore::require_state(
      std::as_const(*storage_), state);
  for (const auto &pair : storage_->pair_indices) {
    for (const auto &component : record.frame.validity(pair.evolved_entry)) {
      const bool valid = region == ScratchStateRegion::interior
                             ? component.interior
                         : region == ScratchStateRegion::outer
                             ? component.outer
                             : component.ghosts;
      if (!valid)
        return false;
    }
  }
  return true;
}

void ScratchStateTransaction::set_state_valid(
    ScratchStateToken &state, const ScratchStateRegion region,
    const bool valid) {
  auto &record = ScratchStateTransactionCore::require_state(*storage_, state);
  for (const auto &pair : storage_->pair_indices) {
    auto &components = record.frame.mutable_validity(pair.evolved_entry);
    for (auto &component : components) {
      if (region == ScratchStateRegion::interior)
        component.interior = valid;
      else if (region == ScratchStateRegion::outer)
        component.outer = valid;
      else
        component.ghosts = valid;
    }
  }
}

void ScratchStateTransaction::rollback_live_evolved(
    const ScratchStateToken &state) {
  // Preserve the distinction between a pre-existing terminal state and a
  // rollback failure: discarded/faulted transactions are rejected without
  // reviving or otherwise changing them.
  ScratchStateTransactionCore::require_available(*storage_);
  try {
    const auto &record =
        ScratchStateTransactionCore::require_state(*storage_, state);
    if (record.kind != ScratchStateKind::evolved)
      throw std::invalid_argument(
          "live rollback requires an evolved state token");

    // Validate every live identity and every source layout before the first
    // in-place write. The callbacks remain factory-owned, so the public API
    // exposes no driver representation.
    const auto snapshots =
        ScratchStateTransactionCore::read_and_validate_live(*storage_);
    if (record.frame.entry_count() != snapshots.size())
      throw std::runtime_error("live rollback state schema changed");
    for (std::size_t entry = 0; entry < snapshots.size(); ++entry) {
      const auto &snapshot = snapshots[entry];
      if (!snapshot.restore ||
          !same_key(record.frame.key(entry), snapshot.key) ||
          !same_layout(record.frame.multifab(entry), *snapshot.multifab) ||
          record.frame.validity(entry).size() != snapshot.validity.size())
        throw std::runtime_error("live rollback state layout changed");
    }

    for (std::size_t entry = 0; entry < snapshots.size(); ++entry)
      snapshots[entry].restore(record.frame.multifab(entry),
                               record.frame.validity(entry));

    // Rollback is terminal for this attempted primary step. In particular,
    // no PostStep or RHS schedule is rerun after live TL0 is restored.
    discard();
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

void ScratchStateTransaction::restore_state(const ScratchStateToken &state) {
  const auto &record =
      ScratchStateTransactionCore::require_state(*storage_, state);
  if (record.kind != ScratchStateKind::evolved)
    throw std::invalid_argument("only evolved state has a full-frame restore");
  storage_->working_frame->restore_from_transaction(record.frame);
}

void ScratchStateTransaction::restore_left(
    const ScratchStateToken &evolved, const ScratchStateToken &raw_rhs) {
  const auto &evolved_record =
      ScratchStateTransactionCore::require_state(*storage_, evolved);
  const auto &rhs_record =
      ScratchStateTransactionCore::require_state(*storage_, raw_rhs);
  if (evolved_record.kind != ScratchStateKind::evolved ||
      rhs_record.kind != ScratchStateKind::raw_rhs)
    throw std::invalid_argument("scratch left restore state kinds differ");
  storage_->working_frame->restore_from_transaction(evolved_record.frame);
  for (const auto &pair : storage_->pair_indices) {
    storage_->working_frame->copy_interior_for_transaction(
        pair.rhs_entry, rhs_record.frame, pair.evolved_entry);
    storage_->working_frame->copy_interior_validity_for_transaction(
        pair.rhs_entry, rhs_record.frame, pair.evolved_entry, false);
  }
}

void ScratchStateTransaction::linear_combination(
    ScratchStateToken &destination, double destination_scale,
    const std::vector<ScratchLinearTerm> &terms) {
  auto &destination_record =
      ScratchStateTransactionCore::require_state(*storage_, destination);
  if (!std::isfinite(destination_scale))
    throw std::invalid_argument(
        "scratch combination destination scale is not finite");

  std::vector<double> weights;
  std::vector<const ScratchLevelFrame *> sources;
  weights.reserve(terms.size());
  sources.reserve(terms.size());
  for (const auto &term : terms) {
    if (!std::isfinite(term.coefficient) || term.state == nullptr)
      throw std::invalid_argument("scratch combination term is invalid");
    const auto &source =
        ScratchStateTransactionCore::require_state(*storage_, *term.state);
    if (destination_record.kind == ScratchStateKind::raw_rhs &&
        source.kind != ScratchStateKind::raw_rhs)
      throw std::invalid_argument(
          "raw RHS combination requires raw RHS sources");
    if (term.coefficient == 0.0)
      continue;
    if (term.state->state_ == destination.state_) {
      destination_scale += term.coefficient;
      if (!std::isfinite(destination_scale))
        throw std::invalid_argument(
            "scratch combination alias scale is not finite");
      continue;
    }
    weights.push_back(term.coefficient);
    sources.push_back(&source.frame);
  }

  std::vector<std::size_t> source_entries(sources.size());
  for (const auto &pair : storage_->pair_indices) {
    std::fill(source_entries.begin(), source_entries.end(),
              pair.evolved_entry);
    destination_record.frame.linear_combination_interior_for_transaction(
        pair.evolved_entry, destination_scale, weights, sources,
        source_entries);
    auto &destination_validity =
        destination_record.frame.mutable_validity(pair.evolved_entry);
    for (std::size_t component = 0;
         component < destination_validity.size(); ++component) {
      bool interior_valid =
          destination_scale == 0.0
              ? true
              : destination_validity[component].interior;
      for (std::size_t source = 0; source < sources.size(); ++source)
        if (weights[source] != 0.0)
          interior_valid =
              interior_valid &&
              sources[source]
                  ->validity(pair.evolved_entry)[component]
                  .interior;
      destination_validity[component].interior = interior_valid;
    }
  }
}

void ScratchStateTransaction::post_step_after_update(
    const StepContext &context, const StagePoint &stage_point) {
  try {
    const auto coordinates = ScratchStateTransactionCore::stage_coordinates(
        *storage_, context, stage_point);
    if (current_scratch_state_transaction() != this)
      throw std::logic_error(
          "active StepContext belongs to another scratch transaction");
    for (const auto entry : storage_->invalidated_entries) {
      auto &components = storage_->working_frame->mutable_validity(entry);
      for (auto &component : components)
        component = {false, false, false};
    }
#ifdef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
    if (!storage_->test_executor)
      throw std::logic_error("scratch test executor is missing");
    const auto receipt = storage_->test_executor(
        ScratchSchedulePhase::post_step, context, coordinates);
#else
    const auto receipt = storage_->executor->execute_post_step(context,
                                                                coordinates);
#endif
    if (receipt.phase != ScratchSchedulePhase::post_step)
      throw std::runtime_error("scratch PostStep receipt phase differs");
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

void ScratchStateTransaction::evaluate_rhs(const StepContext &context,
                                           const StagePoint &stage_point) {
  try {
    const auto coordinates = ScratchStateTransactionCore::stage_coordinates(
        *storage_, context, stage_point);
    if (current_scratch_state_transaction() != this)
      throw std::logic_error(
          "active StepContext belongs to another scratch transaction");
#ifdef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
    if (!storage_->test_executor)
      throw std::logic_error("scratch test executor is missing");
    const auto receipt = storage_->test_executor(ScratchSchedulePhase::rhs,
                                                  context, coordinates);
#else
    const auto receipt = storage_->executor->execute_rhs(context, coordinates);
#endif
    if (receipt.phase != ScratchSchedulePhase::rhs)
      throw std::runtime_error("scratch RHS receipt phase differs");
    ++storage_->rhs_evaluations;
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

std::size_t ScratchStateTransaction::rhs_evaluation_count() const noexcept {
  return storage_->rhs_evaluations;
}

void ScratchStateTransaction::commit_dense(
    const StepContext &context, const DenseOutputProvider &provider,
    const DenseIntervalId &interval,
    const std::vector<ScratchDenseSampleRef> &samples) {
  ScratchStateTransactionCore::require_available(*storage_);
  if (storage_->dense_committed_once || storage_->committed_dense != nullptr)
    throw std::logic_error("scratch dense interval is already committed");
  try {
    if (current_step_context() != &context ||
        current_scratch_state_transaction() != this)
      throw std::logic_error(
          "scratch dense commit requires its active StepContext scope");
    const auto &provider_capability = provider.capability();
    const auto level_delta_clock =
        storage_->base_delta_clock / storage_->time_refinement_factor;
    const auto level_delta_time =
        storage_->base_delta_time /
        static_cast<double>(storage_->time_refinement_factor);
    if (context.level != storage_->level ||
        interval.level != context.level ||
        interval.begin_clock != context.begin_clock ||
        interval.end_clock != context.end_clock ||
        interval.begin_time != context.begin_time ||
        interval.end_time != context.end_time ||
        interval.method != context.method ||
        provider_capability.method != context.method ||
        interval.tableau_fingerprint !=
            provider_capability.tableau_fingerprint ||
        interval.end_clock - interval.begin_clock != level_delta_clock ||
        !close_time(interval.end_time - interval.begin_time, level_delta_time,
                    interval.begin_time, interval.end_time))
      throw std::invalid_argument("scratch dense interval clock differs");
    const auto &constraints =
        reference_dense_stencil(interval.method).specification().constraints;
    if (samples.size() != constraints.size())
      throw std::invalid_argument("scratch dense sample count differs");

    std::unordered_set<std::uint64_t> unique_states;
    std::vector<std::unique_ptr<OwnedMultiFabDenseState>> owned_samples;
    std::vector<DenseMFabRawSampleView> sample_views;
    owned_samples.reserve(samples.size());
    sample_views.reserve(samples.size());
    for (std::size_t sample = 0; sample < samples.size(); ++sample) {
      const auto &reference = samples[sample];
      if (reference.state == nullptr ||
          !same_constraint(sample_constraint(reference), constraints[sample]))
        throw std::invalid_argument("scratch dense sample metadata differs");
      const auto &record = ScratchStateTransactionCore::require_state(
          *storage_, *reference.state);
      const auto required_kind =
          reference.kind == ScratchDenseSampleKind::value
              ? ScratchStateKind::evolved
              : ScratchStateKind::raw_rhs;
      if (record.kind != required_kind ||
          !unique_states.insert(reference.state->state_).second)
        throw std::invalid_argument("scratch dense sample token differs");
      for (const auto &pair : storage_->pair_indices)
        for (const auto &component :
             record.frame.validity(pair.evolved_entry))
          if (!component.interior)
            throw std::invalid_argument(
                "scratch dense sample interior is invalid");

      std::vector<DenseMFabView> views;
      views.reserve(storage_->pair_indices.size());
      for (const auto &pair : storage_->pair_indices) {
        const auto &key = record.frame.key(pair.evolved_entry);
        views.push_back({{key.hierarchy_epoch, key.patch, key.level,
                          pair.groups.evolved_group},
                         &record.frame.multifab(pair.evolved_entry)});
      }
      owned_samples.push_back(OwnedMultiFabDenseState::copy_of(views));
      sample_views.push_back({constraints[sample], owned_samples.back().get()});
    }

    auto controls = build_reference_dense_controls(
        interval.method, level_delta_time, sample_views);
    auto builder = provider.begin_interval(interval);
    for (auto &control : std::move(controls).release_controls())
      builder->add_control(std::move(control));
    auto sealed = builder->seal();
    storage_->committed_dense = sealed;
    storage_->dense_committed_once = true;
  } catch (...) {
    storage_->fault_and_discard();
    throw;
  }
}

std::shared_ptr<const DenseInterval>
ScratchStateTransaction::take_committed_dense() noexcept {
  try {
    if (storage_->discarded || storage_->faulted)
      return nullptr;
    if (!storage_->epoch_reader ||
        storage_->epoch_reader() != storage_->hierarchy_epoch) {
      storage_->fault_and_discard();
      return nullptr;
    }
    auto result = std::move(storage_->committed_dense);
    storage_->committed_dense.reset();
    return result;
  } catch (...) {
    storage_->fault_and_discard();
    return nullptr;
  }
}

bool ScratchStateTransaction::faulted() const noexcept {
  return storage_->faulted;
}

bool ScratchStateTransaction::discarded() const noexcept {
  return storage_->discarded;
}

void ScratchStateTransaction::discard() noexcept {
  if (storage_->discarded)
    return;
  storage_->discarded = true;
  storage_->states.clear();
  storage_->committed_dense.reset();
  storage_->release_execution_storage();
}

#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
std::unique_ptr<ScratchStateTransaction>
ScratchStateTransactionFactory::create_native(
    const GHExt &ghext, CertifiedLocalScheduleRegistry &registry,
    CertifiedScratchLevelFrame &&frame,
    ScratchStateTransactionFactoryMetadata metadata) {
  DECLARE_CCTK_PARAMETERS;
  const int patch = frame.patch();
  const int level = frame.frame().level();
  if (ghext.patchdata.size() != 1 || patch != 0 ||
      ghext.patchdata.front().leveldata.empty() ||
      ghext.patchdata.front().leveldata.size() > 2 || level < 0 ||
      level > 1 ||
      static_cast<std::size_t>(level) >=
          ghext.patchdata.front().leveldata.size())
    throw std::invalid_argument(
        "native scratch transaction requires patch zero and level zero or one");
  if (!metadata.live_entry_readers.empty())
    throw std::invalid_argument(
        "native scratch transaction constructs its own live readers");
  const auto &coarse_level = ghext.patchdata.front().leveldata.front();
  const auto &native_level = ghext.patchdata.front().leveldata.at(
      static_cast<std::size_t>(level));
  const auto *const native_gh = native_level.get_patch_cctkGH();
  const int expected_time_refinement_factor = 1 << level;
  const auto coarse_delta_clock =
      step_clock_t(coarse_level.delta_iteration.num,
                   coarse_level.delta_iteration.den);
  const auto native_delta_clock =
      step_clock_t(native_level.delta_iteration.num,
                   native_level.delta_iteration.den);
  if (native_gh == nullptr || native_level.patch != patch ||
      native_level.level != level || native_gh->cctk_iteration < 0 ||
      native_gh->cctk_timefac != expected_time_refinement_factor ||
      !(step_clock_t(0) < coarse_delta_clock) ||
      native_delta_clock !=
          coarse_delta_clock / expected_time_refinement_factor ||
      !finite_positive(native_gh->cctk_delta_time))
    throw std::invalid_argument(
        "native scratch transaction time envelope is unsupported");
  metadata.hierarchy_epoch =
      static_cast<std::int64_t>(CarpetX_GetEpoch());
  metadata.patch_count = static_cast<int>(ghext.patchdata.size());
  metadata.level_count =
      static_cast<int>(ghext.patchdata.front().leveldata.size());
  metadata.level_iteration = native_gh->cctk_iteration;
  metadata.base_delta_clock = coarse_delta_clock;
  metadata.base_delta_time = native_gh->cctk_delta_time;
  metadata.time_refinement_factor = native_gh->cctk_timefac;
  metadata.poison_undefined_values = poison_undefined_values;
  metadata.epoch_reader = [] {
    return static_cast<std::int64_t>(CarpetX_GetEpoch());
  };
  metadata.live_entry_readers.reserve(frame.frame().entry_count());
  for (std::size_t entry = 0; entry < frame.frame().entry_count(); ++entry) {
    const int group_index = frame.frame().key(entry).group_index;
    const auto epoch = metadata.hierarchy_epoch;
    metadata.live_entry_readers.push_back(
        [&ghext, group_index, epoch, level] {
          const auto &leveldata = ghext.patchdata.at(0).leveldata.at(
              static_cast<std::size_t>(level));
          if (group_index < 0 ||
              static_cast<std::size_t>(group_index) >=
                  leveldata.groupdata.size() ||
              leveldata.groupdata[static_cast<std::size_t>(group_index)] ==
                  nullptr)
            throw std::runtime_error(
                "native scratch live group identity disappeared");
          const auto &groupdata = *leveldata.groupdata[
              static_cast<std::size_t>(group_index)];
          if (groupdata.mfab.empty() || groupdata.valid.empty())
            throw std::runtime_error(
                "native scratch live TL0 identity disappeared");
          const auto *const multifab = groupdata.mfab.at(0).get();
          const auto &why_valid = groupdata.valid.at(0);
          const auto *const level_identity = &leveldata;
          const auto *const group_identity = &groupdata;
          const auto *const tl0_identity = &groupdata.mfab.at(0);
          const auto *const validity_identity = &why_valid;
          std::vector<ScratchValidity> validity;
          validity.reserve(why_valid.size());
          for (const auto &why : why_valid) {
            const auto &bits = why.get();
            validity.push_back({bits.valid_int, bits.valid_outer,
                                bits.valid_ghosts});
          }
          cGroup cactus_group{};
          const bool grid_function_real =
              CCTK_GroupData(group_index, &cactus_group) == 0 &&
              cactus_group.grouptype == CCTK_GF &&
              cactus_group.vartype == CCTK_VARIABLE_REAL;
          ScratchLiveEntryRestorer restore =
              [&ghext, group_index, epoch, level, level_identity,
               group_identity, tl0_identity, multifab, validity_identity](
                  const amrex::MultiFab &state,
                  const std::vector<ScratchValidity> &restored_validity) {
                if (static_cast<std::int64_t>(CarpetX_GetEpoch()) != epoch)
                  throw std::runtime_error(
                      "native live rollback hierarchy epoch changed");
                auto &mutable_ghext = const_cast<GHExt &>(ghext);
                auto &current_level =
                    mutable_ghext.patchdata.at(0).leveldata.at(
                        static_cast<std::size_t>(level));
                if (group_index < 0 ||
                    static_cast<std::size_t>(group_index) >=
                        current_level.groupdata.size() ||
                    current_level.groupdata[
                        static_cast<std::size_t>(group_index)] == nullptr)
                  throw std::runtime_error(
                      "native live rollback group identity disappeared");
                auto &current_group = *current_level.groupdata[
                    static_cast<std::size_t>(group_index)];
                if (current_group.mfab.empty() ||
                    current_group.valid.empty())
                  throw std::runtime_error(
                      "native live rollback TL0 identity disappeared");
                auto *const destination = current_group.mfab.at(0).get();
                auto &current_validity = current_group.valid.at(0);
                if (&current_level != level_identity ||
                    &current_group != group_identity ||
                    &current_group.mfab.at(0) != tl0_identity ||
                    destination != multifab ||
                    &current_validity != validity_identity ||
                    destination == nullptr || !destination->isDefined() ||
                    !same_layout(*destination, state) ||
                    current_validity.size() != restored_validity.size())
                  throw std::runtime_error(
                      "native live rollback identity or layout changed");

                amrex::MultiFab::Copy(*destination, state, 0, 0,
                                      state.nComp(), state.nGrowVect());
                for (std::size_t component = 0;
                     component < restored_validity.size(); ++component) {
                  valid_t bits;
                  bits.valid_int = restored_validity[component].interior;
                  bits.valid_outer = restored_validity[component].outer;
                  bits.valid_ghosts = restored_validity[component].ghosts;
                  current_validity[component].set_all(bits, [] {
                    return "ScratchStateTransaction primary rollback";
                  });
                }
              };
          return ScratchLiveEntrySnapshot{
              {epoch, level, 0, group_index},
              level_identity,
              group_identity,
              tl0_identity,
              multifab,
              validity_identity,
              multifab,
              std::move(validity),
              grid_function_real,
              std::move(restore)};
        });
  }
  auto storage = std::make_unique<ScratchStateTransaction::Storage>();
  ScratchStateTransactionCore::validate_factory_metadata(
      frame.frame(), metadata, storage->pair_indices,
      storage->invalidated_entries, storage->expected_live);
  storage->hierarchy_epoch = metadata.hierarchy_epoch;
  storage->group_pairs = metadata.group_pairs;
  storage->dependent_groups = metadata.dependent_groups;
  std::sort(storage->dependent_groups.begin(),
            storage->dependent_groups.end());
  storage->level = level;
  storage->level_iteration = metadata.level_iteration;
  storage->base_delta_clock = metadata.base_delta_clock;
  storage->base_delta_time = metadata.base_delta_time;
  storage->time_refinement_factor = metadata.time_refinement_factor;
  storage->epoch_reader = std::move(metadata.epoch_reader);
  storage->live_readers = std::move(metadata.live_entry_readers);
  storage->binding = ScratchLocalLevelBinding::bind(ghext, std::move(frame));
  storage->working_frame = &storage->binding->frame_for_transaction();
  storage->executor =
      CertifiedScratchStageExecutor::create(registry, *storage->binding);
  return std::unique_ptr<ScratchStateTransaction>(
      new ScratchStateTransaction(std::move(storage)));
}
#else
std::unique_ptr<ScratchStateTransaction>
ScratchStateTransactionFactory::create_native(
    const GHExt &, CertifiedLocalScheduleRegistry &,
    CertifiedScratchLevelFrame &&,
    ScratchStateTransactionFactoryMetadata) {
  throw std::logic_error("native scratch transaction is not in the unit gate");
}

std::unique_ptr<ScratchStateTransaction>
ScratchStateTransactionFactory::create_for_test(
    ScratchLevelFrame working_frame,
    ScratchStateTransactionFactoryMetadata metadata,
    TestScheduleExecutor executor) {
  if (!executor)
    throw std::invalid_argument("scratch test executor is missing");
  auto storage = std::make_unique<ScratchStateTransaction::Storage>();
  ScratchStateTransactionCore::validate_factory_metadata(
      working_frame, metadata, storage->pair_indices,
      storage->invalidated_entries, storage->expected_live);
  storage->hierarchy_epoch = metadata.hierarchy_epoch;
  storage->group_pairs = metadata.group_pairs;
  storage->dependent_groups = metadata.dependent_groups;
  std::sort(storage->dependent_groups.begin(),
            storage->dependent_groups.end());
  storage->level = working_frame.level();
  storage->level_iteration = metadata.level_iteration;
  storage->base_delta_clock = metadata.base_delta_clock;
  storage->base_delta_time = metadata.base_delta_time;
  storage->time_refinement_factor = metadata.time_refinement_factor;
  storage->epoch_reader = std::move(metadata.epoch_reader);
  storage->live_readers = std::move(metadata.live_entry_readers);
  storage->test_executor = std::move(executor);
  storage->test_frame.emplace(std::move(working_frame));
  storage->working_frame = &*storage->test_frame;
  return std::unique_ptr<ScratchStateTransaction>(
      new ScratchStateTransaction(std::move(storage)));
}

const void *ScratchStateTransactionFactory::working_entry_address_for_test(
    const ScratchStateTransaction &transaction, const std::size_t entry) {
  ScratchStateTransactionCore::require_available(*transaction.storage_);
  return transaction.storage_->working_frame->entry_address_for_transaction(
      entry);
}

const amrex::MultiFab &
ScratchStateTransactionFactory::working_multifab_for_test(
    const ScratchStateTransaction &transaction, const std::size_t entry) {
  ScratchStateTransactionCore::require_available(*transaction.storage_);
  return transaction.storage_->working_frame->multifab(entry);
}

amrex::MultiFab &ScratchStateTransactionFactory::mutable_working_multifab_for_test(
    ScratchStateTransaction &transaction, const std::size_t entry) {
  ScratchStateTransactionCore::require_available(*transaction.storage_);
  return transaction.storage_->working_frame->mutable_multifab(entry);
}

const std::vector<ScratchValidity> &
ScratchStateTransactionFactory::working_validity_for_test(
    const ScratchStateTransaction &transaction, const std::size_t entry) {
  ScratchStateTransactionCore::require_available(*transaction.storage_);
  return transaction.storage_->working_frame->validity(entry);
}

std::vector<ScratchValidity> &
ScratchStateTransactionFactory::mutable_working_validity_for_test(
    ScratchStateTransaction &transaction, const std::size_t entry) {
  ScratchStateTransactionCore::require_available(*transaction.storage_);
  return transaction.storage_->working_frame->mutable_validity(entry);
}

ScratchStateToken ScratchStateTransactionFactory::stale_epoch_token_for_test(
    const ScratchStateTransaction &transaction,
    const ScratchStateToken &source) {
  auto token = ScratchStateToken(source.owner_, source.state_, source.epoch_ + 1,
                                 source.schema_, source.kind_);
  if (token.epoch_ == transaction.storage_->hierarchy_epoch)
    ++token.epoch_;
  return token;
}

ScratchStateToken
ScratchStateTransactionFactory::incompatible_schema_token_for_test(
    const ScratchStateTransaction &, const ScratchStateToken &source) {
  return ScratchStateToken(source.owner_, source.state_, source.epoch_,
                           source.schema_ + 1, source.kind_);
}
#endif

} // namespace CarpetX
