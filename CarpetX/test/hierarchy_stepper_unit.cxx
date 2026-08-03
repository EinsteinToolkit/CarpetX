#include "hierarchy_stepper.hxx"

#include <cassert>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

using CarpetX::DenseCapability;
using CarpetX::DenseInterval;
using CarpetX::DenseIntervalBuilder;
using CarpetX::DenseIntervalId;
using CarpetX::DenseOutputProvider;
using CarpetX::DenseOutputRegistry;
using CarpetX::DenseStateVector;
using CarpetX::HierarchyAdvanceResult;
using CarpetX::HierarchyEvolutionAdapter;
using CarpetX::HierarchyStepper;
using CarpetX::HierarchyStepperConfig;
using CarpetX::LevelStepSession;
using CarpetX::LevelAdvanceConfig;
using CarpetX::LevelAdvanceResult;
using CarpetX::ScratchStateTransaction;
using CarpetX::StepContext;
using CarpetX::SubcyclingODEMethod;
using CarpetX::TableauFingerprint;
using CarpetX::hierarchy_time_t;

std::string time_string(const hierarchy_time_t time) {
  return std::to_string(time.num) + "/" + std::to_string(time.den);
}

TableauFingerprint fingerprint(const std::uint8_t marker) {
  TableauFingerprint result{};
  result.fill(marker);
  return result;
}

DenseCapability capability(const SubcyclingODEMethod method,
                           const TableauFingerprint &tableau_fingerprint) {
  switch (method) {
  case SubcyclingODEMethod::rk4:
    return DenseCapability{method, tableau_fingerprint, 4, 4, 4, 4, 5,
                           true, true};
  case SubcyclingODEMethod::rkf78_order7:
    return DenseCapability{method, tableau_fingerprint, 7, 7, 11, 23, 8,
                           true, true};
  case SubcyclingODEMethod::dp87_order8:
    return DenseCapability{method, tableau_fingerprint, 8, 8, 13, 39, 9,
                           true, true};
  }
  throw std::logic_error("test requested an unknown method");
}

class ScalarState final : public DenseStateVector {
public:
  explicit ScalarState(const double value) : value(value) {}

  bool compatible(const DenseStateVector &other) const noexcept override {
    return dynamic_cast<const ScalarState *>(&other) != nullptr;
  }

  void copy_from(const DenseStateVector &other) override {
    value = dynamic_cast<const ScalarState &>(other).value;
  }

  void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const DenseStateVector *> &sources) override {
    assert(weights.size() == sources.size());
    value = 0.0;
    for (std::size_t index = 0; index < weights.size(); ++index)
      value += weights[index] *
               dynamic_cast<const ScalarState &>(*sources[index]).value;
  }

  double value;
};

std::shared_ptr<const DenseInterval>
make_linear_time_interval(const StepContext &context,
                          const DenseCapability &dense_capability) {
  DenseIntervalBuilder builder(
      dense_capability,
      DenseIntervalId{context.level, context.begin_clock, context.end_clock,
                      context.begin_time, context.end_time, context.method,
                      dense_capability.tableau_fingerprint});
  const int degree = dense_capability.persistent_vector_count - 1;
  for (int control = 0; control <= degree; ++control) {
    const double fraction = static_cast<double>(control) / degree;
    const double value = context.begin_time +
                         fraction * (context.end_time - context.begin_time);
    builder.add_control(std::make_unique<ScalarState>(value));
  }
  return builder.seal();
}

struct StageObservation {
  int level;
  hierarchy_time_t begin_clock;
  hierarchy_time_t end_clock;
  double stage_time;
  int parent_level;
  double parent_value;
};

struct DenseLifetimeRecord {
  int level;
  hierarchy_time_t begin_clock;
  hierarchy_time_t end_clock;
  std::weak_ptr<const DenseInterval> interval;
};

class RecordingAdapter : public HierarchyEvolutionAdapter {
public:
  explicit RecordingAdapter(std::vector<DenseCapability> capabilities,
                            const hierarchy_time_t initial_clock = 0,
                            const double initial_physical_time = 0.0,
                            const double coarse_dt = 1.0)
      : capabilities_(std::move(capabilities)),
        transaction_slots_(capabilities_.size()),
        initial_clock_(initial_clock),
        initial_physical_time_(initial_physical_time), coarse_dt_(coarse_dt) {}

  std::vector<std::string> events;
  std::vector<StageObservation> stages;
  std::vector<DenseLifetimeRecord> dense_lifetimes;
  int throw_on_level{-1};
  bool omit_required_dense{false};
  bool round_child_endpoint_up{false};
  bool throw_during_sync{false};
  bool observers_saw_all_dense_released{false};
  bool observer_stop_requested{false};
  bool throw_during_observer{false};
  std::function<void()> sync_observer_hook;
  std::size_t begin_level_step_calls{0};
  std::size_t committed_level_sessions{0};
  std::size_t discarded_level_sessions{0};
  std::vector<std::pair<std::string, ScratchStateTransaction *>>
      stage_transactions;

  class Session final : public LevelStepSession {
  public:
    Session(RecordingAdapter &adapter, StepContext context,
            const bool require_dense)
        : adapter_(adapter), context_(std::move(context)),
          require_dense_(require_dense) {}

    ~Session() override {
      if (!committed_)
        ++adapter_.discarded_level_sessions;
    }

    ScratchStateTransaction *transaction() noexcept override {
      return adapter_.transaction_for_level(context_.level);
    }

    LevelAdvanceResult advance() override {
      return adapter_.advance_level_impl(context_, require_dense_);
    }

    void commit() override {
      assert(!committed_);
      assert(CarpetX::current_step_context() == nullptr);
      assert(CarpetX::current_scratch_state_transaction() == nullptr);
      committed_ = true;
      ++adapter_.committed_level_sessions;
    }

  private:
    RecordingAdapter &adapter_;
    StepContext context_;
    bool require_dense_;
    bool committed_{false};
  };

  std::unique_ptr<LevelStepSession>
  begin_level_step(const StepContext &context,
                   const bool require_dense) override {
    ++begin_level_step_calls;
    assert(CarpetX::current_step_context() == nullptr);
    assert(CarpetX::current_scratch_state_transaction() == nullptr);
    assert(context.require_dense_output == require_dense);
    assert(context.endpoint_accepted_step > 0);
    return std::make_unique<Session>(*this, context, require_dense);
  }

  virtual LevelAdvanceResult advance_level_impl(const StepContext &context,
                                                const bool require_dense) {
    assert(CarpetX::current_scratch_state_transaction() ==
           transaction_for_level(context.level));
    assert_covering_ancestor_dense_is_live(context);
    events.push_back("advance L" + std::to_string(context.level) + " [" +
                     time_string(context.begin_clock) + "," +
                     time_string(context.end_clock) + "] dense=" +
                     (require_dense ? "required" : "none"));
    if (context.level == throw_on_level)
      throw std::runtime_error("requested test advance failure");

    next_stage_kind_ = "primary";
    CarpetX::prepare_stage(
        {CarpetX::StageKind::primary, 1, 3, hierarchy_time_t(0),
         context.begin_time},
        context.method);
    next_stage_kind_ = "auxiliary";
    CarpetX::prepare_stage(
        {CarpetX::StageKind::fractional, 2, 3,
         hierarchy_time_t(1, 2),
         0.5 * (context.begin_time + context.end_time)},
        context.method);
    const double endpoint =
        round_child_endpoint_up && context.level > 0
            ? std::nextafter(context.end_time,
                             std::numeric_limits<double>::infinity())
            : context.end_time;
    next_stage_kind_ = "primary";
    CarpetX::prepare_stage(
        {CarpetX::StageKind::primary, 3, 3, hierarchy_time_t(1), endpoint},
        context.method);

    if (!require_dense || omit_required_dense)
      return {};

    auto interval = make_linear_time_interval(
        context, capabilities_.at(static_cast<std::size_t>(context.level)));
    dense_lifetimes.push_back(DenseLifetimeRecord{
        context.level, context.begin_clock, context.end_clock, interval});
    return LevelAdvanceResult{std::move(interval)};
  }

  void prepare_stage(const StepContext &context,
                      const CarpetX::StagePoint &stage_point,
                      const DenseInterval *parent_dense) override {
    if (next_stage_kind_ == "auxiliary")
      assert(stage_point.kind == CarpetX::StageKind::fractional);
    else
      assert(stage_point.kind == CarpetX::StageKind::primary);
    stage_transactions.emplace_back(
        next_stage_kind_, CarpetX::current_scratch_state_transaction());
    assert(stage_transactions.back().second ==
           transaction_for_level(context.level));
    if (context.level == 0) {
      assert(parent_dense == nullptr);
      stages.push_back(StageObservation{context.level, context.begin_clock,
                                        context.end_clock,
                                        stage_point.stage_time, -1,
                                        stage_point.stage_time});
      return;
    }

    assert(parent_dense != nullptr);
    assert(parent_dense->id().level == context.level - 1);
    assert(parent_dense->id().begin_clock <= context.begin_clock);
    assert(parent_dense->id().end_clock >= context.end_clock);
    const auto exact_stage_clock =
        context.begin_clock +
        stage_point.stage_fraction * (context.end_clock - context.begin_clock);
    const auto exact_theta =
        (exact_stage_clock - parent_dense->id().begin_clock) /
        (parent_dense->id().end_clock - parent_dense->id().begin_clock);
    const double theta = static_cast<double>(exact_theta);
    ScalarState parent_value(-1.0);
    parent_dense->evaluate(theta, parent_value);
    assert(std::abs(parent_value.value - stage_point.stage_time) < 1.0e-13);
    stages.push_back(StageObservation{
        context.level, context.begin_clock, context.end_clock,
        stage_point.stage_time,
        parent_dense->id().level, parent_value.value});
  }

  void begin_synchronization(const int coarse_level, const int fine_level,
                             const hierarchy_time_t time) override {
    assert_dense_live(coarse_level, time);
    events.push_back("sync_begin L" + std::to_string(coarse_level) +
                     "<-L" + std::to_string(fine_level) + " at " +
                     time_string(time));
  }

  void synchronize_levels(const int coarse_level, const int fine_level,
                          const hierarchy_time_t time) override {
    assert_dense_live(coarse_level, time);
    events.push_back("synchronize L" + std::to_string(coarse_level) +
                     "<-L" + std::to_string(fine_level) + " at " +
                     time_string(time));
    if (throw_during_sync)
      throw std::runtime_error("requested synchronization failure");
  }

  void end_synchronization(const int coarse_level, const int fine_level,
                           const hierarchy_time_t time) noexcept override {
    assert_dense_live(coarse_level, time);
    events.push_back("sync_end L" + std::to_string(coarse_level) + "<-L" +
                     std::to_string(fine_level) + " at " +
                     time_string(time));
  }

  void run_sync_observers(const hierarchy_time_t time,
                          const double physical_time,
                          const std::uint64_t completed_epoch,
                          const bool stop_requested) override {
    observer_stop_requested = stop_requested;
    observers_saw_all_dense_released = true;
    for (const auto &record : dense_lifetimes)
      observers_saw_all_dense_released &= record.interval.expired();
    const double expected_physical_time =
        initial_physical_time_ +
        static_cast<double>(time - initial_clock_) * coarse_dt_;
    assert(physical_time == expected_physical_time);
    events.push_back("sync_observers epoch=" +
                     std::to_string(completed_epoch) + " at " +
                     time_string(time));
    if (sync_observer_hook)
      sync_observer_hook();
    if (throw_during_observer)
      throw std::runtime_error("requested sync-observer failure");
  }

private:
  ScratchStateTransaction *transaction_for_level(const int level) noexcept {
    return reinterpret_cast<ScratchStateTransaction *>(
        &transaction_slots_.at(static_cast<std::size_t>(level)));
  }

  void assert_covering_ancestor_dense_is_live(
      const StepContext &context) const {
    for (const auto &record : dense_lifetimes) {
      const bool covers = record.level < context.level &&
                          record.begin_clock <= context.begin_clock &&
                          record.end_clock >= context.end_clock;
      if (covers)
        assert(!record.interval.expired());
    }
  }

  void assert_dense_live(const int level, const hierarchy_time_t end) const {
    for (auto record = dense_lifetimes.rbegin();
         record != dense_lifetimes.rend(); ++record) {
      if (record->level == level && record->end_clock == end) {
        assert(!record->interval.expired());
        return;
      }
    }
    assert(false && "synchronization has no retained dense interval");
  }

  std::vector<DenseCapability> capabilities_;
  std::vector<std::uintptr_t> transaction_slots_;
  std::string next_stage_kind_;
  hierarchy_time_t initial_clock_{0};
  double initial_physical_time_{0.0};
  double coarse_dt_{1.0};
};

class ConcurrentStopAdapter final : public RecordingAdapter {
public:
  using RecordingAdapter::RecordingAdapter;

  std::atomic<bool> root_entered{false};
  std::atomic<bool> release_root{false};

  LevelAdvanceResult advance_level_impl(const StepContext &context,
                                        const bool require_dense) override {
    if (context.level == 0) {
      root_entered.store(true, std::memory_order_relaxed);
      while (!release_root.load(std::memory_order_relaxed))
        std::this_thread::yield();
    }
    return RecordingAdapter::advance_level_impl(context, require_dense);
  }
};

struct Fixture {
  DenseOutputRegistry registry;
  std::vector<DenseCapability> capabilities;
  HierarchyStepperConfig config;

  explicit Fixture(const int level_count) {
    assert(level_count > 0);
    capabilities.reserve(static_cast<std::size_t>(level_count));
    config.levels.reserve(static_cast<std::size_t>(level_count));
    for (int level = 0; level < level_count; ++level) {
      const auto method = level % 3 == 0
                              ? SubcyclingODEMethod::rk4
                              : level % 3 == 1
                                    ? SubcyclingODEMethod::rkf78_order7
                                    : SubcyclingODEMethod::dp87_order8;
      const auto table = fingerprint(static_cast<std::uint8_t>(0x31 + level));
      capabilities.push_back(capability(method, table));
      config.levels.push_back(LevelAdvanceConfig{method, table});
      registry.register_provider(
          std::make_shared<DenseOutputProvider>(capabilities.back()));
    }
  }
};

void test_two_levels_use_parent_dense_at_every_exact_child_stage() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);

  const HierarchyAdvanceResult result = stepper.advance_one_epoch(adapter);

  assert(result.synchronized_time == hierarchy_time_t(1));
  assert(result.synchronized_physical_time == 1.0);
  assert(result.epoch == 1);
  assert(!result.stop_requested);
  assert(adapter.events == std::vector<std::string>({
                               "advance L0 [0/1,1/1] dense=required",
                               "advance L1 [0/1,1/2] dense=none",
                               "advance L1 [1/2,1/1] dense=none",
                               "sync_begin L0<-L1 at 1/1",
                               "synchronize L0<-L1 at 1/1",
                               "sync_end L0<-L1 at 1/1",
                               "sync_observers epoch=1 at 1/1",
                           }));

  std::vector<double> child_stage_times;
  for (const auto &stage : adapter.stages) {
    if (stage.level == 1) {
      assert(stage.parent_level == 0);
      assert(std::abs(stage.parent_value - stage.stage_time) < 1.0e-13);
      child_stage_times.push_back(stage.stage_time);
    }
  }
  assert(child_stage_times ==
         std::vector<double>({0.0, 0.25, 0.5, 0.5, 0.75, 1.0}));
  assert(adapter.observers_saw_all_dense_released);
  assert(stepper.hierarchy_synchronized());
  assert(stepper.clocks()[0].accepted_steps == 1);
  assert(stepper.clocks()[1].accepted_steps == 2);
  assert(adapter.committed_level_sessions == 3);
  assert(adapter.discarded_level_sessions == 0);
}

void test_three_level_trace_retains_each_parent_until_child_catchup() {
  Fixture fixture(3);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);

  stepper.advance_one_epoch(adapter);

  assert(adapter.events == std::vector<std::string>({
                               "advance L0 [0/1,1/1] dense=required",
                               "advance L1 [0/1,1/2] dense=required",
                               "advance L2 [0/1,1/4] dense=none",
                               "advance L2 [1/4,1/2] dense=none",
                               "sync_begin L1<-L2 at 1/2",
                               "synchronize L1<-L2 at 1/2",
                               "sync_end L1<-L2 at 1/2",
                               "advance L1 [1/2,1/1] dense=required",
                               "advance L2 [1/2,3/4] dense=none",
                               "advance L2 [3/4,1/1] dense=none",
                               "sync_begin L1<-L2 at 1/1",
                               "synchronize L1<-L2 at 1/1",
                               "sync_end L1<-L2 at 1/1",
                               "sync_begin L0<-L1 at 1/1",
                               "synchronize L0<-L1 at 1/1",
                               "sync_end L0<-L1 at 1/1",
                               "sync_observers epoch=1 at 1/1",
                           }));
  assert(adapter.observers_saw_all_dense_released);
  assert(stepper.clocks()[0].time == hierarchy_time_t(1));
  assert(stepper.clocks()[1].time == hierarchy_time_t(1));
  assert(stepper.clocks()[2].time == hierarchy_time_t(1));
  assert(stepper.clocks()[0].accepted_steps == 1);
  assert(stepper.clocks()[1].accepted_steps == 2);
  assert(stepper.clocks()[2].accepted_steps == 4);
}

void test_level_session_installs_adapter_transaction_for_primary_and_auxiliary() {
  Fixture fixture(3);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);

  stepper.advance_one_epoch(adapter);

  bool saw_primary = false;
  bool saw_auxiliary = false;
  for (const auto &observation : adapter.stage_transactions) {
    assert(observation.second != nullptr);
    saw_primary |= observation.first == "primary";
    saw_auxiliary |= observation.first == "auxiliary";
  }
  assert(saw_primary);
  assert(saw_auxiliary);
  assert(CarpetX::current_scratch_state_transaction() == nullptr);
}

void test_preflight_rejects_unsupported_method_or_dense_capability() {
  Fixture fixture(2);

  auto unsupported = fixture.config;
  unsupported.levels[0].method = static_cast<SubcyclingODEMethod>(99);
  bool rejected_method = false;
  try {
    HierarchyStepper stepper(unsupported, fixture.registry);
  } catch (const std::invalid_argument &) {
    rejected_method = true;
  }
  assert(rejected_method);

  auto missing = fixture.config;
  missing.levels[0].tableau_fingerprint = fingerprint(0x7f);
  bool rejected_capability = false;
  try {
    HierarchyStepper stepper(missing, fixture.registry);
  } catch (const std::invalid_argument &) {
    rejected_capability = true;
  }
  assert(rejected_capability);

  auto missing_leaf = fixture.config;
  missing_leaf.levels[1].tableau_fingerprint = fingerprint(0x6e);
  bool rejected_leaf_capability = false;
  try {
    HierarchyStepper stepper(missing_leaf, fixture.registry);
  } catch (const std::invalid_argument &) {
    rejected_leaf_capability = true;
  }
  assert(rejected_leaf_capability);

  auto invalid_ratio = fixture.config;
  invalid_ratio.refinement_ratio = 3;
  bool rejected_ratio = false;
  try {
    HierarchyStepper stepper(invalid_ratio, fixture.registry);
  } catch (const std::invalid_argument &) {
    rejected_ratio = true;
  }
  assert(rejected_ratio);
}

void test_preflight_rejects_noncanonical_method_capability_shape() {
  const std::vector<SubcyclingODEMethod> methods{
      SubcyclingODEMethod::rk4, SubcyclingODEMethod::rkf78_order7,
      SubcyclingODEMethod::dp87_order8};

  for (std::size_t method_index = 0; method_index < methods.size();
       ++method_index) {
    const auto table =
        fingerprint(static_cast<std::uint8_t>(0x60 + method_index));
    const auto canonical = capability(methods[method_index], table);
    std::vector<DenseCapability> malformed;

    auto wrong_endpoint = canonical;
    --wrong_endpoint.endpoint_order;
    malformed.push_back(wrong_endpoint);

    auto wrong_dense_order = canonical;
    ++wrong_dense_order.dense_uniform_order;
    ++wrong_dense_order.persistent_vector_count;
    malformed.push_back(wrong_dense_order);

    auto wrong_stage_count = canonical;
    ++wrong_stage_count.stage_count;
    malformed.push_back(wrong_stage_count);

    auto wrong_rhs_cost = canonical;
    ++wrong_rhs_cost.extra_rhs_evaluations;
    malformed.push_back(wrong_rhs_cost);

    for (const auto &candidate : malformed) {
      DenseOutputRegistry registry;
      registry.register_provider(
          std::make_shared<DenseOutputProvider>(candidate));
      HierarchyStepperConfig config;
      config.levels.push_back(LevelAdvanceConfig{candidate.method, table});

      bool rejected = false;
      try {
        HierarchyStepper stepper(config, registry);
      } catch (const std::invalid_argument &) {
        rejected = true;
      }
      assert(rejected);
    }
  }
}

void test_missing_dense_fails_before_parent_clock_commit() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  adapter.omit_required_dense = true;

  bool rejected = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::logic_error &) {
    rejected = true;
  }
  assert(rejected);
  assert(stepper.clocks()[0].time == hierarchy_time_t(0));
  assert(stepper.clocks()[0].accepted_steps == 0);
  assert(adapter.committed_level_sessions == 0);
  assert(adapter.discarded_level_sessions == 1);
  assert(!stepper.hierarchy_synchronized());
}

void test_clock_overflow_fails_before_adapter_session_creation() {
  Fixture fixture(1);
  fixture.config.initial_clock = hierarchy_time_t(
      std::numeric_limits<std::int64_t>::max());
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities, fixture.config.initial_clock);

  bool rejected = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::overflow_error &) {
    rejected = true;
  } catch (...) {
  }

  assert(rejected);
  assert(adapter.begin_level_step_calls == 0);
  assert(stepper.clocks()[0].time == fixture.config.initial_clock);
  assert(stepper.clocks()[0].accepted_steps == 0);
}

void test_epoch_wide_step_counter_overflow_fails_before_adapter_call() {
  Fixture fixture(2);
  const auto maximum = std::numeric_limits<std::uint64_t>::max();
  fixture.config.initial_epoch = maximum / 2;
  fixture.config.initial_accepted_steps = {maximum / 2, maximum - 1};
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);

  bool rejected = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::overflow_error &) {
    rejected = true;
  }

  assert(rejected);
  assert(adapter.begin_level_step_calls == 0);
  assert(stepper.clocks()[0].accepted_steps == maximum / 2);
  assert(stepper.clocks()[1].accepted_steps == maximum - 1);
}

void test_child_failure_poisoning_keeps_failing_clock_uncommitted() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  adapter.throw_on_level = 1;

  bool threw = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  assert(threw);
  assert(stepper.clocks()[0].time == hierarchy_time_t(1));
  assert(stepper.clocks()[0].accepted_steps == 1);
  assert(stepper.clocks()[1].time == hierarchy_time_t(0));
  assert(stepper.clocks()[1].accepted_steps == 0);
  assert(adapter.committed_level_sessions == 1);
  assert(adapter.discarded_level_sessions == 1);
  assert(!stepper.hierarchy_synchronized());
  for (const auto &record : adapter.dense_lifetimes)
    assert(record.interval.expired());
}

void test_stop_request_is_consumed_only_after_sync_observers() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  stepper.request_stop_at_next_sync();
  adapter.sync_observer_hook = [&] { assert(stepper.stop_pending()); };

  const auto stopped = stepper.advance_one_epoch(adapter);
  assert(stopped.stop_requested);
  assert(adapter.observer_stop_requested);
  assert(!stepper.stop_pending());
  assert(adapter.events.back() == "sync_observers epoch=1 at 1/1");

  RecordingAdapter resumed_adapter(fixture.capabilities);
  const auto resumed = stepper.advance_one_epoch(resumed_adapter);
  assert(!resumed.stop_requested);
  assert(!resumed_adapter.observer_stop_requested);
  assert(resumed.synchronized_time == hierarchy_time_t(2));
  assert(resumed.epoch == 2);
}

void test_stop_request_arriving_during_observer_is_not_lost() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter first_adapter(fixture.capabilities);
  stepper.request_stop_at_next_sync();
  first_adapter.sync_observer_hook =
      [&] { stepper.request_stop_at_next_sync(); };

  const auto first = stepper.advance_one_epoch(first_adapter);
  assert(first.stop_requested);
  assert(first_adapter.observer_stop_requested);

  RecordingAdapter second_adapter(fixture.capabilities);
  const auto second = stepper.advance_one_epoch(second_adapter);
  assert(second.stop_requested);
  assert(second_adapter.observer_stop_requested);

  RecordingAdapter third_adapter(fixture.capabilities);
  const auto third = stepper.advance_one_epoch(third_adapter);
  assert(!third.stop_requested);
  assert(!third_adapter.observer_stop_requested);
}

void test_stop_request_is_not_consumed_when_observer_fails() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  adapter.throw_during_observer = true;
  stepper.request_stop_at_next_sync();

  bool threw = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::runtime_error &) {
    threw = true;
  }

  assert(threw);
  assert(adapter.observer_stop_requested);
  assert(stepper.stop_pending());
  assert(!stepper.hierarchy_synchronized());
}

void test_endpoint_roundoff_is_canonicalized_before_parent_dense_evaluation() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  adapter.round_child_endpoint_up = true;

  const auto result = stepper.advance_one_epoch(adapter);

  assert(result.synchronized_time == hierarchy_time_t(1));
  std::vector<double> child_endpoint_times;
  for (const auto &stage : adapter.stages) {
    if (stage.level == 1 &&
        stage.stage_time == static_cast<double>(stage.end_clock))
      child_endpoint_times.push_back(stage.stage_time);
  }
  assert(child_endpoint_times == std::vector<double>({0.5, 1.0}));
}

void test_synchronization_scope_ends_and_releases_dense_on_failure() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities);
  adapter.throw_during_sync = true;

  bool threw = false;
  try {
    stepper.advance_one_epoch(adapter);
  } catch (const std::runtime_error &) {
    threw = true;
  }

  assert(threw);
  assert(adapter.events.size() >= 3);
  assert(adapter.events[adapter.events.size() - 3] ==
         "sync_begin L0<-L1 at 1/1");
  assert(adapter.events[adapter.events.size() - 2] ==
         "synchronize L0<-L1 at 1/1");
  assert(adapter.events.back() == "sync_end L0<-L1 at 1/1");
  for (const auto &record : adapter.dense_lifetimes)
    assert(record.interval.expired());
}

void test_recovery_metadata_initializes_epoch_clocks_and_step_counts() {
  Fixture fixture(2);
  fixture.config.initial_clock = hierarchy_time_t(3, 2);
  fixture.config.initial_physical_time = 10.0;
  fixture.config.coarse_dt = 0.25;
  fixture.config.initial_epoch = 7;
  fixture.config.initial_accepted_steps = {7, 14};
  HierarchyStepper stepper(fixture.config, fixture.registry);
  RecordingAdapter adapter(fixture.capabilities, fixture.config.initial_clock,
                           fixture.config.initial_physical_time,
                           fixture.config.coarse_dt);

  assert(stepper.epoch() == 7);
  assert(stepper.clocks()[0].time == hierarchy_time_t(3, 2));
  assert(stepper.clocks()[1].time == hierarchy_time_t(3, 2));
  assert(stepper.clocks()[0].accepted_steps == 7);
  assert(stepper.clocks()[1].accepted_steps == 14);

  const auto result = stepper.advance_one_epoch(adapter);
  assert(result.synchronized_time == hierarchy_time_t(5, 2));
  assert(result.synchronized_physical_time == 10.25);
  assert(result.epoch == 8);
  assert(stepper.clocks()[0].accepted_steps == 8);
  assert(stepper.clocks()[1].accepted_steps == 16);
}

void test_stop_request_is_safe_from_an_asynchronous_thread() {
  Fixture fixture(2);
  HierarchyStepper stepper(fixture.config, fixture.registry);
  ConcurrentStopAdapter adapter(fixture.capabilities);
  HierarchyAdvanceResult result;
  std::exception_ptr failure;

  std::thread evolution([&] {
    try {
      result = stepper.advance_one_epoch(adapter);
    } catch (...) {
      failure = std::current_exception();
    }
  });
  while (!adapter.root_entered.load(std::memory_order_relaxed))
    std::this_thread::yield();
  stepper.request_stop_at_next_sync();
  adapter.release_root.store(true, std::memory_order_relaxed);
  evolution.join();

  assert(failure == nullptr);
  assert(result.stop_requested);
  assert(result.epoch == 1);
}

} // namespace

int main() {
  test_two_levels_use_parent_dense_at_every_exact_child_stage();
  test_three_level_trace_retains_each_parent_until_child_catchup();
  test_level_session_installs_adapter_transaction_for_primary_and_auxiliary();
  test_preflight_rejects_unsupported_method_or_dense_capability();
  test_preflight_rejects_noncanonical_method_capability_shape();
  test_missing_dense_fails_before_parent_clock_commit();
  test_clock_overflow_fails_before_adapter_session_creation();
  test_epoch_wide_step_counter_overflow_fails_before_adapter_call();
  test_child_failure_poisoning_keeps_failing_clock_uncommitted();
  test_stop_request_is_consumed_only_after_sync_observers();
  test_stop_request_arriving_during_observer_is_not_lost();
  test_stop_request_is_not_consumed_when_observer_fails();
  test_endpoint_roundoff_is_canonicalized_before_parent_dense_evaluation();
  test_synchronization_scope_ends_and_releases_dense_on_failure();
  test_recovery_metadata_initializes_epoch_clocks_and_step_counts();
  test_stop_request_is_safe_from_an_asynchronous_thread();
  std::cout << "HierarchyStepper dense recursion tests passed\n";
  return 0;
}
