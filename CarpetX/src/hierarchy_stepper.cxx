#include "hierarchy_stepper.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace CarpetX {
namespace {

bool known_method(const SubcyclingODEMethod method) noexcept {
  switch (method) {
  case SubcyclingODEMethod::rk4:
  case SubcyclingODEMethod::rkf78_order7:
  case SubcyclingODEMethod::dp87_order8:
    return true;
  }
  return false;
}

bool zero_fingerprint(const TableauFingerprint &fingerprint) noexcept {
  return std::all_of(fingerprint.begin(), fingerprint.end(),
                     [](const std::uint8_t value) { return value == 0; });
}

DenseCapability
canonical_capability(const SubcyclingODEMethod method,
                     const TableauFingerprint &tableau_fingerprint) {
  switch (method) {
  case SubcyclingODEMethod::rk4:
    return DenseCapability{method, tableau_fingerprint, 4, 4, 4, 4, 5, true,
                           true};
  case SubcyclingODEMethod::rkf78_order7:
    return DenseCapability{method, tableau_fingerprint, 7, 7, 11, 23, 8, true,
                           true};
  case SubcyclingODEMethod::dp87_order8:
    return DenseCapability{method, tableau_fingerprint, 8, 8, 13, 39, 9, true,
                           true};
  }
  throw std::invalid_argument("unsupported hierarchy ODE method");
}

bool same_capability(const DenseCapability &left,
                     const DenseCapability &right) noexcept {
  return left.method == right.method &&
         left.tableau_fingerprint == right.tableau_fingerprint &&
         left.endpoint_order == right.endpoint_order &&
         left.dense_uniform_order == right.dense_uniform_order &&
         left.stage_count == right.stage_count &&
         left.extra_rhs_evaluations == right.extra_rhs_evaluations &&
         left.persistent_vector_count == right.persistent_vector_count &&
         left.arbitrary_theta == right.arbitrary_theta &&
         left.verified == right.verified;
}

std::int64_t checked_multiply_by_positive(const std::int64_t value,
                                          const std::int64_t factor) {
  if (factor <= 0)
    throw std::logic_error("clock denominator factor must be positive");
  if (value > std::numeric_limits<std::int64_t>::max() / factor ||
      value < std::numeric_limits<std::int64_t>::min() / factor)
    throw std::overflow_error("exact hierarchy clock multiplication overflow");
  return value * factor;
}

std::int64_t checked_add_integer(const std::int64_t left,
                                 const std::int64_t right) {
  if ((right > 0 &&
       left > std::numeric_limits<std::int64_t>::max() - right) ||
      (right < 0 &&
       left < std::numeric_limits<std::int64_t>::min() - right))
    throw std::overflow_error("exact hierarchy clock addition overflow");
  return left + right;
}

hierarchy_time_t checked_clock_add(const hierarchy_time_t left,
                                   const hierarchy_time_t right) {
  if (left.den <= 0 || right.den <= 0)
    throw std::invalid_argument(
        "exact hierarchy clocks require positive denominators");
  const auto common = std::gcd(left.den, right.den);
  const auto left_factor = right.den / common;
  const auto right_factor = left.den / common;
  const auto left_term = checked_multiply_by_positive(left.num, left_factor);
  const auto right_term = checked_multiply_by_positive(right.num, right_factor);
  const auto numerator = checked_add_integer(left_term, right_term);
  const auto denominator =
      checked_multiply_by_positive(left.den, left_factor);
  return hierarchy_time_t(numerator, denominator);
}

class HierarchyStagePreparer final : public StagePreparer {
public:
  HierarchyStagePreparer(HierarchyEvolutionAdapter &adapter,
                         const StepContext &expected_context,
                         std::shared_ptr<const DenseInterval> parent_dense)
      : adapter_(adapter), expected_context_(expected_context),
        parent_dense_(std::move(parent_dense)) {
    if (expected_context_.level == 0 && parent_dense_ != nullptr)
      throw std::logic_error("root level cannot have parent dense state");
    if (expected_context_.level > 0 && parent_dense_ == nullptr)
      throw std::logic_error("fine level requires parent dense state");
  }

  void prepare_stage(const StepContext &context,
                     const StagePoint &stage_point) override {
    if (&context != &expected_context_)
      throw std::logic_error("stage context does not belong to this level step");

    const auto exact_stage_clock =
        context.begin_clock +
        stage_point.stage_fraction * (context.end_clock - context.begin_clock);
    StagePoint canonical_stage_point = stage_point;
    canonical_stage_point.stage_time =
        context.begin_time + static_cast<double>(stage_point.stage_fraction) *
                                 (context.end_time - context.begin_time);
    const double time_scale =
        std::max({1.0, std::abs(context.begin_time),
                  std::abs(context.end_time),
                  std::abs(canonical_stage_point.stage_time)});
    const double tolerance =
        16.0 * std::numeric_limits<double>::epsilon() * time_scale;

    if (parent_dense_ != nullptr) {
      const auto &id = parent_dense_->id();
      if (id.level != context.level - 1 ||
          id.begin_clock > context.begin_clock ||
          id.end_clock < context.end_clock)
        throw std::logic_error(
            "parent dense interval does not cover the child step exactly");
      if (exact_stage_clock < id.begin_clock ||
          exact_stage_clock > id.end_clock ||
          canonical_stage_point.stage_time < id.begin_time - tolerance ||
          canonical_stage_point.stage_time > id.end_time + tolerance)
        throw std::out_of_range(
            "parent dense interval does not cover the requested stage time");
    }

    adapter_.prepare_stage(context, canonical_stage_point,
                           parent_dense_.get());
  }

private:
  HierarchyEvolutionAdapter &adapter_;
  const StepContext &expected_context_;
  std::shared_ptr<const DenseInterval> parent_dense_;
};

class SynchronizationScope final {
public:
  SynchronizationScope(HierarchyEvolutionAdapter &adapter,
                       const int coarse_level, const int fine_level,
                       const hierarchy_time_t time)
      : adapter_(adapter), coarse_level_(coarse_level),
        fine_level_(fine_level), time_(time) {
    adapter_.begin_synchronization(coarse_level_, fine_level_, time_);
  }

  ~SynchronizationScope() {
    adapter_.end_synchronization(coarse_level_, fine_level_, time_);
  }

  SynchronizationScope(const SynchronizationScope &) = delete;
  SynchronizationScope &operator=(const SynchronizationScope &) = delete;

private:
  HierarchyEvolutionAdapter &adapter_;
  int coarse_level_;
  int fine_level_;
  hierarchy_time_t time_;
};

} // namespace

HierarchyStepper::HierarchyStepper(
    HierarchyStepperConfig config,
    const DenseOutputRegistry &dense_registry)
    : level_configs_(std::move(config.levels)),
      initial_clock_(config.initial_clock),
      initial_physical_time_(config.initial_physical_time),
      coarse_dt_(config.coarse_dt), epoch_(config.initial_epoch) {
  if (level_configs_.empty())
    throw std::invalid_argument("hierarchy must contain at least one level");
  if (config.refinement_ratio != 2)
    throw std::invalid_argument("only factor-two time refinement is supported");
  if (level_configs_.size() > 62)
    throw std::invalid_argument("hierarchy depth exceeds exact clock range");
  if (!std::isfinite(initial_physical_time_) || !std::isfinite(coarse_dt_) ||
      !(coarse_dt_ > 0.0))
    throw std::invalid_argument(
        "initial physical time and coarse dt must be finite, with dt positive");
  if (config.initial_accepted_steps.empty()) {
    if (config.initial_epoch != 0)
      throw std::invalid_argument(
          "a recovered hierarchy epoch requires accepted-step counts");
  } else {
    if (config.initial_accepted_steps.size() != level_configs_.size())
      throw std::invalid_argument(
          "accepted-step count must be supplied for every level");
    if (config.initial_accepted_steps.front() != config.initial_epoch)
      throw std::invalid_argument(
          "coarse accepted-step count must equal the hierarchy epoch");
    for (std::size_t level = 1; level < level_configs_.size(); ++level) {
      const auto coarse_steps = config.initial_accepted_steps[level - 1];
      if (coarse_steps > std::numeric_limits<std::uint64_t>::max() / 2 ||
          config.initial_accepted_steps[level] != 2 * coarse_steps)
        throw std::invalid_argument(
            "recovered accepted-step counts must follow factor-two refinement");
    }
  }

  expected_dense_capabilities_.reserve(level_configs_.size());
  active_dense_.resize(level_configs_.size());
  clocks_.reserve(level_configs_.size());

  hierarchy_time_t level_dt{1};
  for (std::size_t level = 0; level < level_configs_.size(); ++level) {
    const auto &level_config = level_configs_[level];
    if (!known_method(level_config.method))
      throw std::invalid_argument("hierarchy level uses an unsupported method");
    if (zero_fingerprint(level_config.tableau_fingerprint))
      throw std::invalid_argument(
          "hierarchy level tableau fingerprint must be non-zero");

    try {
      const auto provider = dense_registry.require(
          level_config.method, level_config.tableau_fingerprint);
      const auto expected = canonical_capability(
          level_config.method, level_config.tableau_fingerprint);
      if (!same_capability(provider->capability(), expected))
        throw std::invalid_argument(
            "dense capability does not match the fixed method contract");
      expected_dense_capabilities_.push_back(expected);
    } catch (const std::out_of_range &) {
      throw std::invalid_argument(
          "no verified dense capability matches a hierarchy level");
    }

    const auto accepted_steps = config.initial_accepted_steps.empty()
                                    ? std::uint64_t{0}
                                    : config.initial_accepted_steps[level];
    clocks_.push_back(LevelClock{initial_clock_, level_dt, accepted_steps});
    level_dt /= 2;
  }
}

HierarchyAdvanceResult
HierarchyStepper::advance_one_epoch(HierarchyEvolutionAdapter &adapter) {
  if (in_progress_ || faulted_)
    throw std::logic_error("cannot advance an active or faulted hierarchy");
  if (!clocks_aligned() || !no_dense_in_flight())
    throw std::logic_error("cannot advance a desynchronized hierarchy");
  if (epoch_ == std::numeric_limits<std::uint64_t>::max())
    throw std::overflow_error("hierarchy epoch counter exhausted");
  preflight_epoch_capacity();

  in_progress_ = true;
  try {
    advance_level(0, adapter);
    if (!clocks_aligned() || !no_dense_in_flight())
      throw std::logic_error("hierarchy did not synchronize after an epoch");

    const auto synchronized_time = clocks_.front().time;
    const auto synchronized_physical_time = physical_time(synchronized_time);
    const auto completed_epoch = epoch_ + 1;
    const auto stop_requested =
        stop_requested_.exchange(false, std::memory_order_acq_rel);
    stop_snapshot_in_flight_.store(stop_requested, std::memory_order_release);
    try {
      adapter.run_sync_observers(synchronized_time,
                                 synchronized_physical_time, completed_epoch,
                                 stop_requested);
    } catch (...) {
      if (stop_requested)
        stop_requested_.store(true, std::memory_order_release);
      stop_snapshot_in_flight_.store(false, std::memory_order_release);
      throw;
    }
    stop_snapshot_in_flight_.store(false, std::memory_order_release);
    epoch_ = completed_epoch;

    in_progress_ = false;
    return HierarchyAdvanceResult{synchronized_time,
                                  synchronized_physical_time, epoch_,
                                  stop_requested};
  } catch (...) {
    for (auto &interval : active_dense_)
      interval.reset();
    faulted_ = true;
    in_progress_ = false;
    throw;
  }
}

void HierarchyStepper::request_stop_at_next_sync() noexcept {
  stop_requested_.store(true, std::memory_order_release);
}

bool HierarchyStepper::stop_pending() const noexcept {
  return stop_requested_.load(std::memory_order_acquire) ||
         stop_snapshot_in_flight_.load(std::memory_order_acquire);
}

bool HierarchyStepper::hierarchy_synchronized() const noexcept {
  return !in_progress_ && !faulted_ && clocks_aligned() &&
         no_dense_in_flight();
}

std::uint64_t HierarchyStepper::epoch() const noexcept { return epoch_; }

const std::vector<LevelClock> &HierarchyStepper::clocks() const noexcept {
  return clocks_;
}

bool HierarchyStepper::clocks_aligned() const noexcept {
  if (clocks_.empty())
    return false;
  const auto synchronized_time = clocks_.front().time;
  return std::all_of(clocks_.begin(), clocks_.end(),
                     [&](const LevelClock &clock) {
                       return clock.time == synchronized_time;
                     });
}

bool HierarchyStepper::no_dense_in_flight() const noexcept {
  return std::all_of(active_dense_.begin(), active_dense_.end(),
                     [](const auto &interval) { return interval == nullptr; });
}

void HierarchyStepper::preflight_epoch_capacity() const {
  std::uint64_t accepted_step_increment = 1;
  const auto coarse_dt = clocks_.front().dt;
  for (std::size_t level = 0; level < clocks_.size(); ++level) {
    const auto &clock = clocks_[level];
    if (clock.accepted_steps >
        std::numeric_limits<std::uint64_t>::max() - accepted_step_increment)
      throw std::overflow_error(
          "level accepted-step counter cannot cover the next epoch");

    static_cast<void>(checked_clock_add(clock.time, clock.dt));
    static_cast<void>(checked_clock_add(clock.time, coarse_dt));

    if (level + 1 < clocks_.size()) {
      if (accepted_step_increment >
          std::numeric_limits<std::uint64_t>::max() / 2)
        throw std::overflow_error(
            "hierarchy depth exceeds accepted-step counter capacity");
      accepted_step_increment *= 2;
    }
  }
}

double HierarchyStepper::physical_time(const hierarchy_time_t clock) const {
  const double result =
      initial_physical_time_ +
      static_cast<double>(clock - initial_clock_) * coarse_dt_;
  if (!std::isfinite(result))
    throw std::overflow_error("hierarchy physical time is not finite");
  return result;
}

void HierarchyStepper::validate_dense_interval(
    const int level, const StepContext &context,
    const DenseInterval &interval) const {
  const auto index = static_cast<std::size_t>(level);
  const auto &expected_capability = expected_dense_capabilities_[index];
  const auto &id = interval.id();
  if (id.level != level || id.begin_clock != context.begin_clock ||
      id.end_clock != context.end_clock ||
      id.begin_time != context.begin_time || id.end_time != context.end_time ||
      id.method != context.method ||
      id.tableau_fingerprint !=
          level_configs_[index].tableau_fingerprint)
    throw std::logic_error(
        "published dense interval does not match the accepted level step");
  if (!same_capability(interval.capability(), expected_capability))
    throw std::logic_error(
        "published dense interval capability differs from preflight");
}

void HierarchyStepper::advance_level(const int level,
                                     HierarchyEvolutionAdapter &adapter) {
  if (level < 0 || static_cast<std::size_t>(level) >= clocks_.size())
    throw std::logic_error("invalid hierarchy level");

  const auto index = static_cast<std::size_t>(level);
  auto &clock = clocks_[index];
  const auto begin = clock.time;
  const auto end = checked_clock_add(begin, clock.dt);
  const bool require_dense = index + 1 < clocks_.size();

  if (require_dense && clocks_[index + 1].time != begin)
    throw std::logic_error("child clock is not aligned with its parent");

  std::shared_ptr<const DenseInterval> parent_dense;
  if (level > 0) {
    parent_dense = active_dense_[index - 1];
    if (parent_dense == nullptr)
      throw std::logic_error("fine level has no retained parent dense interval");
  }

  const StepContext context{level,
                            begin,
                            end,
                            physical_time(begin),
                            physical_time(end),
                            level_configs_[index].method,
                            require_dense,
                            clock.accepted_steps + 1};
  HierarchyStagePreparer stage_preparer(adapter, context,
                                        std::move(parent_dense));

  auto session = adapter.begin_level_step(context, require_dense);
  auto *const transaction =
      session == nullptr ? nullptr : session->transaction();
  if (transaction == nullptr)
    throw std::logic_error(
        "level-step session must own a scratch-state transaction");

  LevelAdvanceResult result;
  {
    ScopedStepContext scope(context, stage_preparer, transaction);
    result = session->advance();
  }

  if (require_dense) {
    if (result.dense_interval == nullptr)
      throw std::logic_error("parent level did not publish dense state");
    validate_dense_interval(level, context, *result.dense_interval);
  } else if (result.dense_interval != nullptr) {
    throw std::logic_error("leaf level unexpectedly published dense state");
  }
  session->commit();

  clock.time = end;
  if (clock.accepted_steps == std::numeric_limits<std::uint64_t>::max())
    throw std::overflow_error("level accepted-step counter exhausted");
  ++clock.accepted_steps;
  session.reset();

  if (!require_dense)
    return;

  active_dense_[index] = std::move(result.dense_interval);
  try {
    advance_level(level + 1, adapter);
    advance_level(level + 1, adapter);
    if (clocks_[index + 1].time != end)
      throw std::logic_error("child clock did not catch up with its parent");

    const SynchronizationScope synchronization_scope(adapter, level, level + 1,
                                                     end);
    adapter.synchronize_levels(level, level + 1, end);
  } catch (...) {
    active_dense_[index].reset();
    throw;
  }
  active_dense_[index].reset();
}

} // namespace CarpetX
