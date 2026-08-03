#include "explicit_rk.hxx"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace allocation_probe {
inline bool enabled{false};
inline std::size_t count{0};
}

void *operator new(const std::size_t size) {
  if (allocation_probe::enabled)
    ++allocation_probe::count;
  if (void *const memory = std::malloc(size))
    return memory;
  throw std::bad_alloc();
}

void *operator new[](const std::size_t size) { return ::operator new(size); }
void operator delete(void *const memory) noexcept { std::free(memory); }
void operator delete[](void *const memory) noexcept { std::free(memory); }
void operator delete(void *const memory, std::size_t) noexcept {
  std::free(memory);
}
void operator delete[](void *const memory, std::size_t) noexcept {
  std::free(memory);
}

namespace {

using ODESolvers::ExplicitRKMethod;
using ODESolvers::ExplicitRKTableau;
using ODESolvers::InitialRHSMode;
using ODESolvers::LoadedRHSProvenance;

std::uint64_t bits(const double value) {
  std::uint64_t result = 0;
  static_assert(sizeof result == sizeof value);
  std::memcpy(&result, &value, sizeof result);
  return result;
}

std::string encoded(const double value) {
  std::ostringstream stream;
  stream << std::hex << bits(value);
  return stream.str();
}

struct VectorState {
  std::vector<double> values;
  int identity{-1};
};

struct RecordedStagePoint {
  int kind;
  int stage_index;
  int stage_count;
  std::int64_t numerator;
  std::int64_t denominator;
};

struct RecordingOps {
  using scalar_type = double;
  using state_type = VectorState;

  VectorState state_value;
  VectorState rhs_value;
  double time{0.0};
  int next_identity{10};
  int rhs_calls{0};
  int validations{0};
  int mutations{0};
  std::uint64_t generation{0};
  bool loaded_rhs{false};
  std::vector<std::string> trace;
  std::vector<RecordedStagePoint> stage_points;

  explicit RecordingOps(std::vector<double> values)
      : state_value{std::move(values), 1},
        rhs_value{std::vector<double>(state_value.values.size()), 2} {}

  const VectorState &state() const noexcept { return state_value; }
  const VectorState &rhs() const noexcept { return rhs_value; }

  template <class StagePoint> void set_stage_point(const StagePoint &point) {
    stage_points.push_back(
        {static_cast<int>(point.kind), point.stage_index, point.stage_count,
         point.parent_fraction.numerator, point.parent_fraction.denominator});
  }

  VectorState snapshot_state() {
    trace.push_back("copy-y:" + std::to_string(next_identity));
    return VectorState{state_value.values, next_identity++};
  }

  VectorState snapshot_rhs() {
    trace.push_back("copy-f:" + std::to_string(next_identity));
    return VectorState{rhs_value.values, next_identity++};
  }

  void prepare_initial(const double stage_time) {
    ++mutations;
    time = stage_time;
    trace.push_back("prepare-initial:" + encoded(stage_time));
  }

  void update_state(const int update_index, const double stage_time,
                    const double destination_scale,
                    const ODESolvers::LinearCombinationView<double, VectorState>
                        combination) {
    ++mutations;
    ++generation;
    loaded_rhs = false;
    std::ostringstream entry;
    entry << "update:" << update_index << ':' << encoded(stage_time) << ':'
          << encoded(destination_scale);
    std::vector<double> updated(state_value.values.size(), 0.0);
    for (std::size_t component = 0; component < updated.size(); ++component)
      updated[component] = destination_scale * state_value.values[component];
    for (std::size_t source = 0; source < combination.size; ++source) {
      entry << ':' << encoded(combination.factors[source]) << '@'
            << combination.sources[source]->identity;
      for (std::size_t component = 0; component < updated.size(); ++component)
        updated[component] +=
            combination.factors[source] *
            combination.sources[source]->values[component];
    }
    state_value.values = std::move(updated);
    time = stage_time;
    trace.push_back(entry.str());
    trace.push_back("prepare-stage:" + encoded(stage_time));
    trace.push_back("poststep:" + encoded(stage_time));
  }

  void update_state(const int update_index, const double stage_time,
                    const double destination_scale,
                    const std::vector<double> &factors,
                    const std::vector<const VectorState *> &sources) {
    assert(factors.size() == sources.size());
    update_state(update_index, stage_time, destination_scale,
                 {factors.data(), sources.data(), factors.size()});
  }

  void evaluate_rhs(const int stage_index) {
    ++rhs_calls;
    trace.push_back("rhs:" + std::to_string(stage_index) + ':' +
                    encoded(time));
    for (std::size_t component = 0; component < state_value.values.size();
         ++component) {
      const double y = state_value.values[component];
      rhs_value.values[component] =
          static_cast<double>(component + 1) * std::sin(y) +
          0.125 * time * y * y + 0.01 * static_cast<double>(component + 1);
    }
  }

  void validate_rhs(const int stage_index) {
    ++validations;
    trace.push_back("valid:" + std::to_string(stage_index));
    for (const double value : rhs_value.values)
      if (!std::isfinite(value))
        throw std::runtime_error("non-finite RHS");
  }

  void accumulate_rk4(VectorState &accumulator, const double factor,
                      const VectorState &increment) {
    ++mutations;
    trace.push_back("rk4-accum:" + encoded(factor) + '@' +
                    std::to_string(increment.identity));
    for (std::size_t component = 0; component < accumulator.values.size();
         ++component)
      accumulator.values[component] += factor * increment.values[component];
  }

  void preload_exact_left_rhs() {
    const int saved_calls = rhs_calls;
    const auto saved_trace = trace;
    evaluate_rhs(1);
    rhs_calls = saved_calls;
    trace = saved_trace;
    loaded_rhs = true;
  }

  std::uint64_t state_generation() const noexcept { return generation; }
  LoadedRHSProvenance<double>
  loaded_rhs_provenance(const double left_time) const {
    if (!loaded_rhs || time != left_time)
      throw std::logic_error("no exact loaded left RHS");
    return {generation, left_time};
  }
  void validate_loaded_rhs_provenance(
      const LoadedRHSProvenance<double> &provenance) const {
    if (!loaded_rhs || provenance.state_generation != generation ||
        provenance.left_time != time)
      throw std::logic_error("loaded RHS provenance is stale");
  }
  void consume_loaded_rhs(
      const LoadedRHSProvenance<double> &provenance) {
    validate_loaded_rhs_provenance(provenance);
    loaded_rhs = false;
  }
  void invalidate_loaded_rhs_generation() noexcept { ++generation; }
};

struct RecordingObserver {
  std::vector<std::string> events;

  void initial_state(const double time, const VectorState &) {
    events.push_back("initial-state:" + encoded(time));
  }
  void initial_rhs(const double time, const VectorState &) {
    events.push_back("initial-rhs:" + encoded(time));
  }
  void stage_rhs(const int stage, const double time, const VectorState &) {
    events.push_back("stage-rhs:" + std::to_string(stage) + ':' +
                     encoded(time));
  }
  void accepted_endpoint(const double time, const VectorState &) {
    events.push_back("accepted:" + encoded(time));
  }
};

struct LegacyTableau {
  std::vector<std::vector<double>> a;
  std::vector<double> b;
  std::vector<double> c;
};

LegacyTableau legacy_tableau(const ExplicitRKMethod method) {
  const auto R = [](const double n, const double d) { return n / d; };
  if (method == ExplicitRKMethod::rkf78) {
    return {{{},
             {R(2, 27)},
             {R(1, 36), R(3, 36)},
             {R(1, 24), 0, R(3, 24)},
             {R(20, 48), 0, R(-75, 48), R(75, 48)},
             {R(1, 20), 0, 0, R(5, 20), R(4, 20)},
             {R(-25, 108), 0, 0, R(125, 108), R(-260, 108), R(250, 108)},
             {R(31, 300), 0, 0, 0, R(61, 225), R(-2, 9), R(13, 900)},
             {2, 0, 0, R(-53, 6), R(704, 45), R(-107, 9), R(67, 90), 3},
             {R(-91, 108), 0, 0, R(23, 108), R(-976, 135), R(311, 54),
              R(-19, 60), R(17, 6), R(-1, 12)},
             {R(2383, 4100), 0, 0, R(-341, 164), R(4496, 1025), R(-301, 82),
              R(2133, 4100), R(45, 82), R(45, 164), R(18, 41)}},
            {R(41, 840), 0, 0, 0, 0, R(34, 105), R(9, 35), R(9, 35),
             R(9, 280), R(9, 280), R(41, 840)},
            {0, R(2, 27), R(1, 9), R(1, 6), R(5, 12), R(1, 2), R(5, 6),
             R(1, 6), R(2, 3), R(1, 3), 1}};
  }
  if (method == ExplicitRKMethod::dp87) {
    return {{{},
             {R(1, 18)},
             {R(1, 48), R(1, 16)},
             {R(1, 32), 0, R(3, 32)},
             {R(5, 16), 0, -R(75, 64), R(75, 64)},
             {R(3, 80), 0, 0, R(3, 16), R(3, 20)},
             {R(29443841, 614563906), 0, 0, R(77736538, 692538347),
              -R(28693883, 1125000000), R(23124283, 1800000000)},
             {R(16016141, 946692911), 0, 0, R(61564180, 158732637),
              R(22789713, 633445777), R(545815736, 2771057229),
              -R(180193667, 1043307555)},
             {R(39632708, 573591083), 0, 0, -R(433636366, 683701615),
              -R(421739975, 2616292301), R(100302831, 723423059),
              R(790204164, 839813087), R(800635310, 3783071287)},
             {R(246121993, 1340847787), 0, 0, -R(37695042795, 15268766246),
              -R(309121744, 1061227803), -R(12992083, 490766935),
              R(6005943493, 2108947869), R(393006217, 1396673457),
              R(123872331, 1001029789)},
             {-R(1028468189, 846180014), 0, 0, R(8478235783, 508512852),
              R(1311729495, 1432422823), -R(10304129995, 1701304382),
              -R(48777925059, 3047939560), R(15336726248, 1032824649),
              -R(45442868181, 3398467696), R(3065993473, 597172653)},
             {R(185892177, 718116043), 0, 0, -R(3185094517, 667107341),
              -R(477755414, 1098053517), -R(703635378, 230739211),
              R(5731566787, 1027545527), R(5232866602, 850066563),
              -R(4093664535, 808688257), R(3962137247, 1805957418),
              R(65686358, 487910083)},
             {R(403863854, 491063109), 0, 0, -R(5068492393, 434740067),
              -R(411421997, 543043805), R(652783627, 914296604),
              R(11173962825, 925320556), -R(13158990841, 6184727034),
              R(3936647629, 1978049680), -R(160528059, 685178525),
              R(248638103, 1413531060), 0}},
            {R(14005451, 335480064), 0, 0, 0, 0,
             -R(59238493, 1068277825), R(181606767, 758867731),
             R(561292985, 797845732), -R(1041891430, 1371343529),
             R(760417239, 1151165299), R(118820643, 751138087),
             -R(528747749, 2220607170), R(1, 4)},
            {}};
  }
  throw std::invalid_argument("legacy table requested for non-table method");
}

void legacy_advance(const ExplicitRKMethod method, const double begin,
                    const double dt, RecordingOps &ops,
                    RecordingObserver &observer) {
  ops.prepare_initial(begin);
  const auto old = ops.snapshot_state();
  observer.initial_state(begin, ops.state());
  auto rhs_stage = [&](const int stage, const double time) {
    ops.evaluate_rhs(stage);
    ops.validate_rhs(stage);
    if (stage == 1)
      observer.initial_rhs(time, ops.rhs());
    observer.stage_rhs(stage, time, ops.rhs());
  };
  if (method == ExplicitRKMethod::rk4) {
    rhs_stage(1, begin);
    auto accum = ops.snapshot_rhs();
    ops.update_state(1, begin + dt / 2, 1.0, {dt / 2}, {&accum});
    rhs_stage(2, begin + dt / 2);
    ops.accumulate_rk4(accum, 2.0, ops.rhs());
    ops.update_state(2, begin + dt / 2, 0.0, {1.0, dt / 2},
                     {&old, &ops.rhs()});
    rhs_stage(3, begin + dt / 2);
    ops.accumulate_rk4(accum, 2.0, ops.rhs());
    ops.update_state(3, begin + dt, 0.0, {1.0, dt}, {&old, &ops.rhs()});
    rhs_stage(4, begin + dt);
    ops.update_state(4, begin + dt, 0.0, {1.0, dt / 6, dt / 6},
                     {&old, &accum, &ops.rhs()});
  } else {
    const auto table = legacy_tableau(method);
    std::vector<VectorState> stages;
    stages.reserve(table.a.size());
    for (std::size_t step = 0; step < table.a.size(); ++step) {
      double c = 0.0;
      if (method == ExplicitRKMethod::rkf78) {
        c = table.c[step];
      } else {
        for (const double coefficient : table.a[step])
          c += coefficient;
      }
      if (step > 0) {
        std::vector<double> factors{1.0};
        std::vector<const VectorState *> sources{&old};
        for (std::size_t i = 0; i < table.a[step].size(); ++i)
          if (table.a[step][i] != 0.0) {
            factors.push_back(table.a[step][i] * dt);
            sources.push_back(&stages[i]);
          }
        ops.update_state(static_cast<int>(step), begin + c * dt, 0.0,
                         factors, sources);
      }
      rhs_stage(static_cast<int>(step + 1), begin + c * dt);
      stages.push_back(ops.snapshot_rhs());
    }
    std::vector<double> factors{1.0};
    std::vector<const VectorState *> sources{&old};
    for (std::size_t i = 0; i < table.b.size(); ++i)
      if (table.b[i] != 0.0) {
        factors.push_back(table.b[i] * dt);
        sources.push_back(&stages[i]);
      }
    ops.update_state(static_cast<int>(table.a.size()), begin + dt, 0.0,
                     factors, sources);
  }
  observer.accepted_endpoint(begin + dt, ops.state());
}

void assert_same_bits(const std::vector<double> &left,
                      const std::vector<double> &right) {
  assert(left.size() == right.size());
  for (std::size_t i = 0; i < left.size(); ++i)
    assert(bits(left[i]) == bits(right[i]));
}

void test_matches_frozen_legacy_arithmetic_and_complete_trace() {
  for (const auto method : {ExplicitRKMethod::rk4, ExplicitRKMethod::rkf78,
                            ExplicitRKMethod::dp87}) {
    for (const auto &initial : {std::vector<double>{0.1},
                                std::vector<double>{-0.4, 0.25},
                                std::vector<double>{0.9, -0.2, 0.03}}) {
      RecordingOps expected(initial);
      RecordingObserver expected_observer;
      legacy_advance(method, 0.375, 0.0625, expected, expected_observer);

      RecordingOps actual(initial);
      RecordingObserver actual_observer;
      ODESolvers::advance_explicit_rk(method, 0.375, 0.0625,
                                      InitialRHSMode::calculate, actual,
                                      actual_observer);

      assert_same_bits(actual.state().values, expected.state().values);
      assert_same_bits(actual.rhs().values, expected.rhs().values);
      assert(actual.trace == expected.trace);
      assert(actual_observer.events == expected_observer.events);
    }
  }
}

void test_stage_counts_ordering_endpoint_and_rhs_reuse() {
  for (const auto method : {ExplicitRKMethod::rk4, ExplicitRKMethod::rkf78,
                            ExplicitRKMethod::dp87}) {
    const int stages = method == ExplicitRKMethod::rk4
                           ? 4
                           : method == ExplicitRKMethod::rkf78 ? 11 : 13;
    RecordingOps calculated({0.2, -0.1});
    RecordingObserver calculated_observer;
    ODESolvers::advance_explicit_rk(method, 0.0, 0.125,
                                    InitialRHSMode::calculate, calculated,
                                    calculated_observer);
    assert(calculated.rhs_calls == stages);
    assert(calculated.validations == stages);
    assert(calculated_observer.events.front().find("initial-state:") == 0);
    assert(calculated_observer.events.back().find("accepted:") == 0);
    assert(calculated_observer.events[calculated_observer.events.size() - 2]
               .find("stage-rhs:") == 0);

    RecordingOps reused({0.2, -0.1});
    reused.prepare_initial(0.0);
    reused.preload_exact_left_rhs();
    auto token = ODESolvers::make_loaded_rhs_token(reused, 0.0);
    reused.trace.clear();
    reused.mutations = 0;
    RecordingObserver reused_observer;
    ODESolvers::advance_explicit_rk(method, 0.0, 0.125,
                                    InitialRHSMode::reuse_loaded, reused,
                                    reused_observer, token);
    assert(reused.rhs_calls == stages - 1);
    assert(reused.validations == stages);
    assert_same_bits(reused.state().values, calculated.state().values);
    assert_same_bits(reused.rhs().values, calculated.rhs().values);
  }
}

struct ScalarState {
  double value{0.0};
};

struct ExponentialOps {
  using scalar_type = double;
  using state_type = ScalarState;
  ScalarState y;
  ScalarState f;
  double time{0.0};
  int mutations{0};
  std::uint64_t generation{0};

  const ScalarState &state() const noexcept { return y; }
  const ScalarState &rhs() const noexcept { return f; }
  ScalarState snapshot_state() { return y; }
  ScalarState snapshot_rhs() { return f; }
  void prepare_initial(const double t) {
    ++mutations;
    time = t;
  }
  void update_state(const int, const double t, const double scale,
                    const ODESolvers::LinearCombinationView<double, ScalarState>
                        combination) {
    ++mutations;
    ++generation;
    double next = scale * y.value;
    for (std::size_t i = 0; i < combination.size; ++i)
      next += combination.factors[i] * combination.sources[i]->value;
    y.value = next;
    time = t;
  }
  void evaluate_rhs(const int) { f.value = y.value; }
  void validate_rhs(const int) {
    if (!std::isfinite(f.value))
      throw std::runtime_error("bad RHS");
  }
  void accumulate_rk4(ScalarState &accumulator, const double factor,
                      const ScalarState &increment) {
    ++mutations;
    accumulator.value += factor * increment.value;
  }
  std::uint64_t state_generation() const noexcept { return generation; }
  LoadedRHSProvenance<double>
  loaded_rhs_provenance(const double) const {
    throw std::logic_error("ExponentialOps has no loaded RHS");
  }
  void validate_loaded_rhs_provenance(
      const LoadedRHSProvenance<double> &) const {
    throw std::logic_error("ExponentialOps has no loaded RHS");
  }
  void consume_loaded_rhs(const LoadedRHSProvenance<double> &) {
    throw std::logic_error("ExponentialOps has no loaded RHS");
  }
};

// The ordinary calculate path must retain the original Ops contract. In
// particular, it must not instantiate any loaded-RHS provenance API.
struct CalculateOnlyOps {
  using scalar_type = double;
  using state_type = ScalarState;

  ScalarState y{1.0};
  ScalarState f{0.0};
  double time{0.0};

  const ScalarState &state() const noexcept { return y; }
  const ScalarState &rhs() const noexcept { return f; }
  ScalarState snapshot_state() const { return y; }
  ScalarState snapshot_rhs() const { return f; }
  void prepare_initial(const double t) { time = t; }
  void update_state(const int, const double t, const double scale,
                    const ODESolvers::LinearCombinationView<double, ScalarState>
                        combination) {
    double next = scale * y.value;
    for (std::size_t i = 0; i < combination.size; ++i)
      next += combination.factors[i] * combination.sources[i]->value;
    y.value = next;
    time = t;
  }
  void evaluate_rhs(const int) { f.value = y.value + time; }
  void validate_rhs(const int) const {
    if (!std::isfinite(f.value))
      throw std::runtime_error("bad RHS");
  }
  void accumulate_rk4(ScalarState &accumulator, const double factor,
                      const ScalarState &increment) const {
    accumulator.value += factor * increment.value;
  }
};

void test_calculate_only_ops_does_not_require_reuse_api() {
  CalculateOnlyOps ops;
  ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                  InitialRHSMode::calculate, ops);
  assert(std::isfinite(ops.y.value));
}

double integrate_exp(const ExplicitRKMethod method, const int steps) {
  ExponentialOps ops{{1.0}, {0.0}, 0.0, 0, 0};
  for (int step = 0; step < steps; ++step)
    ODESolvers::advance_explicit_rk(
        method, static_cast<double>(step) / steps, 1.0 / steps,
        InitialRHSMode::calculate, ops);
  return ops.y.value;
}

void test_observed_fixed_interval_endpoint_orders() {
  for (const auto &item : {std::pair{ExplicitRKMethod::rk4, 4},
                           std::pair{ExplicitRKMethod::rkf78, 7},
                           std::pair{ExplicitRKMethod::dp87, 8}}) {
    const double coarse = std::abs(integrate_exp(item.first, 2) - std::exp(1.0));
    const double fine = std::abs(integrate_exp(item.first, 4) - std::exp(1.0));
    const double observed = std::log2(coarse / fine);
    assert(observed > static_cast<double>(item.second) - 0.65);
  }
}

template <typename Function> void expect_failure(Function &&function);

void test_reuse_loaded_requires_fresh_one_shot_provenance() {
  {
    RecordingOps ops({0.2});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                      InitialRHSMode::reuse_loaded, ops);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    RecordingOps ops({0.2});
    expect_failure([&] { (void)ODESolvers::make_loaded_rhs_token(ops, 0.0); });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    RecordingOps ops({0.2});
    ops.prepare_initial(0.0);
    ops.preload_exact_left_rhs();
    auto token = ODESolvers::make_loaded_rhs_token(ops, 0.0);
    ops.invalidate_loaded_rhs_generation();
    ops.mutations = 0;
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                      InitialRHSMode::reuse_loaded, ops,
                                      token);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    RecordingOps ops({0.2});
    ops.prepare_initial(0.0);
    ops.preload_exact_left_rhs();
    auto token = ODESolvers::make_loaded_rhs_token(ops, 0.0);
    ops.mutations = 0;
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.125, 0.1,
                                      InitialRHSMode::reuse_loaded, ops,
                                      token);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    RecordingOps ops({0.2});
    ops.prepare_initial(0.0);
    ops.preload_exact_left_rhs();
    auto token = ODESolvers::make_loaded_rhs_token(ops, 0.0);
    ops.mutations = 0;
    ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                    InitialRHSMode::reuse_loaded, ops, token);
    const int mutations = ops.mutations;
    const int rhs_calls = ops.rhs_calls;
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                      InitialRHSMode::reuse_loaded, ops,
                                      token);
    });
    assert(ops.mutations == mutations);
    assert(ops.rhs_calls == rhs_calls);
  }
}

void test_zero_kernel_owned_allocations_after_canonical_warmup() {
  for (const auto method : {ExplicitRKMethod::rk4, ExplicitRKMethod::rkf78,
                            ExplicitRKMethod::dp87}) {
    ExponentialOps warmup{{1.0}, {0.0}, 0.0, 0, 0};
    ODESolvers::advance_explicit_rk(method, 0.0, 0.1,
                                    InitialRHSMode::calculate, warmup);

    ExponentialOps measured{{1.0}, {0.0}, 0.0, 0, 0};
    allocation_probe::count = 0;
    allocation_probe::enabled = true;
    ODESolvers::advance_explicit_rk(method, 0.0, 0.1,
                                    InitialRHSMode::calculate, measured);
    allocation_probe::enabled = false;
    assert(allocation_probe::count == 0);
  }
}

template <typename Function> void expect_failure(Function &&function) {
  bool failed = false;
  try {
    function();
  } catch (...) {
    failed = true;
  }
  assert(failed);
}

template <typename Function>
void expect_failure_containing(Function &&function,
                               const std::string &expected_text) {
  bool failed = false;
  try {
    function();
  } catch (const std::exception &error) {
    failed = true;
    assert(std::string(error.what()).find(expected_text) != std::string::npos);
  } catch (...) {
    assert(false && "expected a standard exception");
  }
  assert(failed);
}

struct ThrowingObserver : ODESolvers::NullExplicitRKObserver {
  void stage_rhs(const int, const double, const VectorState &) {
    throw std::runtime_error("observer failure");
  }
};

void test_fail_closed_before_first_mutation_and_no_retry() {
  for (const auto bad_method : {static_cast<ExplicitRKMethod>(99)}) {
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(bad_method, 0.0, 0.1,
                                      InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  for (const double bad_dt : {0.0,
                              std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::quiet_NaN()}) {
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, bad_dt,
                                      InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
  }
  {
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(
          ExplicitRKMethod::rk4, std::numeric_limits<double>::quiet_NaN(), 0.1,
          InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
  }
  {
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(
          ExplicitRKMethod::rk4, 0.0, 0.1,
          static_cast<InitialRHSMode>(99), ops);
    });
    assert(ops.mutations == 0);
  }
  {
    ExplicitRKTableau malformed =
        ODESolvers::explicit_rk_tableau(ExplicitRKMethod::rkf78);
    malformed.a[1][0].denominator = 0;
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(malformed, 0.0, 0.1,
                                      InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
  }
  {
    ExplicitRKTableau malformed =
        ODESolvers::explicit_rk_tableau(ExplicitRKMethod::rkf78);
    malformed.a[1][0].numerator = std::numeric_limits<std::int64_t>::min();
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(malformed, 0.0, 0.1,
                                      InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    ExplicitRKTableau equivalent_spelling =
        ODESolvers::explicit_rk_tableau(ExplicitRKMethod::rk4);
    equivalent_spelling.a[1][0] = {2, 4};
    RecordingOps ops({0.1});
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(equivalent_spelling, 0.0, 0.1,
                                      InitialRHSMode::calculate, ops);
    });
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    ExplicitRKTableau near_overflow =
        ODESolvers::explicit_rk_tableau(ExplicitRKMethod::rk4);
    near_overflow.b[0] = {std::numeric_limits<std::int64_t>::max(), 2};
    near_overflow.b[1] = {std::numeric_limits<std::int64_t>::max(), 3};
    RecordingOps ops({0.1});
    expect_failure_containing(
        [&] {
          ODESolvers::advance_explicit_rk(near_overflow, 0.0, 0.1,
                                          InitialRHSMode::calculate, ops);
        },
        "overflow");
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    ExplicitRKTableau derived_minimum =
        ODESolvers::explicit_rk_tableau(ExplicitRKMethod::rk4);
    derived_minimum.b[0] = {
        std::numeric_limits<std::int64_t>::min() + 1, 1};
    derived_minimum.b[1] = {-1, 1};
    RecordingOps ops({0.1});
    expect_failure_containing(
        [&] {
          ODESolvers::advance_explicit_rk(derived_minimum, 0.0, 0.1,
                                          InitialRHSMode::calculate, ops);
        },
        "overflow");
    assert(ops.mutations == 0);
    assert(ops.rhs_calls == 0);
  }
  {
    RecordingOps ops({0.1});
    ThrowingObserver observer;
    expect_failure([&] {
      ODESolvers::advance_explicit_rk(ExplicitRKMethod::rk4, 0.0, 0.1,
                                      InitialRHSMode::calculate, ops,
                                      observer);
    });
    assert(ops.rhs_calls == 1);
  }
}

void test_descriptor_shape_orders_and_exact_rational_contract() {
  for (const auto &expected : {
           std::pair{ExplicitRKMethod::rk4, std::pair{4, 4}},
           std::pair{ExplicitRKMethod::rkf78, std::pair{11, 7}},
           std::pair{ExplicitRKMethod::dp87, std::pair{13, 8}}}) {
    const auto &table = ODESolvers::explicit_rk_tableau(expected.first);
    assert(static_cast<int>(table.a.size()) == expected.second.first);
    assert(static_cast<int>(table.b.size()) == expected.second.first);
    assert(table.endpoint_order == expected.second.second);
    ODESolvers::validate_explicit_rk_tableau(table);
    for (const auto &row : table.a)
      for (const auto coefficient : row)
        assert(coefficient.denominator > 0);
  }
}

void test_exact_primary_stage_protocol_points() {
  struct Expected {
    ExplicitRKMethod method;
    std::vector<std::pair<std::int64_t, std::int64_t>> fractions;
  };
  const std::array<Expected, 3> expected{{
      {ExplicitRKMethod::rk4, {{0, 1}, {1, 2}, {1, 2}, {1, 1}}},
      {ExplicitRKMethod::rkf78,
       {{0, 1}, {2, 27}, {1, 9}, {1, 6}, {5, 12}, {1, 2},
        {5, 6}, {1, 6}, {2, 3}, {1, 3}, {1, 1}}},
      {ExplicitRKMethod::dp87,
       {{0, 1}, {1, 18}, {1, 12}, {1, 8}, {5, 16}, {3, 8},
        {59, 400}, {93, 200}, {5490023248LL, 9719169821LL},
        {13, 20}, {1201146811LL, 1299019798LL}, {1, 1}, {1, 1}}},
  }};

  for (const auto &entry : expected) {
    RecordingOps ops({0.2});
    ODESolvers::advance_explicit_rk(entry.method, 0.0, 0.125,
                                    InitialRHSMode::calculate, ops);
    assert(ops.stage_points.size() == entry.fractions.size() + 1);
    for (std::size_t stage = 0; stage < entry.fractions.size(); ++stage) {
      const auto &point = ops.stage_points[stage];
      assert(point.kind == 0);
      assert(point.stage_index == static_cast<int>(stage + 1));
      assert(point.stage_count == static_cast<int>(entry.fractions.size()));
      assert(point.numerator == entry.fractions[stage].first);
      assert(point.denominator == entry.fractions[stage].second);
    }
    const auto &endpoint = ops.stage_points.back();
    assert(endpoint.kind == 0);
    assert(endpoint.stage_index == endpoint.stage_count);
    assert(endpoint.numerator == 1);
    assert(endpoint.denominator == 1);
  }
}

} // namespace

int main() {
  test_descriptor_shape_orders_and_exact_rational_contract();
  test_exact_primary_stage_protocol_points();
  test_matches_frozen_legacy_arithmetic_and_complete_trace();
  test_stage_counts_ordering_endpoint_and_rhs_reuse();
  test_observed_fixed_interval_endpoint_orders();
  test_calculate_only_ops_does_not_require_reuse_api();
  test_reuse_loaded_requires_fresh_one_shot_provenance();
  test_zero_kernel_owned_allocations_after_canonical_warmup();
  test_fail_closed_before_first_mutation_and_no_retry();
  std::cout << "Explicit RK kernel tests passed\n";
  return 0;
}
