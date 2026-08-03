#include "explicit_rk_dense_provider.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using ODESolvers::ExplicitRKMethod;
using ODESolvers::InitialRHSMode;
using ODESolvers::LoadedRHSProvenance;

void require(const bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

template <class Exception, class Function>
void require_throws(Function &&function, const std::string &message) {
  bool threw = false;
  try {
    function();
  } catch (const Exception &) {
    threw = true;
  }
  require(threw, message);
}

std::string fingerprint_hex(
    const ODESolvers::ExplicitRKTableauFingerprint &fingerprint) {
  std::ostringstream stream;
  stream << std::hex << std::setfill('0');
  for (const auto byte : fingerprint)
    stream << std::setw(2) << static_cast<unsigned int>(byte);
  return stream.str();
}

struct TestState {
  static int live_count;
  static int move_attempts;
  static int fail_move_countdown;
  std::vector<long double> values;

  explicit TestState(std::vector<long double> values_)
      : values(std::move(values_)) {
    ++live_count;
  }
  ~TestState() { --live_count; }

  TestState(const TestState &) = delete;
  TestState &operator=(const TestState &) = delete;
  TestState(TestState &&other) : values() {
    ++move_attempts;
    if (fail_move_countdown == 0) {
      fail_move_countdown = -1;
      throw std::bad_alloc();
    }
    if (fail_move_countdown > 0)
      --fail_move_countdown;
    values = std::move(other.values);
    ++live_count;
  }
  TestState &operator=(TestState &&other) noexcept {
    values = std::move(other.values);
    return *this;
  }
};

int TestState::live_count = 0;
int TestState::move_attempts = 0;
int TestState::fail_move_countdown = -1;

TestState clone_state(const TestState &source) {
  return TestState(source.values);
}

bool same_state(const TestState &left, const TestState &right) {
  return left.values == right.values;
}

TestState square_state(const TestState &state) {
  std::vector<long double> result(state.values.size());
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = state.values[i] * state.values[i];
  return TestState(std::move(result));
}

struct ScratchOps {
  using scalar_type = long double;
  using state_type = TestState;

  TestState state_value;
  TestState rhs_value;
  long double time{0.0L};
  std::uint64_t generation{0};
  bool loaded_rhs{false};
  std::size_t rhs_calls{0};
  std::size_t token_consumptions{0};
  int fail_event{0};
  int event_count{0};
  bool record_events{true};
  std::vector<std::string> events;
  std::vector<long double> restored_left_values;
  bool require_parent_fill{false};
  bool parent_fill_ready{false};
  std::optional<ODESolvers::ExplicitRKStagePoint> pending_stage_point;
  std::size_t parent_fill_count{0};
  std::function<void(const ODESolvers::ExplicitRKStagePoint &, long double)>
      driver_prepare;

  explicit ScratchOps(const std::size_t components = 1,
                      const int fail_event_ = 0)
      : state_value(std::vector<long double>(components, 0.0L)),
        rhs_value(std::vector<long double>(components, 0.0L)),
        fail_event(fail_event_) {}

  void event(const char *const name) {
    ++event_count;
    if (record_events)
      events.emplace_back(name);
    if (fail_event == event_count)
      throw std::runtime_error(std::string("injected failure at ") + name);
  }

  const TestState &state() const noexcept { return state_value; }
  const TestState &rhs() const noexcept { return rhs_value; }

  TestState snapshot_state() {
    event("snapshot-state");
    return clone_state(state_value);
  }
  TestState snapshot_rhs() {
    event("snapshot-rhs");
    return clone_state(rhs_value);
  }

  void restore_left(const TestState &left, const TestState &left_rhs,
                    const long double left_time) {
    event("restore-left");
    require(left.values.size() == left_rhs.values.size(),
            "restore_left shape mismatch");
    state_value = clone_state(left);
    rhs_value = clone_state(left_rhs);
    time = left_time;
    ++generation;
    loaded_rhs = true;
    restored_left_values.push_back(left.values.front());
  }

  void restore_state(const TestState &state, const long double state_time) {
    event("restore-state");
    state_value = clone_state(state);
    time = state_time;
    ++generation;
    loaded_rhs = false;
  }

  void set_stage_point(const ODESolvers::ExplicitRKStagePoint &point) {
    pending_stage_point = point;
    parent_fill_ready = false;
  }

  void prepare_current_stage() {
    if (!require_parent_fill)
      return;
    require(pending_stage_point.has_value(),
            "L1 auxiliary stage reached state/time without an exact stage point");
    require(static_cast<bool>(driver_prepare),
            "L1 auxiliary stage has no driver preparation callback");
    driver_prepare(*pending_stage_point, time);
    pending_stage_point.reset();
    parent_fill_ready = true;
    ++parent_fill_count;
  }

  void prepare_initial(const long double stage_time) {
    time = stage_time;
    prepare_current_stage();
  }

  void update_state(
      const int, const long double stage_time,
      const long double destination_scale,
      const ODESolvers::LinearCombinationView<long double, TestState>
          combination) {
    event("stage-update");
    std::vector<long double> updated(state_value.values.size(), 0.0L);
    for (std::size_t component = 0; component < updated.size(); ++component)
      updated[component] = destination_scale * state_value.values[component];
    for (std::size_t source = 0; source < combination.size; ++source)
      for (std::size_t component = 0; component < updated.size(); ++component)
        updated[component] += combination.factors[source] *
                              combination.sources[source]->values[component];
    state_value.values = std::move(updated);
    time = stage_time;
    prepare_current_stage();
    ++generation;
    loaded_rhs = false;
  }

  void evaluate_rhs(const int) {
    event("stage-rhs");
    if (require_parent_fill)
      require(parent_fill_ready,
              "L1 auxiliary RHS ran before the L0 parent fill");
    ++rhs_calls;
    for (std::size_t component = 0; component < state_value.values.size();
         ++component)
      rhs_value.values[component] =
          state_value.values[component] * state_value.values[component];
  }

  void validate_rhs(const int) {
    event("validate-rhs");
    for (const long double value : rhs_value.values)
      if (!std::isfinite(value))
        throw std::runtime_error("non-finite test RHS");
  }

  void accumulate_rk4(TestState &accumulator, const long double factor,
                      const TestState &increment) {
    event("rk4-accumulate");
    for (std::size_t component = 0; component < accumulator.values.size();
         ++component)
      accumulator.values[component] += factor * increment.values[component];
  }

  std::uint64_t state_generation() const noexcept { return generation; }
  LoadedRHSProvenance<long double>
  loaded_rhs_provenance(const long double left_time) const {
    if (!loaded_rhs || left_time != time)
      throw std::logic_error("no loaded left RHS");
    return {generation, left_time};
  }
  void validate_loaded_rhs_provenance(
      const LoadedRHSProvenance<long double> &provenance) const {
    if (!loaded_rhs || provenance.state_generation != generation ||
        provenance.left_time != time)
      throw std::logic_error("stale loaded left RHS");
  }
  void consume_loaded_rhs(
      const LoadedRHSProvenance<long double> &provenance) {
    event("consume-loaded-rhs");
    validate_loaded_rhs_provenance(provenance);
    loaded_rhs = false;
    ++token_consumptions;
  }

  TestState probe_endpoint_rhs(
      const long double probe_time,
      const ODESolvers::ExplicitRKStagePoint &stage_point) {
    event("endpoint-probe");
    if (probe_time != time)
      throw std::logic_error("endpoint probe time mismatch");
    set_stage_point(stage_point);
    prepare_current_stage();
    evaluate_rhs(0);
    validate_rhs(0);
    return snapshot_rhs();
  }

  std::size_t rhs_evaluation_count() const noexcept { return rhs_calls; }
};

class ThreeLevelParentFill final : public CarpetX::StagePreparer {
public:
  std::array<std::size_t, 3> fill_generation{{1, 0, 0}};
  std::vector<CarpetX::StagePoint> points;

  void prepare_stage(const CarpetX::StepContext &context,
                     const CarpetX::StagePoint &point) override {
    require(context.level == 1,
            "synthetic auxiliary solve did not target level 1");
    fill_generation[1] = fill_generation[0] + 1;
    points.push_back(point);
  }
};

static_assert(!std::is_copy_constructible<TestState>::value,
              "fixture state must remain move-only");
static_assert(
    !std::is_copy_constructible<
        ODESolvers::ReferenceDenseSampleBatch<TestState>>::value,
    "sample batch must remain move-only");

struct PrimaryResult {
  TestState endpoint;
  TestState final_legacy_rhs;
};

PrimaryResult primary_step(const ExplicitRKMethod method,
                           const long double t0, const long double dt,
                           const TestState &left,
                           const TestState &left_rhs) {
  ScratchOps ops(left.values.size());
  ops.restore_left(left, left_rhs, t0);
  auto token = ODESolvers::make_loaded_rhs_token(ops, t0);
  ODESolvers::advance_explicit_rk(method, t0, dt,
                                  InitialRHSMode::reuse_loaded, ops, token);
  return {ops.snapshot_state(), ops.snapshot_rhs()};
}

struct MethodExpectation {
  ExplicitRKMethod method;
  CarpetX::SubcyclingODEMethod dense_method;
  const char *fingerprint;
  int order;
  int stages;
  int extra_rhs;
  int controls;
  int fractional_solves;
};

const std::array<MethodExpectation, 3> expectations{{
    {ExplicitRKMethod::rk4, CarpetX::SubcyclingODEMethod::rk4,
     "5a1c84e15a045292e8478588e49ea15a8c7d1dd792ba3e13a31b9bf99556f32d",
     4, 4, 4, 5, 1},
    {ExplicitRKMethod::rkf78, CarpetX::SubcyclingODEMethod::rkf78_order7,
     "3d600f326c614c8ba6f35ac84d4097bb820d7cf4e37cd3c43ae899ef5e9833a5",
     7, 11, 23, 8, 2},
    {ExplicitRKMethod::dp87, CarpetX::SubcyclingODEMethod::dp87_order8,
     "a3faabaef46ae521e6aa524e6f0307f0b999dcb9414ecd50d30eb884e72fdc92",
     8, 13, 39, 9, 3}}};

void test_fingerprints_capabilities_and_registry() {
  CarpetX::DenseOutputRegistry registry;
  for (const auto &expected : expectations) {
    const auto fingerprint =
        ODESolvers::explicit_rk_tableau_fingerprint(expected.method);
    require(fingerprint_hex(fingerprint) == expected.fingerprint,
            "tableau fingerprint changed");
    const auto capability =
        ODESolvers::reference_dense_capability(expected.method);
    require(capability.method == expected.dense_method,
            "dense method mapping changed");
    require(capability.tableau_fingerprint == fingerprint,
            "capability fingerprint mismatch");
    require(capability.endpoint_order == expected.order &&
                capability.dense_uniform_order == expected.order &&
                capability.stage_count == expected.stages &&
                capability.extra_rhs_evaluations == expected.extra_rhs &&
                capability.persistent_vector_count == expected.controls &&
                capability.arbitrary_theta && capability.verified,
            "dense capability fields changed");
    registry.register_provider(
        ODESolvers::make_reference_dense_provider(expected.method));
    require(registry.require(expected.dense_method, fingerprint)->capability()
                .tableau_fingerprint == fingerprint,
            "registered provider lookup failed");
    auto wrong = fingerprint;
    wrong.front() ^= 0x80U;
    require_throws<std::out_of_range>(
        [&] { (void)registry.require(expected.dense_method, wrong); },
        "fingerprint mismatch did not fail closed");
  }
  require_throws<std::invalid_argument>(
      [] {
        (void)ODESolvers::reference_dense_capability(
            static_cast<ExplicitRKMethod>(99));
      },
      "unsupported provider method did not fail before construction");

  CarpetX::DenseOutputRegistry mismatch_registry;
  auto mismatched =
      ODESolvers::reference_dense_capability(ExplicitRKMethod::rk4);
  const auto canonical_rk4_fingerprint = mismatched.tableau_fingerprint;
  mismatched.tableau_fingerprint =
      ODESolvers::explicit_rk_tableau_fingerprint(ExplicitRKMethod::rkf78);
  mismatch_registry.register_provider(
      std::make_shared<const CarpetX::DenseOutputProvider>(mismatched));
  require_throws<std::out_of_range>(
      [&] {
        (void)mismatch_registry.require(CarpetX::SubcyclingODEMethod::rk4,
                                        canonical_rk4_fingerprint);
      },
      "coefficient/method fingerprint mismatch did not fail closed");
}

std::vector<double> numeric_controls(
    const ExplicitRKMethod method, const long double parent_dt,
    const ODESolvers::ReferenceDenseSampleBatch<TestState> &batch,
    const std::size_t component) {
  const auto dense_method = ODESolvers::subcycling_method(method);
  const auto &stencil = CarpetX::reference_dense_stencil(dense_method);
  require(batch.size() == stencil.sample_count(),
          "numeric sample count mismatch");
  std::vector<double> samples(batch.size(), 0.0);
  for (std::size_t i = 0; i < batch.size(); ++i) {
    const auto &sample = batch.samples()[i];
    long double value = sample.payload.values.at(component);
    if (sample.kind == CarpetX::DenseSampleKind::scaled_derivative)
      value *= parent_dt;
    samples[i] = static_cast<double>(value);
  }
  return stencil.make_controls(samples);
}

class TestDenseState final : public CarpetX::DenseStateVector {
public:
  static int live_count;
  explicit TestDenseState(std::vector<double> values_)
      : values(std::move(values_)) {
    ++live_count;
  }
  ~TestDenseState() override { --live_count; }

  bool compatible(const CarpetX::DenseStateVector &other) const noexcept
      override {
    const auto *const typed = dynamic_cast<const TestDenseState *>(&other);
    return typed != nullptr && typed->values.size() == values.size();
  }
  void copy_from(const CarpetX::DenseStateVector &other) override {
    const auto &typed = dynamic_cast<const TestDenseState &>(other);
    values = typed.values;
  }
  void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const CarpetX::DenseStateVector *> &sources) override {
    require(weights.size() == sources.size(),
            "dense state linear combination shape mismatch");
    std::fill(values.begin(), values.end(), 0.0);
    for (std::size_t source = 0; source < sources.size(); ++source) {
      const auto &typed =
          dynamic_cast<const TestDenseState &>(*sources[source]);
      for (std::size_t component = 0; component < values.size(); ++component)
        values[component] += weights[source] * typed.values[component];
    }
  }

  std::vector<double> values;
};

int TestDenseState::live_count = 0;

std::shared_ptr<const CarpetX::DenseInterval> build_test_interval(
    const MethodExpectation &expected, const long double parent_dt,
    const ODESolvers::ReferenceDenseSampleBatch<TestState> &batch,
    const int fail_control = -1) {
  const auto provider =
      ODESolvers::make_reference_dense_provider(expected.method);
  const auto capability = provider->capability();
  const CarpetX::DenseIntervalId id{
      0,
      CarpetX::step_clock_t(0),
      CarpetX::step_clock_t(1),
      0.0,
      static_cast<double>(parent_dt),
      capability.method,
      capability.tableau_fingerprint};
  auto builder = provider->begin_interval(id);
  std::vector<std::vector<double>> component_controls;
  for (std::size_t component = 0;
       component < batch.samples().front().payload.values.size(); ++component)
    component_controls.push_back(
        numeric_controls(expected.method, parent_dt, batch, component));
  for (int control = 0; control < expected.controls; ++control) {
    if (control == fail_control)
      throw std::bad_alloc();
    std::vector<double> values(component_controls.size(), 0.0);
    for (std::size_t component = 0; component < values.size(); ++component)
      values[component] =
          component_controls[component][static_cast<std::size_t>(control)];
    builder->add_control(
        std::make_unique<TestDenseState>(std::move(values)));
  }
  return builder->seal();
}

void test_collection_and_interval(const MethodExpectation &expected) {
  TestState left({0.2L, 0.35L});
  TestState left_rhs = square_state(left);
  auto primary = primary_step(expected.method, 0.125L, 0.0625L, left,
                              left_rhs);
  const auto left_before = left.values;
  const auto rhs_before = left_rhs.values;
  const auto endpoint_before = primary.endpoint.values;
  const auto legacy_before = primary.final_legacy_rhs.values;

  ScratchOps scratch(left.values.size());
  auto batch = ODESolvers::collect_reference_dense_samples(
      expected.method, 0.125L, 0.0625L, left, left_rhs, primary.endpoint,
      scratch);
  require(batch.size() == static_cast<std::size_t>(expected.controls),
          "collector sample count changed");
  require(scratch.rhs_calls == static_cast<std::size_t>(expected.extra_rhs),
          "collector extra RHS count changed");
  require(scratch.token_consumptions ==
              static_cast<std::size_t>(expected.fractional_solves),
          "fractional solves did not use fresh one-shot tokens");
  require(scratch.restored_left_values.size() ==
              static_cast<std::size_t>(expected.fractional_solves + 2),
          "collector did not restore each fractional solve independently");
  require(left.values == left_before && left_rhs.values == rhs_before &&
              primary.endpoint.values == endpoint_before &&
              primary.final_legacy_rhs.values == legacy_before,
          "collector mutated a primary state or final legacy RHS");
  for (const long double restored : scratch.restored_left_values)
    require(restored == left.values.front(),
            "fractional solve was not restored from the common left state");

  const auto &constraints = CarpetX::reference_dense_stencil(
                                expected.dense_method)
                                .specification()
                                .constraints;
  for (std::size_t i = 0; i < constraints.size(); ++i) {
    require(batch.samples()[i].theta == constraints[i].theta &&
                batch.samples()[i].kind == constraints[i].kind,
            "collector sample order differs from Phase 4");
    if (batch.samples()[i].kind ==
        CarpetX::DenseSampleKind::scaled_derivative) {
      const auto &derivative = batch.samples()[i].payload.values;
      const auto &value = batch.samples()[i - 1].payload.values;
      for (std::size_t component = 0; component < derivative.size();
           ++component)
        require(derivative[component] == value[component] * value[component],
                "collector scaled a raw derivative payload");
    }
    if (batch.samples()[i].theta == 1.0 &&
        batch.samples()[i].kind == CarpetX::DenseSampleKind::value)
      require(same_state(batch.samples()[i].payload, primary.endpoint),
              "accepted endpoint sample lost bit identity");
  }

  {
    auto interval = build_test_interval(expected, 0.0625L, batch);
    require(interval->control_count() ==
                static_cast<std::size_t>(expected.controls),
            "interval control count changed");
    TestDenseState begin(std::vector<double>(left.values.size(), -1.0));
    TestDenseState end(std::vector<double>(left.values.size(), -1.0));
    TestDenseState interior(std::vector<double>(left.values.size(), -1.0));
    interval->evaluate(0.0, begin);
    interval->evaluate(1.0, end);
    interval->evaluate(0.37, interior);
    for (std::size_t component = 0; component < left.values.size();
         ++component) {
      require(begin.values[component] == static_cast<double>(left.values[component]),
              "theta=0 endpoint is not bit exact");
      require(end.values[component] ==
                  static_cast<double>(primary.endpoint.values[component]),
              "theta=1 endpoint is not bit exact");
      require(std::isfinite(interior.values[component]),
              "arbitrary theta evaluation is not finite");
    }
  }
  require(TestDenseState::live_count == 0,
          "sealed interval leaked test controls");

  for (int control = 0; control < expected.controls; ++control) {
    require_throws<std::bad_alloc>(
        [&] { (void)build_test_interval(expected, 0.0625L, batch, control); },
        "control allocation failure did not roll back");
    require(TestDenseState::live_count == 0,
            "failed interval publication leaked controls");
  }
}

void test_collector_failure_rollback(const MethodExpectation &expected) {
  TestState left({0.24L});
  TestState left_rhs = square_state(left);
  auto primary = primary_step(expected.method, 0.0L, 0.1L, left, left_rhs);
  const auto left_before = left.values;
  const auto rhs_before = left_rhs.values;
  const auto endpoint_before = primary.endpoint.values;
  const auto legacy_before = primary.final_legacy_rhs.values;

  ScratchOps baseline(1);
  {
    auto batch = ODESolvers::collect_reference_dense_samples(
        expected.method, 0.0L, 0.1L, left, left_rhs, primary.endpoint,
        baseline);
    require(!batch.empty(), "baseline batch unexpectedly empty");
  }
  const int critical_event_count = baseline.event_count;
  require(critical_event_count > expected.extra_rhs,
          "failure fixture did not observe collector operations");

  for (int failure = 1; failure <= critical_event_count; ++failure) {
    const int live_before = TestState::live_count;
    {
      ScratchOps scratch(1, failure);
      require_throws<std::runtime_error>(
          [&] {
            (void)ODESolvers::collect_reference_dense_samples(
                expected.method, 0.0L, 0.1L, left, left_rhs,
                primary.endpoint, scratch);
          },
          "injected collector operation failure was not propagated");
    }
    require(TestState::live_count == live_before,
            "failed collector transaction leaked sample states");
    require(left.values == left_before && left_rhs.values == rhs_before &&
                primary.endpoint.values == endpoint_before &&
                primary.final_legacy_rhs.values == legacy_before,
            "failed collector transaction mutated primary data");
  }

  TestState::move_attempts = 0;
  {
    ScratchOps move_baseline(1);
    move_baseline.record_events = false;
    auto batch = ODESolvers::collect_reference_dense_samples(
        expected.method, 0.0L, 0.1L, left, left_rhs, primary.endpoint,
        move_baseline);
    require(!batch.empty(), "move-allocation baseline batch is empty");
  }
  const int move_positions = TestState::move_attempts;
  require(move_positions >= expected.controls,
          "fixture did not observe all sample/state move positions");
  for (int failure = 0; failure < move_positions; ++failure) {
    const int live_before = TestState::live_count;
    {
      ScratchOps allocation_scratch(1);
      allocation_scratch.record_events = false;
      TestState::fail_move_countdown = failure;
      require_throws<std::bad_alloc>(
          [&] {
            (void)ODESolvers::collect_reference_dense_samples(
                expected.method, 0.0L, 0.1L, left, left_rhs,
                primary.endpoint, allocation_scratch);
          },
          "fixture move/allocation failure was not propagated");
      TestState::fail_move_countdown = -1;
    }
    require(TestState::live_count == live_before,
            "fixture move/allocation failure leaked state");
  }
}

long double exact_solution(const long double y0, const long double time) {
  return y0 / (1.0L - y0 * time);
}

struct ConvergenceError {
  long double dense;
  long double endpoint;
};

constexpr long double convergence_initial = 0.8L;
constexpr long double convergence_final_time = 0.85L;
constexpr std::array<long double, 4> convergence_sample_thetas{{
    0.13L, 0.37L, 0.61L, 0.89L}};

ConvergenceError convergence_error(const MethodExpectation &expected,
                                   const int steps) {
  const long double dt =
      convergence_final_time / static_cast<long double>(steps);
  long double current = convergence_initial;
  long double dense_error_squared = 0.0L;
  std::size_t dense_sample_count = 0;
  for (int step = 0; step < steps; ++step) {
    const long double t0 = static_cast<long double>(step) * dt;
    TestState left({current});
    TestState left_rhs({current * current});
    auto primary = primary_step(expected.method, t0, dt, left, left_rhs);
    ScratchOps scratch(1);
    auto batch = ODESolvers::collect_reference_dense_samples(
        expected.method, t0, dt, left, left_rhs, primary.endpoint, scratch);
    auto interval = build_test_interval(expected, dt, batch);
    for (const long double theta : convergence_sample_thetas) {
      TestDenseState evaluated({0.0});
      interval->evaluate(static_cast<double>(theta), evaluated);
      const long double computed =
          static_cast<long double>(evaluated.values.front());
      const long double exact =
          exact_solution(convergence_initial, t0 + theta * dt);
      const long double error = computed - exact;
      dense_error_squared += error * error;
      ++dense_sample_count;
    }
    current = primary.endpoint.values.front();
  }
  return {std::sqrt(dense_error_squared /
                    static_cast<long double>(dense_sample_count)),
          std::abs(current - exact_solution(convergence_initial,
                                            convergence_final_time))};
}

void test_observed_orders(const MethodExpectation &expected) {
  const auto coarse = convergence_error(expected, 8);
  const auto medium = convergence_error(expected, 16);
  const auto fine = convergence_error(expected, 32);
  require(coarse.dense > medium.dense && medium.dense > fine.dense,
          "dense convergence errors do not decrease");
  require(coarse.endpoint > medium.endpoint &&
              medium.endpoint > fine.endpoint,
          "endpoint convergence errors do not decrease");
  const auto observed = [](const long double left, const long double right) {
    return std::log(left / right) / std::log(2.0L);
  };
  const long double threshold =
      static_cast<long double>(expected.order) - 0.25L;
  const long double dense_order_1 = observed(coarse.dense, medium.dense);
  const long double dense_order_2 = observed(medium.dense, fine.dense);
  const long double endpoint_order_1 =
      observed(coarse.endpoint, medium.endpoint);
  const long double endpoint_order_2 =
      observed(medium.endpoint, fine.endpoint);
  std::cout << std::setprecision(12) << std::scientific
            << "ORDER_TABLE p=" << expected.order
            << " h=" << static_cast<double>(convergence_final_time / 8.0L)
            << ',' << static_cast<double>(convergence_final_time / 16.0L)
            << ',' << static_cast<double>(convergence_final_time / 32.0L)
            << " dense_error=" << static_cast<double>(coarse.dense) << ','
            << static_cast<double>(medium.dense) << ','
            << static_cast<double>(fine.dense)
            << " endpoint_error=" << static_cast<double>(coarse.endpoint)
            << ',' << static_cast<double>(medium.endpoint) << ','
            << static_cast<double>(fine.endpoint)
            << " dense_order=" << static_cast<double>(dense_order_1) << ','
            << static_cast<double>(dense_order_2)
            << " endpoint_order=" << static_cast<double>(endpoint_order_1)
            << ',' << static_cast<double>(endpoint_order_2) << '\n';
  std::ostringstream detail;
  detail << "method order " << expected.order << ": dense="
         << static_cast<double>(dense_order_1) << ','
         << static_cast<double>(dense_order_2) << " endpoint="
         << static_cast<double>(endpoint_order_1) << ','
         << static_cast<double>(endpoint_order_2);
  require(dense_order_1 >= threshold && dense_order_2 >= threshold,
          "dense observed order is below the endpoint order gate: " +
              detail.str());
  require(endpoint_order_1 >= threshold && endpoint_order_2 >= threshold,
          "primary endpoint observed order regressed: " + detail.str());
}

void test_three_level_auxiliary_dense_prepares_l0_fill_before_l1_rhs() {
  const long double t0 = 1000000000000.0L;
  const long double dt = 0.5L;
  TestState left({0.25L});
  TestState left_rhs({0.0625L});
  auto primary = primary_step(ExplicitRKMethod::rk4, t0, dt, left, left_rhs);

  ScratchOps scratch;
  scratch.require_parent_fill = true;
  ThreeLevelParentFill preparer;
  const CarpetX::StepContext context{
      1, CarpetX::step_clock_t(17, 8), CarpetX::step_clock_t(19, 8),
      static_cast<double>(t0), static_cast<double>(t0 + dt),
      CarpetX::SubcyclingODEMethod::rk4};
  auto *const transaction = reinterpret_cast<CarpetX::ScratchStateTransaction *>(
      static_cast<std::uintptr_t>(0x10));
  CarpetX::ScopedStepContext scope(context, preparer, transaction);
  scratch.driver_prepare = [&](const ODESolvers::ExplicitRKStagePoint &point,
                               const long double stage_time) {
    const auto kind =
        point.kind == ODESolvers::ExplicitRKStageKind::fractional
            ? CarpetX::StageKind::fractional
            : CarpetX::StageKind::endpoint_probe;
    CarpetX::prepare_stage(
        {kind, point.stage_index, point.stage_count,
         CarpetX::step_clock_t(point.parent_fraction.numerator,
                               point.parent_fraction.denominator),
         static_cast<double>(stage_time)},
        context.method);
  };

  const auto samples = ODESolvers::collect_reference_dense_samples(
      ExplicitRKMethod::rk4, t0, dt, left, left_rhs, primary.endpoint,
      scratch);
  require(samples.size() == 5,
          "three-level auxiliary dense sample batch is incomplete");
  require(scratch.parent_fill_count > 0,
          "L1 auxiliary dense construction never invoked the L0 fill");
  require(preparer.fill_generation[1] > preparer.fill_generation[0],
          "L0 parent fill did not reach L1");
  require(preparer.fill_generation[2] == 0,
          "L1 auxiliary construction unexpectedly advanced level 2");
  require(std::any_of(preparer.points.begin(), preparer.points.end(),
                      [](const auto &point) {
                        return point.kind == CarpetX::StageKind::fractional;
                      }),
          "fractional stages were not distinguished");
  require(std::any_of(preparer.points.begin(), preparer.points.end(),
                      [](const auto &point) {
                        return point.kind == CarpetX::StageKind::endpoint_probe;
                      }),
          "endpoint probe was not distinguished");
}

} // namespace

int main() {
  test_fingerprints_capabilities_and_registry();
  for (const auto &expected : expectations) {
    test_collection_and_interval(expected);
    test_collector_failure_rollback(expected);
    test_observed_orders(expected);
  }
  test_three_level_auxiliary_dense_prepares_l0_fill_before_l1_rhs();
  require(TestState::live_count == 0, "test states leaked at process exit");
  require(TestDenseState::live_count == 0,
          "dense test states leaked at process exit");
  std::cout << "Explicit RK dense provider tests passed\n";
  return 0;
}
