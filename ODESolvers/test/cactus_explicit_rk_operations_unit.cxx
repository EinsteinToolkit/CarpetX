#include "cactus_explicit_rk_operations.hxx"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct State {
  using scalar_type = double;
  double value{0.0};
  State snapshot_state() const { return *this; }
  State snapshot_rhs() const { return *this; }
};

void test_live_callbacks_remain_the_default_path() {
  State state{2.0};
  State rhs{3.0};
  int prepared = 0;
  int evaluated = 0;
  int validated = 0;
  int updated = 0;
  int accumulated = 0;

  ODESolvers::CactusExplicitRKOperations operations{
      state,
      rhs,
      [&](double) { ++prepared; },
      [&](int) { ++evaluated; },
      [&](int) { ++validated; },
      [&](int, double, double,
          ODESolvers::LinearCombinationView<double, State>) { ++updated; },
      [&](State &, double, const State &) { ++accumulated; }};

  operations.prepare_initial(1.0);
  operations.evaluate_rhs(1);
  operations.validate_rhs(1);
  operations.update_state(1, 1.5, 0.0, {nullptr, nullptr, 0});
  operations.accumulate_rk4(state, 1.0, rhs);
  assert(prepared == 1 && evaluated == 1 && validated == 1 && updated == 1 &&
         accumulated == 1);
  assert(operations.state_generation() == 1);
}

void test_stage_preparation_precedes_poststep_and_rhs() {
  State state{2.0};
  State rhs{3.0};
  std::vector<std::string> trace;
  ODESolvers::CactusExplicitRKOperations operations{
      state,
      rhs,
      [&](double) { trace.emplace_back("set-initial-time"); },
      [&](int) { trace.emplace_back("rhs"); },
      [](int) {},
      [&](int, double, double,
          ODESolvers::LinearCombinationView<double, State>) {
        trace.emplace_back("set-state-time");
      },
      [](State &, double, const State &) {}};
  operations.stage_preparation_callback =
      [&](const ODESolvers::ExplicitRKStagePoint &point, double) {
        assert(point.kind == ODESolvers::ExplicitRKStageKind::primary);
        trace.emplace_back("parent-fill-exchange-prolong");
      };
  operations.live_post_step_callback =
      [&](double) { trace.emplace_back("poststep"); };

  operations.set_stage_point(
      {ODESolvers::ExplicitRKStageKind::primary, 1, 4, {0, 1}});
  operations.prepare_initial(1.0);
  operations.evaluate_rhs(1);
  operations.set_stage_point(
      {ODESolvers::ExplicitRKStageKind::primary, 2, 4, {1, 2}});
  operations.update_state(1, 1.5, 0.0, {nullptr, nullptr, 0});
  operations.evaluate_rhs(2);

  assert((trace == std::vector<std::string>{
                       "set-initial-time", "parent-fill-exchange-prolong",
                       "rhs", "set-state-time",
                       "parent-fill-exchange-prolong", "poststep", "rhs"}));
}

void test_nullable_scratch_hooks_dispatch_without_live_global_state() {
  State state{2.0};
  State rhs{3.0};
  int restore_left = 0;
  int restore_state = 0;
  int probe = 0;
  std::size_t rhs_count = 7;

  using Hooks = ODESolvers::CactusExplicitRKScratchHooks<double, State>;
  Hooks hooks;
  hooks.restore_left = [&](State &dst_state, State &dst_rhs,
                           const State &left_state, const State &left_rhs,
                           double) {
    ++restore_left;
    dst_state = left_state;
    dst_rhs = left_rhs;
  };
  hooks.restore_state = [&](State &dst, const State &source, double) {
    ++restore_state;
    dst = source;
  };
  hooks.probe_endpoint_rhs = [&](State &, State &dst_rhs,
                                 const ODESolvers::ExplicitRKStagePoint &,
                                 double) {
    ++probe;
    dst_rhs.value = 19.0;
    return dst_rhs.snapshot_rhs();
  };
  hooks.rhs_evaluation_count = [&] { return rhs_count; };

  ODESolvers::CactusExplicitRKOperations operations{
      state,
      rhs,
      [](double) {},
      [](int) {},
      [](int) {},
      [](int, double, double,
         ODESolvers::LinearCombinationView<double, State>) {},
      [](State &, double, const State &) {}};
  operations.scratch_hooks = &hooks;
  std::vector<std::string> endpoint_trace;
  int stage_phase = 0;
  operations.stage_materialization_callback =
      [&](State &, const ODESolvers::ExplicitRKStagePoint &point, double) {
        assert(stage_phase == 0);
        stage_phase = 1;
        endpoint_trace.emplace_back(
            point.kind == ODESolvers::ExplicitRKStageKind::fractional
                ? "materialize-fractional"
                : "materialize-endpoint");
      };
  operations.stage_preparation_callback =
      [&](const ODESolvers::ExplicitRKStagePoint &point, double) {
        assert(stage_phase == 1);
        stage_phase = 2;
        endpoint_trace.emplace_back(
            point.kind == ODESolvers::ExplicitRKStageKind::fractional
                ? "prepare-fractional"
                : "prepare-endpoint");
      };
  hooks.post_step_after_update =
      [&](State &, const ODESolvers::ExplicitRKStagePoint &point, double) {
        assert(stage_phase == 2);
        stage_phase = 3;
        endpoint_trace.emplace_back(
            point.kind == ODESolvers::ExplicitRKStageKind::fractional
                ? "scratch-poststep-fractional"
                : "scratch-poststep-endpoint");
      };
  hooks.evaluate_rhs =
      [&](State &, State &, int,
          const ODESolvers::ExplicitRKStagePoint &point, double) {
        assert(stage_phase == 2 || stage_phase == 3);
        stage_phase = 0;
        endpoint_trace.emplace_back(
            point.kind == ODESolvers::ExplicitRKStageKind::fractional
                ? "scratch-rhs-fractional"
                : "scratch-rhs-endpoint");
      };
  hooks.probe_endpoint_rhs = [&](State &, State &dst_rhs,
                                 const ODESolvers::ExplicitRKStagePoint &,
                                 double) {
    assert(stage_phase == 3);
    stage_phase = 0;
    endpoint_trace.emplace_back("scratch-rhs-endpoint");
    ++probe;
    dst_rhs.value = 19.0;
    return dst_rhs.snapshot_rhs();
  };

  const State left{5.0};
  const State left_rhs{7.0};
  operations.restore_left(left, left_rhs, 1.0);
  auto token = ODESolvers::make_loaded_rhs_token(operations, 1.0);
  operations.validate_loaded_rhs_provenance(token.provenance());
  assert(state.value == 5.0 && rhs.value == 7.0 && restore_left == 1);

  operations.set_stage_point(
      {ODESolvers::ExplicitRKStageKind::fractional, 1, 4, {0, 1}});
  operations.prepare_initial(1.0);
  operations.evaluate_rhs(1);
  operations.set_stage_point(
      {ODESolvers::ExplicitRKStageKind::fractional, 2, 4, {1, 4}});
  operations.update_state(1, 1.25, 0.0, {nullptr, nullptr, 0});
  operations.evaluate_rhs(2);

  const State endpoint{11.0};
  operations.restore_state(endpoint, 2.0);
  assert(state.value == 11.0 && restore_state == 1);
  assert(operations
             .probe_endpoint_rhs(
                 2.0, {ODESolvers::ExplicitRKStageKind::endpoint_probe, 1, 1,
                       {1, 1}})
             .value == 19.0 &&
         probe == 1);
  assert((endpoint_trace == std::vector<std::string>{
                                "materialize-fractional",
                                "prepare-fractional",
                                "scratch-rhs-fractional",
                                "materialize-fractional",
                                "prepare-fractional",
                                "scratch-poststep-fractional",
                                "scratch-rhs-fractional",
                                "materialize-endpoint", "prepare-endpoint",
                                "scratch-poststep-endpoint",
                                "scratch-rhs-endpoint"}));
  assert(stage_phase == 0);
  assert(operations.rhs_evaluation_count() == rhs_count);
}

void test_missing_scratch_hook_fails_closed() {
  State state{};
  State rhs{};
  ODESolvers::CactusExplicitRKOperations operations{
      state,
      rhs,
      [](double) {},
      [](int) {},
      [](int) {},
      [](int, double, double,
         ODESolvers::LinearCombinationView<double, State>) {},
      [](State &, double, const State &) {}};
  bool threw = false;
  try {
    operations.restore_state(state, 1.0);
  } catch (const std::logic_error &) {
    threw = true;
  }
  assert(threw);
}

void test_scratch_stage_preparation_requires_materialization() {
  State state{};
  State rhs{};
  using Hooks = ODESolvers::CactusExplicitRKScratchHooks<double, State>;
  Hooks hooks;
  ODESolvers::CactusExplicitRKOperations operations{
      state,
      rhs,
      [](double) {},
      [](int) {},
      [](int) {},
      [](int, double, double,
         ODESolvers::LinearCombinationView<double, State>) {},
      [](State &, double, const State &) {}};
  operations.scratch_hooks = &hooks;
  operations.stage_preparation_callback =
      [](const ODESolvers::ExplicitRKStagePoint &, double) {};
  operations.set_stage_point(
      {ODESolvers::ExplicitRKStageKind::fractional, 1, 4, {0, 1}});
  bool threw = false;
  try {
    operations.prepare_initial(1.0);
  } catch (const std::logic_error &) {
    threw = true;
  }
  assert(threw);
}

} // namespace

int main() {
  test_live_callbacks_remain_the_default_path();
  test_stage_preparation_precedes_poststep_and_rhs();
  test_nullable_scratch_hooks_dispatch_without_live_global_state();
  test_missing_scratch_hook_fails_closed();
  test_scratch_stage_preparation_requires_materialization();
}
