#include "../src/subcycling_runtime_clock.hxx"
#include "../src/subcycling_static_v1_contract.hxx"

#include <array>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using CarpetX::RuntimeClockMetadata;
using CarpetX::StaticV1EvolvePath;
using CarpetX::StaticV1LevelOneCreationEnvelope;
using CarpetX::StaticV1PolicyEnvelope;
using CarpetX::StaticV1Rational;
using CarpetX::StaticV1SyncObserverAction;
using CarpetX::StepContext;
using CarpetX::SubcyclingODEMethod;
using CarpetX::candidate_runtime_clock;
using CarpetX::cycle_then_capture_static_v1_level_state;
using CarpetX::has_exact_static_v1_test_odesolvers2_active_thorns;
using CarpetX::permits_static_v1_level_one_creation;
using CarpetX::select_static_v1_evolve_path;
using CarpetX::static_v1_parent_refinement_ratio_index;
using CarpetX::should_stop_static_v1_initial_regrid_loop;
using CarpetX::static_v1_sync_observer_order;
using CarpetX::step_clock_t;
using CarpetX::validate_static_v1_policy_envelope;

template <class Function> void rejects(Function &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  assert(rejected);
}

StaticV1PolicyEnvelope valid_policy() {
  return {1, 2, 2, 0, false, false, true, false, false,
          true, true, true, true, true, true};
}

StaticV1LevelOneCreationEnvelope valid_level_one_creation() {
  return {true, 1, 0, 1, 2, 1, 0, 0, 0,
          StaticV1Rational{0, 1}, StaticV1Rational{1, 1}, false,
          {2, 2, 2}, 0.0};
}

void test_legacy_selection_is_a_strict_false_branch() {
  assert(select_static_v1_evolve_path(false) ==
         StaticV1EvolvePath::legacy);
  assert(select_static_v1_evolve_path(true) ==
         StaticV1EvolvePath::static_v1);
}

void test_static_policy_fails_closed_before_runtime_mutation() {
  validate_static_v1_policy_envelope(valid_policy());

  auto policy = valid_policy();
  policy.patch_count = 2;
  rejects([&] { validate_static_v1_policy_envelope(policy); });

  policy = valid_policy();
  policy.recovering = true;
  rejects([&] { validate_static_v1_policy_envelope(policy); });

  policy = valid_policy();
  policy.restrict_during_sync = true;
  rejects([&] { validate_static_v1_policy_envelope(policy); });

  policy = valid_policy();
  policy.complete_method_schema = false;
  rejects([&] { validate_static_v1_policy_envelope(policy); });

  policy = valid_policy();
  policy.bounded_test_odesolvers2_configuration = false;
  rejects([&] { validate_static_v1_policy_envelope(policy); });
}

void test_static_v1_level_one_creation_envelope_is_exact() {
  assert(permits_static_v1_level_one_creation(valid_level_one_creation()));

  auto envelope = valid_level_one_creation();
  envelope.subcycling_enabled = false;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.configured_patch_count = 2;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.callback_patch = 1;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.callback_level = 0;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.configured_level_count = 3;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.existing_level_count = 0;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.regrid_every = 1;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_patch = 1;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_level = 1;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_iteration = {1, 1};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_delta_iteration = {2, 1};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_iteration = {1, 2};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_delta_iteration = {3, 2};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.coarse_is_subcycling_level = true;
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.spatial_refinement_ratio = {1, 2, 2};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.spatial_refinement_ratio = {2, 1, 2};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.spatial_refinement_ratio = {2, 2, 1};
  assert(!permits_static_v1_level_one_creation(envelope));

  envelope = valid_level_one_creation();
  envelope.callback_time = 0.5;
  assert(!permits_static_v1_level_one_creation(envelope));
}

void test_static_v1_level_one_uses_the_parent_refinement_ratio() {
  assert(static_v1_parent_refinement_ratio_index(1) == 0);
}

void test_static_v1_stops_initial_regrid_after_one_l1_creation() {
  assert(should_stop_static_v1_initial_regrid_loop(true, 1, 2, 2));
  assert(!should_stop_static_v1_initial_regrid_loop(false, 1, 2, 2));
  assert(!should_stop_static_v1_initial_regrid_loop(true, 2, 2, 2));
  assert(!should_stop_static_v1_initial_regrid_loop(true, 1, 3, 2));
  assert(!should_stop_static_v1_initial_regrid_loop(true, 1, 2, 1));
}

void test_active_thorn_gate_is_an_exact_whitelist() {
  std::vector<std::string> active{
      "zlib",          "AMReX",       "TestODESolvers2", "CarpetX",
      "CarpetXRegrid", "Cactus",      "Arith",            "BoxInBox",
      "IOUtil",        "Loop",        "MPI",              "NSIMD",
      "ODESolvers",    "yaml_cpp"};
  assert(has_exact_static_v1_test_odesolvers2_active_thorns(active));
  active.push_back("ADMBase");
  assert(!has_exact_static_v1_test_odesolvers2_active_thorns(active));
  active.pop_back();
  active.pop_back();
  assert(!has_exact_static_v1_test_odesolvers2_active_thorns(active));
}

void test_active_level_history_rotates_once_before_tl0_capture() {
  std::vector<int> events;
  const int captured = cycle_then_capture_static_v1_level_state(
      [&] { events.push_back(1); },
      [&] {
        events.push_back(2);
        return 17;
      });
  assert(captured == 17);
  assert((events == std::vector<int>{1, 2}));
}

void test_candidate_mapping_uses_endpoint_iterations_without_preincrement() {
  const StepContext coarse{0, step_clock_t(0), step_clock_t(1), 4.0, 4.5,
                           SubcyclingODEMethod::rk4, true, 1};
  const StepContext fine_first{1, step_clock_t(0), step_clock_t(1, 2), 4.0,
                               4.25, SubcyclingODEMethod::rk4, false, 1};
  const StepContext fine_second{1, step_clock_t(1, 2), step_clock_t(1), 4.25,
                                4.5, SubcyclingODEMethod::rk4, false, 2};

  const RuntimeClockMetadata coarse_clock =
      candidate_runtime_clock(coarse, 0.5);
  const RuntimeClockMetadata fine_first_clock =
      candidate_runtime_clock(fine_first, 0.5);
  const RuntimeClockMetadata fine_second_clock =
      candidate_runtime_clock(fine_second, 0.5);

  assert(coarse_clock.iteration == 1 && coarse_clock.timefac == 1);
  assert(fine_first_clock.iteration == 1 && fine_first_clock.timefac == 2);
  assert(fine_second_clock.iteration == 2 && fine_second_clock.timefac == 2);
}

void test_sync_observer_order_is_exact() {
  constexpr auto order = static_v1_sync_observer_order();
  static_assert(order[0] ==
                StaticV1SyncObserverAction::cycle_synchronized_globals);
  static_assert(order[1] == StaticV1SyncObserverAction::postrestrict);
  static_assert(order[2] == StaticV1SyncObserverAction::poststep);
  static_assert(order[3] == StaticV1SyncObserverAction::analysis);
  static_assert(order[4] == StaticV1SyncObserverAction::output);
  static_assert(order[5] == StaticV1SyncObserverAction::checkpoint);
}

} // namespace

int main() {
  test_legacy_selection_is_a_strict_false_branch();
  test_static_policy_fails_closed_before_runtime_mutation();
  test_static_v1_level_one_creation_envelope_is_exact();
  test_static_v1_level_one_uses_the_parent_refinement_ratio();
  test_static_v1_stops_initial_regrid_after_one_l1_creation();
  test_active_thorn_gate_is_an_exact_whitelist();
  test_active_level_history_rotates_once_before_tl0_capture();
  test_candidate_mapping_uses_endpoint_iterations_without_preincrement();
  test_sync_observer_order_is_exact();
  std::cout << "subcycling_static_v1_top_adapter_unit: PASS\n";
}
