#include "subcycling_runtime_clock.hxx"

#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

using CarpetX::RuntimeClockMetadata;
using CarpetX::StaticV1ClockEnvelope;
using CarpetX::StaticV1RecoverySeedEnvelope;
using CarpetX::StaticV1StepperSeed;
using CarpetX::StepContext;
using CarpetX::SubcyclingODEMethod;
using CarpetX::candidate_runtime_clock;
using CarpetX::full_sync_root_runtime_clock;
using CarpetX::static_v1_base_delta_time;
using CarpetX::static_v1_recovery_seed;
using CarpetX::step_clock_t;
using CarpetX::validate_static_v1_clock_envelope;

template <typename Exception, typename Function>
void expect_throw(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const Exception &) {
    threw = true;
  }
  assert(threw);
}

void assert_metadata(const RuntimeClockMetadata &actual, const int iteration,
                     const double end_time, const double base_delta_time,
                     const int timefac) {
  assert(actual.iteration == iteration);
  assert(actual.end_time == end_time);
  assert(actual.base_delta_time == base_delta_time);
  assert(actual.timefac == timefac);
}

StepContext context(const int level, const step_clock_t begin_clock,
                    const step_clock_t end_clock, const double begin_time,
                    const double end_time,
                    const std::uint64_t endpoint_accepted_step) {
  return StepContext{level,
                     begin_clock,
                     end_clock,
                     begin_time,
                     end_time,
                     SubcyclingODEMethod::rk4,
                     level == 0,
                     endpoint_accepted_step};
}

void test_base_dt_and_candidate_metadata() {
  const double base_dt = static_v1_base_delta_time(0.125, 2);
  assert(base_dt == 0.25);

  assert_metadata(candidate_runtime_clock(
                      context(0, step_clock_t(0), step_clock_t(1), 10.0,
                              10.25, 1),
                      base_dt),
                  1, 10.25, 0.25, 1);
  assert_metadata(candidate_runtime_clock(
                      context(1, step_clock_t(0), step_clock_t(1, 2), 10.0,
                              10.125, 1),
                      base_dt),
                  1, 10.125, 0.25, 2);
  assert_metadata(candidate_runtime_clock(
                      context(1, step_clock_t(1, 2), step_clock_t(1),
                              10.125, 10.25, 2),
                      base_dt),
                  2, 10.25, 0.25, 2);
}

void test_full_sync_and_static_v1_envelope() {
  assert_metadata(full_sync_root_runtime_clock(7, 42.5, 0.25), 7, 42.5,
                  0.25, 1);

  validate_static_v1_clock_envelope(StaticV1ClockEnvelope{
      2,
      2,
      false,
      false,
      step_clock_t(7),
      {step_clock_t(7), step_clock_t(7)},
      {step_clock_t(1), step_clock_t(1, 2)},
      7,
      {7, 14}});
}

void test_base_dt_rejects_invalid_static_v1_inputs() {
  expect_throw<std::invalid_argument>(
      [] { static_cast<void>(static_v1_base_delta_time(0.125, 3)); });
  expect_throw<std::invalid_argument>(
      [] { static_cast<void>(static_v1_base_delta_time(0.0, 2)); });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(static_v1_base_delta_time(
        std::numeric_limits<double>::quiet_NaN(), 2));
  });
  expect_throw<std::overflow_error>([] {
    static_cast<void>(static_v1_base_delta_time(
        std::numeric_limits<double>::max(), 2));
  });
}

void test_candidate_rejects_invalid_level_clock_time_and_iteration() {
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(2, step_clock_t(0), step_clock_t(1, 4), 10.0, 10.0625, 1),
        0.25));
  });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(1, step_clock_t(0), step_clock_t(1), 10.0, 10.125, 1),
        0.25));
  });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(1, step_clock_t(0), step_clock_t(1, 2), 10.0, 10.2, 1),
        0.25));
  });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(1, step_clock_t(0), step_clock_t(1, 2), 10.0, 10.125, 2),
        0.25));
  });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(0, step_clock_t(0), step_clock_t(1), 10.0,
                std::numeric_limits<double>::infinity(), 1),
        0.25));
  });
  expect_throw<std::invalid_argument>([] {
    static_cast<void>(candidate_runtime_clock(
        context(0, step_clock_t(0), step_clock_t(1), 10.0, 10.25, 0),
        0.25));
  });
  expect_throw<std::overflow_error>([] {
    static_cast<void>(candidate_runtime_clock(
        context(0, step_clock_t(0), step_clock_t(1), 10.0, 10.25,
                static_cast<std::uint64_t>(INT_MAX) + 1),
        0.25));
  });
}

StaticV1ClockEnvelope valid_envelope() {
  return StaticV1ClockEnvelope{2,
                               2,
                               false,
                               false,
                               step_clock_t(7),
                               {step_clock_t(7), step_clock_t(7)},
                               {step_clock_t(1), step_clock_t(1, 2)},
                               7,
                               {7, 14}};
}

StaticV1RecoverySeedEnvelope valid_recovery_seed_envelope() {
  return StaticV1RecoverySeedEnvelope{
      2,
      2,
      true,
      false,
      2,
      1,
      0.025,
      0.0125,
      {step_clock_t(0), step_clock_t(0)},
      {step_clock_t(1), step_clock_t(1, 2)}};
}

void test_strict_full_sync_recovery_seed_at_epoch_two() {
  const StaticV1StepperSeed seed =
      static_v1_recovery_seed(valid_recovery_seed_envelope());

  assert(seed.initial_clock == step_clock_t(2));
  assert(seed.initial_physical_time == 0.025);
  assert(seed.initial_epoch == 2);
  assert((seed.initial_accepted_steps ==
          std::array<std::uint64_t, 2>{2, 4}));
}

void test_recovery_seed_rejects_desynchronized_rebuilt_state() {
  auto envelope = valid_recovery_seed_envelope();
  envelope.root_timefac = 2;
  expect_throw<std::invalid_argument>(
      [&] { static_cast<void>(static_v1_recovery_seed(envelope)); });

  envelope = valid_recovery_seed_envelope();
  envelope.rebuilt_level_clocks[1] = step_clock_t(1, 2);
  expect_throw<std::invalid_argument>(
      [&] { static_cast<void>(static_v1_recovery_seed(envelope)); });
}

void test_recovery_seed_rejects_malformed_or_nonfinite_state() {
  auto rejects = [](StaticV1RecoverySeedEnvelope envelope) {
    expect_throw<std::invalid_argument>(
        [&] { static_cast<void>(static_v1_recovery_seed(envelope)); });
  };

  auto envelope = valid_recovery_seed_envelope();
  envelope.level_count = 3;
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.refinement_ratio = 4;
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.strict_recovery = false;
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.dynamic_regrid = true;
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.root_iteration = -1;
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.level_delta_clocks[1] = step_clock_t(1, 4);
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.root_time = std::numeric_limits<double>::quiet_NaN();
  rejects(envelope);

  envelope = valid_recovery_seed_envelope();
  envelope.base_delta_time =
      std::numeric_limits<double>::infinity();
  rejects(envelope);
}

void test_static_v1_envelope_fails_closed() {
  auto envelope = valid_envelope();
  envelope.refinement_ratio = 4;
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.recovering = true;
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.dynamic_regrid = true;
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.synchronized_clock = step_clock_t(3, 2);
  envelope.level_clocks = {step_clock_t(3, 2), step_clock_t(3, 2)};
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.level_clocks[1] = step_clock_t(1);
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.level_delta_clocks[1] = step_clock_t(1, 4);
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.accepted_steps[0] = 6;
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.accepted_steps[1] = 13;
  expect_throw<std::invalid_argument>(
      [&] { validate_static_v1_clock_envelope(envelope); });

  envelope = valid_envelope();
  envelope.accepted_steps = {static_cast<std::uint64_t>(INT_MAX),
                             static_cast<std::uint64_t>(INT_MAX) * 2};
  envelope.completed_epoch = static_cast<std::uint64_t>(INT_MAX);
  expect_throw<std::overflow_error>(
      [&] { validate_static_v1_clock_envelope(envelope); });
}

} // namespace

int main() {
  test_base_dt_and_candidate_metadata();
  test_full_sync_and_static_v1_envelope();
  test_base_dt_rejects_invalid_static_v1_inputs();
  test_candidate_rejects_invalid_level_clock_time_and_iteration();
  test_strict_full_sync_recovery_seed_at_epoch_two();
  test_recovery_seed_rejects_desynchronized_rebuilt_state();
  test_recovery_seed_rejects_malformed_or_nonfinite_state();
  test_static_v1_envelope_fails_closed();
  std::cout << "subcycling_runtime_clock_unit: PASS\n";
}
