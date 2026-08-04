#include "subcycling_runtime_clock.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CarpetX {
namespace {

void require_finite_positive_dt(const double delta_time) {
  if (!std::isfinite(delta_time) || !(delta_time > 0.0))
    throw std::invalid_argument("base delta time must be finite and positive");
}

void require_valid_clock(const step_clock_t clock) {
  if ((clock.den != 1 && clock.den != 2) || clock.num < 0 ||
      clock.num > std::numeric_limits<int>::max())
    throw std::invalid_argument(
        "static-v1 exact clock must be a bounded integer or half-integer");
}

bool equal_clock(const step_clock_t lhs, const step_clock_t rhs) {
  require_valid_clock(lhs);
  require_valid_clock(rhs);
  return lhs == rhs;
}

bool clock_equals_integer(const step_clock_t clock,
                          const std::uint64_t value) {
  require_valid_clock(clock);
  return clock.den == 1 && clock.num == static_cast<std::int64_t>(value);
}

bool clock_difference_is(const step_clock_t end,
                         const step_clock_t begin, const int denominator) {
  require_valid_clock(begin);
  require_valid_clock(end);
  if (denominator != 1 && denominator != 2)
    throw std::logic_error("unsupported static-v1 time refinement factor");
  return end - begin == step_clock_t(1, denominator);
}

void require_physical_step(const StepContext &context,
                           const double expected_delta_time) {
  if (!std::isfinite(context.begin_time) ||
      !std::isfinite(context.end_time))
    throw std::invalid_argument("step endpoint times must be finite");

  const long double begin = context.begin_time;
  const long double end = context.end_time;
  const long double expected = expected_delta_time;
  const long double actual = end - begin;
  const long double scale =
      std::max({1.0L, std::fabs(begin), std::fabs(end), std::fabs(expected)});
  const long double tolerance =
      32.0L * std::numeric_limits<double>::epsilon() * scale;
  if (!(actual > 0.0L) || std::fabs(actual - expected) > tolerance)
    throw std::invalid_argument(
        "physical step does not match base delta time and timefac");
}

void require_int_iteration(const std::uint64_t iteration) {
  if (iteration >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::overflow_error("runtime iteration exceeds CCTK_INT range");
}

} // namespace

double static_v1_base_delta_time(const double legacy_finest_delta_time,
                                 const int level_count) {
  if (level_count != 2)
    throw std::invalid_argument("static-v1 requires exactly two levels");
  if (!std::isfinite(legacy_finest_delta_time) ||
      !(legacy_finest_delta_time > 0.0))
    throw std::invalid_argument(
        "legacy finest delta time must be finite and positive");
  if (legacy_finest_delta_time >
      std::numeric_limits<double>::max() / 2.0)
    throw std::overflow_error("base delta time is not representable");
  const double base_delta_time = 2.0 * legacy_finest_delta_time;
  if (!std::isfinite(base_delta_time))
    throw std::overflow_error("base delta time is not finite");
  return base_delta_time;
}

RuntimeClockMetadata candidate_runtime_clock(const StepContext &context,
                                             const double base_delta_time) {
  require_finite_positive_dt(base_delta_time);
  if (context.level != 0 && context.level != 1)
    throw std::invalid_argument("static-v1 supports only levels zero and one");

  const int timefac = context.level == 0 ? 1 : 2;
  if (context.endpoint_accepted_step == 0)
    throw std::invalid_argument("accepted endpoint step must be positive");
  require_int_iteration(context.endpoint_accepted_step);
  require_valid_clock(context.begin_clock);
  require_valid_clock(context.end_clock);
  const step_clock_t expected_begin(
      static_cast<std::int64_t>(context.endpoint_accepted_step - 1), timefac);
  const step_clock_t expected_end(
      static_cast<std::int64_t>(context.endpoint_accepted_step), timefac);
  if (context.begin_clock != expected_begin || context.end_clock != expected_end)
    throw std::invalid_argument(
        "step exact clocks do not match the accepted-step count");
  if (!clock_difference_is(context.end_clock, context.begin_clock, timefac))
    throw std::invalid_argument(
        "step exact-clock delta does not match level timefac");

  const double physical_delta_time = base_delta_time / timefac;
  if (!std::isfinite(physical_delta_time) ||
      !(physical_delta_time > 0.0))
    throw std::overflow_error(
        "level physical delta time is not representable");
  require_physical_step(context, physical_delta_time);

  return RuntimeClockMetadata{
      static_cast<int>(context.endpoint_accepted_step), context.end_time,
      base_delta_time, timefac};
}

RuntimeClockMetadata
full_sync_root_runtime_clock(const std::uint64_t completed_epoch,
                             const double synchronized_time,
                             const double base_delta_time) {
  require_finite_positive_dt(base_delta_time);
  if (!std::isfinite(synchronized_time))
    throw std::invalid_argument("synchronized physical time must be finite");
  require_int_iteration(completed_epoch);
  return RuntimeClockMetadata{static_cast<int>(completed_epoch),
                              synchronized_time, base_delta_time, 1};
}

StaticV1StepperSeed static_v1_recovery_seed(
    const StaticV1RecoverySeedEnvelope &envelope) {
  if (envelope.level_count != 2 || envelope.refinement_ratio != 2)
    throw std::invalid_argument(
        "static-v1 recovery requires two factor-two levels");
  if (!envelope.strict_recovery)
    throw std::invalid_argument(
        "static-v1 recovery requires Cactus strict recovery mode");
  if (envelope.dynamic_regrid)
    throw std::invalid_argument(
        "static-v1 recovery does not support dynamic regridding");
  if (envelope.root_iteration < 0)
    throw std::invalid_argument(
        "static-v1 recovered root iteration must be non-negative");
  if (envelope.root_timefac != 1)
    throw std::invalid_argument(
        "static-v1 recovery requires synchronized root time metadata");
  if (!std::isfinite(envelope.root_time))
    throw std::invalid_argument(
        "static-v1 recovered root time must be finite");
  require_finite_positive_dt(envelope.base_delta_time);

  for (const auto clock : envelope.rebuilt_level_clocks)
    if (!equal_clock(clock, step_clock_t(0)))
      throw std::invalid_argument(
          "static-v1 recovery requires freshly rebuilt zero level clocks");
  if (!clock_difference_is(envelope.level_delta_clocks[0], step_clock_t(0),
                           1) ||
      !clock_difference_is(envelope.level_delta_clocks[1], step_clock_t(0),
                           2))
    throw std::invalid_argument(
        "static-v1 recovered level deltas must be one and one-half");

  const auto epoch = static_cast<std::uint64_t>(envelope.root_iteration);
  return StaticV1StepperSeed{step_clock_t(envelope.root_iteration),
                             envelope.root_time,
                             epoch,
                             {epoch, 2 * epoch}};
}

void validate_static_v1_clock_envelope(
    const StaticV1ClockEnvelope &envelope) {
  if (envelope.level_count != 2 || envelope.refinement_ratio != 2)
    throw std::invalid_argument(
        "static-v1 requires two levels with factor-two refinement");
  if (envelope.recovering)
    throw std::invalid_argument("static-v1 does not support recovery");
  if (envelope.dynamic_regrid)
    throw std::invalid_argument("static-v1 does not support dynamic regridding");

  require_int_iteration(envelope.completed_epoch);
  for (const auto accepted_step : envelope.accepted_steps)
    require_int_iteration(accepted_step);
  if (envelope.accepted_steps[0] != envelope.completed_epoch)
    throw std::invalid_argument(
        "coarse accepted-step count must equal the completed epoch");
  if (envelope.accepted_steps[0] >
          std::numeric_limits<std::uint64_t>::max() / 2 ||
      envelope.accepted_steps[1] != 2 * envelope.accepted_steps[0])
    throw std::invalid_argument(
        "fine accepted-step count must be twice the coarse count");

  if (!clock_equals_integer(envelope.synchronized_clock,
                            envelope.completed_epoch))
    throw std::invalid_argument(
        "synchronized exact clock must equal the completed epoch");
  if (!equal_clock(envelope.level_clocks[0], envelope.synchronized_clock) ||
      !equal_clock(envelope.level_clocks[1], envelope.synchronized_clock))
    throw std::invalid_argument(
        "both level clocks must equal the synchronized clock");
  if (!clock_difference_is(envelope.level_delta_clocks[0], step_clock_t(0),
                           1) ||
      !clock_difference_is(envelope.level_delta_clocks[1], step_clock_t(0),
                           2))
    throw std::invalid_argument(
        "static-v1 level deltas must be exactly one and one-half");
}

} // namespace CarpetX
