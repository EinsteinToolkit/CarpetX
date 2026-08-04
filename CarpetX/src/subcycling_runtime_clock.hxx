#ifndef CARPETX_SUBCYCLING_RUNTIME_CLOCK_HXX
#define CARPETX_SUBCYCLING_RUNTIME_CLOCK_HXX

#include "subcycling_step_context.hxx"

#include <array>
#include <cstdint>

namespace CarpetX {

struct RuntimeClockMetadata {
  int iteration;
  double end_time;
  double base_delta_time;
  int timefac;
};

struct StaticV1ClockEnvelope {
  int level_count;
  int refinement_ratio;
  bool recovering;
  bool dynamic_regrid;
  step_clock_t synchronized_clock;
  std::array<step_clock_t, 2> level_clocks;
  std::array<step_clock_t, 2> level_delta_clocks;
  std::uint64_t completed_epoch;
  std::array<std::uint64_t, 2> accepted_steps;
};

struct StaticV1RecoverySeedEnvelope {
  int level_count;
  int refinement_ratio;
  bool strict_recovery;
  bool dynamic_regrid;
  int root_iteration;
  int root_timefac;
  double root_time;
  double base_delta_time;
  std::array<step_clock_t, 2> rebuilt_level_clocks;
  std::array<step_clock_t, 2> level_delta_clocks;
};

struct StaticV1StepperSeed {
  step_clock_t initial_clock;
  double initial_physical_time;
  std::uint64_t initial_epoch;
  std::array<std::uint64_t, 2> initial_accepted_steps;
};

double static_v1_base_delta_time(double legacy_finest_delta_time,
                                 int level_count);

RuntimeClockMetadata candidate_runtime_clock(const StepContext &context,
                                             double base_delta_time);

RuntimeClockMetadata
full_sync_root_runtime_clock(std::uint64_t completed_epoch,
                             double synchronized_time,
                             double base_delta_time);

StaticV1StepperSeed static_v1_recovery_seed(
    const StaticV1RecoverySeedEnvelope &envelope);

void validate_static_v1_clock_envelope(
    const StaticV1ClockEnvelope &envelope);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_RUNTIME_CLOCK_HXX
