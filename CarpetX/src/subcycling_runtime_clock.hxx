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

double static_v1_base_delta_time(double legacy_finest_delta_time,
                                 int level_count);

RuntimeClockMetadata candidate_runtime_clock(const StepContext &context,
                                             double base_delta_time);

RuntimeClockMetadata
full_sync_root_runtime_clock(std::uint64_t completed_epoch,
                             double synchronized_time,
                             double base_delta_time);

void validate_static_v1_clock_envelope(
    const StaticV1ClockEnvelope &envelope);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_RUNTIME_CLOCK_HXX
