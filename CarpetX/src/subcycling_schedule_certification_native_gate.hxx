#ifndef CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_NATIVE_GATE_HXX
#define CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_NATIVE_GATE_HXX

#include "subcycling_schedule_certification.hxx"

#include <optional>

namespace CarpetX {

// Test-only inventory seam for the same-binary CPU-native gate. It returns
// canonical data only and can never publish executable schedule handles.
struct NativeGateScheduleObservationResult {
  std::optional<ScheduleCertificationExpectation> expectation;
  std::optional<ScheduleCertificationFailure> failure;
  explicit operator bool() const noexcept { return expectation.has_value(); }
};

[[nodiscard]] NativeGateScheduleObservationResult
observe_local_subcycling_schedules_for_native_gate(
    const ScheduleBuildProvenance &observed_provenance);

} // namespace CarpetX

#endif
