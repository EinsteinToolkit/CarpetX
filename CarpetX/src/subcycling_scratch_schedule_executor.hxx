#ifndef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_HXX
#define CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_HXX

#include "subcycling_step_context.hxx"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

namespace CarpetX {

enum class ScratchSchedulePhase : std::uint8_t { post_step, rhs };

struct ScratchStageCoordinates {
  int level_iteration;
  double stage_time;
  step_clock_t base_delta_clock;
  double base_delta_time;
  int time_refinement_factor;
};

struct ScratchScheduleExecutionReceipt {
  ScratchSchedulePhase phase;
  std::size_t routine_count;
  std::size_t local_call_count;
};

enum class ScratchScheduleExecutionErrorCode : std::uint8_t {
  invalid_step_context,
  invalid_coordinates,
  unsupported_configuration,
  stale_binding,
  binding_busy,
  binding_faulted,
  reentrant_execution,
  unresolved_access,
  invalid_read,
  call_failed,
  already_faulted,
  internal_failure,
};

class ScratchScheduleExecutionError final : public std::runtime_error {
public:
  ScratchScheduleExecutionErrorCode code() const noexcept;
  std::optional<ScratchSchedulePhase> phase() const noexcept;
  std::optional<std::size_t> routine_ordinal() const noexcept;
  std::optional<std::size_t> context_ordinal() const noexcept;

private:
  ScratchScheduleExecutionError(
      ScratchScheduleExecutionErrorCode code,
      std::optional<ScratchSchedulePhase> phase,
      std::optional<std::size_t> routine_ordinal,
      std::optional<std::size_t> context_ordinal,
      const std::string &message);
  ScratchScheduleExecutionErrorCode code_;
  std::optional<ScratchSchedulePhase> phase_;
  std::optional<std::size_t> routine_ordinal_;
  std::optional<std::size_t> context_ordinal_;
  friend class CertifiedScratchStageExecutor;
};

class CertifiedLocalScheduleRegistry;
class ScratchLocalLevelBinding;

class CertifiedScratchStageExecutor final {
public:
  [[nodiscard]] static std::unique_ptr<CertifiedScratchStageExecutor>
  create(CertifiedLocalScheduleRegistry &registry,
         ScratchLocalLevelBinding &binding);

  ~CertifiedScratchStageExecutor();
  CertifiedScratchStageExecutor(const CertifiedScratchStageExecutor &) =
      delete;
  CertifiedScratchStageExecutor &
  operator=(const CertifiedScratchStageExecutor &) = delete;
  CertifiedScratchStageExecutor(CertifiedScratchStageExecutor &&) = delete;
  CertifiedScratchStageExecutor &
  operator=(CertifiedScratchStageExecutor &&) = delete;

  ScratchScheduleExecutionReceipt
  execute_post_step(const StepContext &step,
                    const ScratchStageCoordinates &coordinates);
  ScratchScheduleExecutionReceipt
  execute_rhs(const StepContext &step,
              const ScratchStageCoordinates &coordinates);
  bool faulted() const noexcept;

private:
  struct Storage;
  CertifiedScratchStageExecutor(CertifiedLocalScheduleRegistry &registry,
                                ScratchLocalLevelBinding &binding,
                                std::unique_ptr<Storage> storage);
  [[noreturn]] static void
  raise(ScratchScheduleExecutionErrorCode code,
        std::optional<ScratchSchedulePhase> phase,
        std::optional<std::size_t> routine_ordinal,
        std::optional<std::size_t> context_ordinal,
        const std::string &message);
  ScratchScheduleExecutionReceipt
  execute(ScratchSchedulePhase phase, const StepContext &step,
          const ScratchStageCoordinates &coordinates);

  CertifiedLocalScheduleRegistry &registry_;
  ScratchLocalLevelBinding &binding_;
  std::unique_ptr<Storage> storage_;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_HXX
