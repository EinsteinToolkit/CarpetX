#ifndef CARPETX_SUBCYCLING_STEP_CONTEXT_HXX
#define CARPETX_SUBCYCLING_STEP_CONTEXT_HXX

#include <rational.hxx>

#include <cstdint>

namespace CarpetX {

class ScratchStateTransaction;

using step_clock_t = Arith::rational<std::int64_t>;

enum class SubcyclingODEMethod { rk4, rkf78_order7, dp87_order8 };

enum class StageKind : std::uint8_t { primary, fractional, endpoint_probe };

struct StagePoint {
  StageKind kind;
  int stage_index;
  int stage_count;
  step_clock_t stage_fraction;
  double stage_time;
};

struct StepContext {
  int level;
  step_clock_t begin_clock;
  step_clock_t end_clock;
  double begin_time;
  double end_time;
  SubcyclingODEMethod method;
  bool require_dense_output{false};
  std::uint64_t endpoint_accepted_step{0};
};

class StagePreparer {
public:
  virtual ~StagePreparer() = default;
  virtual void prepare_stage(const StepContext &context,
                             const StagePoint &stage_point) = 0;
};

class ScopedStepContext {
public:
  ScopedStepContext(const StepContext &context, StagePreparer &preparer);
  ScopedStepContext(const StepContext &context, StagePreparer &preparer,
                    ScratchStateTransaction *transaction);
  ~ScopedStepContext();
  ScopedStepContext(const ScopedStepContext &) = delete;
  ScopedStepContext &operator=(const ScopedStepContext &) = delete;

private:
  ScratchStateTransaction *prior_transaction_{nullptr};
};

const StepContext *current_step_context() noexcept;
bool step_context_active() noexcept;
ScratchStateTransaction *current_scratch_state_transaction() noexcept;
void prepare_stage(const StagePoint &stage_point, SubcyclingODEMethod method);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_STEP_CONTEXT_HXX
