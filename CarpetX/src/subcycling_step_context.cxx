#include "subcycling_step_context.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace CarpetX {
namespace {

struct ActiveStepContext {
  const StepContext *context{nullptr};
  StagePreparer *preparer{nullptr};
};

thread_local ActiveStepContext active_step_context;
thread_local ScratchStateTransaction *active_scratch_transaction{nullptr};

} // namespace

ScopedStepContext::ScopedStepContext(const StepContext &context,
                                     StagePreparer &preparer)
    : ScopedStepContext(context, preparer, nullptr) {}

ScopedStepContext::ScopedStepContext(const StepContext &context,
                                     StagePreparer &preparer,
                                     ScratchStateTransaction *transaction)
    : prior_transaction_(active_scratch_transaction) {
  if (active_step_context.context != nullptr)
    throw std::logic_error("a StepContext is already active on this thread");
  if (context.level < 0)
    throw std::invalid_argument("StepContext level must be non-negative");
  if (!(context.begin_clock < context.end_clock))
    throw std::invalid_argument("StepContext clocks must increase exactly");
  if (!std::isfinite(context.begin_time) ||
      !std::isfinite(context.end_time) ||
      !(context.begin_time < context.end_time))
    throw std::invalid_argument("StepContext times must be finite and increase");

  active_step_context = ActiveStepContext{&context, &preparer};
  active_scratch_transaction = transaction;
}

ScopedStepContext::~ScopedStepContext() {
  active_step_context = ActiveStepContext{};
  active_scratch_transaction = prior_transaction_;
}

const StepContext *current_step_context() noexcept {
  return active_step_context.context;
}

bool step_context_active() noexcept {
  return active_step_context.context != nullptr;
}

ScratchStateTransaction *current_scratch_state_transaction() noexcept {
  return active_scratch_transaction;
}

void prepare_stage(const StagePoint &stage_point,
                   const SubcyclingODEMethod method) {
  const auto context = active_step_context.context;
  if (context == nullptr)
    return;
  if (method != context->method)
    throw std::invalid_argument("ODE method does not match active StepContext");
  if (stage_point.stage_index <= 0 || stage_point.stage_count <= 0 ||
      stage_point.stage_index > stage_point.stage_count)
    throw std::invalid_argument("stage index/count are invalid");
  switch (stage_point.kind) {
  case StageKind::primary:
  case StageKind::fractional:
    break;
  case StageKind::endpoint_probe:
    if (stage_point.stage_index != 1 || stage_point.stage_count != 1)
      throw std::invalid_argument(
          "endpoint probe must use the singleton stage index");
    break;
  default:
    throw std::invalid_argument("stage kind is invalid");
  }
  if (stage_point.stage_fraction < step_clock_t(0) ||
      stage_point.stage_fraction > step_clock_t(1))
    throw std::out_of_range("exact stage fraction is outside [0,1]");
  if (!std::isfinite(stage_point.stage_time))
    throw std::invalid_argument("stage time must be finite");

  if (stage_point.kind != StageKind::primary &&
      active_scratch_transaction == nullptr)
    throw std::logic_error(
        "auxiliary stage requires the active scratch transaction");

  const auto time_scale =
      std::max({1.0, std::abs(context->begin_time),
                std::abs(context->end_time),
                std::abs(stage_point.stage_time)});
  const auto tolerance =
      32.0 * std::numeric_limits<double>::epsilon() * time_scale;
  if (stage_point.stage_time < context->begin_time - tolerance ||
      stage_point.stage_time > context->end_time + tolerance)
    throw std::out_of_range("stage time is outside the active StepContext");

  const auto fraction = static_cast<double>(stage_point.stage_fraction);
  const auto canonical_stage_time =
      context->begin_time +
      fraction * (context->end_time - context->begin_time);
  if (!std::isfinite(canonical_stage_time) ||
      std::abs(stage_point.stage_time - canonical_stage_time) > tolerance)
    throw std::invalid_argument(
        "stage time does not match the exact stage fraction");

  active_step_context.preparer->prepare_stage(*context, stage_point);
}

} // namespace CarpetX
