#ifndef ODESOLVERS_CACTUS_EXPLICIT_RK_OPERATIONS_HXX
#define ODESOLVERS_CACTUS_EXPLICIT_RK_OPERATIONS_HXX

#include "explicit_rk.hxx"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace ODESolvers {

template <typename Scalar, typename State>
struct CactusExplicitRKScratchHooks {
  std::function<void(State &, State &, const State &, const State &, Scalar)>
      restore_left;
  std::function<void(State &, const State &, Scalar)> restore_state;
  std::function<void(State &, const ExplicitRKStagePoint &, Scalar)>
      post_step_after_update;
  std::function<void(State &, State &, int, const ExplicitRKStagePoint &,
                     Scalar)>
      evaluate_rhs;
  std::function<State(State &, State &, const ExplicitRKStagePoint &, Scalar)>
      probe_endpoint_rhs;
  std::function<std::size_t()> rhs_evaluation_count;
};

template <typename State, typename PrepareInitial, typename EvaluateRHS,
          typename ValidateRHS, typename UpdateState, typename AccumulateRK4>
struct CactusExplicitRKOperations {
  using scalar_type = typename State::scalar_type;
  using state_type = State;
  using scratch_hooks_type = CactusExplicitRKScratchHooks<scalar_type, State>;

  State &var;
  State &rhs_component;
  PrepareInitial prepare_initial_callback;
  EvaluateRHS evaluate_rhs_callback;
  ValidateRHS validate_rhs_callback;
  UpdateState update_state_callback;
  AccumulateRK4 accumulate_rk4_callback;
  std::function<void(State &, const ExplicitRKStagePoint &, scalar_type)>
      stage_materialization_callback;
  std::function<void(const ExplicitRKStagePoint &, scalar_type)>
      stage_preparation_callback;
  std::function<void(scalar_type)> live_post_step_callback;
  scratch_hooks_type *scratch_hooks = nullptr;
  std::uint64_t state_generation_value = 0;
  bool loaded_rhs_available = false;
  scalar_type loaded_rhs_time =
      std::numeric_limits<scalar_type>::quiet_NaN();
  scalar_type current_stage_time =
      std::numeric_limits<scalar_type>::quiet_NaN();
  std::optional<ExplicitRKStagePoint> current_stage_point;
  bool current_stage_prepared = false;

  CactusExplicitRKOperations(State &var_, State &rhs_component_,
                             PrepareInitial prepare_initial_callback_,
                             EvaluateRHS evaluate_rhs_callback_,
                             ValidateRHS validate_rhs_callback_,
                             UpdateState update_state_callback_,
                             AccumulateRK4 accumulate_rk4_callback_)
      : var(var_), rhs_component(rhs_component_),
        prepare_initial_callback(std::move(prepare_initial_callback_)),
        evaluate_rhs_callback(std::move(evaluate_rhs_callback_)),
        validate_rhs_callback(std::move(validate_rhs_callback_)),
        update_state_callback(std::move(update_state_callback_)),
        accumulate_rk4_callback(std::move(accumulate_rk4_callback_)) {}

  const State &state() const noexcept { return var; }
  const State &rhs() const noexcept { return rhs_component; }

  State snapshot_state() const { return var.snapshot_state(); }
  State snapshot_rhs() const { return rhs_component.snapshot_rhs(); }

  void set_stage_point(const ExplicitRKStagePoint &point) {
    if (point.stage_index <= 0 || point.stage_count <= 0 ||
        point.stage_index > point.stage_count ||
        point.parent_fraction.denominator <= 0 ||
        point.parent_fraction.numerator < 0 ||
        point.parent_fraction.numerator > point.parent_fraction.denominator)
      throw std::invalid_argument("explicit RK stage point is invalid");
    switch (point.kind) {
    case ExplicitRKStageKind::primary:
    case ExplicitRKStageKind::fractional:
      break;
    case ExplicitRKStageKind::endpoint_probe:
      if (point.stage_index != 1 || point.stage_count != 1)
        throw std::invalid_argument(
            "endpoint probe must use the singleton stage index");
      break;
    default:
      throw std::invalid_argument("explicit RK stage kind is invalid");
    }
    current_stage_point = point;
    current_stage_prepared = false;
  }

  void prepare_initial(const scalar_type stage_time) {
    current_stage_time = stage_time;
    prepare_initial_callback(stage_time);
    prepare_active_stage();
  }

  void evaluate_rhs(const int stage) {
    if ((stage_materialization_callback || stage_preparation_callback) &&
        !current_stage_prepared)
      throw std::logic_error(
          "RHS evaluation reached an unprepared explicit RK stage");
    if (scratch_hooks != nullptr) {
      const auto &stage_point = require_current_stage_point();
      require_hook(static_cast<bool>(scratch_hooks->evaluate_rhs),
                   "scratch evaluate_rhs hook is missing");
      scratch_hooks->evaluate_rhs(var, rhs_component, stage, stage_point,
                                  current_stage_time);
    } else {
      evaluate_rhs_callback(stage);
    }
  }

  void validate_rhs(const int stage) { validate_rhs_callback(stage); }

  void update_state(
      const int update_index, const scalar_type stage_time,
      const scalar_type destination_scale,
      const LinearCombinationView<scalar_type, State> combination) {
    update_state_callback(update_index, stage_time, destination_scale,
                          combination);
    current_stage_time = stage_time;
    prepare_active_stage();
    if (scratch_hooks != nullptr) {
      const auto &stage_point = require_current_stage_point();
      require_hook(static_cast<bool>(scratch_hooks->post_step_after_update),
                   "scratch PostStep hook is missing");
      scratch_hooks->post_step_after_update(var, stage_point, stage_time);
    } else if (live_post_step_callback) {
      live_post_step_callback(stage_time);
    }
    ++state_generation_value;
    loaded_rhs_available = false;
  }

  void accumulate_rk4(State &accumulator, const scalar_type factor,
                      const State &increment) {
    accumulate_rk4_callback(accumulator, factor, increment);
  }

  void restore_left(const State &left_state, const State &left_rhs,
                    const scalar_type left_time) {
    require_scratch_hook_set();
    require_hook(static_cast<bool>(scratch_hooks->restore_left),
                 "scratch restore_left hook is missing");
    scratch_hooks->restore_left(var, rhs_component, left_state, left_rhs,
                                left_time);
    current_stage_time = left_time;
    ++state_generation_value;
    loaded_rhs_available = true;
    loaded_rhs_time = left_time;
    current_stage_point.reset();
    current_stage_prepared = false;
  }

  void restore_state(const State &state, const scalar_type state_time) {
    require_scratch_hook_set();
    require_hook(static_cast<bool>(scratch_hooks->restore_state),
                 "scratch restore_state hook is missing");
    scratch_hooks->restore_state(var, state, state_time);
    current_stage_time = state_time;
    ++state_generation_value;
    loaded_rhs_available = false;
    current_stage_point.reset();
    current_stage_prepared = false;
  }

  State probe_endpoint_rhs(const scalar_type state_time,
                           const ExplicitRKStagePoint &stage_point) {
    require_scratch_hook_set();
    set_stage_point(stage_point);
    current_stage_time = state_time;
    prepare_active_stage();
    require_hook(static_cast<bool>(scratch_hooks->post_step_after_update),
                 "scratch PostStep hook is missing");
    scratch_hooks->post_step_after_update(var, stage_point, state_time);
    require_hook(static_cast<bool>(scratch_hooks->probe_endpoint_rhs),
                 "scratch endpoint RHS hook is missing");
    return scratch_hooks->probe_endpoint_rhs(var, rhs_component, stage_point,
                                             state_time);
  }

  std::size_t rhs_evaluation_count() const {
    require_scratch_hook_set();
    require_hook(static_cast<bool>(scratch_hooks->rhs_evaluation_count),
                 "scratch RHS count hook is missing");
    return scratch_hooks->rhs_evaluation_count();
  }

  std::uint64_t state_generation() const noexcept {
    return state_generation_value;
  }

  LoadedRHSProvenance<scalar_type>
  loaded_rhs_provenance(const scalar_type left_time) const {
    if (!loaded_rhs_available || left_time != loaded_rhs_time)
      throw std::logic_error("no exact loaded left RHS is available");
    return {state_generation_value, left_time};
  }

  void validate_loaded_rhs_provenance(
      const LoadedRHSProvenance<scalar_type> &provenance) const {
    if (!loaded_rhs_available ||
        provenance.state_generation != state_generation_value ||
        provenance.left_time != loaded_rhs_time)
      throw std::logic_error("loaded RHS provenance is stale");
  }

  void consume_loaded_rhs(
      const LoadedRHSProvenance<scalar_type> &provenance) {
    validate_loaded_rhs_provenance(provenance);
    loaded_rhs_available = false;
  }

private:
  const ExplicitRKStagePoint &require_current_stage_point() const {
    if (!current_stage_point.has_value())
      throw std::logic_error("scratch operation has no exact stage point");
    return *current_stage_point;
  }

  void prepare_active_stage() {
    if (!stage_materialization_callback && !stage_preparation_callback)
      return;
    if (!current_stage_point.has_value())
      throw std::logic_error(
          "driver stage preparation has no exact stage point");
    if (scratch_hooks != nullptr && stage_preparation_callback &&
        !stage_materialization_callback)
      throw std::logic_error(
          "scratch driver stage preparation has no materialization callback");
    if (stage_materialization_callback)
      stage_materialization_callback(var, *current_stage_point,
                                     current_stage_time);
    if (stage_preparation_callback)
      stage_preparation_callback(*current_stage_point, current_stage_time);
    current_stage_prepared = true;
  }

  static void require_hook(const bool present, const char *const message) {
    if (!present)
      throw std::logic_error(message);
  }

  void require_scratch_hook_set() const {
    if (scratch_hooks == nullptr)
      throw std::logic_error("scratch operation requested without hooks");
  }
};

template <typename State, typename PrepareInitial, typename EvaluateRHS,
          typename ValidateRHS, typename UpdateState, typename AccumulateRK4>
CactusExplicitRKOperations(State &, State &, PrepareInitial, EvaluateRHS,
                           ValidateRHS, UpdateState, AccumulateRK4)
    -> CactusExplicitRKOperations<State, PrepareInitial, EvaluateRHS,
                                  ValidateRHS, UpdateState, AccumulateRK4>;

} // namespace ODESolvers

#endif // ODESOLVERS_CACTUS_EXPLICIT_RK_OPERATIONS_HXX
