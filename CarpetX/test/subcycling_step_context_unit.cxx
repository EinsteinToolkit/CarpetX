#include "subcycling_step_context.hxx"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {

using CarpetX::ScopedStepContext;
using CarpetX::StagePreparer;
using CarpetX::StageKind;
using CarpetX::StagePoint;
using CarpetX::StepContext;
using CarpetX::SubcyclingODEMethod;
using CarpetX::current_step_context;
using CarpetX::prepare_stage;
using CarpetX::step_clock_t;
using CarpetX::step_context_active;

static_assert(
    std::is_same_v<decltype(current_step_context()), const StepContext *>);

class RecordingPreparer final : public StagePreparer {
public:
  std::vector<const StepContext *> contexts;
  std::vector<StagePoint> stage_points;
  bool throw_on_prepare{false};

  void prepare_stage(const StepContext &context,
                     const StagePoint &stage_point) override {
    contexts.push_back(&context);
    stage_points.push_back(stage_point);
    if (throw_on_prepare)
      throw std::runtime_error("requested prepare failure");
  }
};

StepContext valid_context() {
  return StepContext{2, step_clock_t(5, 8), step_clock_t(7, 8), 1.25, 2.5,
                     SubcyclingODEMethod::rk4};
}

StagePoint point(const StageKind kind, const int stage_index,
                 const int stage_count, const step_clock_t fraction,
                 const double time) {
  return StagePoint{kind, stage_index, stage_count, fraction, time};
}

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

void test_absent_context_is_a_noop() {
  assert(!step_context_active());
  assert(current_step_context() == nullptr);
  prepare_stage(point(StageKind::endpoint_probe, 1, 1, step_clock_t(1),
                      std::numeric_limits<double>::quiet_NaN()),
                SubcyclingODEMethod::dp87_order8);
  assert(!step_context_active());
}

void test_scope_exposes_const_exact_context_and_records_any_in_range_order() {
  RecordingPreparer preparer;
  const auto context = valid_context();

  {
    ScopedStepContext scope(context, preparer);
    assert(step_context_active());
    assert(current_step_context() == &context);
    assert(current_step_context()->level == 2);
    assert(current_step_context()->begin_clock == step_clock_t(5, 8));
    assert(current_step_context()->end_clock == step_clock_t(7, 8));

    prepare_stage(point(StageKind::primary, 1, 4, step_clock_t(0), 1.25),
                  SubcyclingODEMethod::rk4);
    prepare_stage(point(StageKind::primary, 4, 4, step_clock_t(1), 2.5),
                  SubcyclingODEMethod::rk4);
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(2, 5), 1.75),
                  SubcyclingODEMethod::rk4);
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(2, 5), 1.75),
                  SubcyclingODEMethod::rk4);
  }

  assert(!step_context_active());
  assert(current_step_context() == nullptr);
  assert(preparer.contexts ==
         std::vector<const StepContext *>({&context, &context, &context,
                                           &context}));
  assert(preparer.stage_points.size() == 4);
  assert(preparer.stage_points[0].stage_fraction == step_clock_t(0));
  assert(preparer.stage_points[1].stage_fraction == step_clock_t(1));
  assert(preparer.stage_points[2].stage_fraction == step_clock_t(2, 5));
  assert(preparer.stage_points[3].stage_fraction == step_clock_t(2, 5));
}

void test_active_prepare_rejects_method_and_range_violations() {
  RecordingPreparer preparer;
  const auto context = valid_context();
  ScopedStepContext scope(context, preparer);

  expect_throw<std::invalid_argument>([&] {
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(1, 5), 1.5),
                  SubcyclingODEMethod::rkf78_order7);
  });
  expect_throw<std::invalid_argument>([&] {
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(1, 5),
                        std::numeric_limits<double>::quiet_NaN()),
                  SubcyclingODEMethod::rk4);
  });
  expect_throw<std::out_of_range>(
      [&] {
        prepare_stage(point(StageKind::primary, 1, 4, step_clock_t(0), 1.0),
                      SubcyclingODEMethod::rk4);
      });
  expect_throw<std::out_of_range>(
      [&] {
        prepare_stage(point(StageKind::primary, 4, 4, step_clock_t(1), 3.0),
                      SubcyclingODEMethod::rk4);
      });
  expect_throw<std::invalid_argument>([&] {
    prepare_stage(point(StageKind::primary, 0, 4, step_clock_t(0), 1.25),
                  SubcyclingODEMethod::rk4);
  });
  expect_throw<std::invalid_argument>([&] {
    prepare_stage(point(StageKind::primary, 2, 1, step_clock_t(1, 2), 1.875),
                  SubcyclingODEMethod::rk4);
  });
  expect_throw<std::out_of_range>([&] {
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(-1, 8), 1.25),
                  SubcyclingODEMethod::rk4);
  });
  expect_throw<std::invalid_argument>([&] {
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(1, 2), 1.5),
                  SubcyclingODEMethod::rk4);
  });
  assert(preparer.stage_points.empty());
}

void test_installation_rejects_invalid_contexts_and_nesting() {
  RecordingPreparer preparer;

  auto check_invalid = [&](const StepContext &context) {
    expect_throw<std::invalid_argument>(
        [&] { ScopedStepContext scope(context, preparer); });
    assert(!step_context_active());
  };

  auto context = valid_context();
  context.level = -1;
  check_invalid(context);

  context = valid_context();
  context.end_clock = context.begin_clock;
  check_invalid(context);

  context = valid_context();
  context.end_clock = step_clock_t(1, 2);
  check_invalid(context);

  context = valid_context();
  context.begin_time = std::numeric_limits<double>::quiet_NaN();
  check_invalid(context);

  context = valid_context();
  context.end_time = std::numeric_limits<double>::infinity();
  check_invalid(context);

  context = valid_context();
  context.end_time = context.begin_time;
  check_invalid(context);

  context = valid_context();
  context.end_time = 1.0;
  check_invalid(context);

  const auto outer_context = valid_context();
  const auto inner_context = valid_context();
  {
    ScopedStepContext outer(outer_context, preparer);
    expect_throw<std::logic_error>(
        [&] { ScopedStepContext inner(inner_context, preparer); });
    assert(current_step_context() == &outer_context);
  }
  assert(!step_context_active());
}

void test_unwinding_restores_absence_and_propagates_callback_exceptions() {
  RecordingPreparer preparer;
  const auto context = valid_context();

  expect_throw<std::runtime_error>([&] {
    ScopedStepContext scope(context, preparer);
    throw std::runtime_error("requested scope unwind");
  });
  assert(!step_context_active());

  preparer.throw_on_prepare = true;
  expect_throw<std::runtime_error>([&] {
    ScopedStepContext scope(context, preparer);
    prepare_stage(point(StageKind::primary, 2, 4, step_clock_t(1, 5), 1.5),
                  SubcyclingODEMethod::rk4);
  });
  assert(!step_context_active());
  assert(current_step_context() == nullptr);
}

void test_auxiliary_stage_requires_and_preserves_transaction_binding() {
  RecordingPreparer preparer;
  const auto context = valid_context();
  auto *const transaction = reinterpret_cast<CarpetX::ScratchStateTransaction *>(
      static_cast<std::uintptr_t>(0x10));
  {
    ScopedStepContext scope(context, preparer, transaction);
    prepare_stage(
        point(StageKind::fractional, 2, 4, step_clock_t(1, 4), 1.5625),
        SubcyclingODEMethod::rk4);
    prepare_stage(
        point(StageKind::endpoint_probe, 1, 1, step_clock_t(3, 4), 2.1875),
        SubcyclingODEMethod::rk4);
  }
  assert(preparer.stage_points.size() == 2);

  ScopedStepContext scope(context, preparer);
  expect_throw<std::logic_error>([&] {
    prepare_stage(
        point(StageKind::fractional, 2, 4, step_clock_t(1, 4), 1.5625),
        SubcyclingODEMethod::rk4);
  });
}

} // namespace

int main() {
  test_absent_context_is_a_noop();
  test_scope_exposes_const_exact_context_and_records_any_in_range_order();
  test_active_prepare_rejects_method_and_range_violations();
  test_installation_rejects_invalid_contexts_and_nesting();
  test_unwinding_restores_absence_and_propagates_callback_exceptions();
  test_auxiliary_stage_requires_and_preserves_transaction_binding();
  std::cout << "Subcycling StepContext tests passed\n";
  return 0;
}
