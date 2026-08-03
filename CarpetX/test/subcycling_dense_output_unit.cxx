#include "subcycling_dense_output.hxx"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using CarpetX::DenseCapability;
using CarpetX::DenseInterval;
using CarpetX::DenseIntervalBuilder;
using CarpetX::DenseIntervalId;
using CarpetX::DenseOutputProvider;
using CarpetX::DenseOutputRegistry;
using CarpetX::DenseStateVector;
using CarpetX::SubcyclingODEMethod;
using CarpetX::TableauFingerprint;
using CarpetX::step_clock_t;

static_assert(!std::is_copy_constructible_v<DenseIntervalBuilder>);
static_assert(!std::is_copy_assignable_v<DenseIntervalBuilder>);
static_assert(std::is_move_constructible_v<DenseIntervalBuilder>);
static_assert(std::is_move_assignable_v<DenseIntervalBuilder>);
static_assert(std::is_final_v<DenseOutputProvider>);

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

TableauFingerprint fingerprint(const std::uint8_t seed = 1) {
  TableauFingerprint result{};
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = static_cast<std::uint8_t>(seed + i);
  return result;
}

DenseCapability capability(const int degree = 4) {
  return DenseCapability{SubcyclingODEMethod::rk4,
                         fingerprint(),
                         degree,
                         degree,
                         4,
                         1,
                         degree + 1,
                         true,
                         true};
}

DenseIntervalId interval_id() {
  return DenseIntervalId{2,
                         step_clock_t(5, 8),
                         step_clock_t(7, 8),
                         1.25,
                         2.5,
                         SubcyclingODEMethod::rk4,
                         fingerprint()};
}

class ScalarState final : public DenseStateVector {
public:
  explicit ScalarState(const double value = 0.0, const int shape = 1,
                       int *destructor_count = nullptr)
      : value(value), shape(shape), destructor_count(destructor_count) {}

  ~ScalarState() override {
    if (destructor_count != nullptr)
      ++*destructor_count;
  }

  bool compatible(const DenseStateVector &other) const noexcept override {
    const auto *scalar = dynamic_cast<const ScalarState *>(&other);
    return scalar != nullptr && scalar->shape == shape;
  }

  void copy_from(const DenseStateVector &other) override {
    const auto *scalar = dynamic_cast<const ScalarState *>(&other);
    if (scalar == nullptr || scalar->shape != shape)
      throw std::invalid_argument("incompatible scalar copy");
    value = scalar->value;
    ++copy_calls;
  }

  void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const DenseStateVector *> &sources) override {
    if (weights.size() != sources.size())
      throw std::invalid_argument("weight/source mismatch");
    double result = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i) {
      const auto *scalar = dynamic_cast<const ScalarState *>(sources[i]);
      if (scalar == nullptr || scalar->shape != shape)
        throw std::invalid_argument("incompatible scalar combination");
      result += weights[i] * scalar->value;
    }
    value = result;
    ++linear_combination_calls;
  }

  double value;
  int shape;
  int *destructor_count;
  int copy_calls{0};
  int linear_combination_calls{0};
};

std::unique_ptr<DenseStateVector> scalar_control(const double value,
                                                const int shape = 1,
                                                int *destroyed = nullptr) {
  return std::make_unique<ScalarState>(value, shape, destroyed);
}

double binomial(const int n, const int k) {
  if (k < 0 || k > n)
    return 0.0;
  double result = 1.0;
  for (int i = 1; i <= k; ++i)
    result *= static_cast<double>(n - k + i) / static_cast<double>(i);
  return result;
}

void test_empty_registry_and_capability_rejections() {
  DenseOutputRegistry registry;
  expect_throw<std::out_of_range>([&] {
    registry.require(SubcyclingODEMethod::rk4, fingerprint());
  });

  auto invalid = capability();
  invalid.verified = false;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.dense_uniform_order = invalid.endpoint_order - 1;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.endpoint_order = 0;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.arbitrary_theta = false;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.tableau_fingerprint = TableauFingerprint{};
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.stage_count = 0;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.extra_rhs_evaluations = -1;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.persistent_vector_count = 0;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  invalid = capability();
  invalid.persistent_vector_count = invalid.dense_uniform_order;
  expect_throw<std::invalid_argument>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(invalid));
  });

  expect_throw<std::invalid_argument>(
      [&] { registry.register_provider(nullptr); });
}

void test_exact_lookup_and_duplicate_method_rejection() {
  DenseOutputRegistry registry;
  auto caller_capability = capability();
  auto provider =
      std::make_shared<DenseOutputProvider>(caller_capability);
  registry.register_provider(provider);
  assert(registry.require(SubcyclingODEMethod::rk4, fingerprint()) == provider);

  caller_capability.method = SubcyclingODEMethod::dp87_order8;
  caller_capability.tableau_fingerprint = fingerprint(9);
  caller_capability.verified = false;
  assert(provider->capability().method == SubcyclingODEMethod::rk4);
  assert(provider->capability().tableau_fingerprint == fingerprint());
  assert(provider->capability().verified);
  assert(provider->begin_interval(interval_id()) != nullptr);

  expect_throw<std::out_of_range>([&] {
    registry.require(SubcyclingODEMethod::rk4, fingerprint(2));
  });
  expect_throw<std::out_of_range>([&] {
    registry.require(SubcyclingODEMethod::dp87_order8, fingerprint());
  });

  auto duplicate = capability(5);
  duplicate.tableau_fingerprint = fingerprint(7);
  expect_throw<std::logic_error>([&] {
    registry.register_provider(
        std::make_shared<DenseOutputProvider>(duplicate));
  });
}

void test_interval_id_validation() {
  const auto provider = DenseOutputProvider(capability());
  auto invalid = interval_id();

  invalid.level = -1;
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.end_clock = invalid.begin_clock;
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.end_clock = step_clock_t(1, 2);
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.begin_time = std::numeric_limits<double>::quiet_NaN();
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.end_time = std::numeric_limits<double>::infinity();
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.end_time = invalid.begin_time;
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.method = SubcyclingODEMethod::rkf78_order7;
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
  invalid = interval_id();
  invalid.tableau_fingerprint = fingerprint(2);
  expect_throw<std::invalid_argument>([&] { provider.begin_interval(invalid); });
}

void test_builder_state_machine_and_owned_lifetime() {
  const auto provider = DenseOutputProvider(capability());
  auto builder = provider.begin_interval(interval_id());

  expect_throw<std::invalid_argument>(
      [&] { builder->add_control(nullptr); });
  builder->add_control(scalar_control(0.0));
  expect_throw<std::invalid_argument>(
      [&] { builder->add_control(scalar_control(1.0, 2)); });
  expect_throw<std::logic_error>([&] { builder->seal(); });

  builder->add_control(scalar_control(0.25));
  builder->add_control(scalar_control(0.5));
  builder->add_control(scalar_control(0.75));
  builder->add_control(scalar_control(1.0));
  expect_throw<std::length_error>(
      [&] { builder->add_control(scalar_control(2.0)); });

  auto interval = builder->seal();
  static_assert(std::is_same_v<decltype(interval),
                               std::shared_ptr<const CarpetX::DenseInterval>>);
  expect_throw<std::logic_error>(
      [&] { builder->add_control(scalar_control(3.0)); });
  expect_throw<std::logic_error>([&] { builder->seal(); });

  auto movable = provider.begin_interval(interval_id());
  DenseIntervalBuilder moved = std::move(*movable);
  expect_throw<std::logic_error>(
      [&] { movable->add_control(scalar_control(0.0)); });
  moved.add_control(scalar_control(0.0));

  int destroyed = 0;
  std::shared_ptr<const CarpetX::DenseInterval> retained;
  {
    auto lifetime_builder = provider.begin_interval(interval_id());
    for (int i = 0; i < 5; ++i)
      lifetime_builder->add_control(scalar_control(i, 1, &destroyed));
    retained = lifetime_builder->seal();
    assert(destroyed == 0);
    auto second_reference = retained;
    retained.reset();
    assert(destroyed == 0);
    retained = std::move(second_reference);
  }
  assert(destroyed == 0);
  retained.reset();
  assert(destroyed == 5);
}

void test_move_assignment_rejects_inactive_or_nonempty_builders() {
  const auto provider = DenseOutputProvider(capability());

  auto sealed_target = provider.begin_interval(interval_id());
  for (int i = 0; i < 5; ++i)
    sealed_target->add_control(scalar_control(i));
  const auto retained_target = sealed_target->seal();
  auto active_source = provider.begin_interval(interval_id());
  expect_throw<std::logic_error>(
      [&] { *sealed_target = std::move(*active_source); });
  active_source->add_control(scalar_control(0.0));
  assert(retained_target->control_count() == 5);

  auto active_target = provider.begin_interval(interval_id());
  auto sealed_source = provider.begin_interval(interval_id());
  for (int i = 0; i < 5; ++i)
    sealed_source->add_control(scalar_control(i));
  const auto retained_source = sealed_source->seal();
  expect_throw<std::logic_error>(
      [&] { *active_target = std::move(*sealed_source); });
  active_target->add_control(scalar_control(0.0));
  assert(retained_source->control_count() == 5);

  auto moved_target = provider.begin_interval(interval_id());
  DenseIntervalBuilder target_owner = std::move(*moved_target);
  auto second_active_source = provider.begin_interval(interval_id());
  expect_throw<std::logic_error>(
      [&] { *moved_target = std::move(*second_active_source); });
  second_active_source->add_control(scalar_control(0.0));
  target_owner.add_control(scalar_control(0.0));

  auto second_active_target = provider.begin_interval(interval_id());
  auto moved_source = provider.begin_interval(interval_id());
  DenseIntervalBuilder source_owner = std::move(*moved_source);
  expect_throw<std::logic_error>(
      [&] { *second_active_target = std::move(*moved_source); });
  second_active_target->add_control(scalar_control(0.0));
  source_owner.add_control(scalar_control(0.0));

  auto nonempty_target = provider.begin_interval(interval_id());
  nonempty_target->add_control(scalar_control(0.0));
  auto third_active_source = provider.begin_interval(interval_id());
  expect_throw<std::logic_error>(
      [&] { *nonempty_target = std::move(*third_active_source); });
  third_active_source->add_control(scalar_control(0.0));

  auto empty_target = provider.begin_interval(interval_id());
  auto populated_source = provider.begin_interval(interval_id());
  populated_source->add_control(scalar_control(0.0));
  *empty_target = std::move(*populated_source);
  empty_target->add_control(scalar_control(0.25));
  expect_throw<std::logic_error>(
      [&] { populated_source->add_control(scalar_control(1.0)); });
}

void test_bernstein_reproduces_monomials() {
  for (const int degree : {4, 7, 8}) {
    auto cap = capability(degree);
    DenseOutputProvider provider(cap);
    for (int power = 0; power <= degree; ++power) {
      auto builder = provider.begin_interval(interval_id());
      for (int k = 0; k <= degree; ++k) {
        const auto coefficient =
            k < power ? 0.0 : binomial(k, power) / binomial(degree, power);
        builder->add_control(scalar_control(coefficient));
      }
      const auto interval = builder->seal();
      for (const double theta : {0.0, 0.13, 0.5, 0.87, 1.0}) {
        ScalarState result;
        interval->evaluate(theta, result);
        assert(std::abs(result.value - std::pow(theta, power)) < 2.0e-12);
      }
    }
  }
}

void test_endpoint_and_interior_dispatch_and_theta_validation() {
  DenseOutputProvider provider(capability());
  auto builder = provider.begin_interval(interval_id());
  for (int i = 0; i < 5; ++i)
    builder->add_control(scalar_control(static_cast<double>(i)));
  const auto interval = builder->seal();

  ScalarState destination;
  interval->evaluate(0.0, destination);
  assert(destination.value == 0.0);
  assert(destination.copy_calls == 1);
  assert(destination.linear_combination_calls == 0);
  interval->evaluate(1.0, destination);
  assert(destination.value == 4.0);
  assert(destination.copy_calls == 2);
  assert(destination.linear_combination_calls == 0);
  interval->evaluate(0.5, destination);
  assert(destination.value == 2.0);
  assert(destination.copy_calls == 2);
  assert(destination.linear_combination_calls == 1);

  expect_throw<std::invalid_argument>([&] {
    interval->evaluate(std::numeric_limits<double>::quiet_NaN(), destination);
  });
  expect_throw<std::out_of_range>([&] { interval->evaluate(-0.01, destination); });
  expect_throw<std::out_of_range>([&] { interval->evaluate(1.01, destination); });

  ScalarState incompatible(0.0, 2);
  expect_throw<std::invalid_argument>(
      [&] { interval->evaluate(0.5, incompatible); });
}

void test_exact_rational_theta_overload_and_validation() {
  const auto exact_overload = static_cast<void (DenseInterval::*)(
      step_clock_t, DenseStateVector &) const>(&DenseInterval::evaluate);
  assert(exact_overload != nullptr);

  DenseOutputProvider provider(capability());
  auto builder = provider.begin_interval(interval_id());
  for (int i = 0; i < 5; ++i)
    builder->add_control(scalar_control(static_cast<double>(i)));
  const auto interval = builder->seal();

  ScalarState destination;
  interval->evaluate(step_clock_t(0), destination);
  assert(destination.value == 0.0);
  interval->evaluate(step_clock_t(1), destination);
  assert(destination.value == 4.0);
  interval->evaluate(step_clock_t(1, 4), destination);
  assert(std::abs(destination.value - 1.0) < 1.0e-12);

  expect_throw<std::out_of_range>(
      [&] { interval->evaluate(step_clock_t(-1, 4), destination); });
  expect_throw<std::out_of_range>(
      [&] { interval->evaluate(step_clock_t(5, 4), destination); });
  const step_clock_t zero_denominator(
      1, 0, step_clock_t::no_normalize{});
  expect_throw<std::invalid_argument>(
      [&] { interval->evaluate(zero_denominator, destination); });
}

} // namespace

int main() {
  test_empty_registry_and_capability_rejections();
  test_exact_lookup_and_duplicate_method_rejection();
  test_interval_id_validation();
  test_builder_state_machine_and_owned_lifetime();
  test_move_assignment_rejects_inactive_or_nonempty_builders();
  test_bernstein_reproduces_monomials();
  test_endpoint_and_interior_dispatch_and_theta_validation();
  test_exact_rational_theta_overload_and_validation();
  std::cout << "Subcycling dense-output contract tests passed\n";
  return 0;
}
