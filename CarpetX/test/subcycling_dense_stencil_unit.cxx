#include "subcycling_dense_stencil.hxx"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using CarpetX::DenseSampleConstraint;
using CarpetX::DenseSampleKind;
using CarpetX::DenseStencil;
using CarpetX::DenseStencilSpecification;
using CarpetX::SubcyclingODEMethod;

static_assert(!std::is_copy_assignable_v<DenseStencil>);
static_assert(!std::is_move_assignable_v<DenseStencil>);

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

double binomial(const int n, const int k) {
  if (k < 0 || k > n)
    return 0.0;
  const int smaller = k < n - k ? k : n - k;
  double result = 1.0;
  for (int j = 1; j <= smaller; ++j)
    result *= static_cast<double>(n - smaller + j) / static_cast<double>(j);
  return result;
}

double evaluate_bernstein(const std::vector<double> &controls,
                          const double theta) {
  const int degree = static_cast<int>(controls.size()) - 1;
  double result = 0.0;
  for (int i = 0; i <= degree; ++i)
    result += controls[static_cast<std::size_t>(i)] * binomial(degree, i) *
              std::pow(theta, i) * std::pow(1.0 - theta, degree - i);
  return result;
}

double evaluate_bernstein_derivative(const std::vector<double> &controls,
                                     const double theta) {
  const int degree = static_cast<int>(controls.size()) - 1;
  if (degree == 0)
    return 0.0;
  double result = 0.0;
  for (int i = 0; i < degree; ++i)
    result += static_cast<double>(degree) *
              (controls[static_cast<std::size_t>(i + 1)] -
               controls[static_cast<std::size_t>(i)]) *
              binomial(degree - 1, i) * std::pow(theta, i) *
              std::pow(1.0 - theta, degree - 1 - i);
  return result;
}

DenseStencilSpecification valid_rk4_specification() {
  return DenseStencilSpecification{
      SubcyclingODEMethod::rk4,
      4,
      4,
      4,
      4,
      {{0.0, DenseSampleKind::value},
       {0.0, DenseSampleKind::scaled_derivative},
       {0.5, DenseSampleKind::value},
       {1.0, DenseSampleKind::value},
       {1.0, DenseSampleKind::scaled_derivative}}};
}

void assert_constraint(const DenseSampleConstraint &constraint,
                       const double theta, const DenseSampleKind kind) {
  assert(constraint.theta == theta);
  assert(constraint.kind == kind);
}

void test_exact_reference_method_specifications() {
  const auto &rk4 =
      CarpetX::reference_dense_stencil(SubcyclingODEMethod::rk4);
  const auto &rk4_spec = rk4.specification();
  assert(rk4_spec.endpoint_order == 4);
  assert(rk4_spec.dense_uniform_order == 4);
  assert(rk4_spec.stage_count == 4);
  assert(rk4_spec.extra_rhs_evaluations == 4);
  assert(rk4_spec.constraints.size() == 5);
  assert_constraint(rk4_spec.constraints[0], 0.0, DenseSampleKind::value);
  assert_constraint(rk4_spec.constraints[1], 0.0,
                    DenseSampleKind::scaled_derivative);
  assert_constraint(rk4_spec.constraints[2], 0.5, DenseSampleKind::value);
  assert_constraint(rk4_spec.constraints[3], 1.0, DenseSampleKind::value);
  assert_constraint(rk4_spec.constraints[4], 1.0,
                    DenseSampleKind::scaled_derivative);

  const auto &rkf78 =
      CarpetX::reference_dense_stencil(SubcyclingODEMethod::rkf78_order7);
  const auto &rkf78_spec = rkf78.specification();
  assert(rkf78_spec.endpoint_order == 7);
  assert(rkf78_spec.dense_uniform_order == 7);
  assert(rkf78_spec.stage_count == 11);
  assert(rkf78_spec.extra_rhs_evaluations == 23);
  assert(rkf78_spec.constraints.size() == 8);
  for (std::size_t i = 0; i < 4; ++i) {
    const double theta = static_cast<double>(i) / 3.0;
    assert_constraint(rkf78_spec.constraints[2 * i], theta,
                      DenseSampleKind::value);
    assert_constraint(rkf78_spec.constraints[2 * i + 1], theta,
                      DenseSampleKind::scaled_derivative);
  }

  const auto &dp87 =
      CarpetX::reference_dense_stencil(SubcyclingODEMethod::dp87_order8);
  const auto &dp87_spec = dp87.specification();
  assert(dp87_spec.endpoint_order == 8);
  assert(dp87_spec.dense_uniform_order == 8);
  assert(dp87_spec.stage_count == 13);
  assert(dp87_spec.extra_rhs_evaluations == 39);
  assert(dp87_spec.constraints.size() == 9);
  assert_constraint(dp87_spec.constraints[0], 0.0, DenseSampleKind::value);
  assert_constraint(dp87_spec.constraints[1], 0.0,
                    DenseSampleKind::scaled_derivative);
  assert_constraint(dp87_spec.constraints[2], 0.25, DenseSampleKind::value);
  assert_constraint(dp87_spec.constraints[3], 0.25,
                    DenseSampleKind::scaled_derivative);
  assert_constraint(dp87_spec.constraints[4], 0.5, DenseSampleKind::value);
  assert_constraint(dp87_spec.constraints[5], 0.5,
                    DenseSampleKind::scaled_derivative);
  assert_constraint(dp87_spec.constraints[6], 0.75, DenseSampleKind::value);
  assert_constraint(dp87_spec.constraints[7], 1.0, DenseSampleKind::value);
  assert_constraint(dp87_spec.constraints[8], 1.0,
                    DenseSampleKind::scaled_derivative);
  for (const auto &constraint : dp87_spec.constraints)
    assert(!(constraint.theta == 0.75 &&
             constraint.kind == DenseSampleKind::scaled_derivative));

  assert(&rk4 == &CarpetX::reference_dense_stencil(SubcyclingODEMethod::rk4));
  assert(&rkf78 ==
         &CarpetX::reference_dense_stencil(SubcyclingODEMethod::rkf78_order7));
  assert(&dp87 ==
         &CarpetX::reference_dense_stencil(SubcyclingODEMethod::dp87_order8));
}

void test_invalid_specifications_fail_closed() {
  auto spec = valid_rk4_specification();
  spec.constraints[2] = {0.0, DenseSampleKind::value};
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  spec.constraints.insert(spec.constraints.begin() + 2,
                          {0.0, DenseSampleKind::scaled_derivative});
  spec.constraints.pop_back();
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  std::swap(spec.constraints[2], spec.constraints[3]);
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  std::swap(spec.constraints[0], spec.constraints[1]);
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  spec.constraints[2].theta = std::numeric_limits<double>::quiet_NaN();
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  for (const double invalid_theta : {-0.01, 1.01}) {
    spec = valid_rk4_specification();
    spec.constraints[2].theta = invalid_theta;
    expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  }

  spec = valid_rk4_specification();
  spec.endpoint_order = 0;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.dense_uniform_order = 0;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.dense_uniform_order = 3;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.dense_uniform_order = std::numeric_limits<int>::max();
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.constraints.pop_back();
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  spec.constraints[0].theta = 0.25;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.constraints[3].theta = 0.75;
  spec.constraints[4].theta = 0.75;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });

  spec = valid_rk4_specification();
  spec.stage_count = 0;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.extra_rhs_evaluations = -1;
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.method = static_cast<SubcyclingODEMethod>(99);
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
  spec = valid_rk4_specification();
  spec.constraints[2].kind = static_cast<DenseSampleKind>(99);
  expect_throw<std::invalid_argument>([&] { DenseStencil stencil(spec); });
}

std::vector<double> monomial_samples(const DenseStencil &stencil,
                                     const int power) {
  std::vector<double> samples;
  for (const auto &constraint : stencil.specification().constraints) {
    if (constraint.kind == DenseSampleKind::value) {
      samples.push_back(std::pow(constraint.theta, power));
    } else if (power == 0) {
      samples.push_back(0.0);
    } else {
      samples.push_back(static_cast<double>(power) *
                        std::pow(constraint.theta, power - 1));
    }
  }
  return samples;
}

void test_monomial_values_and_derivatives_are_reproduced() {
  for (const auto method : {SubcyclingODEMethod::rk4,
                            SubcyclingODEMethod::rkf78_order7,
                            SubcyclingODEMethod::dp87_order8}) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    const int degree = stencil.specification().dense_uniform_order;
    for (int power = 0; power <= degree; ++power) {
      const auto controls = stencil.make_controls(monomial_samples(stencil, power));
      for (const double theta : {0.0, 0.07, 0.23, 0.5, 0.81, 1.0})
        assert(std::abs(evaluate_bernstein(controls, theta) -
                        std::pow(theta, power)) < 2.0e-11);
      for (const auto &constraint : stencil.specification().constraints) {
        const double expected =
            power == 0
                ? 0.0
                : static_cast<double>(power) *
                      std::pow(constraint.theta, power - 1);
        assert(std::abs(evaluate_bernstein_derivative(
                            controls, constraint.theta) -
                        expected) < 2.0e-10);
      }
    }
  }
}

void test_weight_matrix_is_finite_immutable_and_endpoint_exact() {
  for (const auto method : {SubcyclingODEMethod::rk4,
                            SubcyclingODEMethod::rkf78_order7,
                            SubcyclingODEMethod::dp87_order8}) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    const std::size_t count = stencil.sample_count();
    assert(stencil.control_count() == count);
    assert(stencil.weights().size() == count * count);
    for (const double weight : stencil.weights())
      assert(std::isfinite(weight));

    const auto &constraints = stencil.specification().constraints;
    std::size_t begin_value = count;
    std::size_t end_value = count;
    for (std::size_t j = 0; j < count; ++j) {
      if (constraints[j].theta == 0.0 &&
          constraints[j].kind == DenseSampleKind::value)
        begin_value = j;
      if (constraints[j].theta == 1.0 &&
          constraints[j].kind == DenseSampleKind::value)
        end_value = j;
    }
    assert(begin_value < count);
    assert(end_value < count);
    for (std::size_t j = 0; j < count; ++j) {
      assert(stencil.weight(0, j) == (j == begin_value ? 1.0 : 0.0));
      assert(stencil.weight(count - 1, j) ==
             (j == end_value ? 1.0 : 0.0));
    }

    const auto first_weights = stencil.weights();
    const auto second_controls =
        CarpetX::reference_dense_stencil(method).make_controls(
            std::vector<double>(count, 1.25));
    assert(stencil.weights() == first_weights);
    assert(second_controls.size() == count);
    expect_throw<std::out_of_range>([&] { stencil.weight(count, 0); });
    expect_throw<std::out_of_range>([&] { stencil.weight(0, count); });
    expect_throw<std::invalid_argument>(
        [&] { stencil.make_controls(std::vector<double>(count - 1, 0.0)); });
  }
}

} // namespace

int main() {
  test_exact_reference_method_specifications();
  test_invalid_specifications_fail_closed();
  test_monomial_values_and_derivatives_are_reproduced();
  test_weight_matrix_is_finite_immutable_and_endpoint_exact();
  std::cout << "Subcycling dense-stencil tests passed\n";
  return 0;
}
