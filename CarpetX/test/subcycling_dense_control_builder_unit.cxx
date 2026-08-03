#include "subcycling_dense_control_builder.hxx"

#include <AMReX_MultiFab.H>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using CarpetX::DenseMFabControlSet;
using CarpetX::DenseMFabKey;
using CarpetX::DenseMFabRawSampleView;
using CarpetX::DenseMFabView;
using CarpetX::DenseSampleConstraint;
using CarpetX::DenseSampleKind;
using CarpetX::OwnedMultiFabDenseState;
using CarpetX::SubcyclingODEMethod;

static_assert(!std::is_copy_constructible_v<DenseMFabControlSet>);
static_assert(!std::is_copy_assignable_v<DenseMFabControlSet>);
static_assert(std::is_nothrow_move_constructible_v<DenseMFabControlSet>);
static_assert(std::is_nothrow_move_assignable_v<DenseMFabControlSet>);

constexpr SubcyclingODEMethod methods[] = {
    SubcyclingODEMethod::rk4, SubcyclingODEMethod::rkf78_order7,
    SubcyclingODEMethod::dp87_order8};

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

DenseMFabKey key(const std::int64_t epoch = 19, const int patch = 0,
                 const int level = 1, const int group = 7) {
  return DenseMFabKey{epoch, patch, level, group};
}

std::unique_ptr<amrex::MultiFab>
borrowed_mfab(const DenseMFabKey &, const double value,
              const int box_array = 11, const int distribution = 13,
              const int components = 1, const int grow = 2) {
  auto result = std::make_unique<amrex::MultiFab>(
      amrex::BoxArray(box_array), amrex::DistributionMapping(distribution),
      components, amrex::IntVect(grow));
  result->setVal(static_cast<amrex::Real>(value), 0, components, 0);
  for (int component = 0; component < components; ++component) {
    result->test_set_valid(
        component, 1,
        static_cast<amrex::Real>(value + 0.125 * (component + 1)));
    result->test_set_ghost(
        component, static_cast<amrex::Real>(900.0 + value + component));
  }
  return result;
}

std::unique_ptr<OwnedMultiFabDenseState>
state(const double value, const DenseMFabKey &entry_key = key(),
      const int box_array = 11, const int distribution = 13,
      const int components = 1, const int grow = 2) {
  auto borrowed = borrowed_mfab(entry_key, value, box_array, distribution,
                                components, grow);
  return OwnedMultiFabDenseState::copy_of(
      std::vector<DenseMFabView>{{entry_key, borrowed.get()}});
}

double first_value(const OwnedMultiFabDenseState &source) {
  return static_cast<double>(source.multifab(0).test_valid(0, 0));
}

double first_ghost(const OwnedMultiFabDenseState &source) {
  return static_cast<double>(source.multifab(0).test_ghost(0));
}

struct Samples {
  std::vector<std::unique_ptr<OwnedMultiFabDenseState>> owned;
  std::vector<DenseMFabRawSampleView> views;
};

Samples samples_from_values(const SubcyclingODEMethod method,
                            const std::vector<double> &values) {
  const auto &constraints =
      CarpetX::reference_dense_stencil(method).specification().constraints;
  assert(values.size() == constraints.size());

  Samples result;
  result.owned.reserve(values.size());
  result.views.reserve(values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    result.owned.push_back(state(values[index]));
    result.views.push_back(
        DenseMFabRawSampleView{constraints[index], result.owned.back().get()});
  }
  return result;
}

Samples monomial_samples(const SubcyclingODEMethod method, const double dt,
                         const int power) {
  const auto &constraints =
      CarpetX::reference_dense_stencil(method).specification().constraints;
  std::vector<double> values;
  values.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    if (constraint.kind == DenseSampleKind::value) {
      values.push_back(std::pow(constraint.theta, power));
    } else if (power == 0) {
      values.push_back(0.0);
    } else {
      values.push_back(static_cast<double>(power) *
                       std::pow(constraint.theta, power - 1) / dt);
    }
  }
  return samples_from_values(method, values);
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

double evaluate_bernstein(const DenseMFabControlSet &controls,
                          const double theta) {
  const int degree = static_cast<int>(controls.control_count()) - 1;
  double value = 0.0;
  for (int i = 0; i <= degree; ++i) {
    value += first_value(controls.control(static_cast<std::size_t>(i))) *
             binomial(degree, i) * std::pow(theta, i) *
             std::pow(1.0 - theta, degree - i);
  }
  return value;
}

void assert_near(const double actual, const double expected,
                 const double tolerance = 2.0e-5) {
  assert(std::abs(actual - expected) <= tolerance);
}

void test_builds_owned_controls_and_reproduces_reference_polynomials() {
  constexpr double dt = 0.375;
  for (const auto method : methods) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    const int degree = stencil.specification().dense_uniform_order;
    for (int power = 0; power <= degree; ++power) {
      auto samples = monomial_samples(method, dt, power);
      std::vector<double> scalar_samples;
      scalar_samples.reserve(samples.views.size());
      for (std::size_t index = 0; index < samples.views.size(); ++index) {
        const auto raw = first_value(*samples.views[index].state);
        scalar_samples.push_back(
            samples.views[index].constraint.kind == DenseSampleKind::value
                ? raw
                : dt * raw);
      }
      const auto scalar_controls = stencil.make_controls(scalar_samples);

      const auto controls =
          CarpetX::build_reference_dense_controls(method, dt, samples.views);
      assert(controls.method() == method);
      assert(controls.parent_dt() == dt);
      assert(controls.control_count() == stencil.control_count());
      for (std::size_t index = 0; index < controls.control_count(); ++index) {
        assert_near(first_value(controls.control(index)),
                    scalar_controls[index]);
      }
      for (const double theta : {0.0, 0.07, 0.31, 0.5, 0.83, 1.0}) {
        assert_near(evaluate_bernstein(controls, theta),
                    std::pow(theta, power), 7.0e-5);
      }
    }
  }
}

void test_only_derivative_columns_receive_parent_dt() {
  constexpr double dt = 3.25;
  for (const auto method : methods) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    std::vector<double> raw_samples(stencil.sample_count());
    std::vector<double> scaled_samples(stencil.sample_count());
    for (std::size_t index = 0; index < raw_samples.size(); ++index) {
      raw_samples[index] = 0.5 + static_cast<double>(index);
      scaled_samples[index] =
          stencil.specification().constraints[index].kind ==
                  DenseSampleKind::scaled_derivative
              ? dt * raw_samples[index]
              : raw_samples[index];
    }
    auto samples = samples_from_values(method, raw_samples);
    const auto expected = stencil.make_controls(scaled_samples);
    const auto controls =
        CarpetX::build_reference_dense_controls(method, dt, samples.views);
    for (std::size_t index = 0; index < controls.control_count(); ++index)
      assert_near(first_value(controls.control(index)), expected[index]);
  }
}

void test_endpoint_controls_are_exact_copy_operations() {
  auto samples = samples_from_values(SubcyclingODEMethod::rk4,
                                    {-0.0, 3.0, 4.0, -0.0, 5.0});
  amrex::MultiFab::reset_operation_log();

  const auto controls = CarpetX::build_reference_dense_controls(
      SubcyclingODEMethod::rk4, 0.5, samples.views);

  assert(std::signbit(first_value(controls.control(0))));
  assert(std::signbit(first_value(controls.control(4))));
  assert(amrex::MultiFab::copy_count() == 2);
  assert(amrex::MultiFab::all_operations_interior_only());
}

std::size_t expected_interior_saxpy_count(const SubcyclingODEMethod method,
                                          const double dt) {
  const auto &stencil = CarpetX::reference_dense_stencil(method);
  std::size_t result = 0;
  for (std::size_t control = 1; control + 1 < stencil.control_count();
       ++control) {
    for (std::size_t sample = 0; sample < stencil.sample_count(); ++sample) {
      double weight = stencil.weight(control, sample);
      if (stencil.specification().constraints[sample].kind ==
          DenseSampleKind::scaled_derivative)
        weight *= dt;
      const auto device_weight = static_cast<amrex::Real>(weight);
      if (device_weight != amrex::Real(0))
        ++result;
    }
  }
  return result;
}

void test_zero_coefficients_are_removed_before_saxpy() {
  constexpr double dt = 0.75;
  for (const auto method : methods) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    std::vector<double> values(stencil.sample_count(), 1.0);
    auto samples = samples_from_values(method, values);
    amrex::MultiFab::reset_operation_log();

    const auto controls =
        CarpetX::build_reference_dense_controls(method, dt, samples.views);

    assert(controls.control_count() == stencil.control_count());
    assert(amrex::MultiFab::saxpy_count() ==
           expected_interior_saxpy_count(method, dt));
    assert(amrex::MultiFab::setval_count() ==
           static_cast<int>(stencil.control_count() - 2));
    assert(amrex::MultiFab::copy_count() == 2);
  }
}

template <typename Function>
void expect_preallocation_failure(Function &&function) {
  amrex::MultiFab::reset_construction_count();
  expect_throw<std::invalid_argument>(std::forward<Function>(function));
  assert(amrex::MultiFab::construction_count() == 0);
}

void test_invalid_method_dt_count_and_metadata_fail_before_allocation() {
  auto samples = samples_from_values(SubcyclingODEMethod::rk4,
                                    {1.0, 2.0, 3.0, 4.0, 5.0});
  expect_preallocation_failure([&] {
    CarpetX::build_reference_dense_controls(
        static_cast<SubcyclingODEMethod>(99), 1.0, samples.views);
  });
  for (const double invalid_dt :
       {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()}) {
    expect_preallocation_failure([&] {
      CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4,
                                               invalid_dt, samples.views);
    });
  }

  auto wrong_count = samples.views;
  wrong_count.pop_back();
  expect_preallocation_failure([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                             wrong_count);
  });

  for (std::size_t index = 0; index < samples.views.size(); ++index) {
    auto wrong_theta = samples.views;
    wrong_theta[index].constraint.theta =
        std::nextafter(wrong_theta[index].constraint.theta, 2.0);
    expect_preallocation_failure([&] {
      CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                               wrong_theta);
    });

    auto wrong_kind = samples.views;
    wrong_kind[index].constraint.kind =
        wrong_kind[index].constraint.kind == DenseSampleKind::value
            ? DenseSampleKind::scaled_derivative
            : DenseSampleKind::value;
    expect_preallocation_failure([&] {
      CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                               wrong_kind);
    });
  }

  auto null_sample = samples.views;
  null_sample[2].state = nullptr;
  expect_preallocation_failure([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                             null_sample);
  });

  auto duplicate_sample = samples.views;
  duplicate_sample[2].state = duplicate_sample[0].state;
  expect_preallocation_failure([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                             duplicate_sample);
  });
}

void test_all_samples_are_globally_prevalidated_even_when_a_weight_is_zero() {
  auto samples = samples_from_values(SubcyclingODEMethod::rk4,
                                    {1.0, 2.0, 3.0, 4.0, 5.0});
  const auto &stencil =
      CarpetX::reference_dense_stencil(SubcyclingODEMethod::rk4);
  std::size_t zero_weight_sample = stencil.sample_count();
  for (std::size_t sample = 0; sample < stencil.sample_count(); ++sample) {
    if (stencil.weight(1, sample) == 0.0) {
      zero_weight_sample = sample;
      break;
    }
  }
  assert(zero_weight_sample < stencil.sample_count());

  auto incompatible_epoch = state(8.0, key(20, 0, 1, 7));
  auto epoch_views = samples.views;
  epoch_views[zero_weight_sample].state = incompatible_epoch.get();
  expect_preallocation_failure([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                             epoch_views);
  });

  const std::vector<std::unique_ptr<OwnedMultiFabDenseState>> mismatches = [&] {
    std::vector<std::unique_ptr<OwnedMultiFabDenseState>> result;
    result.push_back(state(8.0, key(19, 1, 1, 7)));
    result.push_back(state(8.0, key(19, 0, 2, 7)));
    result.push_back(state(8.0, key(19, 0, 1, 8)));
    result.push_back(state(8.0, key(), 17, 13, 1, 2));
    result.push_back(state(8.0, key(), 11, 19, 1, 2));
    result.push_back(state(8.0, key(), 11, 13, 2, 2));
    result.push_back(state(8.0, key(), 11, 13, 1, 3));
    return result;
  }();
  for (const auto &mismatch : mismatches) {
    auto views = samples.views;
    views[zero_weight_sample].state = mismatch.get();
    expect_preallocation_failure([&] {
      CarpetX::build_reference_dense_controls(SubcyclingODEMethod::rk4, 1.0,
                                               views);
    });
  }
}

void test_nonfinite_effective_coefficients_fail_before_allocation() {
  bool found_double_overflow_case = false;
  for (const auto method : methods) {
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    for (std::size_t control = 1; control + 1 < stencil.control_count();
         ++control) {
      for (std::size_t sample = 0; sample < stencil.sample_count(); ++sample) {
        if (stencil.specification().constraints[sample].kind !=
            DenseSampleKind::scaled_derivative)
          continue;
        const double weight = std::abs(stencil.weight(control, sample));
        if (!(weight > 1.0))
          continue;
        const double dt =
            std::nextafter(std::numeric_limits<double>::max() / weight,
                           std::numeric_limits<double>::max());
        if (!std::isfinite(dt) ||
            std::isfinite(dt * stencil.weight(control, sample)))
          continue;
        std::vector<double> values(stencil.sample_count(), 1.0);
        auto samples = samples_from_values(method, values);
        expect_preallocation_failure([&] {
          CarpetX::build_reference_dense_controls(method, dt, samples.views);
        });
        found_double_overflow_case = true;
        break;
      }
      if (found_double_overflow_case)
        break;
    }
    if (found_double_overflow_case)
      break;
  }
  assert(found_double_overflow_case);

  if constexpr (std::numeric_limits<amrex::Real>::max() <
                std::numeric_limits<double>::max()) {
    const auto method = SubcyclingODEMethod::rk4;
    const auto &stencil = CarpetX::reference_dense_stencil(method);
    std::vector<double> values(stencil.sample_count(), 1.0);
    auto samples = samples_from_values(method, values);
    const double dt =
        static_cast<double>(std::numeric_limits<amrex::Real>::max()) * 8.0;
    assert(std::isfinite(dt));
    bool has_finite_out_of_range_coefficient = false;
    for (std::size_t control = 0; control < stencil.control_count();
         ++control) {
      for (std::size_t sample = 0; sample < stencil.sample_count(); ++sample) {
        if (stencil.specification().constraints[sample].kind !=
            DenseSampleKind::scaled_derivative)
          continue;
        const double effective = dt * stencil.weight(control, sample);
        if (std::isfinite(effective) &&
            (effective < static_cast<double>(
                             std::numeric_limits<amrex::Real>::lowest()) ||
             effective > static_cast<double>(
                             std::numeric_limits<amrex::Real>::max()))) {
          has_finite_out_of_range_coefficient = true;
          break;
        }
      }
      if (has_finite_out_of_range_coefficient)
        break;
    }
    assert(has_finite_out_of_range_coefficient);
    expect_preallocation_failure([&] {
      CarpetX::build_reference_dense_controls(method, dt, samples.views);
    });
  }
}

void test_sources_and_ghosts_are_unchanged() {
  auto samples = samples_from_values(SubcyclingODEMethod::rk4,
                                    {1.0, 2.0, 3.0, 4.0, 5.0});
  std::vector<double> interiors;
  std::vector<double> ghosts;
  for (const auto &sample : samples.views) {
    interiors.push_back(first_value(*sample.state));
    ghosts.push_back(first_ghost(*sample.state));
  }

  const auto controls = CarpetX::build_reference_dense_controls(
      SubcyclingODEMethod::rk4, 0.5, samples.views);

  for (std::size_t index = 0; index < samples.views.size(); ++index) {
    assert(first_value(*samples.views[index].state) == interiors[index]);
    assert(first_ghost(*samples.views[index].state) == ghosts[index]);
  }
  for (std::size_t index = 0; index < controls.control_count(); ++index)
    assert(first_ghost(controls.control(index)) ==
           static_cast<double>(amrex::MultiFab::uninitialized_value()));
}

void test_allocation_and_operation_failures_release_all_local_controls() {
  auto samples = samples_from_values(SubcyclingODEMethod::dp87_order8,
                                    {1.0, 2.0, 3.0, 4.0, 5.0,
                                     6.0, 7.0, 8.0, 9.0});
  const int baseline = amrex::MultiFab::live_count();

  amrex::MultiFab::reset_failures();
  amrex::MultiFab::fail_construction_after(2);
  expect_throw<std::runtime_error>([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::dp87_order8,
                                             0.5, samples.views);
  });
  assert(amrex::MultiFab::live_count() == baseline);

  amrex::MultiFab::reset_failures();
  amrex::MultiFab::fail_setval_after(0);
  expect_throw<std::runtime_error>([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::dp87_order8,
                                             0.5, samples.views);
  });
  assert(amrex::MultiFab::live_count() == baseline);

  amrex::MultiFab::reset_failures();
  amrex::MultiFab::fail_saxpy_after(5);
  expect_throw<std::runtime_error>([&] {
    CarpetX::build_reference_dense_controls(SubcyclingODEMethod::dp87_order8,
                                             0.5, samples.views);
  });
  assert(amrex::MultiFab::live_count() == baseline);
  amrex::MultiFab::reset_failures();
}

void test_move_release_and_bounds_contract() {
  auto samples = samples_from_values(SubcyclingODEMethod::rk4,
                                    {1.0, 2.0, 3.0, 4.0, 5.0});
  auto controls = CarpetX::build_reference_dense_controls(
      SubcyclingODEMethod::rk4, 0.5, samples.views);
  expect_throw<std::out_of_range>([&] { controls.control(5); });

  DenseMFabControlSet moved(std::move(controls));
  assert(moved.control_count() == 5);
  auto released = std::move(moved).release_controls();
  assert(released.size() == 5);
  for (const auto &control : released)
    assert(control != nullptr);
}

} // namespace

int main() {
  test_builds_owned_controls_and_reproduces_reference_polynomials();
  test_only_derivative_columns_receive_parent_dt();
  test_endpoint_controls_are_exact_copy_operations();
  test_zero_coefficients_are_removed_before_saxpy();
  test_invalid_method_dt_count_and_metadata_fail_before_allocation();
  test_all_samples_are_globally_prevalidated_even_when_a_weight_is_zero();
  test_nonfinite_effective_coefficients_fail_before_allocation();
  test_sources_and_ghosts_are_unchanged();
  test_allocation_and_operation_failures_release_all_local_controls();
  test_move_release_and_bounds_contract();
  std::cout << "Dense MultiFab control-builder tests passed\n";
  return 0;
}
