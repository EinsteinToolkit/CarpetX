#include "subcycling_dense_stencil.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace CarpetX {
namespace {

using long_matrix = std::vector<std::vector<long double>>;

bool known_method(const SubcyclingODEMethod method) noexcept {
  switch (method) {
  case SubcyclingODEMethod::rk4:
  case SubcyclingODEMethod::rkf78_order7:
  case SubcyclingODEMethod::dp87_order8:
    return true;
  }
  return false;
}

void validate_specification(const DenseStencilSpecification &specification) {
  if (!known_method(specification.method))
    throw std::invalid_argument("dense stencil has an unsupported ODE method");
  if (specification.endpoint_order <= 0 ||
      specification.dense_uniform_order <= 0 ||
      specification.dense_uniform_order < specification.endpoint_order)
    throw std::invalid_argument("dense stencil orders are invalid");
  if (specification.stage_count <= 0 ||
      specification.extra_rhs_evaluations < 0)
    throw std::invalid_argument("dense stencil stage cost is invalid");

  const auto expected_count =
      static_cast<std::size_t>(specification.dense_uniform_order) + 1U;
  if (specification.constraints.size() != expected_count)
    throw std::invalid_argument(
        "dense stencil needs exactly dense_uniform_order + 1 samples");

  bool has_begin_value = false;
  bool has_end_value = false;
  double previous_theta = 0.0;
  DenseSampleKind previous_kind = DenseSampleKind::value;
  bool have_previous = false;
  for (const auto &constraint : specification.constraints) {
    if (!std::isfinite(constraint.theta) || constraint.theta < 0.0 ||
        constraint.theta > 1.0)
      throw std::invalid_argument(
          "dense stencil sample theta must be finite and in [0, 1]");
    if (constraint.kind != DenseSampleKind::value &&
        constraint.kind != DenseSampleKind::scaled_derivative)
      throw std::invalid_argument("dense stencil has an unknown sample kind");

    if (!have_previous) {
      if (constraint.kind != DenseSampleKind::value)
        throw std::invalid_argument(
            "dense stencil derivative cannot precede its value");
    } else if (constraint.theta < previous_theta) {
      throw std::invalid_argument("dense stencil sample nodes must be sorted");
    } else if (constraint.theta == previous_theta) {
      if (previous_kind != DenseSampleKind::value ||
          constraint.kind != DenseSampleKind::scaled_derivative)
        throw std::invalid_argument(
            "dense stencil sample constraints are duplicated or misordered");
    } else if (constraint.kind != DenseSampleKind::value) {
      throw std::invalid_argument(
          "dense stencil derivative cannot precede its value");
    }

    if (constraint.theta == 0.0 &&
        constraint.kind == DenseSampleKind::value)
      has_begin_value = true;
    if (constraint.theta == 1.0 &&
        constraint.kind == DenseSampleKind::value)
      has_end_value = true;
    previous_theta = constraint.theta;
    previous_kind = constraint.kind;
    have_previous = true;
  }
  if (!has_begin_value || !has_end_value)
    throw std::invalid_argument(
        "dense stencil requires endpoint value constraints at 0 and 1");
}

long double binomial(const int n, const int k) {
  if (k < 0 || k > n)
    return 0.0L;
  const int smaller = std::min(k, n - k);
  long double result = 1.0L;
  for (int j = 1; j <= smaller; ++j)
    result *= static_cast<long double>(n - smaller + j) /
              static_cast<long double>(j);
  return result;
}

long_matrix inverse_confluent_matrix(
    const DenseStencilSpecification &specification) {
  const int degree = specification.dense_uniform_order;
  const std::size_t count = specification.constraints.size();
  long_matrix augmented(count, std::vector<long double>(2 * count, 0.0L));
  long double matrix_scale = 0.0L;

  for (std::size_t row = 0; row < count; ++row) {
    const auto &constraint = specification.constraints[row];
    const long double theta = static_cast<long double>(constraint.theta);
    if (constraint.kind == DenseSampleKind::value) {
      long double theta_power = 1.0L;
      for (int column = 0; column <= degree; ++column) {
        augmented[row][static_cast<std::size_t>(column)] = theta_power;
        matrix_scale = std::max(matrix_scale, std::abs(theta_power));
        theta_power *= theta;
      }
    } else {
      // k=0 is exactly zero.  Do not form theta^(-1) at theta=0.
      augmented[row][0] = 0.0L;
      long double theta_power = 1.0L;
      for (int column = 1; column <= degree; ++column) {
        const long double entry =
            static_cast<long double>(column) * theta_power;
        augmented[row][static_cast<std::size_t>(column)] = entry;
        matrix_scale = std::max(matrix_scale, std::abs(entry));
        theta_power *= theta;
      }
    }
    augmented[row][count + row] = 1.0L;
  }

  const long double singularity_threshold =
      std::numeric_limits<long double>::epsilon() *
      std::max(1.0L, matrix_scale) * static_cast<long double>(count) * 128.0L;

  for (std::size_t column = 0; column < count; ++column) {
    std::size_t pivot_row = column;
    long double pivot_magnitude = std::abs(augmented[column][column]);
    for (std::size_t row = column + 1; row < count; ++row) {
      const long double candidate = std::abs(augmented[row][column]);
      if (candidate > pivot_magnitude) {
        pivot_magnitude = candidate;
        pivot_row = row;
      }
    }
    if (pivot_magnitude <= singularity_threshold)
      throw std::invalid_argument(
          "dense stencil constraints form a singular interpolation matrix");
    if (pivot_row != column)
      std::swap(augmented[pivot_row], augmented[column]);

    const long double pivot = augmented[column][column];
    for (long double &entry : augmented[column])
      entry /= pivot;

    for (std::size_t row = 0; row < count; ++row) {
      if (row == column)
        continue;
      const long double factor = augmented[row][column];
      if (factor == 0.0L)
        continue;
      for (std::size_t entry = 0; entry < 2 * count; ++entry)
        augmented[row][entry] -= factor * augmented[column][entry];
    }
  }

  long_matrix inverse(count, std::vector<long double>(count, 0.0L));
  for (std::size_t row = 0; row < count; ++row)
    for (std::size_t column = 0; column < count; ++column)
      inverse[row][column] = augmented[row][count + column];
  return inverse;
}

DenseStencilSpecification rk4_specification() {
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

DenseStencilSpecification rkf78_specification() {
  return DenseStencilSpecification{
      SubcyclingODEMethod::rkf78_order7,
      7,
      7,
      11,
      23,
      {{0.0, DenseSampleKind::value},
       {0.0, DenseSampleKind::scaled_derivative},
       {1.0 / 3.0, DenseSampleKind::value},
       {1.0 / 3.0, DenseSampleKind::scaled_derivative},
       {2.0 / 3.0, DenseSampleKind::value},
       {2.0 / 3.0, DenseSampleKind::scaled_derivative},
       {1.0, DenseSampleKind::value},
       {1.0, DenseSampleKind::scaled_derivative}}};
}

DenseStencilSpecification dp87_specification() {
  return DenseStencilSpecification{
      SubcyclingODEMethod::dp87_order8,
      8,
      8,
      13,
      39,
      {{0.0, DenseSampleKind::value},
       {0.0, DenseSampleKind::scaled_derivative},
       {0.25, DenseSampleKind::value},
       {0.25, DenseSampleKind::scaled_derivative},
       {0.5, DenseSampleKind::value},
       {0.5, DenseSampleKind::scaled_derivative},
       {0.75, DenseSampleKind::value},
       {1.0, DenseSampleKind::value},
       {1.0, DenseSampleKind::scaled_derivative}}};
}

} // namespace

DenseStencil::DenseStencil(DenseStencilSpecification specification)
    : specification_(std::move(specification)),
      weights_(build_weights(specification_)) {}

const DenseStencilSpecification &DenseStencil::specification() const noexcept {
  return specification_;
}

std::size_t DenseStencil::control_count() const noexcept {
  return specification_.constraints.size();
}

std::size_t DenseStencil::sample_count() const noexcept {
  return specification_.constraints.size();
}

const std::vector<double> &DenseStencil::weights() const noexcept {
  return weights_;
}

double DenseStencil::weight(const std::size_t control,
                            const std::size_t sample) const {
  const std::size_t count = sample_count();
  if (control >= count || sample >= count)
    throw std::out_of_range("dense stencil weight index is out of range");
  return weights_[control * count + sample];
}

std::vector<double>
DenseStencil::make_controls(const std::vector<double> &samples) const {
  const std::size_t count = sample_count();
  if (samples.size() != count)
    throw std::invalid_argument("dense stencil sample count does not match");
  std::vector<double> controls(count, 0.0);
  for (std::size_t control = 0; control < count; ++control)
    for (std::size_t sample = 0; sample < count; ++sample)
      controls[control] += weight(control, sample) * samples[sample];
  return controls;
}

std::vector<double> DenseStencil::build_weights(
    const DenseStencilSpecification &specification) {
  validate_specification(specification);
  const int degree = specification.dense_uniform_order;
  const std::size_t count = specification.constraints.size();
  const auto inverse = inverse_confluent_matrix(specification);
  long_matrix long_weights(count, std::vector<long double>(count, 0.0L));

  for (int control = 0; control <= degree; ++control) {
    for (std::size_t sample = 0; sample < count; ++sample) {
      long double result = 0.0L;
      for (int power = 0; power <= control; ++power) {
        const long double coefficient =
            binomial(control, power) / binomial(degree, power);
        result += coefficient *
                  inverse[static_cast<std::size_t>(power)][sample];
      }
      long_weights[static_cast<std::size_t>(control)][sample] = result;
    }
  }

  std::size_t begin_value = count;
  std::size_t end_value = count;
  for (std::size_t sample = 0; sample < count; ++sample) {
    const auto &constraint = specification.constraints[sample];
    if (constraint.theta == 0.0 &&
        constraint.kind == DenseSampleKind::value)
      begin_value = sample;
    if (constraint.theta == 1.0 &&
        constraint.kind == DenseSampleKind::value)
      end_value = sample;
  }

  const long double endpoint_tolerance = 2.0e-12L;
  for (std::size_t sample = 0; sample < count; ++sample) {
    const long double expected_begin = sample == begin_value ? 1.0L : 0.0L;
    const long double expected_end = sample == end_value ? 1.0L : 0.0L;
    if (std::abs(long_weights.front()[sample] - expected_begin) >
            endpoint_tolerance ||
        std::abs(long_weights.back()[sample] - expected_end) >
            endpoint_tolerance)
      throw std::runtime_error(
          "dense stencil endpoint rows failed identity verification");
  }

  std::vector<double> result(count * count, 0.0);
  for (std::size_t control = 0; control < count; ++control) {
    for (std::size_t sample = 0; sample < count; ++sample) {
      const long double entry = long_weights[control][sample];
      if (!std::isfinite(entry))
        throw std::runtime_error("dense stencil produced a non-finite weight");
      const double converted = static_cast<double>(entry);
      if (!std::isfinite(converted))
        throw std::runtime_error("dense stencil weight does not fit in double");
      result[control * count + sample] = converted;
    }
  }
  for (std::size_t sample = 0; sample < count; ++sample) {
    result[sample] = sample == begin_value ? 1.0 : 0.0;
    result[(count - 1) * count + sample] = sample == end_value ? 1.0 : 0.0;
  }
  return result;
}

const DenseStencil &reference_dense_stencil(const SubcyclingODEMethod method) {
  static const DenseStencil rk4(rk4_specification());
  static const DenseStencil rkf78(rkf78_specification());
  static const DenseStencil dp87(dp87_specification());
  switch (method) {
  case SubcyclingODEMethod::rk4:
    return rk4;
  case SubcyclingODEMethod::rkf78_order7:
    return rkf78;
  case SubcyclingODEMethod::dp87_order8:
    return dp87;
  }
  throw std::invalid_argument("no reference dense stencil for ODE method");
}

} // namespace CarpetX
