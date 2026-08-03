#ifndef ODESOLVERS_EXPLICIT_RK_DENSE_PROVIDER_HXX
#define ODESOLVERS_EXPLICIT_RK_DENSE_PROVIDER_HXX

#include "explicit_rk.hxx"

#include <subcycling_dense_output.hxx>
#include <subcycling_dense_stencil.hxx>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace ODESolvers {

template <class State> struct ReferenceDenseSample {
  double theta;
  CarpetX::DenseSampleKind kind;
  State payload;

  ReferenceDenseSample(const double theta_,
                       const CarpetX::DenseSampleKind kind_, State payload_)
      : theta(theta_), kind(kind_), payload(std::move(payload_)) {}

  ReferenceDenseSample(const ReferenceDenseSample &) = delete;
  ReferenceDenseSample &operator=(const ReferenceDenseSample &) = delete;
  ReferenceDenseSample(ReferenceDenseSample &&) noexcept(
      std::is_nothrow_move_constructible<State>::value) = default;
  ReferenceDenseSample &operator=(ReferenceDenseSample &&) noexcept(
      std::is_nothrow_move_assignable<State>::value) = default;
};

template <class State> class ReferenceDenseSampleBatch {
public:
  using sample_type = ReferenceDenseSample<State>;

  ReferenceDenseSampleBatch(const ReferenceDenseSampleBatch &) = delete;
  ReferenceDenseSampleBatch &operator=(const ReferenceDenseSampleBatch &) =
      delete;
  ReferenceDenseSampleBatch(ReferenceDenseSampleBatch &&) noexcept = default;
  ReferenceDenseSampleBatch &
  operator=(ReferenceDenseSampleBatch &&) noexcept = default;

  const std::vector<sample_type> &samples() const noexcept { return samples_; }
  std::size_t size() const noexcept { return samples_.size(); }
  bool empty() const noexcept { return samples_.empty(); }

private:
  explicit ReferenceDenseSampleBatch(std::vector<sample_type> samples)
      : samples_(std::move(samples)) {}

  std::vector<sample_type> samples_;

  template <class ScratchOps>
  friend ReferenceDenseSampleBatch<typename ScratchOps::state_type>
  collect_reference_dense_samples(
      ExplicitRKMethod, typename ScratchOps::scalar_type,
      typename ScratchOps::scalar_type,
      const typename ScratchOps::state_type &,
      const typename ScratchOps::state_type &,
      const typename ScratchOps::state_type &, ScratchOps &);
};

inline CarpetX::SubcyclingODEMethod
subcycling_method(const ExplicitRKMethod method) {
  switch (method) {
  case ExplicitRKMethod::rk4:
    return CarpetX::SubcyclingODEMethod::rk4;
  case ExplicitRKMethod::rkf78:
    return CarpetX::SubcyclingODEMethod::rkf78_order7;
  case ExplicitRKMethod::dp87:
    return CarpetX::SubcyclingODEMethod::dp87_order8;
  }
  throw std::invalid_argument("unsupported explicit RK dense method");
}

inline CarpetX::DenseCapability
reference_dense_capability(const ExplicitRKMethod method) {
  const auto dense_method = subcycling_method(method);
  const auto &specification =
      CarpetX::reference_dense_stencil(dense_method).specification();
  const auto fingerprint = explicit_rk_tableau_fingerprint(method);
  return CarpetX::DenseCapability{
      dense_method,
      fingerprint,
      specification.endpoint_order,
      specification.dense_uniform_order,
      specification.stage_count,
      specification.extra_rhs_evaluations,
      specification.dense_uniform_order + 1,
      true,
      true};
}

inline std::shared_ptr<const CarpetX::DenseOutputProvider>
make_reference_dense_provider(const ExplicitRKMethod method) {
  return std::make_shared<const CarpetX::DenseOutputProvider>(
      reference_dense_capability(method));
}

inline RationalCoefficient
reference_dense_sample_fraction(const ExplicitRKMethod method,
                                const double theta) {
  const auto match = [theta](const double expected,
                             const RationalCoefficient exact) {
    if (theta == expected)
      return exact;
    return RationalCoefficient{std::numeric_limits<std::int64_t>::min(), 1};
  };
  for (const auto candidate : {
           match(0.0, {0, 1}), match(1.0, {1, 1}),
           method == ExplicitRKMethod::rk4
               ? match(0.5, {1, 2})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1},
           method == ExplicitRKMethod::rkf78
               ? match(1.0 / 3.0, {1, 3})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1},
           method == ExplicitRKMethod::rkf78
               ? match(2.0 / 3.0, {2, 3})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1},
           method == ExplicitRKMethod::dp87
               ? match(0.25, {1, 4})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1},
           method == ExplicitRKMethod::dp87
               ? match(0.5, {1, 2})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1},
           method == ExplicitRKMethod::dp87
               ? match(0.75, {3, 4})
               : RationalCoefficient{
                     std::numeric_limits<std::int64_t>::min(), 1}}) {
    if (candidate.numerator != std::numeric_limits<std::int64_t>::min())
      return candidate;
  }
  throw std::logic_error(
      "reference dense sample has no audited exact fraction");
}

template <class ScratchOps>
ReferenceDenseSampleBatch<typename ScratchOps::state_type>
collect_reference_dense_samples(
    const ExplicitRKMethod method, const typename ScratchOps::scalar_type t0,
    const typename ScratchOps::scalar_type parent_dt,
    const typename ScratchOps::state_type &left_state,
    const typename ScratchOps::state_type &left_rhs,
    const typename ScratchOps::state_type &accepted_endpoint,
    ScratchOps &scratch) {
  using Scalar = typename ScratchOps::scalar_type;
  using State = typename ScratchOps::state_type;
  using Sample = ReferenceDenseSample<State>;

  using std::isfinite;
  if (!isfinite(t0) || !isfinite(parent_dt) || !(parent_dt > Scalar(0)) ||
      !isfinite(t0 + parent_dt))
    throw std::invalid_argument(
        "reference dense interval must be finite with positive dt");

  const auto dense_method = subcycling_method(method);
  const auto &stencil = CarpetX::reference_dense_stencil(dense_method);
  const auto &specification = stencil.specification();
  const std::size_t initial_rhs_count = scratch.rhs_evaluation_count();
  std::vector<Sample> candidate;
  candidate.reserve(specification.constraints.size());

  const auto append = [&](const CarpetX::DenseSampleConstraint constraint,
                          State payload) {
    candidate.emplace_back(constraint.theta, constraint.kind,
                           std::move(payload));
  };

  std::size_t index = 0;
  while (index < specification.constraints.size()) {
    const auto value_constraint = specification.constraints[index];
    if (value_constraint.kind != CarpetX::DenseSampleKind::value)
      throw std::logic_error(
          "reference dense stencil value/derivative order is invalid");

    const Scalar theta = static_cast<Scalar>(value_constraint.theta);
    const Scalar sample_time = t0 + theta * parent_dt;
    const auto sample_fraction =
        reference_dense_sample_fraction(method, value_constraint.theta);
    if (value_constraint.theta == 0.0) {
      scratch.restore_left(left_state, left_rhs, t0);
      append(value_constraint, scratch.snapshot_state());
    } else if (value_constraint.theta == 1.0) {
      scratch.restore_state(accepted_endpoint, t0 + parent_dt);
      append(value_constraint, scratch.snapshot_state());
    } else {
      scratch.restore_left(left_state, left_rhs, t0);
      auto token = make_loaded_rhs_token(scratch, t0);
      advance_explicit_rk(method, t0, theta * parent_dt,
                          InitialRHSMode::reuse_loaded,
                          {ExplicitRKStageKind::fractional, {0, 1},
                           sample_fraction},
                          scratch, token);
      append(value_constraint, scratch.snapshot_state());
    }
    ++index;

    if (index < specification.constraints.size()) {
      const auto derivative_constraint = specification.constraints[index];
      if (derivative_constraint.theta == value_constraint.theta) {
        if (derivative_constraint.kind !=
            CarpetX::DenseSampleKind::scaled_derivative)
          throw std::logic_error(
              "reference dense stencil has duplicate value samples");
        if (value_constraint.theta == 0.0) {
          scratch.restore_left(left_state, left_rhs, t0);
          append(derivative_constraint, scratch.snapshot_rhs());
        } else {
          scratch.restore_state(candidate.back().payload, sample_time);
          append(derivative_constraint,
                 scratch.probe_endpoint_rhs(
                     sample_time,
                     {ExplicitRKStageKind::endpoint_probe, 1, 1,
                      sample_fraction}));
        }
        ++index;
      }
    }
  }

  const std::size_t final_rhs_count = scratch.rhs_evaluation_count();
  if (final_rhs_count < initial_rhs_count ||
      final_rhs_count - initial_rhs_count !=
          static_cast<std::size_t>(specification.extra_rhs_evaluations))
    throw std::runtime_error(
        "reference dense provider observed an unexpected extra RHS count");
  if (candidate.size() != stencil.sample_count())
    throw std::runtime_error(
        "reference dense provider produced an incomplete sample batch");
  for (std::size_t sample = 0; sample < candidate.size(); ++sample) {
    const auto expected = specification.constraints[sample];
    if (candidate[sample].theta != expected.theta ||
        candidate[sample].kind != expected.kind)
      throw std::runtime_error(
          "reference dense provider sample order changed unexpectedly");
  }

  return ReferenceDenseSampleBatch<State>(std::move(candidate));
}

} // namespace ODESolvers

#endif // ODESOLVERS_EXPLICIT_RK_DENSE_PROVIDER_HXX
