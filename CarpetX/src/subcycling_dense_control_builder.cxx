#include "subcycling_dense_control_builder.hxx"

#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace CarpetX {
namespace {

struct ValidatedControlMetadata {
  const DenseStencil *stencil;
  std::size_t begin_value_sample;
  std::size_t end_value_sample;
  std::vector<std::vector<double>> effective_weights;
};

bool same_constraint(const DenseSampleConstraint &left,
                     const DenseSampleConstraint &right) noexcept {
  return left.theta == right.theta && left.kind == right.kind;
}

ValidatedControlMetadata validate_and_build_metadata(
    const SubcyclingODEMethod method, const double parent_dt,
    const std::vector<DenseMFabRawSampleView> &samples) {
  const auto &stencil = reference_dense_stencil(method);
  const auto &constraints = stencil.specification().constraints;

  if (!std::isfinite(parent_dt) || !(parent_dt > 0.0))
    throw std::invalid_argument(
        "dense MultiFab control parent dt must be finite and positive");
  if (samples.size() != constraints.size())
    throw std::invalid_argument(
        "dense MultiFab control sample count does not match stencil");

  std::unordered_set<const OwnedMultiFabDenseState *> unique_states;
  unique_states.reserve(samples.size());
  for (std::size_t index = 0; index < samples.size(); ++index) {
    if (!same_constraint(samples[index].constraint, constraints[index]))
      throw std::invalid_argument(
          "dense MultiFab sample metadata does not match stencil");
    if (samples[index].state == nullptr)
      throw std::invalid_argument("dense MultiFab sample state must not be null");
    if (!unique_states.insert(samples[index].state).second)
      throw std::invalid_argument("dense MultiFab sample states must be unique");
  }

  const auto &reference = *samples.front().state;
  for (const auto &sample : samples) {
    if (!reference.compatible(*sample.state) ||
        !sample.state->compatible(reference))
      throw std::invalid_argument(
          "dense MultiFab samples must be mutually compatible");
  }

  std::size_t begin_value_sample = samples.size();
  std::size_t end_value_sample = samples.size();
  for (std::size_t index = 0; index < constraints.size(); ++index) {
    if (constraints[index].kind != DenseSampleKind::value)
      continue;
    if (constraints[index].theta == 0.0)
      begin_value_sample = index;
    if (constraints[index].theta == 1.0)
      end_value_sample = index;
  }
  if (begin_value_sample == samples.size() ||
      end_value_sample == samples.size())
    throw std::logic_error("reference dense stencil lacks endpoint values");

  std::vector<std::vector<double>> effective_weights(
      stencil.control_count(), std::vector<double>(stencil.sample_count()));
  for (std::size_t control = 0; control < stencil.control_count(); ++control) {
    for (std::size_t sample = 0; sample < stencil.sample_count(); ++sample) {
      double effective = stencil.weight(control, sample);
      if (constraints[sample].kind == DenseSampleKind::scaled_derivative)
        effective *= parent_dt;
      if (!std::isfinite(effective))
        throw std::invalid_argument(
            "dense MultiFab effective coefficient is not finite");

      const double representable_lowest = static_cast<double>(
          std::numeric_limits<amrex::Real>::lowest());
      const double representable_max =
          static_cast<double>(std::numeric_limits<amrex::Real>::max());
      if (effective < representable_lowest || effective > representable_max)
        throw std::invalid_argument(
            "dense MultiFab coefficient is outside amrex Real range");
      const amrex::Real representable = static_cast<amrex::Real>(effective);
      if (!std::isfinite(representable))
        throw std::invalid_argument(
            "dense MultiFab coefficient is not representable by amrex Real");
      effective_weights[control][sample] =
          static_cast<double>(representable);
    }
  }

  return ValidatedControlMetadata{&stencil, begin_value_sample,
                                  end_value_sample,
                                  std::move(effective_weights)};
}

std::unique_ptr<OwnedMultiFabDenseState>
copy_control(const OwnedMultiFabDenseState &source) {
  auto result = OwnedMultiFabDenseState::empty_like(source);
  result->copy_from(source);
  return result;
}

std::unique_ptr<OwnedMultiFabDenseState> combination_control(
    const std::vector<double> &weights,
    const std::vector<DenseMFabRawSampleView> &samples) {
  std::vector<double> nonzero_weights;
  std::vector<const OwnedMultiFabDenseState *> nonzero_sources;
  nonzero_weights.reserve(weights.size());
  nonzero_sources.reserve(weights.size());
  for (std::size_t index = 0; index < weights.size(); ++index) {
    if (weights[index] == 0.0)
      continue;
    nonzero_weights.push_back(weights[index]);
    nonzero_sources.push_back(samples[index].state);
  }
  if (nonzero_sources.empty())
    throw std::logic_error("reference dense control has no nonzero weights");
  return OwnedMultiFabDenseState::linear_combination_of(nonzero_weights,
                                                         nonzero_sources);
}

} // namespace

DenseMFabControlSet::DenseMFabControlSet(
    const SubcyclingODEMethod method, const double parent_dt,
    std::vector<std::unique_ptr<OwnedMultiFabDenseState>> controls)
    : method_(method), parent_dt_(parent_dt), controls_(std::move(controls)) {}

SubcyclingODEMethod DenseMFabControlSet::method() const noexcept {
  return method_;
}

double DenseMFabControlSet::parent_dt() const noexcept { return parent_dt_; }

std::size_t DenseMFabControlSet::control_count() const noexcept {
  return controls_.size();
}

const OwnedMultiFabDenseState &
DenseMFabControlSet::control(const std::size_t index) const {
  return *controls_.at(index);
}

std::vector<std::unique_ptr<OwnedMultiFabDenseState>>
DenseMFabControlSet::release_controls() && {
  return std::move(controls_);
}

DenseMFabControlSet build_reference_dense_controls(
    const SubcyclingODEMethod method, const double parent_dt,
    const std::vector<DenseMFabRawSampleView> &samples) {
  const auto metadata =
      validate_and_build_metadata(method, parent_dt, samples);

  std::vector<std::unique_ptr<OwnedMultiFabDenseState>> controls;
  controls.reserve(metadata.stencil->control_count());
  controls.push_back(copy_control(*samples[metadata.begin_value_sample].state));
  for (std::size_t control = 1;
       control + 1 < metadata.stencil->control_count(); ++control) {
    controls.push_back(
        combination_control(metadata.effective_weights[control], samples));
  }
  controls.push_back(copy_control(*samples[metadata.end_value_sample].state));
  return DenseMFabControlSet(method, parent_dt, std::move(controls));
}

} // namespace CarpetX
