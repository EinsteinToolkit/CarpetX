#include "subcycling_dense_output.hxx"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace CarpetX {
namespace {

bool known_method(const SubcyclingODEMethod method) noexcept {
  switch (method) {
  case SubcyclingODEMethod::rk4:
  case SubcyclingODEMethod::rkf78_order7:
  case SubcyclingODEMethod::dp87_order8:
    return true;
  }
  return false;
}

bool zero_fingerprint(const TableauFingerprint &fingerprint) noexcept {
  return std::all_of(fingerprint.begin(), fingerprint.end(),
                     [](const std::uint8_t value) { return value == 0; });
}

void validate_capability(const DenseCapability &capability) {
  if (!known_method(capability.method))
    throw std::invalid_argument("dense provider has an unknown ODE method");
  if (zero_fingerprint(capability.tableau_fingerprint))
    throw std::invalid_argument("dense provider fingerprint must be non-zero");
  if (capability.endpoint_order <= 0 || capability.dense_uniform_order <= 0)
    throw std::invalid_argument("dense provider orders must be positive");
  if (capability.dense_uniform_order < capability.endpoint_order)
    throw std::invalid_argument(
        "dense uniform order must not be below endpoint order");
  if (capability.stage_count <= 0 || capability.extra_rhs_evaluations < 0)
    throw std::invalid_argument("dense provider cost must be non-negative");
  if (capability.persistent_vector_count <= 0 ||
      static_cast<long long>(capability.persistent_vector_count) !=
          static_cast<long long>(capability.dense_uniform_order) + 1LL)
    throw std::invalid_argument(
        "dense provider must own one Bernstein control per polynomial degree");
  if (!capability.arbitrary_theta)
    throw std::invalid_argument(
        "dense provider must support arbitrary interior theta");
  if (!capability.verified)
    throw std::invalid_argument("dense provider is not verified");
}

void validate_interval_id(const DenseCapability &capability,
                          const DenseIntervalId &id) {
  if (id.level < 0)
    throw std::invalid_argument("dense interval level must be non-negative");
  if (!(id.begin_clock < id.end_clock))
    throw std::invalid_argument("dense interval clocks must increase exactly");
  if (!std::isfinite(id.begin_time) || !std::isfinite(id.end_time) ||
      !(id.begin_time < id.end_time))
    throw std::invalid_argument(
        "dense interval physical times must be finite and increase");
  if (id.method != capability.method)
    throw std::invalid_argument(
        "dense interval method does not match provider capability");
  if (id.tableau_fingerprint != capability.tableau_fingerprint)
    throw std::invalid_argument(
        "dense interval fingerprint does not match provider capability");
}

double binomial(const int n, const int k) {
  if (k < 0 || k > n)
    return 0.0;
  const int smaller = std::min(k, n - k);
  double result = 1.0;
  for (int i = 1; i <= smaller; ++i)
    result *= static_cast<double>(n - smaller + i) / static_cast<double>(i);
  return result;
}

} // namespace

DenseInterval::DenseInterval(
    DenseCapability capability, DenseIntervalId id,
    std::vector<std::unique_ptr<DenseStateVector>> controls)
    : capability_(std::move(capability)), id_(std::move(id)),
      controls_(std::move(controls)) {}

const DenseIntervalId &DenseInterval::id() const noexcept { return id_; }

const DenseCapability &DenseInterval::capability() const noexcept {
  return capability_;
}

std::size_t DenseInterval::control_count() const noexcept {
  return controls_.size();
}

void DenseInterval::evaluate(const double theta,
                             DenseStateVector &destination) const {
  if (!std::isfinite(theta))
    throw std::invalid_argument("dense evaluation theta must be finite");
  if (theta < 0.0 || theta > 1.0)
    throw std::out_of_range("dense evaluation theta must be in [0, 1]");
  if (controls_.empty())
    throw std::logic_error("dense interval has no controls");
  if (!destination.compatible(*controls_.front()) ||
      !controls_.front()->compatible(destination))
    throw std::invalid_argument(
        "dense evaluation destination is incompatible with controls");

  if (theta == 0.0) {
    destination.copy_from(*controls_.front());
    return;
  }
  if (theta == 1.0) {
    destination.copy_from(*controls_.back());
    return;
  }

  const int degree = static_cast<int>(controls_.size()) - 1;
  std::vector<double> weights;
  std::vector<const DenseStateVector *> sources;
  weights.reserve(controls_.size());
  sources.reserve(controls_.size());
  for (int k = 0; k <= degree; ++k) {
    weights.push_back(binomial(degree, k) * std::pow(theta, k) *
                      std::pow(1.0 - theta, degree - k));
    sources.push_back(controls_[static_cast<std::size_t>(k)].get());
  }
  destination.linear_combination(weights, sources);
}

void DenseInterval::evaluate(const step_clock_t theta,
                             DenseStateVector &destination) const {
  if (theta.den <= 0)
    throw std::invalid_argument(
        "exact dense evaluation theta must have a positive denominator");
  if (theta < 0 || theta > 1)
    throw std::out_of_range(
        "exact dense evaluation theta must be in [0, 1]");
  evaluate(static_cast<double>(theta), destination);
}

DenseIntervalBuilder::DenseIntervalBuilder(DenseCapability capability,
                                           DenseIntervalId id)
    : capability_(std::move(capability)), id_(std::move(id)) {
  validate_capability(capability_);
  validate_interval_id(capability_, id_);
  controls_.reserve(
      static_cast<std::size_t>(capability_.persistent_vector_count));
}

DenseIntervalBuilder::DenseIntervalBuilder(
    DenseIntervalBuilder &&other) noexcept
    : capability_(std::move(other.capability_)), id_(std::move(other.id_)),
      controls_(std::move(other.controls_)), state_(other.state_) {
  other.state_ = State::moved_from;
}

DenseIntervalBuilder &
DenseIntervalBuilder::operator=(DenseIntervalBuilder &&other) {
  require_active();
  other.require_active();
  if (this == &other)
    return *this;
  if (!controls_.empty())
    throw std::logic_error(
        "dense interval move-assignment target must be empty");

  capability_ = std::move(other.capability_);
  id_ = std::move(other.id_);
  controls_ = std::move(other.controls_);
  state_ = State::active;
  other.state_ = State::moved_from;
  return *this;
}

void DenseIntervalBuilder::require_active() const {
  if (state_ == State::sealed)
    throw std::logic_error("dense interval builder is already sealed");
  if (state_ == State::moved_from)
    throw std::logic_error("dense interval builder was moved from");
}

void DenseIntervalBuilder::add_control(
    std::unique_ptr<DenseStateVector> control) {
  require_active();
  if (control == nullptr)
    throw std::invalid_argument("dense interval control must be owned");
  if (controls_.size() >=
      static_cast<std::size_t>(capability_.persistent_vector_count))
    throw std::length_error("dense interval has excess controls");
  if (!controls_.empty() &&
      (!controls_.front()->compatible(*control) ||
       !control->compatible(*controls_.front())))
    throw std::invalid_argument("dense interval controls are incompatible");
  controls_.push_back(std::move(control));
}

std::shared_ptr<const DenseInterval> DenseIntervalBuilder::seal() {
  require_active();
  if (controls_.size() !=
      static_cast<std::size_t>(capability_.persistent_vector_count))
    throw std::logic_error("dense interval is missing controls");

  auto result = std::shared_ptr<const DenseInterval>(new DenseInterval(
      capability_, id_, std::move(controls_)));
  state_ = State::sealed;
  return result;
}

DenseOutputProvider::DenseOutputProvider(DenseCapability capability)
    : capability_(std::move(capability)) {
  validate_capability(capability_);
}

const DenseCapability &DenseOutputProvider::capability() const noexcept {
  return capability_;
}

std::unique_ptr<DenseIntervalBuilder> DenseOutputProvider::begin_interval(
    const DenseIntervalId &id) const {
  return std::make_unique<DenseIntervalBuilder>(capability_, id);
}

void DenseOutputRegistry::register_provider(
    std::shared_ptr<const DenseOutputProvider> provider) {
  if (provider == nullptr)
    throw std::invalid_argument("cannot register a null dense provider");
  const auto capability = provider->capability();
  validate_capability(capability);
  const auto duplicate =
      std::find_if(entries_.begin(), entries_.end(), [&](const Entry &entry) {
        return entry.capability.method == capability.method;
      });
  if (duplicate != entries_.end())
    throw std::logic_error("a dense provider already owns this ODE method");
  entries_.push_back(Entry{capability, std::move(provider)});
}

std::shared_ptr<const DenseOutputProvider> DenseOutputRegistry::require(
    const SubcyclingODEMethod method,
    const TableauFingerprint &tableau_fingerprint) const {
  const auto entry =
      std::find_if(entries_.begin(), entries_.end(), [&](const Entry &candidate) {
        return candidate.capability.method == method &&
               candidate.capability.tableau_fingerprint == tableau_fingerprint;
      });
  if (entry == entries_.end())
    throw std::out_of_range(
        "no verified dense provider matches method and fingerprint");
  return entry->provider;
}

} // namespace CarpetX
