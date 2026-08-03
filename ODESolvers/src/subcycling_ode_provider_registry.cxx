#include "subcycling_ode_provider_registry.hxx"

#include "explicit_rk_dense_provider.hxx"

#include <stdexcept>
#include <utility>

namespace ODESolvers {
namespace {

template <ExplicitRKMethod Method>
ExplicitRKStagePoint
exact_stage_metadata(const ExplicitRKAdvanceFrame &frame,
                     const int stage_index) {
  return explicit_rk_stage_point(Method, frame, stage_index);
}

SubcyclingODEProviderCapability
make_capability(const std::string_view parameter_name,
                const ExplicitRKMethod method,
                const ExactStageMetadataProvider metadata_provider) {
  const auto &tableau = explicit_rk_tableau(method);
  auto dense = reference_dense_capability(method);
  const auto expected_stage_count = static_cast<int>(tableau.a.size());

  if (parameter_name.empty() || metadata_provider == nullptr ||
      dense.method != subcycling_method(method) ||
      dense.tableau_fingerprint != explicit_rk_tableau_fingerprint(method) ||
      dense.endpoint_order != tableau.endpoint_order ||
      dense.dense_uniform_order < tableau.endpoint_order ||
      dense.stage_count != expected_stage_count || !dense.arbitrary_theta ||
      !dense.verified)
    throw std::logic_error(
        "inconsistent subcycling ODE provider capability");

  return SubcyclingODEProviderCapability{parameter_name, method,
                                         std::move(dense),
                                         metadata_provider};
}

} // namespace

const SubcyclingODEProviderRegistry &subcycling_ode_provider_registry() {
  static const SubcyclingODEProviderRegistry registry{{
      make_capability("RK4", ExplicitRKMethod::rk4,
                      &exact_stage_metadata<ExplicitRKMethod::rk4>),
      make_capability("RKF78", ExplicitRKMethod::rkf78,
                      &exact_stage_metadata<ExplicitRKMethod::rkf78>),
      make_capability("DP87", ExplicitRKMethod::dp87,
                      &exact_stage_metadata<ExplicitRKMethod::dp87>),
  }};
  return registry;
}

const SubcyclingODEProviderCapability &
require_subcycling_ode_provider(const ExplicitRKMethod method) {
  for (const auto &entry : subcycling_ode_provider_registry())
    if (entry.method == method)
      return entry;
  throw std::invalid_argument("unsupported subcycling ODE method");
}

const SubcyclingODEProviderCapability &
require_subcycling_ode_provider(const std::string_view parameter_name) {
  for (const auto &entry : subcycling_ode_provider_registry())
    if (entry.parameter_name == parameter_name)
      return entry;
  throw std::invalid_argument("unsupported subcycling ODE method");
}

} // namespace ODESolvers
