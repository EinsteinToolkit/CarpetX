#ifndef ODESOLVERS_SUBCYCLING_ODE_PROVIDER_REGISTRY_HXX
#define ODESOLVERS_SUBCYCLING_ODE_PROVIDER_REGISTRY_HXX

#include "explicit_rk.hxx"

#include <subcycling_dense_output.hxx>

#include <array>
#include <string_view>

namespace ODESolvers {

using ExactStageMetadataProvider = ExplicitRKStagePoint (*)(
    const ExplicitRKAdvanceFrame &frame, int stage_index);

struct SubcyclingODEProviderCapability {
  std::string_view parameter_name;
  ExplicitRKMethod method;
  CarpetX::DenseCapability dense;
  ExactStageMetadataProvider exact_stage_metadata;
};

using SubcyclingODEProviderRegistry =
    std::array<SubcyclingODEProviderCapability, 3>;

const SubcyclingODEProviderRegistry &subcycling_ode_provider_registry();

const SubcyclingODEProviderCapability &
require_subcycling_ode_provider(ExplicitRKMethod method);

const SubcyclingODEProviderCapability &
require_subcycling_ode_provider(std::string_view parameter_name);

} // namespace ODESolvers

#endif // ODESOLVERS_SUBCYCLING_ODE_PROVIDER_REGISTRY_HXX
