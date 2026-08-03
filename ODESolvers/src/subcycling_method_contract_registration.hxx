#ifndef ODESOLVERS_SUBCYCLING_METHOD_CONTRACT_REGISTRATION_HXX
#define ODESOLVERS_SUBCYCLING_METHOD_CONTRACT_REGISTRATION_HXX

#include "subcycling_ode_provider_registry.hxx"

#include <subcycling_method_contract.hxx>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace ODESolvers {

inline std::optional<CarpetX::SubcyclingMethodContract>
make_subcycling_method_contract(const bool subcycling_enabled,
                                const std::string_view method) {
  if (!subcycling_enabled)
    return std::nullopt;

  const auto &provider = require_subcycling_ode_provider(method);
  return CarpetX::SubcyclingMethodContract{
      "ODESolvers::method", std::string(provider.parameter_name),
      provider.dense};
}

inline void validate_subcycling_method_steering(
    const CarpetX::SubcyclingMethodContractSnapshot &snapshot,
    const std::string_view new_value) {
  if (snapshot.contract.provider_id != "ODESolvers::method")
    throw std::logic_error(
        "the active subcycling method contract belongs to another provider");
  if (new_value != snapshot.contract.parameter_value)
    throw std::logic_error(
        "ODESolvers::method is frozen while CarpetX subcycling is active");
}

} // namespace ODESolvers

#endif // ODESOLVERS_SUBCYCLING_METHOD_CONTRACT_REGISTRATION_HXX
