#include "subcycling_method_contract_registration.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameter.h>
#include <cctk_Parameters.h>

#include <subcycling_method_contract.hxx>

#include <exception>
#include <stdexcept>

namespace ODESolvers {
namespace {

constexpr const char *steering_guard_name =
    "ODESolvers_SUBCYCLING_METHOD_STEERING_GUARD";
bool steering_guard_registered = false;

bool carpetx_subcycling_enabled() {
  const auto *const enabled = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("use_subcycling_wip", "CarpetX", nullptr));
  if (enabled == nullptr)
    throw std::logic_error(
        "CarpetX::use_subcycling_wip parameter is unavailable");
  return *enabled != 0;
}

void method_steering_callback(void *, const char *, const char *,
                              const char *const new_value) {
  try {
    const auto snapshot =
        CarpetX::require_registered_subcycling_method_contract();
    validate_subcycling_method_steering(snapshot, new_value);
  } catch (const std::exception &error) {
    CCTK_VERROR("Rejected ODESolvers::method steering while CarpetX "
                "subcycling is active: %s",
                error.what());
  } catch (...) {
    CCTK_ERROR("Rejected ODESolvers::method steering while CarpetX "
               "subcycling is active with an unknown contract error");
  }
}

} // namespace
} // namespace ODESolvers

extern "C" void
ODESolvers_RegisterSubcyclingMethodContract(CCTK_ARGUMENTS) {
  static_cast<void>(cctkGH);
  DECLARE_CCTK_PARAMETERS;

  try {
    const bool enabled = ODESolvers::carpetx_subcycling_enabled();
    const auto contract =
        ODESolvers::make_subcycling_method_contract(enabled, method);
    if (!contract)
      return;

    CarpetX::register_subcycling_method_contract(*contract);
    if (!ODESolvers::steering_guard_registered) {
      const int status = CCTK_ParameterSetNotifyRegister(
          ODESolvers::method_steering_callback, nullptr,
          ODESolvers::steering_guard_name, "^ODESolvers$", "^method$");
      if (status != 0)
        throw std::runtime_error(
            "could not register the ODESolvers method steering guard");
      ODESolvers::steering_guard_registered = true;
    }
  } catch (const std::exception &error) {
    CCTK_VERROR("CarpetX subcycling method-contract registration failed for "
                "ODESolvers::method=\"%s\": %s",
                method, error.what());
  } catch (...) {
    CCTK_VERROR("CarpetX subcycling method-contract registration failed for "
                "ODESolvers::method=\"%s\" with an unknown error",
                method);
  }
}

extern "C" void
ODESolvers_UnregisterSubcyclingMethodSteeringGuard(CCTK_ARGUMENTS) {
  static_cast<void>(cctkGH);
  if (!ODESolvers::steering_guard_registered)
    return;
  const int status =
      CCTK_ParameterSetNotifyUnregister(ODESolvers::steering_guard_name);
  if (status != 0)
    CCTK_VWARN(CCTK_WARN_ALERT,
               "Could not unregister the ODESolvers subcycling method "
               "steering guard: %d",
               status);
  ODESolvers::steering_guard_registered = false;
}
