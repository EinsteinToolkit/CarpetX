#include "subcycling_method_contract_registration.hxx"

#include <subcycling_method_contract.hxx>

#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void check(const bool condition, const char *const message) {
  if (!condition)
    throw std::runtime_error(message);
}

template <typename Exception, typename Function>
void expect_throw(Function &&function, const char *const message) {
  try {
    std::forward<Function>(function)();
  } catch (const Exception &) {
    return;
  }
  throw std::runtime_error(message);
}

bool same_dense(const CarpetX::DenseCapability &left,
                const CarpetX::DenseCapability &right) {
  return left.method == right.method &&
         left.tableau_fingerprint == right.tableau_fingerprint &&
         left.endpoint_order == right.endpoint_order &&
         left.dense_uniform_order == right.dense_uniform_order &&
         left.stage_count == right.stage_count &&
         left.extra_rhs_evaluations == right.extra_rhs_evaluations &&
         left.persistent_vector_count == right.persistent_vector_count &&
         left.arbitrary_theta == right.arbitrary_theta &&
         left.verified == right.verified;
}

} // namespace

int main() {
  using CarpetX::SubcyclingGroupSchema;
  using CarpetX::SubcyclingMethodContractRegistry;
  using ODESolvers::make_subcycling_method_contract;
  using ODESolvers::validate_subcycling_method_steering;

  const auto disabled =
      make_subcycling_method_contract(false, "legacy-unsupported-method");
  check(!disabled.has_value(),
        "disabled subcycling must not interpret the legacy method");

  const auto enabled = make_subcycling_method_contract(true, "RK4");
  check(enabled.has_value(), "RK4 subcycling contract was not created");
  check(enabled->provider_id == "ODESolvers::method",
        "contract provider identity is not stable");
  check(enabled->parameter_value == "RK4",
        "contract did not preserve the frozen parameter value");
  check(enabled->dense.method == CarpetX::SubcyclingODEMethod::rk4,
        "contract did not preserve the CarpetX method");
  check(enabled->dense.verified && enabled->dense.arbitrary_theta,
        "contract is not a verified arbitrary-theta dense provider");

  SubcyclingMethodContractRegistry registry;
  check(registry.empty(), "new method-contract registry is not empty");
  expect_throw<std::out_of_range>(
      [&] { static_cast<void>(registry.require_snapshot()); },
      "empty registry returned a method contract");

  const auto epoch = registry.register_contract(*enabled);
  check(epoch == 1, "first method-contract registration epoch is not one");
  check(!registry.empty(), "registered method contract was not retained");

  const auto snapshot = registry.require_snapshot();
  check(snapshot.registration_epoch == epoch,
        "snapshot registration epoch differs from registration");
  check(snapshot.contract.provider_id == enabled->provider_id &&
            snapshot.contract.parameter_value == enabled->parameter_value &&
            same_dense(snapshot.contract.dense, enabled->dense),
        "snapshot differs from the registered contract");
  check(!snapshot.group_schema.has_value(),
        "method-only registration invented a group schema");
  registry.validate_snapshot(snapshot);
  expect_throw<std::out_of_range>(
      [&] { static_cast<void>(registry.require_complete_snapshot()); },
      "incomplete registry returned a transaction-ready snapshot");

  const auto duplicate_epoch = registry.register_contract(*enabled);
  check(duplicate_epoch == epoch,
        "identical duplicate registration changed the epoch");

  auto foreign_owner = *enabled;
  foreign_owner.provider_id = "AnotherODESolver";
  expect_throw<std::logic_error>(
      [&] { static_cast<void>(registry.register_contract(foreign_owner)); },
      "a second provider claimed the frozen contract");

  const auto rkf78 = make_subcycling_method_contract(true, "RKF78");
  check(rkf78.has_value(), "RKF78 subcycling contract was not created");
  expect_throw<std::logic_error>(
      [&] { static_cast<void>(registry.register_contract(*rkf78)); },
      "a mismatched method replaced the frozen contract");

  const SubcyclingGroupSchema schema{
      {{3, 8}, {11, 14}},
      {8, 14, 21},
  };
  const auto schema_epoch =
      registry.register_group_schema("ODESolvers::method", schema);
  check(schema_epoch == epoch + 1,
        "first group-schema registration did not advance the epoch");
  const auto complete = registry.require_complete_snapshot();
  check(complete.registration_epoch == schema_epoch &&
            complete.group_schema.has_value(),
        "complete snapshot omitted the registered group schema");
  check(complete.group_schema->ordered_group_pairs.size() == 2 &&
            complete.group_schema->ordered_group_pairs[0].evolved_group == 3 &&
            complete.group_schema->ordered_group_pairs[0].rhs_group == 8 &&
            complete.group_schema->ordered_group_pairs[1].evolved_group == 11 &&
            complete.group_schema->ordered_group_pairs[1].rhs_group == 14,
        "complete snapshot changed ordered evolved/RHS pairs");
  check(complete.group_schema->dependent_groups ==
            std::vector<int>({8, 14, 21}),
        "complete snapshot changed the sorted dependent group set");
  check(registry.register_group_schema("ODESolvers::method", schema) ==
            schema_epoch,
        "identical duplicate group schema changed the epoch");

  auto changed_schema = schema;
  changed_schema.dependent_groups.push_back(34);
  expect_throw<std::logic_error>(
      [&] {
        static_cast<void>(registry.register_group_schema(
            "ODESolvers::method", changed_schema));
      },
      "mismatched group schema replaced the frozen schema");
  expect_throw<std::logic_error>(
      [&] {
        static_cast<void>(
            registry.register_group_schema("AnotherODESolver", schema));
      },
      "foreign provider registered a group schema");

  validate_subcycling_method_steering(complete, "RK4");
  expect_throw<std::logic_error>(
      [&] { validate_subcycling_method_steering(complete, "RKF78"); },
      "method steering changed a frozen subcycling contract");

  expect_throw<std::logic_error>(
      [&] { registry.validate_snapshot(snapshot); },
      "registry accepted a snapshot made stale by schema registration");

  auto tampered = complete;
  tampered.contract.parameter_value = "DP87";
  expect_throw<std::logic_error>(
      [&] { registry.validate_snapshot(tampered); },
      "registry accepted a modified read-only snapshot");

  std::cout << "Subcycling method-contract handoff tests passed\n";
  return 0;
}
