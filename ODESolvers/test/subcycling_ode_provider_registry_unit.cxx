#include "subcycling_ode_provider_registry.hxx"

#include <array>
#include <cassert>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using ODESolvers::ExplicitRKAdvanceFrame;
using ODESolvers::ExplicitRKMethod;
using ODESolvers::ExplicitRKStageKind;
using ODESolvers::RationalCoefficient;

template <class Exception, class Function>
void require_throws(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Exception &) {
    threw = true;
  }
  assert(threw);
}

bool same_rational(const RationalCoefficient value, const std::int64_t numerator,
                   const std::int64_t denominator) {
  return value.numerator == numerator && value.denominator == denominator;
}

struct MethodExpectation {
  std::string_view name;
  ExplicitRKMethod method;
  CarpetX::SubcyclingODEMethod dense_method;
  int endpoint_order;
  int dense_order;
  std::vector<RationalCoefficient> stage_fractions;
};

const std::array<MethodExpectation, 3> expectations{{
    {"RK4",
     ExplicitRKMethod::rk4,
     CarpetX::SubcyclingODEMethod::rk4,
     4,
     4,
     {{0, 1}, {1, 2}, {1, 2}, {1, 1}}},
    {"RKF78",
     ExplicitRKMethod::rkf78,
     CarpetX::SubcyclingODEMethod::rkf78_order7,
     7,
     7,
     {{0, 1}, {2, 27}, {1, 9}, {1, 6}, {5, 12}, {1, 2},
      {5, 6}, {1, 6}, {2, 3}, {1, 3}, {1, 1}}},
    {"DP87",
     ExplicitRKMethod::dp87,
     CarpetX::SubcyclingODEMethod::dp87_order8,
     8,
     8,
     {{0, 1},
      {1, 18},
      {1, 12},
      {1, 8},
      {5, 16},
      {3, 8},
      {59, 400},
      {93, 200},
      {5490023248LL, 9719169821LL},
      {13, 20},
      {1201146811LL, 1299019798LL},
      {1, 1},
      {1, 1}}},
}};

void test_registry_exposes_only_exact_supported_capabilities() {
  const auto &registry = ODESolvers::subcycling_ode_provider_registry();
  assert(registry.size() == expectations.size());

  for (std::size_t index = 0; index < expectations.size(); ++index) {
    const auto &expected = expectations[index];
    const auto &entry = registry[index];
    assert(entry.parameter_name == expected.name);
    assert(entry.method == expected.method);
    assert(entry.dense.method == expected.dense_method);
    assert(entry.dense.endpoint_order == expected.endpoint_order);
    assert(entry.dense.dense_uniform_order == expected.dense_order);
    assert(entry.dense.stage_count ==
           static_cast<int>(expected.stage_fractions.size()));
    assert(entry.dense.arbitrary_theta);
    assert(entry.dense.verified);
    assert(entry.exact_stage_metadata != nullptr);
    assert(entry.dense.tableau_fingerprint ==
           ODESolvers::explicit_rk_tableau_fingerprint(expected.method));

    assert(&ODESolvers::require_subcycling_ode_provider(expected.name) ==
           &entry);
    assert(&ODESolvers::require_subcycling_ode_provider(expected.method) ==
           &entry);
  }
}

void test_stage_metadata_provider_routes_every_exact_stage_without_cloning() {
  const ExplicitRKAdvanceFrame primary{ExplicitRKStageKind::primary,
                                       {0, 1}, {1, 1}};
  for (const auto &expected : expectations) {
    const auto &entry =
        ODESolvers::require_subcycling_ode_provider(expected.method);
    for (std::size_t stage = 0; stage < expected.stage_fractions.size();
         ++stage) {
      const auto point =
          entry.exact_stage_metadata(primary, static_cast<int>(stage + 1));
      assert(point.kind == ExplicitRKStageKind::primary);
      assert(point.stage_index == static_cast<int>(stage + 1));
      assert(point.stage_count ==
             static_cast<int>(expected.stage_fractions.size()));
      assert(point.parent_fraction.numerator ==
             expected.stage_fractions[stage].numerator);
      assert(point.parent_fraction.denominator ==
             expected.stage_fractions[stage].denominator);
    }

    const ExplicitRKAdvanceFrame fractional{ExplicitRKStageKind::fractional,
                                            {1, 4}, {1, 2}};
    const auto first = entry.exact_stage_metadata(fractional, 1);
    const auto last = entry.exact_stage_metadata(
        fractional, static_cast<int>(expected.stage_fractions.size()));
    assert(first.kind == ExplicitRKStageKind::fractional);
    assert(same_rational(first.parent_fraction, 1, 4));
    assert(same_rational(last.parent_fraction, 3, 4));
  }
}

void test_unsupported_methods_and_invalid_stage_metadata_fail_closed() {
  for (const std::string_view name : {std::string_view{},
                                      std::string_view{"rk4"},
                                      std::string_view{"RK45"},
                                      std::string_view{"SSPRK3"}})
    require_throws<std::invalid_argument>([name] {
      (void)ODESolvers::require_subcycling_ode_provider(name);
    });

  require_throws<std::invalid_argument>([] {
    (void)ODESolvers::require_subcycling_ode_provider(
        static_cast<ExplicitRKMethod>(99));
  });

  const auto &entry =
      ODESolvers::require_subcycling_ode_provider(ExplicitRKMethod::rk4);
  const ExplicitRKAdvanceFrame primary{ExplicitRKStageKind::primary,
                                       {0, 1}, {1, 1}};
  require_throws<std::out_of_range>(
      [&] { (void)entry.exact_stage_metadata(primary, 0); });
  require_throws<std::out_of_range>([&] {
    (void)entry.exact_stage_metadata(primary, entry.dense.stage_count + 1);
  });
}

} // namespace

int main() {
  test_registry_exposes_only_exact_supported_capabilities();
  test_stage_metadata_provider_routes_every_exact_stage_without_cloning();
  test_unsupported_methods_and_invalid_stage_metadata_fail_closed();
}
