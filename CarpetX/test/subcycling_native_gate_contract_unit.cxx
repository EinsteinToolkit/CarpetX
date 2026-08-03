#include "../src/subcycling_native_gate.hxx"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void check(const bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

template <class Function> void rejects(Function &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  check(rejected, "unsupported gate iteration was accepted");
}

} // namespace

int main() {
  const auto rk4 = CarpetX::native_gate_method_contract(1);
  check(rk4.method == CarpetX::SubcyclingODEMethod::rk4, "RK4 method");
  check(std::string(rk4.ode_parameter_value) == "RK4", "RK4 name");
  check(rk4.extra_rhs_evaluations == 4, "RK4 RHS count");
  check(rk4.control_count == 5, "RK4 controls");

  const auto rkf78 = CarpetX::native_gate_method_contract(2);
  check(rkf78.method == CarpetX::SubcyclingODEMethod::rkf78_order7,
        "RKF78 method");
  check(std::string(rkf78.ode_parameter_value) == "RKF78", "RKF78 name");
  check(rkf78.extra_rhs_evaluations == 23, "RKF78 RHS count");
  check(rkf78.control_count == 8, "RKF78 controls");

  const auto dp87 = CarpetX::native_gate_method_contract(3);
  check(dp87.method == CarpetX::SubcyclingODEMethod::dp87_order8,
        "DP87 method");
  check(std::string(dp87.ode_parameter_value) == "DP87", "DP87 name");
  check(dp87.extra_rhs_evaluations == 39, "DP87 RHS count");
  check(dp87.control_count == 9, "DP87 controls");

  rejects([] { (void)CarpetX::native_gate_method_contract(0); });
  rejects([] { (void)CarpetX::native_gate_method_contract(4); });

  const CarpetX::native_gate_detail::PatchContextInput patch{
      0, 1, 0, 1, 2, 1, CarpetX::step_clock_t(2),
      CarpetX::step_clock_t(1), 1.25, 0.125};
  const auto context = CarpetX::native_gate_detail::make_patch_step_context(
      patch, CarpetX::SubcyclingODEMethod::rkf78_order7);
  check(context.level == 0, "patch level");
  check(context.begin_clock == CarpetX::step_clock_t(1),
        "patch begin clock");
  check(context.end_clock == CarpetX::step_clock_t(2), "patch end clock");
  check(context.begin_time == 1.125, "patch begin time");
  check(context.end_time == 1.25, "patch end time");
  check(context.method == CarpetX::SubcyclingODEMethod::rkf78_order7,
        "patch method");

  auto bad_active_levels = patch;
  bad_active_levels.active_max_level = 2;
  rejects([&] {
    (void)CarpetX::native_gate_detail::make_patch_step_context(
        bad_active_levels, CarpetX::SubcyclingODEMethod::rk4);
  });
  auto bad_clock = patch;
  bad_clock.level_iteration = CarpetX::step_clock_t(3);
  rejects([&] {
    (void)CarpetX::native_gate_detail::make_patch_step_context(
        bad_clock, CarpetX::SubcyclingODEMethod::rk4);
  });
  auto bad_timefac = patch;
  bad_timefac.time_refinement_factor = 2;
  rejects([&] {
    (void)CarpetX::native_gate_detail::make_patch_step_context(
        bad_timefac, CarpetX::SubcyclingODEMethod::rk4);
  });
  std::cout << "Phase 8D native gate contract tests passed\n";
  return 0;
}
