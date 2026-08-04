#ifndef CARPETX_SUBCYCLING_STATIC_V1_CONTRACT_HXX
#define CARPETX_SUBCYCLING_STATIC_V1_CONTRACT_HXX

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace CarpetX {

enum class StaticV1EvolvePath { legacy, static_v1 };

constexpr StaticV1EvolvePath
select_static_v1_evolve_path(const bool use_subcycling_wip) noexcept {
  return use_subcycling_wip ? StaticV1EvolvePath::static_v1
                           : StaticV1EvolvePath::legacy;
}

struct StaticV1LevelOneCreationEnvelope {
  bool subcycling_enabled;
  int configured_patch_count;
  int callback_patch;
  int callback_level;
  int configured_level_count;
  int existing_level_count;
  int regrid_every;
  int coarse_patch;
  int coarse_level;
  struct StaticV1Rational {
    std::int64_t numerator;
    std::int64_t denominator;
  };
  StaticV1Rational coarse_iteration;
  StaticV1Rational coarse_delta_iteration;
  bool coarse_is_subcycling_level;
  std::array<int, 3> spatial_refinement_ratio;
  double callback_time;
};

constexpr bool permits_static_v1_level_one_creation(
    const StaticV1LevelOneCreationEnvelope &envelope) noexcept {
  return envelope.subcycling_enabled &&
         envelope.configured_patch_count == 1 &&
         envelope.callback_patch == 0 && envelope.callback_level == 1 &&
         envelope.configured_level_count == 2 &&
         envelope.existing_level_count == 1 && envelope.regrid_every == 0 &&
         envelope.coarse_patch == 0 && envelope.coarse_level == 0 &&
         envelope.coarse_iteration.numerator == 0 &&
         envelope.coarse_iteration.denominator == 1 &&
         envelope.coarse_delta_iteration.numerator == 1 &&
         envelope.coarse_delta_iteration.denominator == 1 &&
         !envelope.coarse_is_subcycling_level &&
         envelope.spatial_refinement_ratio[0] == 2 &&
         envelope.spatial_refinement_ratio[1] == 2 &&
         envelope.spatial_refinement_ratio[2] == 2 &&
         envelope.callback_time == 0.0;
}

using StaticV1Rational = StaticV1LevelOneCreationEnvelope::StaticV1Rational;

constexpr int
static_v1_parent_refinement_ratio_index(const int callback_level) noexcept {
  return callback_level - 1;
}

constexpr bool should_stop_static_v1_initial_regrid_loop(
    const bool subcycling_enabled, const int configured_patch_count,
    const int configured_level_count, const int existing_level_count) noexcept {
  return subcycling_enabled && configured_patch_count == 1 &&
         configured_level_count == 2 && existing_level_count == 2;
}

struct StaticV1PolicyEnvelope {
  int patch_count;
  int level_count;
  int spatial_refinement_ratio;
  int regrid_every;
  bool recovering;
  bool do_reflux;
  bool do_restrict;
  bool restrict_during_sync;
  bool poison_undefined_values;
  bool level_one_subcycles;
  bool root_iteration_zero;
  bool level_clocks_zero;
  bool complete_method_schema;
  bool evolved_schema_supported;
  bool bounded_test_odesolvers2_configuration;
};

inline void
validate_static_v1_policy_envelope(const StaticV1PolicyEnvelope &policy) {
  if (policy.patch_count != 1)
    throw std::invalid_argument(
        "static-v1 requires exactly one CarpetX patch");
  if (policy.level_count != 2)
    throw std::invalid_argument(
        "static-v1 requires exactly two CarpetX levels");
  if (policy.spatial_refinement_ratio != 2)
    throw std::invalid_argument(
        "static-v1 requires factor-two spatial refinement");
  if (policy.regrid_every != 0)
    throw std::invalid_argument(
        "static-v1 requires CarpetX::regrid_every=0");
  if (policy.do_reflux)
    throw std::invalid_argument(
        "static-v1 requires CarpetX::do_reflux=no");
  if (!policy.do_restrict)
    throw std::invalid_argument(
        "static-v1 requires CarpetX::do_restrict=yes");
  if (policy.restrict_during_sync)
    throw std::invalid_argument(
        "static-v1 requires CarpetX::restrict_during_sync=no");
  if (policy.poison_undefined_values)
    throw std::invalid_argument(
        "static-v1 requires CarpetX::poison_undefined_values=no");
  if (!policy.level_one_subcycles)
    throw std::invalid_argument(
        "static-v1 requires level one to be marked for subcycling");
  if (policy.recovering) {
    if (!policy.level_clocks_zero)
      throw std::invalid_argument(
          "static-v1 recovery requires freshly rebuilt zero level clocks");
  } else if (!policy.root_iteration_zero || !policy.level_clocks_zero) {
    throw std::invalid_argument(
        "static-v1 requires zero root and level clocks at a new start");
  }
  if (!policy.complete_method_schema)
    throw std::invalid_argument(
        "static-v1 requires a frozen complete ODE method/group schema");
  if (!policy.evolved_schema_supported)
    throw std::invalid_argument(
        "static-v1 requires GF REAL evolved, RHS, and dependent groups with "
        "one evolved timelevel");
  if (!policy.bounded_test_odesolvers2_configuration)
    throw std::invalid_argument(
        "static-v1 outer PRESTEP/EVOL schedule certification is not yet "
        "implemented; only the exact two-level TestODESolvers2 "
        "state/RHS/dependent schema is supported, and BBH or general thorn "
        "configurations are not supported");
}

inline bool has_exact_static_v1_test_odesolvers2_active_thorns(
    std::vector<std::string> active_thorns) {
  static const std::array<const char *, 14> p1_expected{
      "AMReX",          "Arith", "BoxInBox",    "Cactus", "CarpetX",
      "CarpetXRegrid",  "IOUtil", "Loop",        "MPI",    "NSIMD",
      "ODESolvers",     "TestODESolvers2", "yaml_cpp", "zlib"};
  static const std::array<const char *, 16> p2_checkpoint_expected{
      "AMReX",         "Arith",       "BoxInBox", "Cactus",
      "CarpetX",       "CarpetXRegrid", "HDF5",   "IOUtil",
      "Loop",          "MPI",         "NSIMD",    "ODESolvers",
      "Silo",          "TestODESolvers2", "yaml_cpp", "zlib"};
  std::sort(active_thorns.begin(), active_thorns.end());
  const auto matches = [&](const auto &expected) {
    return active_thorns.size() == expected.size() &&
           std::equal(active_thorns.begin(), active_thorns.end(),
                      expected.begin(), [](const std::string &active,
                                           const char *const required) {
                        return active == required;
                      });
  };
  return matches(p1_expected) || matches(p2_checkpoint_expected);
}

// Accepted active-level GF history must rotate before the canonical TL0
// snapshot is copied. Keeping this sequencing helper Cactus-free makes the
// exact once-and-before contract directly regression-testable.
template <class Cycle, class Capture>
decltype(auto)
cycle_then_capture_static_v1_level_state(Cycle &&cycle, Capture &&capture) {
  std::forward<Cycle>(cycle)();
  return std::forward<Capture>(capture)();
}

enum class StaticV1SyncObserverAction {
  cycle_synchronized_globals,
  postrestrict,
  poststep,
  analysis,
  output,
  checkpoint,
};

constexpr std::array<StaticV1SyncObserverAction, 6>
static_v1_sync_observer_order() noexcept {
  return {StaticV1SyncObserverAction::cycle_synchronized_globals,
          StaticV1SyncObserverAction::postrestrict,
          StaticV1SyncObserverAction::poststep,
          StaticV1SyncObserverAction::analysis,
          StaticV1SyncObserverAction::output,
          StaticV1SyncObserverAction::checkpoint};
}

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_STATIC_V1_CONTRACT_HXX
