#include "subcycling_schedule_state.hxx"

#include <cassert>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace Internal = CarpetX::ScheduleInternal;

struct ActiveRange {
  int first;
  int last;
};

struct FakeGH {
  int cctk_iteration{};
  double cctk_time{};
  double cctk_delta_time{};
  int cctk_timefac{};
};

void assert_metadata(const FakeGH &gh, const int iteration, const double time,
                     const double delta_time, const int timefac) {
  assert(gh.cctk_iteration == iteration);
  assert(gh.cctk_time == time);
  assert(gh.cctk_delta_time == delta_time);
  assert(gh.cctk_timefac == timefac);
}

void test_scoped_active_levels_restore_nested_and_exceptional_state() {
  std::optional<ActiveRange> active(ActiveRange{0, 5});
  {
    Internal::ScopedActiveLevels<ActiveRange> outer(active,
                                                     ActiveRange{1, 4});
    assert(active->first == 1 && active->last == 4);
    {
      Internal::ScopedActiveLevels<ActiveRange> inner(active,
                                                       ActiveRange{2, 3});
      assert(active->first == 2 && active->last == 3);
    }
    assert(active->first == 1 && active->last == 4);
  }
  assert(active->first == 0 && active->last == 5);

  try {
    Internal::ScopedActiveLevels<ActiveRange> scope(active,
                                                     ActiveRange{3, 3});
    throw std::runtime_error("requested unwind");
  } catch (const std::runtime_error &) {
  }
  assert(active->first == 0 && active->last == 5);
}

void test_level_time_metadata_copies_timefac_and_scope_restores_every_field() {
  FakeGH source{7, 1.25, 0.5, 8};
  FakeGH target{2, 0.25, 1.0, 1};
  Internal::copy_level_time_metadata(target, source);
  assert_metadata(target, 7, 1.25, 0.5, 8);

  std::vector<Internal::LevelTimeMetadata<double>> propagated;
  auto propagate = [&propagated](FakeGH &gh) noexcept {
    propagated.push_back(Internal::capture_level_time_metadata(gh));
  };
  const Internal::LevelTimeMetadata<double> scoped{11, 2.5, 0.125, 16};
  {
    Internal::ScopedLevelTimeMetadata<FakeGH, decltype(propagate)> scope(
        target, scoped, propagate);
    assert_metadata(target, 11, 2.5, 0.125, 16);
  }

  assert_metadata(target, 7, 1.25, 0.5, 8);
  assert(propagated.size() == 2);
  assert(propagated.front() == scoped);
  const Internal::LevelTimeMetadata<double> restored{7, 1.25, 0.5, 8};
  assert(propagated.back() == restored);
}

void test_group_domain_router_preserves_legacy_order_and_filters_splits() {
  const auto classify = [](const int group) {
    return group == 0 || group == 2
               ? Internal::TimelevelGroupKind::level_grid_function
               : Internal::TimelevelGroupKind::synchronized_global;
  };
  const auto record = [&](const Internal::TimelevelDomain domain) {
    std::vector<std::string> events;
    Internal::for_each_timelevel_group(
        4, domain, classify,
        [&events](const int group, const Internal::TimelevelGroupKind kind) {
          events.push_back(
              std::string(kind ==
                                  Internal::TimelevelGroupKind::level_grid_function
                              ? "L"
                              : "G") +
              std::to_string(group));
        });
    return events;
  };

  assert(record(Internal::TimelevelDomain::legacy_all) ==
         std::vector<std::string>({"L0", "G1", "L2", "G3"}));
  assert(record(Internal::TimelevelDomain::level_grid_functions) ==
         std::vector<std::string>({"L0", "L2"}));
  assert(record(Internal::TimelevelDomain::synchronized_globals) ==
         std::vector<std::string>({"G1", "G3"}));
}

} // namespace

int main() {
  test_scoped_active_levels_restore_nested_and_exceptional_state();
  test_level_time_metadata_copies_timefac_and_scope_restores_every_field();
  test_group_domain_router_preserves_legacy_order_and_filters_splits();
}
