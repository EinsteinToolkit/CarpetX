#ifndef CARPETX_SUBCYCLING_SCHEDULE_STATE_HXX
#define CARPETX_SUBCYCLING_SCHEDULE_STATE_HXX

#include <optional>
#include <type_traits>
#include <utility>

namespace CarpetX::ScheduleInternal {

template <typename ActiveLevels> class ScopedActiveLevels {
  static_assert(std::is_nothrow_move_constructible_v<ActiveLevels>);
  static_assert(std::is_nothrow_move_assignable_v<ActiveLevels>);

  std::optional<ActiveLevels> &slot_;
  std::optional<ActiveLevels> saved_;

public:
  ScopedActiveLevels(std::optional<ActiveLevels> &slot, ActiveLevels active)
      noexcept
      : slot_(slot), saved_(std::move(slot)) {
    slot_.emplace(std::move(active));
  }

  ScopedActiveLevels(const ScopedActiveLevels &) = delete;
  ScopedActiveLevels &operator=(const ScopedActiveLevels &) = delete;
  ScopedActiveLevels(ScopedActiveLevels &&) = delete;
  ScopedActiveLevels &operator=(ScopedActiveLevels &&) = delete;

  ~ScopedActiveLevels() noexcept { slot_ = std::move(saved_); }
};

template <typename Real> struct LevelTimeMetadata {
  int iteration;
  Real time;
  Real delta_time;
  int timefac;
};

template <typename Real>
constexpr bool operator==(const LevelTimeMetadata<Real> &left,
                          const LevelTimeMetadata<Real> &right) noexcept {
  return left.iteration == right.iteration && left.time == right.time &&
         left.delta_time == right.delta_time && left.timefac == right.timefac;
}

template <typename Real>
constexpr bool operator!=(const LevelTimeMetadata<Real> &left,
                          const LevelTimeMetadata<Real> &right) noexcept {
  return !(left == right);
}

template <typename GH>
constexpr auto capture_level_time_metadata(const GH &gh) noexcept {
  using Real = std::remove_cv_t<
      std::remove_reference_t<decltype(gh.cctk_time)>>;
  return LevelTimeMetadata<Real>{gh.cctk_iteration, gh.cctk_time,
                                 gh.cctk_delta_time, gh.cctk_timefac};
}

template <typename GH, typename Real>
constexpr void apply_level_time_metadata(
    GH &gh, const LevelTimeMetadata<Real> &metadata) noexcept {
  gh.cctk_iteration = metadata.iteration;
  gh.cctk_time = metadata.time;
  gh.cctk_delta_time = metadata.delta_time;
  gh.cctk_timefac = metadata.timefac;
}

template <typename TargetGH, typename SourceGH>
constexpr void copy_level_time_metadata(TargetGH &target,
                                        const SourceGH &source) noexcept {
  apply_level_time_metadata(target, capture_level_time_metadata(source));
}

template <typename GH, typename Propagate> class ScopedLevelTimeMetadata {
  static_assert(std::is_nothrow_invocable_v<Propagate &, GH &>);

  GH &gh_;
  decltype(capture_level_time_metadata(std::declval<const GH &>())) saved_;
  Propagate propagate_;

public:
  using metadata_type = decltype(saved_);

  ScopedLevelTimeMetadata(GH &gh, const metadata_type &metadata,
                          Propagate propagate) noexcept
      : gh_(gh), saved_(capture_level_time_metadata(gh)),
        propagate_(std::move(propagate)) {
    apply_level_time_metadata(gh_, metadata);
    propagate_(gh_);
  }

  ScopedLevelTimeMetadata(const ScopedLevelTimeMetadata &) = delete;
  ScopedLevelTimeMetadata &operator=(const ScopedLevelTimeMetadata &) = delete;
  ScopedLevelTimeMetadata(ScopedLevelTimeMetadata &&) = delete;
  ScopedLevelTimeMetadata &operator=(ScopedLevelTimeMetadata &&) = delete;

  ~ScopedLevelTimeMetadata() noexcept {
    apply_level_time_metadata(gh_, saved_);
    propagate_(gh_);
  }
};

enum class TimelevelGroupKind {
  level_grid_function,
  synchronized_global,
};

enum class TimelevelDomain {
  legacy_all,
  level_grid_functions,
  synchronized_globals,
};

constexpr bool includes_timelevel_group(const TimelevelDomain domain,
                                        const TimelevelGroupKind kind) noexcept {
  switch (domain) {
  case TimelevelDomain::legacy_all:
    return true;
  case TimelevelDomain::level_grid_functions:
    return kind == TimelevelGroupKind::level_grid_function;
  case TimelevelDomain::synchronized_globals:
    return kind == TimelevelGroupKind::synchronized_global;
  }
  return false;
}

template <typename Classify, typename Visit>
void for_each_timelevel_group(const int num_groups,
                              const TimelevelDomain domain,
                              Classify &&classify, Visit &&visit) {
  for (int group = 0; group < num_groups; ++group) {
    const auto kind = classify(group);
    if (includes_timelevel_group(domain, kind))
      visit(group, kind);
  }
}

} // namespace CarpetX::ScheduleInternal

#endif
