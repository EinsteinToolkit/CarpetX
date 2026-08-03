#include "subcycling_scratch_schedule_executor.hxx"

#include "subcycling_schedule_certification.hxx"
#include "subcycling_scratch_local_gh.hxx"

#ifndef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
#include <cctk.h>
#include <cctk_Groups.h>
#include <cctk_Parameters.h>
extern "C" CCTK_INT CarpetX_GetEpoch(void);
#endif

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace CarpetX {

#ifndef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
void synchronize();
#endif

#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
bool subcycling_scratch_executor_poison_enabled_for_test() noexcept;
void subcycling_scratch_executor_synchronize_for_test();
void subcycling_scratch_executor_metadata_snapshotted_for_test();
void subcycling_scratch_executor_after_idle_publish_for_test();
#endif

struct ResolvedScratchScheduleAccess {
  std::size_t frame_entry{};
  std::size_t component{};
  int read_mask{};
  int write_mask{};
  int invalidate_mask{};
};

struct ResolvedScratchScheduleRoutine {
  std::vector<ResolvedScratchScheduleAccess> accesses;
};

namespace {

constexpr int supported_validity_mask =
    CCTK_VALID_GHOSTS | CCTK_VALID_BOUNDARY | CCTK_VALID_INTERIOR;

struct MetadataSnapshot {
  cGH *gh{};
  int iteration{};
  double time{};
  double delta_time{};
  int timefac{};
  const void *scheduled_function{};
};

[[nodiscard]] std::size_t phase_index(const ScratchSchedulePhase phase) {
  return phase == ScratchSchedulePhase::post_step ? 0U : 1U;
}

[[nodiscard]] SubcyclingScheduleTarget
target_for_phase(const ScratchSchedulePhase phase) {
  return phase == ScratchSchedulePhase::post_step
             ? SubcyclingScheduleTarget::post_step
             : SubcyclingScheduleTarget::rhs;
}

[[nodiscard]] bool power_of_two(const int value) noexcept {
  return value > 0 && (value & (value - 1)) == 0;
}

[[nodiscard]] double time_tolerance(const double a, const double b,
                                    const double c = 0.0) noexcept {
  const double scale =
      std::max({1.0, std::abs(a), std::abs(b), std::abs(c)});
  return 16.0 * std::numeric_limits<double>::epsilon() * scale;
}

[[nodiscard]] bool poison_enabled() noexcept {
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
  return subcycling_scratch_executor_poison_enabled_for_test();
#else
  DECLARE_CCTK_PARAMETERS;
  return poison_undefined_values;
#endif
}

void drain_device_work() {
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
  subcycling_scratch_executor_synchronize_for_test();
#else
  synchronize();
#endif
}

[[nodiscard]] bool validity_bit(const ScratchValidity &validity,
                                const int bit) noexcept {
  if (bit == CCTK_VALID_INTERIOR)
    return validity.interior;
  if (bit == CCTK_VALID_BOUNDARY)
    return validity.outer;
  return validity.ghosts;
}

void set_validity_bit(ScratchValidity &validity, const int bit,
                      const bool value) noexcept {
  if (bit == CCTK_VALID_INTERIOR)
    validity.interior = value;
  else if (bit == CCTK_VALID_BOUNDARY)
    validity.outer = value;
  else
    validity.ghosts = value;
}

template <typename F> void for_each_validity_bit(const int mask, F &&function) {
  for (const int bit :
       {CCTK_VALID_GHOSTS, CCTK_VALID_BOUNDARY, CCTK_VALID_INTERIOR})
    if ((mask & bit) != 0)
      function(bit);
}

} // namespace

struct CertifiedScratchStageExecutor::Storage {
  enum class Lifecycle : std::uint8_t { idle, preflight, started, faulted };
  std::array<std::vector<ResolvedScratchScheduleRoutine>, 2> routines;
  std::atomic<Lifecycle> lifecycle{Lifecycle::idle};
};

#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
void ScratchLocalLevelBinding::install_stage_coordinates_for_executor(
    const ScratchStageCoordinates &coordinates) noexcept {
  for (auto &context : contexts_) {
    context->gh.cctk_iteration = coordinates.level_iteration;
    context->gh.cctk_time = coordinates.stage_time;
    context->gh.cctk_delta_time = coordinates.base_delta_time;
    context->gh.cctk_timefac = coordinates.time_refinement_factor;
  }
}
#endif

ScratchScheduleExecutionError::ScratchScheduleExecutionError(
    const ScratchScheduleExecutionErrorCode code,
    const std::optional<ScratchSchedulePhase> phase,
    const std::optional<std::size_t> routine_ordinal,
    const std::optional<std::size_t> context_ordinal,
    const std::string &message)
    : std::runtime_error(message), code_(code), phase_(phase),
      routine_ordinal_(routine_ordinal), context_ordinal_(context_ordinal) {}

ScratchScheduleExecutionErrorCode
ScratchScheduleExecutionError::code() const noexcept {
  return code_;
}

std::optional<ScratchSchedulePhase>
ScratchScheduleExecutionError::phase() const noexcept {
  return phase_;
}

std::optional<std::size_t>
ScratchScheduleExecutionError::routine_ordinal() const noexcept {
  return routine_ordinal_;
}

std::optional<std::size_t>
ScratchScheduleExecutionError::context_ordinal() const noexcept {
  return context_ordinal_;
}

[[noreturn]] void CertifiedScratchStageExecutor::raise(
    const ScratchScheduleExecutionErrorCode code,
    const std::optional<ScratchSchedulePhase> phase,
    const std::optional<std::size_t> routine_ordinal,
    const std::optional<std::size_t> context_ordinal,
    const std::string &message) {
  throw ScratchScheduleExecutionError(code, phase, routine_ordinal,
                                      context_ordinal, message);
}

CertifiedScratchStageExecutor::CertifiedScratchStageExecutor(
    CertifiedLocalScheduleRegistry &registry,
    ScratchLocalLevelBinding &binding, std::unique_ptr<Storage> storage)
    : registry_(registry), binding_(binding), storage_(std::move(storage)) {}

CertifiedScratchStageExecutor::~CertifiedScratchStageExecutor() {
  binding_.release_for_executor();
}

std::unique_ptr<CertifiedScratchStageExecutor>
CertifiedScratchStageExecutor::create(CertifiedLocalScheduleRegistry &registry,
                                      ScratchLocalLevelBinding &binding) {
  if (binding.executor_faulted_.load())
    raise(ScratchScheduleExecutionErrorCode::binding_faulted, std::nullopt,
          std::nullopt, std::nullopt, "scratch binding is faulted");
  if (binding.executor_busy_.load())
    raise(ScratchScheduleExecutionErrorCode::binding_busy, std::nullopt,
          std::nullopt, std::nullopt,
          "scratch binding already has an executor");
  try {
    binding.claim_for_executor();
  } catch (const std::logic_error &) {
    const auto code = binding.executor_faulted_.load()
                          ? ScratchScheduleExecutionErrorCode::binding_faulted
                          : ScratchScheduleExecutionErrorCode::binding_busy;
    raise(code, std::nullopt, std::nullopt, std::nullopt,
          "scratch binding lease could not be acquired");
  }

  bool published = false;
  try {
    auto storage = std::make_unique<Storage>();
    auto &frame = binding.frame_for_executor();
    constexpr std::array<ScratchSchedulePhase, 2> phases{
        ScratchSchedulePhase::post_step, ScratchSchedulePhase::rhs};
    for (const auto phase : phases) {
      const auto target = target_for_phase(phase);
      auto &resolved_phase = storage->routines[phase_index(phase)];
      resolved_phase.reserve(registry.size(target));
      for (std::size_t ordinal = 0; ordinal < registry.size(target);
           ++ordinal) {
        const auto &record = registry.record_for_executor(target, ordinal);
        ResolvedScratchScheduleRoutine resolved;
        resolved.accesses.reserve(record.item.accesses.size());
        for (const auto &canonical : record.item.accesses) {
          if (canonical.timelevel != 0 ||
              ((canonical.read_mask | canonical.write_mask |
                canonical.invalidate_mask) &
               ~supported_validity_mask) != 0)
            raise(ScratchScheduleExecutionErrorCode::unresolved_access, phase,
                  ordinal, std::nullopt,
                  "certified access has an unsupported timelevel or mask");
          const int variable = CCTK_VarIndex(canonical.variable_name.c_str());
          const char *const roundtrip =
              variable >= 0 ? CCTK_FullVarName(variable) : nullptr;
          if (variable < 0 || roundtrip == nullptr ||
              canonical.variable_name != roundtrip ||
              CCTK_VarTypeI(variable) != CCTK_VARIABLE_REAL)
            raise(ScratchScheduleExecutionErrorCode::unresolved_access, phase,
                  ordinal, std::nullopt,
                  "certified variable mapping did not round-trip");
          const int group = CCTK_GroupIndexFromVarI(variable);
          const int first = group >= 0 ? CCTK_FirstVarIndexI(group) : -1;
          cGroup group_data{};
          if (group < 0 || first < 0 ||
              CCTK_GroupData(group, &group_data) != 0 ||
              group_data.grouptype != CCTK_GF || group_data.numvars <= 0 ||
              variable < first || variable >= first + group_data.numvars)
            raise(ScratchScheduleExecutionErrorCode::unresolved_access, phase,
                  ordinal, std::nullopt,
                  "certified variable group mapping is inconsistent");

          std::optional<std::size_t> frame_entry;
          for (std::size_t entry = 0; entry < frame.entry_count(); ++entry) {
            if (frame.key(entry).group_index == group) {
              if (frame_entry.has_value())
                raise(ScratchScheduleExecutionErrorCode::unresolved_access,
                      phase, ordinal, std::nullopt,
                      "scratch frame contains a duplicate group mapping");
              frame_entry = entry;
            }
          }
          const auto component = static_cast<std::size_t>(variable - first);
          if (!frame_entry.has_value() ||
              component >= frame.validity(*frame_entry).size())
            raise(ScratchScheduleExecutionErrorCode::unresolved_access, phase,
                  ordinal, std::nullopt,
                  "scratch frame does not contain the certified variable");
          resolved.accesses.push_back(
              ResolvedScratchScheduleAccess{
                  *frame_entry, component, canonical.read_mask,
                  canonical.write_mask, canonical.invalidate_mask});
        }
        resolved_phase.push_back(std::move(resolved));
      }
    }
    auto result = std::unique_ptr<CertifiedScratchStageExecutor>(
        new CertifiedScratchStageExecutor(registry, binding,
                                           std::move(storage)));
    published = true;
    return result;
  } catch (const ScratchScheduleExecutionError &) {
    if (!published)
      binding.release_for_executor();
    throw;
  } catch (const std::exception &error) {
    if (!published)
      binding.release_for_executor();
    raise(ScratchScheduleExecutionErrorCode::internal_failure, std::nullopt,
          std::nullopt, std::nullopt,
          std::string("executor factory failed: ") + error.what());
  } catch (...) {
    if (!published)
      binding.release_for_executor();
    raise(ScratchScheduleExecutionErrorCode::internal_failure, std::nullopt,
          std::nullopt, std::nullopt,
          "executor factory failed with a non-standard exception");
  }
}

ScratchScheduleExecutionReceipt CertifiedScratchStageExecutor::execute_post_step(
    const StepContext &step, const ScratchStageCoordinates &coordinates) {
  return execute(ScratchSchedulePhase::post_step, step, coordinates);
}

ScratchScheduleExecutionReceipt CertifiedScratchStageExecutor::execute_rhs(
    const StepContext &step, const ScratchStageCoordinates &coordinates) {
  return execute(ScratchSchedulePhase::rhs, step, coordinates);
}

bool CertifiedScratchStageExecutor::faulted() const noexcept {
  return storage_->lifecycle.load() == Storage::Lifecycle::faulted;
}

ScratchScheduleExecutionReceipt CertifiedScratchStageExecutor::execute(
    const ScratchSchedulePhase phase, const StepContext &step,
    const ScratchStageCoordinates &coordinates) {
  using Lifecycle = Storage::Lifecycle;
  for (;;) {
    Lifecycle observed = storage_->lifecycle.load();
    if (observed == Lifecycle::idle) {
      if (storage_->lifecycle.compare_exchange_weak(observed,
                                                    Lifecycle::preflight))
        break;
      continue;
    }
    if (observed == Lifecycle::faulted)
      raise(ScratchScheduleExecutionErrorCode::already_faulted, phase,
            std::nullopt, std::nullopt,
            "scratch executor is already faulted");
    const Lifecycle active_state = observed;
    if (storage_->lifecycle.compare_exchange_strong(observed,
                                                    Lifecycle::faulted)) {
      if (active_state == Lifecycle::started)
        binding_.fault_for_executor();
      raise(ScratchScheduleExecutionErrorCode::reentrant_execution, phase,
            std::nullopt, std::nullopt,
            "scratch executor execution is already active");
    }
  }

  bool execution_started = false;
  try {

  const auto preflight_failure =
      [&](const ScratchScheduleExecutionErrorCode code,
          const char *const message) -> void {
    raise(code, phase, std::nullopt, std::nullopt, message);
  };

  if (current_step_context() != &step || step.level < 0 ||
      binding_.level_for_executor() != step.level ||
      !std::isfinite(step.begin_time) || !std::isfinite(step.end_time) ||
      !(step.begin_time < step.end_time) ||
      !(step.begin_clock < step.end_clock))
    preflight_failure(ScratchScheduleExecutionErrorCode::invalid_step_context,
                      "active StepContext or binding level is invalid");
  if (coordinates.level_iteration < 0 ||
      !std::isfinite(coordinates.stage_time) ||
      !std::isfinite(coordinates.base_delta_time) ||
      !(coordinates.base_delta_time > 0.0) ||
      !power_of_two(coordinates.time_refinement_factor) ||
      !(step_clock_t(0) < coordinates.base_delta_clock))
    preflight_failure(ScratchScheduleExecutionErrorCode::invalid_coordinates,
                      "scratch stage coordinates are invalid");

  const double stage_tolerance =
      time_tolerance(step.begin_time, step.end_time, coordinates.stage_time);
  if (coordinates.stage_time < step.begin_time - stage_tolerance ||
      coordinates.stage_time > step.end_time + stage_tolerance)
    preflight_failure(ScratchScheduleExecutionErrorCode::invalid_coordinates,
                      "scratch stage time is outside the active step");
  const auto level_delta_clock = step.end_clock - step.begin_clock;
  if (level_delta_clock != coordinates.base_delta_clock /
                               coordinates.time_refinement_factor)
    preflight_failure(ScratchScheduleExecutionErrorCode::invalid_coordinates,
                      "base and level clocks do not match exactly");
  const double expected_level_dt =
      coordinates.base_delta_time /
      static_cast<double>(coordinates.time_refinement_factor);
  const double observed_level_dt = step.end_time - step.begin_time;
  if (std::abs(observed_level_dt - expected_level_dt) >
      time_tolerance(observed_level_dt, expected_level_dt,
                     coordinates.base_delta_time))
    preflight_failure(ScratchScheduleExecutionErrorCode::invalid_coordinates,
                      "base and level double time steps do not match");
  if (binding_.hierarchy_epoch() !=
      static_cast<std::int64_t>(CarpetX_GetEpoch()))
    preflight_failure(ScratchScheduleExecutionErrorCode::stale_binding,
                      "scratch binding hierarchy epoch is stale");
  if (poison_enabled())
    preflight_failure(
        ScratchScheduleExecutionErrorCode::unsupported_configuration,
        "poison_undefined_values must be disabled for scratch execution");

  std::vector<MetadataSnapshot> metadata;
  try {
    metadata.reserve(binding_.local_context_count());
    for (std::size_t ordinal = 0;
         ordinal < binding_.local_context_count(); ++ordinal) {
      cGH *const gh = binding_.context_for_executor(ordinal);
      if (gh == nullptr)
        preflight_failure(ScratchScheduleExecutionErrorCode::internal_failure,
                          "scratch binding returned a null context");
      metadata.push_back(MetadataSnapshot{
          gh, gh->cctk_iteration, gh->cctk_time, gh->cctk_delta_time,
          gh->cctk_timefac,
          static_cast<const void *>(gh->current_scheduled_function)});
    }
  } catch (const ScratchScheduleExecutionError &) {
    throw;
  } catch (const std::exception &error) {
    raise(ScratchScheduleExecutionErrorCode::internal_failure, phase,
          std::nullopt, std::nullopt,
          std::string("could not snapshot scratch metadata: ") +
              error.what());
  }
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
  subcycling_scratch_executor_metadata_snapshotted_for_test();
#endif

  std::size_t routine_count = 0;
  std::size_t local_calls = 0;
  {
  const std::size_t context_count = binding_.local_context_count();
  std::vector<int> results;
  std::vector<std::exception_ptr> failures;
  try {
    results.resize(context_count);
    failures.resize(context_count);
  } catch (const std::exception &error) {
    raise(ScratchScheduleExecutionErrorCode::internal_failure, phase,
          std::nullopt, std::nullopt,
          std::string("could not allocate scratch execution state: ") +
              error.what());
  }

  struct MetadataGuard {
    ScratchLocalLevelBinding &binding;
    std::vector<MetadataSnapshot> &snapshots;
    ~MetadataGuard() {
      for (const auto &snapshot : snapshots) {
        snapshot.gh->cctk_iteration = snapshot.iteration;
        snapshot.gh->cctk_time = snapshot.time;
        snapshot.gh->cctk_delta_time = snapshot.delta_time;
        snapshot.gh->cctk_timefac = snapshot.timefac;
        snapshot.gh->current_scheduled_function =
            reinterpret_cast<decltype(snapshot.gh->current_scheduled_function)>(
                const_cast<void *>(snapshot.scheduled_function));
      }
    }
  } metadata_guard{binding_, metadata};
  binding_.install_stage_coordinates_for_executor(coordinates);

  auto &frame = binding_.frame_for_executor();
  const auto &routines = storage_->routines[phase_index(phase)];
  const auto target = target_for_phase(phase);
  for (std::size_t routine_ordinal = 0;
       routine_ordinal < routines.size(); ++routine_ordinal) {
    if (execution_started &&
        storage_->lifecycle.load() != Lifecycle::started)
      raise(ScratchScheduleExecutionErrorCode::reentrant_execution, phase,
            routine_ordinal, std::nullopt,
            "scratch executor lifecycle changed during execution");
    const auto &routine = routines[routine_ordinal];
    for (const auto &resolved : routine.accesses) {
      const auto &validity =
          frame.validity(resolved.frame_entry).at(resolved.component);
      bool reads_valid = true;
      for_each_validity_bit(resolved.read_mask, [&](const int bit) {
        reads_valid = reads_valid && validity_bit(validity, bit);
      });
      if (!reads_valid) {
        raise(ScratchScheduleExecutionErrorCode::invalid_read, phase,
              routine_ordinal, std::nullopt,
              "scratch routine read requires invalid data");
      }
    }

    if (!execution_started) {
      Lifecycle expected = Lifecycle::preflight;
      if (!storage_->lifecycle.compare_exchange_strong(expected,
                                                       Lifecycle::started))
        raise(ScratchScheduleExecutionErrorCode::reentrant_execution, phase,
              routine_ordinal, std::nullopt,
              "scratch executor lifecycle changed before first call");
      execution_started = true;
    }

    std::fill(results.begin(), results.end(), 0);
    std::fill(failures.begin(), failures.end(), nullptr);
#pragma omp parallel for schedule(static)
    for (std::int64_t context_ordinal = 0;
         context_ordinal < static_cast<std::int64_t>(context_count);
         ++context_ordinal) {
      const auto index = static_cast<std::size_t>(context_ordinal);
      try {
        cGH *const gh = binding_.context_for_executor(index);
        if (gh == nullptr)
          throw std::runtime_error("scratch context disappeared");
        results[index] =
            registry_.invoke_for_executor(target, routine_ordinal, gh);
      } catch (...) {
        failures[index] = std::current_exception();
      }
    }

    std::exception_ptr barrier_failure;
    try {
      drain_device_work();
    } catch (...) {
      barrier_failure = std::current_exception();
    }
    local_calls += context_count;

    std::optional<std::size_t> failed_context;
    for (std::size_t index = 0; index < context_count; ++index)
      if (failures[index] != nullptr || results[index] != 0) {
        failed_context = index;
        break;
    }
    if (failed_context.has_value() || barrier_failure != nullptr) {
      raise(ScratchScheduleExecutionErrorCode::call_failed, phase,
            routine_ordinal, failed_context,
            barrier_failure != nullptr
                ? "scratch routine device barrier failed"
                : "scratch routine worker or call bridge failed");
    }
    if (storage_->lifecycle.load() != Lifecycle::started)
      raise(ScratchScheduleExecutionErrorCode::reentrant_execution, phase,
            routine_ordinal, std::nullopt,
            "scratch executor lifecycle changed during routine execution");

    for (const auto &resolved : routine.accesses) {
      auto &validity =
          frame.mutable_validity(resolved.frame_entry).at(resolved.component);
      for_each_validity_bit(resolved.write_mask, [&](const int bit) {
        set_validity_bit(validity, bit, true);
      });
    }
    for (const auto &resolved : routine.accesses) {
      auto &validity =
          frame.mutable_validity(resolved.frame_entry).at(resolved.component);
      for_each_validity_bit(resolved.invalidate_mask, [&](const int bit) {
        set_validity_bit(validity, bit, false);
      });
    }
  }
  routine_count = routines.size();
  }

  Lifecycle expected =
      execution_started ? Lifecycle::started : Lifecycle::preflight;
  if (!storage_->lifecycle.compare_exchange_strong(expected, Lifecycle::idle))
    raise(ScratchScheduleExecutionErrorCode::reentrant_execution, phase,
          std::nullopt, std::nullopt,
          "scratch executor lifecycle changed before receipt publication");
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE
  subcycling_scratch_executor_after_idle_publish_for_test();
#endif
  return ScratchScheduleExecutionReceipt{phase, routine_count, local_calls};
  } catch (const ScratchScheduleExecutionError &) {
    const Lifecycle previous =
        storage_->lifecycle.exchange(Lifecycle::faulted);
    if (execution_started || previous == Lifecycle::started)
      binding_.fault_for_executor();
    throw;
  } catch (const std::exception &error) {
    const Lifecycle previous =
        storage_->lifecycle.exchange(Lifecycle::faulted);
    if (execution_started || previous == Lifecycle::started)
      binding_.fault_for_executor();
    raise(ScratchScheduleExecutionErrorCode::internal_failure, phase,
          std::nullopt, std::nullopt,
          std::string("unexpected scratch executor failure: ") + error.what());
  } catch (...) {
    const Lifecycle previous =
        storage_->lifecycle.exchange(Lifecycle::faulted);
    if (execution_started || previous == Lifecycle::started)
      binding_.fault_for_executor();
    raise(ScratchScheduleExecutionErrorCode::internal_failure, phase,
          std::nullopt, std::nullopt,
          "unknown scratch executor failure");
  }
}

} // namespace CarpetX
