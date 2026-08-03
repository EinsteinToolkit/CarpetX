#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// Keep the production StepContext header but avoid pulling a generated cctk.h
// into this bounded standalone test.
#define CARPETX_ARITH_RATIONAL_HXX
namespace Arith {
template <typename I> struct rational {
  I num{0};
  I den{1};
  constexpr rational() = default;
  constexpr rational(I value) : num(value), den(1) {}
  constexpr rational(I numerator, I denominator)
      : num(numerator), den(denominator) {
    if (den < 0) {
      num = -num;
      den = -den;
    }
    const I divisor = std::gcd(num, den);
    num /= divisor;
    den /= divisor;
  }
  friend constexpr rational operator-(const rational &a,
                                      const rational &b) {
    return rational(a.num * b.den - a.den * b.num, a.den * b.den);
  }
  friend constexpr rational operator/(const rational &a, I divisor) {
    return rational(a.num, a.den * divisor);
  }
  friend constexpr bool operator==(const rational &a, const rational &b) {
    return a.num * b.den == a.den * b.num;
  }
  friend constexpr bool operator!=(const rational &a, const rational &b) {
    return !(a == b);
  }
  friend constexpr bool operator<(const rational &a, const rational &b) {
    return a.num * b.den < a.den * b.num;
  }
  friend constexpr bool operator>(const rational &a, const rational &b) {
    return b < a;
  }
  explicit constexpr operator double() const {
    return static_cast<double>(num) / static_cast<double>(den);
  }
};
} // namespace Arith

// Reuse the already validated Phase 8B2 standalone geometry/MultiFab fixture
// verbatim. Its main is renamed; this test adds only executor behavior.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main phase8b2_imported_main
#include "subcycling_scratch_local_gh_unit.cxx"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#define CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_UNIT
#define CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE

constexpr int CCTK_VARIABLE_REAL = 10;
constexpr int CCTK_VALID_GHOSTS = 1;
constexpr int CCTK_VALID_BOUNDARY = 2;
constexpr int CCTK_VALID_INTERIOR = 4;

namespace {

std::array<std::string, 8> variable_names{
    "Other::a", "Test::u", "Test::v", "Other::b",
    "Other::c", "Test::w", "Test::z", "Other::d"};
bool break_var_roundtrip = false;

} // namespace

int CCTK_VarIndex(const char *name) {
  if (name == nullptr)
    return -1;
  for (std::size_t i = 0; i < variable_names.size(); ++i)
    if (variable_names[i] == name)
      return static_cast<int>(i);
  return -1;
}

const char *CCTK_FullVarName(int vi) {
  if (vi < 0 || static_cast<std::size_t>(vi) >= variable_names.size())
    return nullptr;
  if (break_var_roundtrip && vi == 1)
    return "Test::changed";
  return variable_names[static_cast<std::size_t>(vi)].c_str();
}

int CCTK_GroupIndexFromVarI(int vi) {
  if (vi == 1 || vi == 2)
    return 0;
  if (vi == 5 || vi == 6)
    return 2;
  return -1;
}

int CCTK_VarTypeI(int vi) {
  return CCTK_GroupIndexFromVarI(vi) >= 0 ? CCTK_VARIABLE_REAL : -1;
}

#include "../src/subcycling_schedule_certification.hxx"
#include "../src/subcycling_scratch_schedule_executor.hxx"
#include "../src/subcycling_step_context.cxx"

namespace {

using namespace CarpetX;

struct FakeRoutine {
  CanonicalScheduleRecord record;
  void *handle{};
};

struct CallEvent {
  void *handle{};
  std::size_t context{};
  int epoch{};
  int iteration{};
  double time{};
  double base_dt{};
  int timefac{};
};

struct MetadataObservation {
  int iteration{};
  double time{};
  double delta_time{};
  int timefac{};
  const void *scheduled_function{};
};

struct CallControl {
  std::mutex mutex;
  std::condition_variable worker_cv;
  std::vector<cGH *> contexts;
  std::vector<CallEvent> events;
  std::atomic<int> barrier_count{0};
  std::atomic<bool> metadata_ok{true};
  std::atomic<bool> poison{false};
  std::optional<std::size_t> throw_context;
  std::optional<std::size_t> nonzero_context;
  std::optional<int> throw_barrier;
  int expected_iteration{64};
  double expected_time{0.25};
  double expected_base_dt{1.0};
  int expected_timefac{2};
  CertifiedScratchStageExecutor *reenter_executor{};
  const StepContext *reenter_step{};
  const ScratchStageCoordinates *reenter_coordinates{};
  std::optional<ScratchScheduleExecutionErrorCode> reentry_code;
  bool hold_first_worker{};
  bool first_worker_entered{};
  bool release_first_worker{};
  bool record_metadata_snapshots{};
  std::vector<MetadataObservation> metadata_snapshots;
  bool hold_after_idle_publish{};
  bool after_idle_publish_entered{};
  bool release_after_idle_publish{};

  void reset() {
    std::lock_guard<std::mutex> lock(mutex);
    contexts.clear();
    events.clear();
    barrier_count.store(0);
    metadata_ok.store(true);
    poison.store(false);
    throw_context.reset();
    nonzero_context.reset();
    throw_barrier.reset();
    reenter_executor = nullptr;
    reenter_step = nullptr;
    reenter_coordinates = nullptr;
    reentry_code.reset();
    hold_first_worker = false;
    first_worker_entered = false;
    release_first_worker = false;
    record_metadata_snapshots = false;
    metadata_snapshots.clear();
    hold_after_idle_publish = false;
    after_idle_publish_entered = false;
    release_after_idle_publish = false;
  }
};

CallControl calls;
int post_handle_0 = 10;
int post_handle_1 = 11;
int rhs_handle_0 = 20;
int rhs_handle_1 = 21;

std::size_t context_ordinal(cGH *gh) {
  const auto found = std::find(calls.contexts.begin(), calls.contexts.end(), gh);
  if (found == calls.contexts.end())
    throw std::runtime_error("unknown scratch context");
  return static_cast<std::size_t>(found - calls.contexts.begin());
}

int invoke_fake_handle(void *handle, cGH *gh) {
  const auto ordinal = context_ordinal(gh);
  const int epoch = calls.barrier_count.load();
  if (gh->cctk_iteration != calls.expected_iteration ||
      gh->cctk_time != calls.expected_time ||
      gh->cctk_delta_time != calls.expected_base_dt ||
      gh->cctk_timefac != calls.expected_timefac)
    calls.metadata_ok.store(false);
  gh->current_scheduled_function = handle;
  {
    std::unique_lock<std::mutex> lock(calls.mutex);
    calls.events.push_back(CallEvent{handle, ordinal, epoch,
                                     gh->cctk_iteration, gh->cctk_time,
                                     gh->cctk_delta_time, gh->cctk_timefac});
    if (calls.hold_first_worker && ordinal == 0 &&
        handle == &post_handle_0) {
      calls.first_worker_entered = true;
      calls.worker_cv.notify_all();
      calls.worker_cv.wait(lock, [] { return calls.release_first_worker; });
    }
  }
  if (gh->data != nullptr && gh->data[1] != nullptr && gh->data[1][0] != nullptr)
    *static_cast<double *>(gh->data[1][0]) += 1.0;
  if (calls.reenter_executor != nullptr && ordinal == 0 &&
      handle == &post_handle_0) {
    try {
      (void)calls.reenter_executor->execute_rhs(
          *calls.reenter_step, *calls.reenter_coordinates);
    } catch (const ScratchScheduleExecutionError &error) {
      calls.reentry_code = error.code();
      throw;
    }
  }
  if (calls.throw_context == ordinal)
    throw std::runtime_error("fake worker failure");
  gh->current_scheduled_function = nullptr;
  return calls.nonzero_context == ordinal ? 9 : 0;
}

} // namespace

namespace CarpetX {

struct CertifiedLocalScheduleRegistry::Storage {
  CanonicalScheduleBundle manifest;
  ScheduleBuildProvenance provenance;
  std::array<std::vector<FakeRoutine>, 2> functions;
};

CertifiedLocalScheduleRegistry::CertifiedLocalScheduleRegistry(
    std::unique_ptr<Storage> storage)
    : storage_(std::move(storage)) {}
CertifiedLocalScheduleRegistry::~CertifiedLocalScheduleRegistry() = default;
const CanonicalScheduleBundle &
CertifiedLocalScheduleRegistry::manifest() const noexcept {
  return storage_->manifest;
}
const ScheduleBuildProvenance &
CertifiedLocalScheduleRegistry::provenance() const noexcept {
  return storage_->provenance;
}
std::size_t CertifiedLocalScheduleRegistry::size(
    SubcyclingScheduleTarget target) const noexcept {
  const std::size_t index =
      target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  return storage_->functions[index].size();
}
const CanonicalScheduleRecord &CertifiedLocalScheduleRegistry::record_for_executor(
    SubcyclingScheduleTarget target, std::size_t ordinal) const {
  const std::size_t index =
      target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  return storage_->functions.at(index).at(ordinal).record;
}
int CertifiedLocalScheduleRegistry::invoke_for_executor(
    SubcyclingScheduleTarget target, std::size_t ordinal, cGH *scratch_gh) const {
  const std::size_t index =
      target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  const auto &routine = storage_->functions.at(index).at(ordinal);
  return invoke_fake_handle(routine.handle, scratch_gh);
}

class CertifiedLocalScheduleRegistryExecutorTestAccess {
public:
  static std::unique_ptr<CertifiedLocalScheduleRegistry> make(
      std::vector<FakeRoutine> post_step, std::vector<FakeRoutine> rhs) {
    auto storage = std::make_unique<CertifiedLocalScheduleRegistry::Storage>();
    storage->functions[0] = std::move(rhs);
    storage->functions[1] = std::move(post_step);
    return std::unique_ptr<CertifiedLocalScheduleRegistry>(
        new CertifiedLocalScheduleRegistry(std::move(storage)));
  }
};

bool subcycling_scratch_executor_poison_enabled_for_test() noexcept {
  return calls.poison.load();
}

void subcycling_scratch_executor_synchronize_for_test() {
  const int ordinal = calls.barrier_count.fetch_add(1);
  if (calls.throw_barrier == ordinal)
    throw std::runtime_error("fake device barrier failure");
}

void subcycling_scratch_executor_metadata_snapshotted_for_test() {
  std::lock_guard<std::mutex> lock(calls.mutex);
  if (!calls.record_metadata_snapshots || calls.contexts.empty())
    return;
  const cGH *const gh = calls.contexts.front();
  calls.metadata_snapshots.push_back(MetadataObservation{
      gh->cctk_iteration, gh->cctk_time, gh->cctk_delta_time,
      gh->cctk_timefac, gh->current_scheduled_function});
}

void subcycling_scratch_executor_after_idle_publish_for_test() {
  std::unique_lock<std::mutex> lock(calls.mutex);
  if (!calls.hold_after_idle_publish || calls.after_idle_publish_entered)
    return;
  calls.after_idle_publish_entered = true;
  calls.worker_cv.notify_all();
  calls.worker_cv.wait(lock,
                       [] { return calls.release_after_idle_publish; });
}

} // namespace CarpetX

#include "../src/subcycling_scratch_schedule_executor.cxx"

namespace {

using namespace CarpetX;

CanonicalAccess access(std::string variable, int read, int write,
                       int invalidate) {
  return CanonicalAccess{std::move(variable), 0, read, write, invalidate};
}

FakeRoutine routine(std::uint64_t ordinal, const char *name, void *handle,
                    std::vector<CanonicalAccess> accesses) {
  CanonicalScheduleRecord record;
  record.traversal_ordinal = ordinal;
  record.item_ordinal = ordinal + 100;
  record.item.routine = name;
  record.item.accesses = std::move(accesses);
  return FakeRoutine{std::move(record), handle};
}

std::unique_ptr<CertifiedLocalScheduleRegistry>
make_registry(const bool invalidate_next_read = false) {
  std::vector<FakeRoutine> post;
  post.push_back(routine(
      0, "Post0", &post_handle_0,
      {access("Test::u", CCTK_VALID_INTERIOR, 0, 0),
       access("Test::v", 0, CCTK_VALID_INTERIOR,
              invalidate_next_read ? CCTK_VALID_INTERIOR
                                   : CCTK_VALID_GHOSTS)}));
  post.push_back(routine(
      1, "Post1", &post_handle_1,
      {access("Test::v", CCTK_VALID_INTERIOR, 0, 0),
       access("Test::u", 0, CCTK_VALID_GHOSTS, CCTK_VALID_BOUNDARY)}));
  std::vector<FakeRoutine> rhs;
  rhs.push_back(routine(
      0, "Rhs0", &rhs_handle_0,
      {access("Test::u", CCTK_VALID_GHOSTS, 0, 0),
       access("Test::w", 0, CCTK_VALID_INTERIOR, CCTK_VALID_BOUNDARY)}));
  rhs.push_back(routine(
      1, "Rhs1", &rhs_handle_1,
      {access("Test::w", CCTK_VALID_INTERIOR, 0, 0),
       access("Test::z", 0,
              CCTK_VALID_INTERIOR | CCTK_VALID_GHOSTS,
              CCTK_VALID_GHOSTS)}));
  return CertifiedLocalScheduleRegistryExecutorTestAccess::make(
      std::move(post), std::move(rhs));
}

class NoopPreparer final : public StagePreparer {
public:
  void prepare_stage(const StepContext &, const StagePoint &) override {}
};

struct Fixture {
  GHExt ghext;
  std::unique_ptr<ScratchLocalLevelBinding> binding;
  std::unique_ptr<CertifiedLocalScheduleRegistry> registry;
  std::unique_ptr<CertifiedScratchStageExecutor> executor;
  StepContext step{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                   SubcyclingODEMethod::rk4};
  ScratchStageCoordinates coordinates{64, 0.25, step_clock_t(1), 1.0, 2};
  NoopPreparer preparer;

  explicit Fixture(bool zero_contexts = false,
                   bool invalidate_next_read = false)
      : ghext(make_ghext()) {
    if (zero_contexts) {
      auto &level = ghext.patchdata[0].leveldata[0];
      level.fab->tiles.clear();
      level.local_cctkGHs.clear();
      level.geometry.clear();
    }
    auto certified = make_certified(ghext);
    binding = ScratchLocalLevelBinding::bind(ghext, std::move(certified));
    registry = make_registry(invalidate_next_read);
    calls.reset();
    for (std::size_t i = 0; i < binding->local_context_count(); ++i)
      calls.contexts.push_back(const_cast<cGH *>(
          &ScratchLocalLevelBindingTestAccess::gh(*binding, i)));
    executor = CertifiedScratchStageExecutor::create(*registry, *binding);
  }

  ScratchScheduleExecutionReceipt post() {
    ScopedStepContext scope(step, preparer);
    return executor->execute_post_step(step, coordinates);
  }
  ScratchScheduleExecutionReceipt rhs() {
    ScopedStepContext scope(step, preparer);
    return executor->execute_rhs(step, coordinates);
  }
};

template <typename F>
ScratchScheduleExecutionErrorCode expect_execution_error(F &&call) {
  try {
    call();
  } catch (const ScratchScheduleExecutionError &error) {
    return error.code();
  }
  assert(false && "expected ScratchScheduleExecutionError");
  return ScratchScheduleExecutionErrorCode::internal_failure;
}

std::vector<void *> handles() {
  std::lock_guard<std::mutex> lock(calls.mutex);
  std::vector<void *> result;
  for (const auto &event : calls.events)
    result.push_back(event.handle);
  return result;
}

std::vector<int> epochs() {
  std::lock_guard<std::mutex> lock(calls.mutex);
  std::vector<int> result;
  for (const auto &event : calls.events)
    result.push_back(event.epoch);
  return result;
}

void test_factory_resolves_both_phases_atomically() {
  Fixture fixture;
  assert(fixture.executor != nullptr && !fixture.executor->faulted());
  assert(fixture.registry->size(SubcyclingScheduleTarget::post_step) == 2);
  assert(fixture.registry->size(SubcyclingScheduleTarget::rhs) == 2);
}

template <typename T, typename = void> struct has_public_context : std::false_type {};
template <typename T>
struct has_public_context<
    T, std::void_t<decltype(std::declval<T &>().context_for_executor(0))>>
    : std::true_type {};

void test_handles_and_raw_cgh_never_escape_public_api() {
  static_assert(!has_public_context<ScratchLocalLevelBinding>::value);
  static_assert(!std::is_copy_constructible_v<CertifiedScratchStageExecutor>);
  static_assert(!std::is_move_constructible_v<CertifiedScratchStageExecutor>);
  Fixture fixture;
  (void)fixture.post();
  const auto observed = handles();
  assert(std::count(observed.begin(), observed.end(), &post_handle_0) == 3);
  assert(std::count(observed.begin(), observed.end(), &post_handle_1) == 3);
  assert(std::count(observed.begin(), observed.end(), &rhs_handle_0) == 0);
}

void test_stage_coordinates_install_on_every_scratch_context() {
  Fixture fixture;
  const auto receipt = fixture.post();
  assert(receipt.local_call_count == 6 && calls.metadata_ok.load());
}

void test_level_dt_identity_and_step_context_are_enforced() {
  {
    Fixture fixture;
    fixture.coordinates.base_delta_clock = step_clock_t(3, 2);
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_post_step(fixture.step,
                                                       fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_coordinates);
  }
  {
    Fixture fixture;
    StepContext other = fixture.step;
    ScopedStepContext scope(other, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_post_step(fixture.step,
                                                       fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_step_context);
  }
}

void test_poison_undefined_values_is_rejected_on_every_execute() {
  Fixture fixture;
  (void)fixture.post();
  calls.poison.store(true);
  ScopedStepContext scope(fixture.step, fixture.preparer);
  assert(expect_execution_error([&] {
           (void)fixture.executor->execute_rhs(fixture.step,
                                               fixture.coordinates);
         }) == ScratchScheduleExecutionErrorCode::unsupported_configuration);
}

void test_stage_metadata_restores_exactly_on_every_exit() {
  Fixture fixture;
  std::vector<std::array<std::uintptr_t, 5>> before;
  for (cGH *gh : calls.contexts) {
    gh->cctk_iteration = 777;
    gh->cctk_time = 8.5;
    gh->cctk_delta_time = 9.5;
    gh->cctk_timefac = 13;
    gh->current_scheduled_function = reinterpret_cast<void *>(0x1234);
    before.push_back({static_cast<std::uintptr_t>(gh->cctk_iteration),
                      static_cast<std::uintptr_t>(gh->cctk_time * 2),
                      static_cast<std::uintptr_t>(gh->cctk_delta_time * 2),
                      static_cast<std::uintptr_t>(gh->cctk_timefac),
                      reinterpret_cast<std::uintptr_t>(
                          gh->current_scheduled_function)});
  }
  (void)fixture.post();
  for (std::size_t i = 0; i < calls.contexts.size(); ++i) {
    const cGH *gh = calls.contexts[i];
    assert(gh->cctk_iteration == static_cast<int>(before[i][0]));
    assert(gh->cctk_time == static_cast<double>(before[i][1]) / 2.0);
    assert(gh->cctk_delta_time == static_cast<double>(before[i][2]) / 2.0);
    assert(gh->cctk_timefac == static_cast<int>(before[i][3]));
    assert(reinterpret_cast<std::uintptr_t>(gh->current_scheduled_function) ==
           before[i][4]);
  }
}

void test_binding_lease_rejects_second_executor_and_reentry() {
  Fixture fixture;
  assert(expect_execution_error([&] {
           (void)CertifiedScratchStageExecutor::create(*fixture.registry,
                                                       *fixture.binding);
         }) == ScratchScheduleExecutionErrorCode::binding_busy);
  calls.reenter_executor = fixture.executor.get();
  calls.reenter_step = &fixture.step;
  calls.reenter_coordinates = &fixture.coordinates;
  ScopedStepContext scope(fixture.step, fixture.preparer);
  assert(expect_execution_error([&] {
           (void)fixture.executor->execute_post_step(fixture.step,
                                                     fixture.coordinates);
         }) == ScratchScheduleExecutionErrorCode::call_failed);
  assert(calls.reentry_code ==
         ScratchScheduleExecutionErrorCode::reentrant_execution);
}

void test_two_thread_reentry_faults_active_call_and_binding() {
  Fixture fixture;
  {
    std::lock_guard<std::mutex> lock(calls.mutex);
    calls.hold_first_worker = true;
  }

  bool first_returned_receipt = false;
  std::optional<ScratchScheduleExecutionErrorCode> first_error;
  std::optional<ScratchScheduleExecutionErrorCode> second_error;
  std::thread first([&] {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    try {
      (void)fixture.executor->execute_post_step(fixture.step,
                                                fixture.coordinates);
      first_returned_receipt = true;
    } catch (const ScratchScheduleExecutionError &error) {
      first_error = error.code();
    }
  });

  {
    std::unique_lock<std::mutex> lock(calls.mutex);
    calls.worker_cv.wait(lock, [] { return calls.first_worker_entered; });
  }
  std::thread second([&] {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    try {
      (void)fixture.executor->execute_rhs(fixture.step, fixture.coordinates);
    } catch (const ScratchScheduleExecutionError &error) {
      second_error = error.code();
    }
  });
  second.join();
  {
    std::lock_guard<std::mutex> lock(calls.mutex);
    calls.release_first_worker = true;
  }
  calls.worker_cv.notify_all();
  first.join();

  assert(second_error == ScratchScheduleExecutionErrorCode::reentrant_execution);
  assert(first_error == ScratchScheduleExecutionErrorCode::reentrant_execution);
  assert(!first_returned_receipt && fixture.executor->faulted());
  fixture.executor.reset();
  assert(expect_execution_error([&] {
           (void)CertifiedScratchStageExecutor::create(*fixture.registry,
                                                       *fixture.binding);
         }) == ScratchScheduleExecutionErrorCode::binding_faulted);
}

void test_idle_is_published_only_after_metadata_restore() {
  Fixture fixture;
  for (cGH *gh : calls.contexts) {
    gh->cctk_iteration = 777;
    gh->cctk_time = 8.5;
    gh->cctk_delta_time = 9.5;
    gh->cctk_timefac = 13;
    gh->current_scheduled_function = reinterpret_cast<void *>(0x1234);
  }
  {
    std::lock_guard<std::mutex> lock(calls.mutex);
    calls.record_metadata_snapshots = true;
    calls.hold_after_idle_publish = true;
  }

  bool first_receipt = false;
  bool second_receipt = false;
  std::thread first([&] {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    (void)fixture.executor->execute_post_step(fixture.step,
                                              fixture.coordinates);
    first_receipt = true;
  });
  {
    std::unique_lock<std::mutex> lock(calls.mutex);
    calls.worker_cv.wait(lock,
                         [] { return calls.after_idle_publish_entered; });
  }
  std::thread second([&] {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    (void)fixture.executor->execute_rhs(fixture.step, fixture.coordinates);
    second_receipt = true;
  });
  second.join();

  std::vector<MetadataObservation> observations;
  {
    std::lock_guard<std::mutex> lock(calls.mutex);
    observations = calls.metadata_snapshots;
    calls.release_after_idle_publish = true;
  }
  calls.worker_cv.notify_all();
  first.join();

  assert(first_receipt && second_receipt && observations.size() == 2);
  for (const auto &observation : observations) {
    assert(observation.iteration == 777);
    assert(observation.time == 8.5);
    assert(observation.delta_time == 9.5);
    assert(observation.timefac == 13);
    assert(observation.scheduled_function ==
           reinterpret_cast<void *>(0x1234));
  }
}

void test_poststep_and_rhs_are_separate_certified_phases() {
  Fixture fixture;
  const auto receipt = fixture.rhs();
  assert(receipt.phase == ScratchSchedulePhase::rhs &&
         receipt.routine_count == 2);
  const auto observed = handles();
  assert(std::count(observed.begin(), observed.end(), &rhs_handle_0) == 3);
  assert(std::count(observed.begin(), observed.end(), &rhs_handle_1) == 3);
  assert(std::count(observed.begin(), observed.end(), &post_handle_0) == 0);
}

void test_routine_major_calls_have_one_barrier_per_routine() {
  Fixture fixture;
  (void)fixture.post();
  assert(calls.barrier_count.load() == 2);
  const auto observed_epochs = epochs();
  assert(std::count(observed_epochs.begin(), observed_epochs.end(), 0) == 3);
  assert(std::count(observed_epochs.begin(), observed_epochs.end(), 1) == 3);
}

void test_reads_gate_calls_and_writes_then_invalidates_update_validity() {
  {
    Fixture fixture;
    auto &frame = ScratchLocalLevelBindingTestAccess::frame(*fixture.binding);
    frame.mutable_validity(0)[0].interior = false;
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_post_step(fixture.step,
                                                       fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_read);
    assert(handles().empty() && calls.barrier_count.load() == 0);
  }
  {
    Fixture fixture;
    auto &frame = ScratchLocalLevelBindingTestAccess::frame(*fixture.binding);
    (void)fixture.rhs();
    assert(frame.validity(1)[0].interior);
    assert(!frame.validity(1)[0].outer);
    assert(frame.validity(1)[1].interior);
    assert(!frame.validity(1)[1].ghosts);
  }
}

void test_later_invalid_read_faults_binding_after_first_routine() {
  Fixture fixture(false, true);
  ScopedStepContext scope(fixture.step, fixture.preparer);
  assert(expect_execution_error([&] {
           (void)fixture.executor->execute_post_step(fixture.step,
                                                     fixture.coordinates);
         }) == ScratchScheduleExecutionErrorCode::invalid_read);
  assert(calls.barrier_count.load() == 1);
  const auto observed = handles();
  assert(std::count(observed.begin(), observed.end(), &post_handle_0) == 3);
  fixture.executor.reset();
  assert(expect_execution_error([&] {
           (void)CertifiedScratchStageExecutor::create(*fixture.registry,
                                                       *fixture.binding);
         }) == ScratchScheduleExecutionErrorCode::binding_faulted);
}

void test_missing_or_changed_variable_mapping_rejects_factory() {
  auto ghext = make_ghext();
  auto certified = make_certified(ghext);
  auto binding = ScratchLocalLevelBinding::bind(ghext, std::move(certified));
  auto registry = make_registry();
  break_var_roundtrip = true;
  assert(expect_execution_error([&] {
           (void)CertifiedScratchStageExecutor::create(*registry, *binding);
         }) == ScratchScheduleExecutionErrorCode::unresolved_access);
  break_var_roundtrip = false;
  auto executor = CertifiedScratchStageExecutor::create(*registry, *binding);
  assert(executor != nullptr);
}

void test_inactive_wrong_level_stale_epoch_and_invalid_coordinates_fail_before_calls() {
  {
    Fixture fixture;
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(fixture.step,
                                                 fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_step_context);
    assert(handles().empty());
  }
  {
    Fixture fixture;
    StepContext wrong = fixture.step;
    wrong.level = 1;
    ScopedStepContext scope(wrong, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(wrong, fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_step_context);
    assert(handles().empty());
  }
  {
    Fixture fixture;
    ++current_epoch;
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(fixture.step,
                                                 fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::stale_binding);
    --current_epoch;
    assert(handles().empty());
  }
  {
    Fixture fixture;
    fixture.coordinates.stage_time =
        std::numeric_limits<double>::quiet_NaN();
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(fixture.step,
                                                 fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::invalid_coordinates);
    assert(handles().empty());
  }
}

void test_nonzero_call_result_drains_barrier_and_faults_executor() {
  Fixture fixture;
  calls.nonzero_context = 1;
  ScopedStepContext scope(fixture.step, fixture.preparer);
  assert(expect_execution_error([&] {
           (void)fixture.executor->execute_post_step(fixture.step,
                                                     fixture.coordinates);
         }) == ScratchScheduleExecutionErrorCode::call_failed);
  assert(calls.barrier_count.load() == 1 && fixture.executor->faulted());
}

void test_worker_or_barrier_failure_restores_metadata_and_faults_binding() {
  Fixture fixture;
  calls.throw_context = 2;
  ScopedStepContext scope(fixture.step, fixture.preparer);
  assert(expect_execution_error([&] {
           (void)fixture.executor->execute_post_step(fixture.step,
                                                     fixture.coordinates);
         }) == ScratchScheduleExecutionErrorCode::call_failed);
  assert(calls.barrier_count.load() == 1);
  fixture.executor.reset();
  assert(expect_execution_error([&] {
           (void)CertifiedScratchStageExecutor::create(*fixture.registry,
                                                       *fixture.binding);
         }) == ScratchScheduleExecutionErrorCode::binding_faulted);
}

void test_faulted_executor_cannot_be_reused() {
  Fixture fixture;
  calls.throw_barrier = 0;
  {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(fixture.step,
                                                 fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::call_failed);
  }
  {
    ScopedStepContext scope(fixture.step, fixture.preparer);
    assert(expect_execution_error([&] {
             (void)fixture.executor->execute_rhs(fixture.step,
                                                 fixture.coordinates);
           }) == ScratchScheduleExecutionErrorCode::already_faulted);
  }
}

void test_zero_local_context_rank_is_valid() {
  Fixture fixture(true);
  const auto receipt = fixture.post();
  assert(receipt.routine_count == 2 && receipt.local_call_count == 0);
  assert(calls.barrier_count.load() == 2 && handles().empty());
}

void test_scratch_calls_do_not_mutate_live_tl0_or_live_validity() {
  Fixture fixture;
  auto &group = *fixture.ghext.patchdata[0].leveldata[0].groupdata[0];
  double *live = group.mfab[0]->array(2).ptr(15, -1, -1, 0);
  const double before = *live;
  const auto live_validity = group.valid[0];
  (void)fixture.post();
  assert(*live == before && group.valid[0] == live_validity);
}

} // namespace

int main() {
  test_factory_resolves_both_phases_atomically();
  test_handles_and_raw_cgh_never_escape_public_api();
  test_stage_coordinates_install_on_every_scratch_context();
  test_level_dt_identity_and_step_context_are_enforced();
  test_poison_undefined_values_is_rejected_on_every_execute();
  test_stage_metadata_restores_exactly_on_every_exit();
  test_binding_lease_rejects_second_executor_and_reentry();
  test_two_thread_reentry_faults_active_call_and_binding();
  test_idle_is_published_only_after_metadata_restore();
  test_poststep_and_rhs_are_separate_certified_phases();
  test_routine_major_calls_have_one_barrier_per_routine();
  test_reads_gate_calls_and_writes_then_invalidates_update_validity();
  test_later_invalid_read_faults_binding_after_first_routine();
  test_missing_or_changed_variable_mapping_rejects_factory();
  test_inactive_wrong_level_stale_epoch_and_invalid_coordinates_fail_before_calls();
  test_nonzero_call_result_drains_barrier_and_faults_executor();
  test_worker_or_barrier_failure_restores_metadata_and_faults_binding();
  test_faulted_executor_cannot_be_reused();
  test_zero_local_context_rank_is_valid();
  test_scratch_calls_do_not_mutate_live_tl0_or_live_validity();
  std::cout << "Phase 8B3C certified scratch executor tests passed\n";
}
