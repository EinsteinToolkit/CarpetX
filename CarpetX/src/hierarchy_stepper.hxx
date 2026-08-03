#ifndef CARPETX_HIERARCHY_STEPPER_HXX
#define CARPETX_HIERARCHY_STEPPER_HXX

#include "subcycling_dense_output.hxx"

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

namespace CarpetX {

using hierarchy_time_t = step_clock_t;

struct LevelClock {
  hierarchy_time_t time{0};
  hierarchy_time_t dt{1};
  std::uint64_t accepted_steps{0};
};

struct LevelAdvanceConfig {
  SubcyclingODEMethod method;
  TableauFingerprint tableau_fingerprint;
};

struct HierarchyStepperConfig {
  std::vector<LevelAdvanceConfig> levels;
  hierarchy_time_t initial_clock{0};
  double initial_physical_time{0.0};
  double coarse_dt{1.0};
  int refinement_ratio{2};
  std::uint64_t initial_epoch{0};
  std::vector<std::uint64_t> initial_accepted_steps;
};

struct LevelAdvanceResult {
  std::shared_ptr<const DenseInterval> dense_interval;
};

struct HierarchyAdvanceResult {
  hierarchy_time_t synchronized_time{0};
  double synchronized_physical_time{0.0};
  std::uint64_t epoch{0};
  bool stop_requested{false};
};

class LevelStepSession {
public:
  virtual ~LevelStepSession() = default;
  virtual ScratchStateTransaction *transaction() noexcept = 0;
  virtual LevelAdvanceResult advance() = 0;
  virtual void commit() = 0;
};

class HierarchyEvolutionAdapter {
public:
  virtual ~HierarchyEvolutionAdapter() = default;

  virtual std::unique_ptr<LevelStepSession>
  begin_level_step(const StepContext &context, bool require_dense) = 0;
  virtual void prepare_stage(const StepContext &context,
                             const StagePoint &stage_point,
                             const DenseInterval *parent_dense) = 0;

  virtual void begin_synchronization(int coarse_level, int fine_level,
                                     hierarchy_time_t time) = 0;
  virtual void synchronize_levels(int coarse_level, int fine_level,
                                  hierarchy_time_t time) = 0;
  virtual void end_synchronization(int coarse_level, int fine_level,
                                   hierarchy_time_t time) noexcept = 0;

  virtual void run_sync_observers(hierarchy_time_t time, double physical_time,
                                  std::uint64_t completed_epoch,
                                  bool stop_requested) = 0;
};

class HierarchyStepper {
public:
  HierarchyStepper(HierarchyStepperConfig config,
                   const DenseOutputRegistry &dense_registry);

  HierarchyAdvanceResult advance_one_epoch(HierarchyEvolutionAdapter &adapter);
  void request_stop_at_next_sync() noexcept;
  bool stop_pending() const noexcept;
  bool hierarchy_synchronized() const noexcept;
  std::uint64_t epoch() const noexcept;
  const std::vector<LevelClock> &clocks() const noexcept;

private:
  void advance_level(int level, HierarchyEvolutionAdapter &adapter);
  bool clocks_aligned() const noexcept;
  bool no_dense_in_flight() const noexcept;
  void preflight_epoch_capacity() const;
  double physical_time(hierarchy_time_t clock) const;
  void validate_dense_interval(int level, const StepContext &context,
                               const DenseInterval &interval) const;

  std::vector<LevelAdvanceConfig> level_configs_;
  std::vector<LevelClock> clocks_;
  std::vector<DenseCapability> expected_dense_capabilities_;
  std::vector<std::shared_ptr<const DenseInterval>> active_dense_;
  hierarchy_time_t initial_clock_{0};
  double initial_physical_time_{0.0};
  double coarse_dt_{1.0};
  bool in_progress_{false};
  bool faulted_{false};
  std::atomic_bool stop_requested_{false};
  std::atomic_bool stop_snapshot_in_flight_{false};
  std::uint64_t epoch_{0};
};

} // namespace CarpetX

#endif // CARPETX_HIERARCHY_STEPPER_HXX
