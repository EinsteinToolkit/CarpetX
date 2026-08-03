#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_HXX
#define CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_HXX

#include "subcycling_step_context.hxx"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace CarpetX {

class DenseInterval;
struct DenseIntervalId;
class DenseOutputProvider;
class ScratchStateTransactionFactory;
class ScratchStateTransactionCore;

enum class ScratchStateKind : std::uint8_t { evolved, raw_rhs };
enum class ScratchStateRegion : std::uint8_t { interior, outer, ghosts };

// A token names transaction-owned storage. It intentionally exposes neither
// the storage representation nor any driver object.
class ScratchStateToken final {
public:
  ScratchStateToken() noexcept;
  ~ScratchStateToken() = default;
  ScratchStateToken(const ScratchStateToken &) = delete;
  ScratchStateToken &operator=(const ScratchStateToken &) = delete;
  ScratchStateToken(ScratchStateToken &&other) noexcept;
  ScratchStateToken &operator=(ScratchStateToken &&other) noexcept;

  bool valid() const noexcept;

private:
  ScratchStateToken(std::uint64_t owner, std::uint64_t state,
                    std::int64_t epoch, std::uint64_t schema,
                    ScratchStateKind kind) noexcept;
  void reset() noexcept;

  std::uint64_t owner_{0};
  std::uint64_t state_{0};
  std::int64_t epoch_{-1};
  std::uint64_t schema_{0};
  ScratchStateKind kind_{ScratchStateKind::evolved};

  friend class ScratchStateTransaction;
  friend class ScratchStateTransactionFactory;
  friend class ScratchStateTransactionCore;
};

struct ScratchGroupPair {
  int evolved_group;
  int rhs_group;
};

struct ScratchLinearTerm {
  double coefficient;
  const ScratchStateToken *state;
};

enum class ScratchDenseSampleKind : std::uint8_t {
  value,
  raw_derivative,
};

struct ScratchDenseSampleRef {
  double theta;
  ScratchDenseSampleKind kind;
  const ScratchStateToken *state;
};

class ScratchStateTransaction final {
public:
  ~ScratchStateTransaction();
  ScratchStateTransaction(const ScratchStateTransaction &) = delete;
  ScratchStateTransaction &operator=(const ScratchStateTransaction &) =
      delete;
  ScratchStateTransaction(ScratchStateTransaction &&) = delete;
  ScratchStateTransaction &operator=(ScratchStateTransaction &&) = delete;

  const std::vector<ScratchGroupPair> &group_pairs() const noexcept;
  const std::vector<int> &dependent_groups() const noexcept;

  ScratchStateToken capture_live_evolved();
  ScratchStateToken capture_live_rhs();
  ScratchStateToken capture_scratch_evolved();
  ScratchStateToken capture_scratch_rhs();
  ScratchStateToken clone_state(const ScratchStateToken &state);

  ScratchStateKind state_kind(const ScratchStateToken &state) const;
  bool state_valid(const ScratchStateToken &state,
                   ScratchStateRegion region) const;
  void set_state_valid(ScratchStateToken &state, ScratchStateRegion region,
                       bool valid);
  // Terminal failure recovery for a captured primary-step left state. This
  // restores live TL0 in place and never executes PostStep or RHS schedules.
  void rollback_live_evolved(const ScratchStateToken &state);
  void restore_state(const ScratchStateToken &state);
  void restore_left(const ScratchStateToken &evolved,
                    const ScratchStateToken &raw_rhs);
  void linear_combination(ScratchStateToken &destination,
                          double destination_scale,
                          const std::vector<ScratchLinearTerm> &terms);

  void post_step_after_update(const StepContext &context,
                              const StagePoint &stage_point);
  void evaluate_rhs(const StepContext &context,
                    const StagePoint &stage_point);
  std::size_t rhs_evaluation_count() const noexcept;

  void commit_dense(const StepContext &context,
                    const DenseOutputProvider &provider,
                    const DenseIntervalId &interval,
                    const std::vector<ScratchDenseSampleRef> &samples);
  std::shared_ptr<const DenseInterval> take_committed_dense() noexcept;

  bool faulted() const noexcept;
  bool discarded() const noexcept;
  void discard() noexcept;

private:
  struct Storage;
  explicit ScratchStateTransaction(std::unique_ptr<Storage> storage);
  std::unique_ptr<Storage> storage_;

  friend class ScratchStateTransactionFactory;
  friend class ScratchStateTransactionCore;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_HXX
