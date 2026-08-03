#ifndef CARPETX_TRANSACTION_LEVEL_STEP_SESSION_HXX
#define CARPETX_TRANSACTION_LEVEL_STEP_SESSION_HXX

#include "hierarchy_stepper.hxx"
#include "subcycling_scratch_state_transaction.hxx"

#include <functional>
#include <memory>

namespace CarpetX {

using TransactionLevelEvolution =
    std::function<void(ScratchStateTransaction &)>;

// This callback is optional. A driver adapter must finish validation before
// mutating accepted LevelData metadata and must not throw after mutation.
using ValidatedLevelCommit = std::function<void()>;

class TransactionLevelStepSession final : public LevelStepSession {
public:
  TransactionLevelStepSession(
      std::unique_ptr<ScratchStateTransaction> transaction,
      bool require_dense, TransactionLevelEvolution evolution,
      ValidatedLevelCommit validated_commit = {});
  ~TransactionLevelStepSession() noexcept override;

  TransactionLevelStepSession(const TransactionLevelStepSession &) = delete;
  TransactionLevelStepSession &
  operator=(const TransactionLevelStepSession &) = delete;
  TransactionLevelStepSession(TransactionLevelStepSession &&) = delete;
  TransactionLevelStepSession &operator=(TransactionLevelStepSession &&) =
      delete;

  ScratchStateTransaction *transaction() noexcept override;
  LevelAdvanceResult advance() override;
  void commit() override;

private:
  enum class Lifecycle {
    ready,
    advancing,
    advanced,
    committing,
    committed,
    rolled_back,
  };

  void rollback_or_terminate() noexcept;

  std::unique_ptr<ScratchStateTransaction> transaction_;
  bool require_dense_{false};
  TransactionLevelEvolution evolution_;
  ValidatedLevelCommit validated_commit_;
  ScratchStateToken primary_left_;
  Lifecycle lifecycle_{Lifecycle::ready};
};

} // namespace CarpetX

#endif // CARPETX_TRANSACTION_LEVEL_STEP_SESSION_HXX
