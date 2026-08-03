#include "transaction_level_step_session.hxx"

#include <exception>
#include <stdexcept>
#include <utility>

namespace CarpetX {

TransactionLevelStepSession::TransactionLevelStepSession(
    std::unique_ptr<ScratchStateTransaction> transaction,
    const bool require_dense, TransactionLevelEvolution evolution,
    ValidatedLevelCommit validated_commit)
    : transaction_(std::move(transaction)), require_dense_(require_dense),
      evolution_(std::move(evolution)),
      validated_commit_(std::move(validated_commit)) {
  if (transaction_ == nullptr)
    throw std::invalid_argument(
        "transaction level session requires a transaction");
  if (!evolution_)
    throw std::invalid_argument(
        "transaction level session requires an evolution callback");
  transaction_->arm_live_evolved_rollback();
}

TransactionLevelStepSession::~TransactionLevelStepSession() noexcept {
  if (lifecycle_ != Lifecycle::committed &&
      lifecycle_ != Lifecycle::rolled_back)
    rollback_or_terminate();
}

ScratchStateTransaction *
TransactionLevelStepSession::transaction() noexcept {
  switch (lifecycle_) {
  case Lifecycle::ready:
  case Lifecycle::advancing:
  case Lifecycle::advanced:
  case Lifecycle::committing:
    return transaction_.get();
  case Lifecycle::committed:
  case Lifecycle::rolled_back:
    return nullptr;
  }
  return nullptr;
}

LevelAdvanceResult TransactionLevelStepSession::advance() {
  if (lifecycle_ != Lifecycle::ready)
    throw std::logic_error(
        "transaction level session advance is not available");
  lifecycle_ = Lifecycle::advancing;
  try {
    evolution_(*transaction_);
    if (transaction_->faulted() || transaction_->discarded())
      throw std::logic_error(
          "level evolution left its transaction unavailable");

    auto dense_interval = transaction_->take_committed_dense();
    if (transaction_->faulted() || transaction_->discarded())
      throw std::logic_error(
          "dense publication left its transaction unavailable");
    if (require_dense_ && dense_interval == nullptr)
      throw std::logic_error(
          "transaction level session requires a dense interval");
    if (!require_dense_ && dense_interval != nullptr)
      throw std::logic_error(
          "transaction leaf level published an unexpected dense interval");

    lifecycle_ = Lifecycle::advanced;
    return LevelAdvanceResult{std::move(dense_interval)};
  } catch (...) {
    rollback_or_terminate();
    throw;
  }
}

void TransactionLevelStepSession::commit() {
  if (lifecycle_ != Lifecycle::advanced)
    throw std::logic_error(
        "transaction level session has no accepted endpoint to commit");
  lifecycle_ = Lifecycle::committing;
  try {
    if (validated_commit_)
      validated_commit_();
    if (transaction_->faulted() || transaction_->discarded())
      throw std::logic_error(
          "accepted metadata callback left its transaction unavailable");
    transaction_->disarm_live_evolved_rollback();
    transaction_->discard();
    lifecycle_ = Lifecycle::committed;
    transaction_.reset();
  } catch (...) {
    rollback_or_terminate();
    throw;
  }
}

void TransactionLevelStepSession::rollback_or_terminate() noexcept {
  if (lifecycle_ == Lifecycle::committed ||
      lifecycle_ == Lifecycle::rolled_back)
    return;
  if (transaction_ == nullptr)
    std::terminate();
  try {
    transaction_->rollback_live_evolved();
    lifecycle_ = Lifecycle::rolled_back;
  } catch (...) {
    std::terminate();
  }
}

} // namespace CarpetX
