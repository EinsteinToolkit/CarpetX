#ifndef ODESOLVERS_SUBCYCLING_TRANSACTION_BRIDGE_HXX
#define ODESOLVERS_SUBCYCLING_TRANSACTION_BRIDGE_HXX

#include <subcycling_scratch_state_transaction.hxx>

#include <initializer_list>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ODESolvers {

// Private ODESolvers ownership wrapper around CarpetX's opaque token.  This
// type deliberately knows no cGH, MultiFab, GroupData, or scratch-frame type.
class TransactionStateBackend final {
public:
  struct LinearTerm {
    double coefficient;
    const TransactionStateBackend *state;
  };

  static TransactionStateBackend
  from_token(CarpetX::ScratchStateTransaction &transaction,
             CarpetX::ScratchStateToken token) {
    if (!token.valid())
      throw std::invalid_argument("scratch state token is empty");
    // Validate ownership/schema/epoch now rather than delaying failure until
    // the first arithmetic operation.
    static_cast<void>(transaction.state_kind(token));
    return TransactionStateBackend(transaction, std::move(token));
  }

  TransactionStateBackend(const TransactionStateBackend &) = delete;
  TransactionStateBackend &operator=(const TransactionStateBackend &) =
      delete;
  TransactionStateBackend(TransactionStateBackend &&) noexcept = default;
  TransactionStateBackend &
  operator=(TransactionStateBackend &&) noexcept = default;

  CarpetX::ScratchStateTransaction &transaction() const noexcept {
    return *transaction_;
  }
  CarpetX::ScratchStateToken &token() noexcept { return token_; }
  const CarpetX::ScratchStateToken &token() const noexcept { return token_; }

  CarpetX::ScratchStateKind kind() const {
    return transaction().state_kind(token_);
  }

  bool valid(const CarpetX::ScratchStateRegion region) const {
    return transaction().state_valid(token_, region);
  }

  void set_valid(const CarpetX::ScratchStateRegion region,
                 const bool value) {
    transaction().set_state_valid(token_, region, value);
  }

  TransactionStateBackend clone() const {
    return from_token(transaction(), transaction().clone_state(token_));
  }

  void linear_combination(
      const double destination_scale,
      const std::initializer_list<LinearTerm> terms) {
    linear_combination(destination_scale,
                       std::vector<LinearTerm>(terms.begin(), terms.end()));
  }

  void linear_combination(const double destination_scale,
                          const std::vector<LinearTerm> &terms) {
    std::vector<CarpetX::ScratchLinearTerm> opaque_terms;
    opaque_terms.reserve(terms.size());
    for (const auto &term : terms) {
      if (term.state == nullptr)
        throw std::invalid_argument("scratch linear term is null");
      if (&term.state->transaction() != transaction_)
        throw std::invalid_argument(
            "scratch linear combination mixes transaction owners");
      opaque_terms.push_back({term.coefficient, &term.state->token()});
    }
    transaction().linear_combination(token_, destination_scale, opaque_terms);
  }

private:
  TransactionStateBackend(CarpetX::ScratchStateTransaction &transaction,
                          CarpetX::ScratchStateToken token) noexcept
      : transaction_(&transaction), token_(std::move(token)) {}

  CarpetX::ScratchStateTransaction *transaction_;
  CarpetX::ScratchStateToken token_;
};

struct TransactionPrimaryCapture {
  TransactionStateBackend left_state;
  TransactionStateBackend left_rhs;
  TransactionStateBackend accepted_endpoint;
};

// Observer for the unchanged primary RK advance.  The observer ignores the
// live statecomp_t reference supplied by the kernel and asks the transaction
// to snapshot the certified live identities at the exact observer events.
class TransactionPrimaryObserver final {
public:
  explicit TransactionPrimaryObserver(
      CarpetX::ScratchStateTransaction &transaction) noexcept
      : transaction_(&transaction) {}

  template <typename Scalar, typename State>
  void initial_state(const Scalar, const State &) {
    if (left_state_ || taken_)
      throw std::logic_error("primary left state was captured twice");
    left_state_.emplace(TransactionStateBackend::from_token(
        *transaction_, transaction_->capture_live_evolved()));
  }

  template <typename Scalar, typename State>
  void initial_rhs(const Scalar, const State &) {
    if (!left_state_ || left_rhs_ || taken_)
      throw std::logic_error("primary left RHS capture is out of order");
    left_rhs_.emplace(TransactionStateBackend::from_token(
        *transaction_, transaction_->capture_live_rhs()));
  }

  template <typename Scalar, typename State>
  void stage_rhs(const int, const Scalar, const State &) const noexcept {}

  template <typename Scalar, typename State>
  void accepted_endpoint(const Scalar, const State &) {
    if (!left_state_ || !left_rhs_ || accepted_endpoint_ || taken_)
      throw std::logic_error("primary endpoint capture is out of order");
    accepted_endpoint_.emplace(TransactionStateBackend::from_token(
        *transaction_, transaction_->capture_live_evolved()));
  }

  TransactionPrimaryCapture take_complete() {
    if (taken_ || !left_state_ || !left_rhs_ || !accepted_endpoint_)
      throw std::logic_error("primary transaction capture is incomplete");
    taken_ = true;
    return {std::move(*left_state_), std::move(*left_rhs_),
            std::move(*accepted_endpoint_)};
  }

private:
  CarpetX::ScratchStateTransaction *transaction_;
  std::optional<TransactionStateBackend> left_state_;
  std::optional<TransactionStateBackend> left_rhs_;
  std::optional<TransactionStateBackend> accepted_endpoint_;
  bool taken_{false};
};

} // namespace ODESolvers

#endif // ODESOLVERS_SUBCYCLING_TRANSACTION_BRIDGE_HXX
