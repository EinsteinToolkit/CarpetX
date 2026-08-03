#ifndef ODESOLVERS_EXPLICIT_RK_HXX
#define ODESOLVERS_EXPLICIT_RK_HXX

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace ODESolvers {

enum class ExplicitRKMethod { rk4, rkf78, dp87 };
enum class InitialRHSMode { calculate, reuse_loaded };

template <typename Scalar> struct LoadedRHSProvenance {
  std::uint64_t state_generation;
  Scalar left_time;
};

template <typename Scalar> class LoadedRHSToken {
public:
  LoadedRHSToken(const LoadedRHSToken &) = delete;
  LoadedRHSToken &operator=(const LoadedRHSToken &) = delete;

  LoadedRHSToken(LoadedRHSToken &&other) noexcept
      : provenance_(other.provenance_), consumed_(other.consumed_) {
    other.consumed_ = true;
  }

  LoadedRHSToken &operator=(LoadedRHSToken &&) = delete;

  const LoadedRHSProvenance<Scalar> &provenance() const noexcept {
    return provenance_;
  }
  bool consumed() const noexcept { return consumed_; }
  void consume_once() {
    if (consumed_)
      throw std::invalid_argument("loaded RHS token has already been consumed");
    consumed_ = true;
  }

private:
  explicit LoadedRHSToken(const LoadedRHSProvenance<Scalar> provenance)
      : provenance_(provenance) {}

  LoadedRHSProvenance<Scalar> provenance_;
  bool consumed_ = false;

  template <typename Ops>
  friend LoadedRHSToken<typename Ops::scalar_type>
  make_loaded_rhs_token(const Ops &, typename Ops::scalar_type);
};

template <typename Ops>
LoadedRHSToken<typename Ops::scalar_type>
make_loaded_rhs_token(const Ops &ops,
                      const typename Ops::scalar_type left_time) {
  return LoadedRHSToken<typename Ops::scalar_type>(
      ops.loaded_rhs_provenance(left_time));
}

template <typename Scalar, typename State> struct LinearCombinationView {
  const Scalar *factors;
  const State *const *sources;
  std::size_t size;
};

struct RationalCoefficient {
  std::int64_t numerator;
  std::int64_t denominator;
};

enum class ExplicitRKStageKind : std::uint8_t {
  primary,
  fractional,
  endpoint_probe,
};

struct ExplicitRKAdvanceFrame {
  ExplicitRKStageKind kind;
  RationalCoefficient begin_fraction;
  RationalCoefficient extent_fraction;
};

struct ExplicitRKStagePoint {
  ExplicitRKStageKind kind;
  int stage_index;
  int stage_count;
  RationalCoefficient parent_fraction;
};

struct ExplicitRKTableau {
  ExplicitRKMethod method;
  int endpoint_order;
  std::vector<std::vector<RationalCoefficient>> a;
  std::vector<RationalCoefficient> b;
  // RK4 and RKF78 retain explicit exact stage abscissae. DP87 intentionally
  // leaves this empty and accumulates each converted A row left-to-right.
  std::vector<RationalCoefficient> c;
};

const ExplicitRKTableau &explicit_rk_tableau(ExplicitRKMethod method);
void validate_explicit_rk_tableau(const ExplicitRKTableau &tableau);
RationalCoefficient explicit_rk_stage_fraction(ExplicitRKMethod method,
                                                int stage_index);
ExplicitRKStagePoint
explicit_rk_stage_point(ExplicitRKMethod method,
                        const ExplicitRKAdvanceFrame &frame,
                        int stage_index);

using ExplicitRKTableauFingerprint = std::array<std::uint8_t, 32>;
ExplicitRKTableauFingerprint
explicit_rk_tableau_fingerprint(ExplicitRKMethod method);

struct NullExplicitRKObserver {
  template <typename Scalar, typename State>
  void initial_state(Scalar, const State &) const noexcept {}
  template <typename Scalar, typename State>
  void initial_rhs(Scalar, const State &) const noexcept {}
  template <typename Scalar, typename State>
  void stage_rhs(int, Scalar, const State &) const noexcept {}
  template <typename Scalar, typename State>
  void accepted_endpoint(Scalar, const State &) const noexcept {}
};

namespace detail {

template <typename Ops, typename = void>
struct accepts_stage_point : std::false_type {};

template <typename Ops>
struct accepts_stage_point<
    Ops, std::void_t<decltype(std::declval<Ops &>().set_stage_point(
             std::declval<const ExplicitRKStagePoint &>()))>>
    : std::true_type {};

template <typename Ops>
void publish_stage_point(Ops &ops, const ExplicitRKStagePoint &point) {
  if constexpr (accepts_stage_point<Ops>::value)
    ops.set_stage_point(point);
}

template <typename Scalar>
Scalar convert_coefficient(const RationalCoefficient coefficient) {
  if (coefficient.denominator <= 0)
    throw std::invalid_argument("explicit RK coefficient denominator is not positive");
  const Scalar value = static_cast<Scalar>(coefficient.numerator) /
                       static_cast<Scalar>(coefficient.denominator);
  using std::isfinite;
  if (!isfinite(value))
    throw std::invalid_argument("explicit RK coefficient conversion is not finite");
  return value;
}

template <typename Scalar>
Scalar scale_rational_like_legacy(const Scalar value,
                                  const RationalCoefficient coefficient) {
  if (coefficient.denominator <= 0)
    throw std::invalid_argument("explicit RK coefficient denominator is not positive");
  return value * static_cast<Scalar>(coefficient.numerator) /
         static_cast<Scalar>(coefficient.denominator);
}

template <typename Scalar>
Scalar ratio_of_rationals(const RationalCoefficient numerator,
                          const RationalCoefficient denominator) {
  if (numerator.denominator <= 0 || denominator.denominator <= 0 ||
      denominator.numerator == 0)
    throw std::invalid_argument("invalid explicit RK rational ratio");
  return (static_cast<Scalar>(numerator.numerator) *
          static_cast<Scalar>(denominator.denominator)) /
         (static_cast<Scalar>(numerator.denominator) *
          static_cast<Scalar>(denominator.numerator));
}

template <typename Scalar> struct ConvertedTableau {
  std::vector<std::vector<Scalar>> a;
  std::vector<Scalar> b;
  std::vector<Scalar> c;
};

template <typename Scalar>
ConvertedTableau<Scalar> convert_tableau(const ExplicitRKTableau &tableau) {
  ConvertedTableau<Scalar> result;
  result.a.reserve(tableau.a.size());
  for (const auto &row : tableau.a) {
    std::vector<Scalar> converted;
    converted.reserve(row.size());
    for (const auto coefficient : row)
      converted.push_back(convert_coefficient<Scalar>(coefficient));
    result.a.push_back(std::move(converted));
  }
  result.b.reserve(tableau.b.size());
  for (const auto coefficient : tableau.b)
    result.b.push_back(convert_coefficient<Scalar>(coefficient));
  result.c.reserve(tableau.c.size());
  for (const auto coefficient : tableau.c)
    result.c.push_back(convert_coefficient<Scalar>(coefficient));
  Scalar b_sum = Scalar(0);
  for (const Scalar coefficient : result.b)
    b_sum += coefficient;
  using std::abs;
  if (!std::isfinite(b_sum) ||
      abs(b_sum - Scalar(1)) >
          Scalar(64) * std::numeric_limits<Scalar>::epsilon() *
              Scalar(std::max<std::size_t>(1U, result.b.size())))
    throw std::invalid_argument("explicit RK b coefficients do not sum to one");
  if (tableau.method == ExplicitRKMethod::rkf78) {
    for (std::size_t stage = 0; stage < result.a.size(); ++stage) {
      Scalar row_sum = Scalar(0);
      for (const Scalar coefficient : result.a[stage])
        row_sum += coefficient;
      if (abs(row_sum - result.c.at(stage)) >
          Scalar(64) * std::numeric_limits<Scalar>::epsilon() *
              Scalar(std::max<std::size_t>(1U,
                                           result.a[stage].size())))
        throw std::invalid_argument("RKF78 converted c does not match A row");
    }
  }
  return result;
}

template <typename Scalar>
const ConvertedTableau<Scalar> &
canonical_converted_tableau(const ExplicitRKMethod method) {
  switch (method) {
  case ExplicitRKMethod::rkf78: {
    static const ConvertedTableau<Scalar> converted =
        convert_tableau<Scalar>(explicit_rk_tableau(ExplicitRKMethod::rkf78));
    return converted;
  }
  case ExplicitRKMethod::dp87: {
    static const ConvertedTableau<Scalar> converted =
        convert_tableau<Scalar>(explicit_rk_tableau(ExplicitRKMethod::dp87));
    return converted;
  }
  case ExplicitRKMethod::rk4:
    break;
  }
  throw std::logic_error("RK4 does not require a converted tableau");
}

template <typename Scalar, typename State, std::size_t N>
LinearCombinationView<Scalar, State>
make_linear_combination_view(const std::array<Scalar, N> &factors,
                             const std::array<const State *, N> &sources) {
  return {factors.data(), sources.data(), N};
}

template <typename Scalar>
void validate_interval_and_mode(const Scalar begin_time, const Scalar dt,
                                const InitialRHSMode initial_rhs_mode) {
  using std::isfinite;
  if (!isfinite(begin_time) || !isfinite(dt) || !(dt > Scalar(0)) ||
      !isfinite(begin_time + dt))
    throw std::invalid_argument("explicit RK interval must be finite with positive dt");
  if (initial_rhs_mode != InitialRHSMode::calculate &&
      initial_rhs_mode != InitialRHSMode::reuse_loaded)
    throw std::invalid_argument("invalid explicit RK initial RHS mode");
}

template <typename Ops>
void validate_and_consume_loaded_rhs(
    const typename Ops::scalar_type begin_time,
    Ops &ops,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  if (loaded_rhs_token.consumed())
    throw std::invalid_argument("loaded RHS token has already been consumed");
  const auto &provenance = loaded_rhs_token.provenance();
  using std::isfinite;
  if (!isfinite(provenance.left_time) ||
      provenance.left_time != begin_time ||
      provenance.state_generation != ops.state_generation())
    throw std::invalid_argument(
        "loaded RHS token does not match the left state/time");
  ops.validate_loaded_rhs_provenance(provenance);
  ops.consume_loaded_rhs(provenance);
  loaded_rhs_token.consume_once();
}

template <typename Ops, typename Observer>
void advance_validated_explicit_rk(const ExplicitRKTableau &tableau,
                                   const typename Ops::scalar_type begin_time,
                                   const typename Ops::scalar_type dt,
                                   const bool reuse_loaded_rhs,
                                   const ExplicitRKAdvanceFrame &advance_frame,
                                   Ops &ops,
                                   Observer &observer) {
  using Scalar = typename Ops::scalar_type;
  using State = typename Ops::state_type;

  const auto stage_point = [&](const int stage_index) {
    return explicit_rk_stage_point(tableau.method, advance_frame,
                                   stage_index);
  };
  const int stage_count = static_cast<int>(tableau.a.size());
  publish_stage_point(ops, stage_point(1));
  ops.prepare_initial(begin_time);
  const auto old = ops.snapshot_state();
  observer.initial_state(begin_time, ops.state());

  const auto evaluate_stage = [&](const int stage, const Scalar stage_time) {
    if (stage != 1 || !reuse_loaded_rhs)
      ops.evaluate_rhs(stage);
    ops.validate_rhs(stage);
    if (stage == 1)
      observer.initial_rhs(stage_time, ops.rhs());
    observer.stage_rhs(stage, stage_time, ops.rhs());
  };

  if (tableau.method == ExplicitRKMethod::rk4) {
    evaluate_stage(1, begin_time);
    auto accumulator = ops.snapshot_rhs();

    const Scalar half_dt =
        scale_rational_like_legacy(dt, tableau.a.at(1).at(0));
    const std::array<Scalar, 1> stage2_factors{{half_dt}};
    const std::array<const State *, 1> stage2_sources{{&accumulator}};
    publish_stage_point(ops, stage_point(2));
    ops.update_state(1, begin_time + half_dt, Scalar(1),
                     make_linear_combination_view(stage2_factors,
                                                  stage2_sources));
    evaluate_stage(2, begin_time + half_dt);
    ops.accumulate_rk4(
        accumulator,
        ratio_of_rationals<Scalar>(tableau.b.at(1), tableau.b.at(0)),
        ops.rhs());

    const std::array<Scalar, 2> stage3_factors{{
        Scalar(1),
        scale_rational_like_legacy(dt, tableau.a.at(2).at(1))}};
    const std::array<const State *, 2> stage3_sources{{&old, &ops.rhs()}};
    publish_stage_point(ops, stage_point(3));
    ops.update_state(2, begin_time + half_dt, Scalar(0),
                     make_linear_combination_view(stage3_factors,
                                                  stage3_sources));
    evaluate_stage(3, begin_time + half_dt);
    ops.accumulate_rk4(
        accumulator,
        ratio_of_rationals<Scalar>(tableau.b.at(2), tableau.b.at(0)),
        ops.rhs());

    const std::array<Scalar, 2> stage4_factors{{
        Scalar(1),
        scale_rational_like_legacy(dt, tableau.a.at(3).at(2))}};
    const std::array<const State *, 2> stage4_sources{{&old, &ops.rhs()}};
    publish_stage_point(ops, stage_point(4));
    ops.update_state(3, begin_time + dt, Scalar(0),
                     make_linear_combination_view(stage4_factors,
                                                  stage4_sources));
    evaluate_stage(4, begin_time + dt);
    const std::array<Scalar, 3> endpoint_factors{{
        Scalar(1), scale_rational_like_legacy(dt, tableau.b.at(0)),
        scale_rational_like_legacy(dt, tableau.b.at(3))}};
    const std::array<const State *, 3> endpoint_sources{{
        &old, &accumulator, &ops.rhs()}};
    publish_stage_point(ops, stage_point(stage_count));
    ops.update_state(4, begin_time + dt, Scalar(0),
                     make_linear_combination_view(endpoint_factors,
                                                  endpoint_sources));
  } else {
    const auto &converted = canonical_converted_tableau<Scalar>(tableau.method);
    std::array<std::optional<State>, 13> stages;
    for (std::size_t step = 0; step < converted.a.size(); ++step) {
      Scalar c = Scalar(0);
      if (tableau.method == ExplicitRKMethod::rkf78) {
        c = converted.c.at(step);
      } else {
        // Preserve the DP87 legacy left-to-right execution-scalar sum.
        for (const Scalar coefficient : converted.a.at(step))
          c += coefficient;
      }
      if (step > 0) {
        std::array<Scalar, 14> factors{};
        std::array<const State *, 14> sources{};
        std::size_t count = 1;
        factors[0] = Scalar(1);
        sources[0] = &old;
        for (std::size_t i = 0; i < converted.a.at(step).size(); ++i) {
          const Scalar coefficient = converted.a.at(step).at(i);
          if (coefficient != Scalar(0)) {
            factors[count] = coefficient * dt;
            sources[count] = &stages.at(i).value();
            ++count;
          }
        }
        publish_stage_point(ops, stage_point(static_cast<int>(step + 1)));
        ops.update_state(static_cast<int>(step), begin_time + c * dt,
                         Scalar(0), {factors.data(), sources.data(), count});
      }
      evaluate_stage(static_cast<int>(step + 1), begin_time + c * dt);
      stages[step].emplace(ops.snapshot_rhs());
    }

    std::array<Scalar, 14> factors{};
    std::array<const State *, 14> sources{};
    std::size_t count = 1;
    factors[0] = Scalar(1);
    sources[0] = &old;
    for (std::size_t i = 0; i < converted.b.size(); ++i) {
      const Scalar coefficient = converted.b.at(i);
      if (coefficient != Scalar(0)) {
        factors[count] = coefficient * dt;
        sources[count] = &stages.at(i).value();
        ++count;
      }
    }
    publish_stage_point(ops, stage_point(stage_count));
    ops.update_state(static_cast<int>(converted.a.size()), begin_time + dt,
                     Scalar(0), {factors.data(), sources.data(), count});
  }
  observer.accepted_endpoint(begin_time + dt, ops.state());
}

} // namespace detail

template <typename Ops, typename Observer>
void advance_explicit_rk(const ExplicitRKTableau &tableau,
                         const typename Ops::scalar_type begin_time,
                         const typename Ops::scalar_type dt,
                         const InitialRHSMode initial_rhs_mode, Ops &ops,
                         Observer &observer) {
  validate_explicit_rk_tableau(tableau);
  detail::validate_interval_and_mode(begin_time, dt, initial_rhs_mode);
  if (initial_rhs_mode != InitialRHSMode::calculate)
    throw std::invalid_argument(
        "reuse_loaded requires an explicit loaded RHS token");
  detail::advance_validated_explicit_rk(tableau, begin_time, dt,
                                        false,
                                        {ExplicitRKStageKind::primary,
                                         {0, 1}, {1, 1}},
                                        ops, observer);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(const ExplicitRKTableau &tableau,
                         const typename Ops::scalar_type begin_time,
                         const typename Ops::scalar_type dt,
                         const InitialRHSMode initial_rhs_mode,
                         const ExplicitRKAdvanceFrame &advance_frame, Ops &ops,
                         Observer &observer) {
  validate_explicit_rk_tableau(tableau);
  detail::validate_interval_and_mode(begin_time, dt, initial_rhs_mode);
  if (initial_rhs_mode != InitialRHSMode::calculate)
    throw std::invalid_argument(
        "reuse_loaded requires an explicit loaded RHS token");
  detail::advance_validated_explicit_rk(tableau, begin_time, dt, false,
                                        advance_frame, ops, observer);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(
    const ExplicitRKTableau &tableau,
    const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode, Ops &ops, Observer &observer,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  validate_explicit_rk_tableau(tableau);
  detail::validate_interval_and_mode(begin_time, dt, initial_rhs_mode);
  if (initial_rhs_mode != InitialRHSMode::reuse_loaded)
    throw std::invalid_argument(
        "loaded RHS token is only valid with reuse_loaded");
  detail::validate_and_consume_loaded_rhs(begin_time, ops, loaded_rhs_token);
  detail::advance_validated_explicit_rk(tableau, begin_time, dt,
                                        true,
                                        {ExplicitRKStageKind::primary,
                                         {0, 1}, {1, 1}},
                                        ops, observer);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(
    const ExplicitRKTableau &tableau,
    const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode,
    const ExplicitRKAdvanceFrame &advance_frame, Ops &ops, Observer &observer,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  validate_explicit_rk_tableau(tableau);
  detail::validate_interval_and_mode(begin_time, dt, initial_rhs_mode);
  if (initial_rhs_mode != InitialRHSMode::reuse_loaded)
    throw std::invalid_argument(
        "loaded RHS token is only valid with reuse_loaded");
  detail::validate_and_consume_loaded_rhs(begin_time, ops, loaded_rhs_token);
  detail::advance_validated_explicit_rk(tableau, begin_time, dt, true,
                                        advance_frame, ops, observer);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(const ExplicitRKMethod method,
                         const typename Ops::scalar_type begin_time,
                         const typename Ops::scalar_type dt,
                         const InitialRHSMode initial_rhs_mode, Ops &ops,
                         Observer &observer) {
  const auto &tableau = explicit_rk_tableau(method);
  advance_explicit_rk(tableau, begin_time, dt, initial_rhs_mode, ops, observer);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(
    const ExplicitRKMethod method, const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode, Ops &ops, Observer &observer,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  const auto &tableau = explicit_rk_tableau(method);
  advance_explicit_rk(tableau, begin_time, dt, initial_rhs_mode, ops, observer,
                      loaded_rhs_token);
}

template <typename Ops, typename Observer>
void advance_explicit_rk(
    const ExplicitRKMethod method, const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode,
    const ExplicitRKAdvanceFrame &advance_frame, Ops &ops, Observer &observer,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  const auto &tableau = explicit_rk_tableau(method);
  advance_explicit_rk(tableau, begin_time, dt, initial_rhs_mode,
                      advance_frame, ops, observer, loaded_rhs_token);
}

template <typename Ops>
void advance_explicit_rk(const ExplicitRKTableau &tableau,
                         const typename Ops::scalar_type begin_time,
                         const typename Ops::scalar_type dt,
                         const InitialRHSMode initial_rhs_mode, Ops &ops) {
  NullExplicitRKObserver observer;
  advance_explicit_rk(tableau, begin_time, dt, initial_rhs_mode, ops,
                      observer);
}

template <typename Ops>
void advance_explicit_rk(
    const ExplicitRKTableau &tableau,
    const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode, Ops &ops,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  NullExplicitRKObserver observer;
  advance_explicit_rk(tableau, begin_time, dt, initial_rhs_mode, ops, observer,
                      loaded_rhs_token);
}

template <typename Ops>
void advance_explicit_rk(
    const ExplicitRKMethod method, const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode,
    const ExplicitRKAdvanceFrame &advance_frame, Ops &ops,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  NullExplicitRKObserver observer;
  advance_explicit_rk(method, begin_time, dt, initial_rhs_mode, advance_frame,
                      ops, observer, loaded_rhs_token);
}

template <typename Ops>
void advance_explicit_rk(const ExplicitRKMethod method,
                         const typename Ops::scalar_type begin_time,
                         const typename Ops::scalar_type dt,
                         const InitialRHSMode initial_rhs_mode, Ops &ops) {
  NullExplicitRKObserver observer;
  advance_explicit_rk(method, begin_time, dt, initial_rhs_mode, ops, observer);
}

template <typename Ops>
void advance_explicit_rk(
    const ExplicitRKMethod method, const typename Ops::scalar_type begin_time,
    const typename Ops::scalar_type dt,
    const InitialRHSMode initial_rhs_mode, Ops &ops,
    LoadedRHSToken<typename Ops::scalar_type> &loaded_rhs_token) {
  NullExplicitRKObserver observer;
  advance_explicit_rk(method, begin_time, dt, initial_rhs_mode, ops, observer,
                      loaded_rhs_token);
}

} // namespace ODESolvers

#endif // ODESOLVERS_EXPLICIT_RK_HXX
