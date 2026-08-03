#include "subcycling_stage_spatial_preparer.hxx"
#include "subcycling_stage_spatial_preparer_internal.hxx"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace CarpetX {
namespace {

double fraction_as_double(const step_clock_t value) {
  return static_cast<double>(value);
}

bool close_stage_time(const double actual, const double expected) {
  const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
  return std::abs(actual - expected) <=
         32.0 * std::numeric_limits<double>::epsilon() * scale;
}

StageSpatialTarget target_for(const StageKind kind) {
  switch (kind) {
  case StageKind::primary:
    return StageSpatialTarget::primary_live_tl0;
  case StageKind::fractional:
  case StageKind::endpoint_probe:
    return StageSpatialTarget::transaction_scratch;
  }
  throw std::invalid_argument("unknown subcycling stage kind");
}

std::string target_name(const StageSpatialTarget target) {
  switch (target) {
  case StageSpatialTarget::primary_live_tl0:
    return "live";
  case StageSpatialTarget::transaction_scratch:
    return "scratch";
  }
  throw std::logic_error("unknown stage spatial target");
}

void validate_stage(const StepContext &context,
                    const StagePoint &stage_point,
                    const step_clock_t stage_clock) {
  if (context.level < 0 || !(context.begin_clock < context.end_clock))
    throw std::invalid_argument("stage spatial context is invalid");
  if (!std::isfinite(context.begin_time) ||
      !std::isfinite(context.end_time) ||
      !(context.begin_time < context.end_time))
    throw std::invalid_argument(
        "stage spatial physical interval is invalid");
  if (stage_point.stage_index <= 0 || stage_point.stage_count <= 0 ||
      stage_point.stage_index > stage_point.stage_count)
    throw std::invalid_argument("stage spatial index/count is invalid");
  if (stage_point.kind == StageKind::endpoint_probe &&
      (stage_point.stage_index != 1 || stage_point.stage_count != 1))
    throw std::invalid_argument(
        "endpoint probe must be a singleton stage");
  if (stage_point.stage_fraction.den <= 0 ||
      stage_point.stage_fraction < 0 || stage_point.stage_fraction > 1)
    throw std::invalid_argument(
        "stage spatial fraction must be exact and in [0,1]");
  if (stage_clock < context.begin_clock || stage_clock > context.end_clock)
    throw std::invalid_argument("stage spatial clock escaped its step");
  if (!std::isfinite(stage_point.stage_time))
    throw std::invalid_argument("stage spatial time must be finite");
  const double expected_time =
      context.begin_time + fraction_as_double(stage_point.stage_fraction) *
                               (context.end_time - context.begin_time);
  if (!close_stage_time(stage_point.stage_time, expected_time))
    throw std::invalid_argument(
        "stage spatial exact clock and physical time disagree");
}

void validate_metadata(const StageSpatialMetadata &metadata,
                       const StepContext &context) {
  if (metadata.patch_count != 1)
    throw std::invalid_argument(
        "stage spatial preparation supports exactly one patch");
  if (metadata.level_count <= 0 || metadata.level_count > 2)
    throw std::invalid_argument(
        "stage spatial preparation supports one or two levels");
  if (metadata.level < 0 || metadata.level >= metadata.level_count ||
      metadata.level != context.level)
    throw std::invalid_argument(
        "stage spatial transaction level changed");
  if (metadata.level_count == 2 &&
      metadata.spatial_refinement_ratio != 2)
    throw std::invalid_argument(
        "stage spatial preparation requires refinement ratio two");
  if (metadata.transaction_epoch < 0 ||
      metadata.transaction_epoch != metadata.observed_epoch)
    throw std::runtime_error("stage spatial hierarchy epoch changed");
  if (metadata.global_sync_ownership_conflict)
    throw std::logic_error(
        "global SYNC conflicts with driver-owned stage spatial fill");
  if (metadata.evolved_group_count == 0)
    throw std::invalid_argument(
        "stage spatial preparation has no evolved groups");
}

std::optional<step_clock_t>
validate_parent(const StageSpatialMetadata &metadata,
                const StepContext &context, const step_clock_t stage_clock,
                const DenseInterval *parent_dense) {
  if (metadata.level == 0) {
    if (parent_dense != nullptr)
      throw std::invalid_argument(
          "coarsest stage must not receive a parent dense interval");
    return std::nullopt;
  }
  if (parent_dense == nullptr)
    throw std::invalid_argument(
        "fine stage requires its covering parent dense interval");
  const auto &id = parent_dense->id();
  if (id.level != metadata.level - 1 || id.method != context.method ||
      parent_dense->capability().method != context.method)
    throw std::invalid_argument(
        "fine stage parent dense interval metadata differs");
  if (!(id.begin_clock < id.end_clock) || stage_clock < id.begin_clock ||
      stage_clock > id.end_clock)
    throw std::out_of_range(
        "fine stage is outside its parent dense interval");
  const auto theta =
      (stage_clock - id.begin_clock) / (id.end_clock - id.begin_clock);
  if (theta.den <= 0 || theta < 0 || theta > 1)
    throw std::out_of_range(
        "fine stage parent theta is outside [0,1]");
  return theta;
}

} // namespace

TwoLevelStageSpatialPreparer::TwoLevelStageSpatialPreparer(
    std::unique_ptr<detail::StageSpatialPreparerBackend> backend)
    : backend_(std::move(backend)) {
  if (backend_ == nullptr)
    throw std::invalid_argument("stage spatial backend is missing");
}

TwoLevelStageSpatialPreparer::~TwoLevelStageSpatialPreparer() = default;
TwoLevelStageSpatialPreparer::TwoLevelStageSpatialPreparer(
    TwoLevelStageSpatialPreparer &&) noexcept = default;
TwoLevelStageSpatialPreparer &TwoLevelStageSpatialPreparer::operator=(
    TwoLevelStageSpatialPreparer &&) noexcept = default;

StageSpatialPreparationReceipt TwoLevelStageSpatialPreparer::prepare_impl(
    ScratchStateTransaction *transaction, const StepContext &context,
    const StagePoint &stage_point, const DenseInterval *parent_dense) {
  try {
    const auto metadata = backend_->inspect(transaction);
    validate_metadata(metadata, context);
    const auto target = target_for(stage_point.kind);
    const auto stage_clock =
        context.begin_clock +
        (context.end_clock - context.begin_clock) *
            stage_point.stage_fraction;
    validate_stage(context, stage_point, stage_clock);
    const auto parent_theta =
        validate_parent(metadata, context, stage_clock, parent_dense);

    if (metadata.level == 0) {
      backend_->prepare_level_zero(transaction, target, stage_clock);
    } else {
      backend_->prepare_level_one(transaction, target, stage_clock,
                                  *parent_theta, *parent_dense);
    }
    // Spatial validity is promoted only after every exchange, boundary fill,
    // dense sample, and prolongation has returned successfully.
    backend_->promote(transaction, target);
    return StageSpatialPreparationReceipt{
        target, 0, metadata.level, stage_clock, parent_theta,
        metadata.evolved_group_count};
  } catch (...) {
    backend_->fault(transaction);
    throw;
  }
}

#ifdef CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_UNIT
namespace {

class TestStageSpatialBackend final
    : public detail::StageSpatialPreparerBackend {
public:
  explicit TestStageSpatialBackend(
      TwoLevelStageSpatialPreparer::TestBackend &test)
      : test_(test) {}

  StageSpatialMetadata inspect(ScratchStateTransaction *) override {
    const auto &m = test_.metadata;
    return {m.patch_count,
            m.level_count,
            m.spatial_refinement_ratio,
            m.level,
            m.transaction_epoch,
            m.observed_epoch,
            m.global_sync_ownership_conflict,
            m.evolved_group_count};
  }

  void prepare_level_zero(ScratchStateTransaction *,
                          const StageSpatialTarget target,
                          step_clock_t) override {
    ++test_.prepare_calls;
    test_.events.push_back("prepare L0 " + target_name(target));
    if (test_.throw_during_prepare)
      throw std::runtime_error("injected stage spatial fill failure");
  }

  void prepare_level_one(ScratchStateTransaction *,
                         const StageSpatialTarget target, step_clock_t,
                         const step_clock_t parent_theta,
                         const DenseInterval &) override {
    ++test_.prepare_calls;
    test_.observed_parent_theta = parent_theta;
    test_.events.push_back(
        "prepare L1 " + target_name(target) + " theta=" +
        std::to_string(parent_theta.num) + "/" +
        std::to_string(parent_theta.den));
    if (test_.throw_during_prepare)
      throw std::runtime_error("injected stage spatial fill failure");
  }

  void promote(ScratchStateTransaction *,
               const StageSpatialTarget target) override {
    ++test_.promote_calls;
    test_.events.push_back("promote " + target_name(target));
    if (test_.throw_during_promote)
      throw std::runtime_error("injected stage spatial promotion failure");
  }

  void fault(ScratchStateTransaction *) noexcept override {
    ++test_.fault_calls;
    test_.events.push_back("fault");
  }

private:
  TwoLevelStageSpatialPreparer::TestBackend &test_;
};

} // namespace

TwoLevelStageSpatialPreparer
TwoLevelStageSpatialPreparer::create_for_test(TestBackend &backend) {
  return TwoLevelStageSpatialPreparer(
      std::make_unique<TestStageSpatialBackend>(backend));
}

StageSpatialPreparationReceipt
TwoLevelStageSpatialPreparer::prepare_for_test(
    const StepContext &context, const StagePoint &stage_point,
    const DenseInterval *parent_dense) {
  return prepare_impl(nullptr, context, stage_point, parent_dense);
}
#endif

} // namespace CarpetX
