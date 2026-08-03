#include "subcycling_stage_spatial_preparer.hxx"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

using namespace CarpetX;

template <typename Exception = std::exception, typename Function>
void expect_throw(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const Exception &) {
    threw = true;
  }
  assert(threw);
}

TableauFingerprint fingerprint() {
  TableauFingerprint value{};
  value.fill(9);
  return value;
}

class ScalarState final : public DenseStateVector {
public:
  explicit ScalarState(const double value) : value_(value) {}

  bool compatible(const DenseStateVector &other) const noexcept override {
    return dynamic_cast<const ScalarState *>(&other) != nullptr;
  }

  void copy_from(const DenseStateVector &other) override {
    value_ = dynamic_cast<const ScalarState &>(other).value_;
  }

  void linear_combination(
      const std::vector<double> &weights,
      const std::vector<const DenseStateVector *> &sources) override {
    value_ = 0.0;
    for (std::size_t i = 0; i < weights.size(); ++i)
      value_ += weights[i] *
                dynamic_cast<const ScalarState &>(*sources[i]).value_;
  }

private:
  double value_;
};

std::shared_ptr<const DenseInterval> parent_interval() {
  const auto table = fingerprint();
  const DenseCapability capability{SubcyclingODEMethod::rk4, table, 4, 4,
                                   4, 4, 5, true, true};
  DenseIntervalBuilder builder(
      capability,
      DenseIntervalId{0, step_clock_t(0), step_clock_t(1), 10.0, 14.0,
                      SubcyclingODEMethod::rk4, table});
  for (int control = 0; control < 5; ++control)
    builder.add_control(std::make_unique<ScalarState>(control));
  return builder.seal();
}

StepContext level_context(const int level) {
  return StepContext{level, step_clock_t(0), step_clock_t(1), 10.0, 14.0,
                     SubcyclingODEMethod::rk4};
}

StagePoint stage(const StageKind kind, const step_clock_t fraction,
                 const double time, const int index = 1,
                 const int count = 1) {
  return StagePoint{kind, index, count, fraction, time};
}

void test_primary_level_zero_prepares_live_then_promotes() {
  TwoLevelStageSpatialPreparer::TestBackend backend;
  backend.metadata.level = 0;
  backend.metadata.level_count = 2;
  backend.metadata.evolved_group_count = 3;
  auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);

  const auto receipt = preparer.prepare_for_test(
      level_context(0), stage(StageKind::primary, step_clock_t(1, 2), 12.0),
      nullptr);

  assert(receipt.target == StageSpatialTarget::primary_live_tl0);
  assert(receipt.patch == 0);
  assert(receipt.level == 0);
  assert(receipt.stage_clock == step_clock_t(1, 2));
  assert(!receipt.parent_theta.has_value());
  assert(receipt.evolved_group_count == 3);
  assert((backend.events == std::vector<std::string>{"prepare L0 live",
                                                      "promote live"}));
  assert(backend.fault_calls == 0);
}

void test_fractional_level_one_uses_exact_parent_theta_and_scratch() {
  TwoLevelStageSpatialPreparer::TestBackend backend;
  backend.metadata.level = 1;
  backend.metadata.level_count = 2;
  backend.metadata.evolved_group_count = 2;
  auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);
  const auto parent = parent_interval();

  const auto receipt = preparer.prepare_for_test(
      level_context(1),
      stage(StageKind::fractional, step_clock_t(1, 4), 11.0, 2, 4),
      parent.get());

  assert(receipt.target == StageSpatialTarget::transaction_scratch);
  assert(receipt.level == 1);
  assert(receipt.stage_clock == step_clock_t(1, 4));
  assert(receipt.parent_theta == step_clock_t(1, 4));
  assert(backend.observed_parent_theta == step_clock_t(1, 4));
  assert((backend.events ==
          std::vector<std::string>{"prepare L1 scratch theta=1/4",
                                   "promote scratch"}));
  assert(backend.fault_calls == 0);
}

void test_endpoint_probe_is_transaction_private_scratch() {
  TwoLevelStageSpatialPreparer::TestBackend backend;
  backend.metadata.level = 0;
  auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);

  const auto receipt = preparer.prepare_for_test(
      level_context(0),
      stage(StageKind::endpoint_probe, step_clock_t(1), 14.0), nullptr);

  assert(receipt.target == StageSpatialTarget::transaction_scratch);
  assert((backend.events == std::vector<std::string>{"prepare L0 scratch",
                                                      "promote scratch"}));
}

void test_preflight_failures_do_not_prepare_or_promote_and_fault_once() {
  struct Case {
    int patch_count{1};
    int level_count{2};
    int ratio{2};
    std::int64_t observed_epoch{7};
    bool ownership_conflict{false};
    bool provide_parent{true};
  };
  const std::vector<Case> cases{
      Case{2, 2, 2, 7, false, true},
      Case{1, 3, 2, 7, false, true},
      Case{1, 2, 4, 7, false, true},
      Case{1, 2, 2, 8, false, true},
      Case{1, 2, 2, 7, true, true},
      Case{1, 2, 2, 7, false, false},
  };

  for (const auto &test : cases) {
    TwoLevelStageSpatialPreparer::TestBackend backend;
    backend.metadata.patch_count = test.patch_count;
    backend.metadata.level_count = test.level_count;
    backend.metadata.spatial_refinement_ratio = test.ratio;
    backend.metadata.level = 1;
    backend.metadata.transaction_epoch = 7;
    backend.metadata.observed_epoch = test.observed_epoch;
    backend.metadata.global_sync_ownership_conflict =
        test.ownership_conflict;
    auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);
    const auto parent = parent_interval();

    expect_throw([&] {
      preparer.prepare_for_test(
          level_context(1),
          stage(StageKind::fractional, step_clock_t(1, 4), 11.0, 2, 4),
          test.provide_parent ? parent.get() : nullptr);
    });
    assert(backend.prepare_calls == 0);
    assert(backend.promote_calls == 0);
    assert(backend.fault_calls == 1);
  }
}

void test_fill_failure_never_promotes_validity_and_faults_once() {
  TwoLevelStageSpatialPreparer::TestBackend backend;
  backend.metadata.level = 0;
  backend.throw_during_prepare = true;
  auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);

  expect_throw<std::runtime_error>([&] {
    preparer.prepare_for_test(
        level_context(0),
        stage(StageKind::fractional, step_clock_t(1, 2), 12.0), nullptr);
  });

  assert(backend.prepare_calls == 1);
  assert(backend.promote_calls == 0);
  assert(backend.fault_calls == 1);
  assert((backend.events == std::vector<std::string>{"prepare L0 scratch",
                                                      "fault"}));
}

void test_metadata_and_stage_drift_fail_closed() {
  TwoLevelStageSpatialPreparer::TestBackend backend;
  backend.metadata.level = 0;
  auto preparer = TwoLevelStageSpatialPreparer::create_for_test(backend);

  expect_throw([&] {
    preparer.prepare_for_test(
        level_context(1), stage(StageKind::primary, step_clock_t(0), 10.0),
        nullptr);
  });
  assert(backend.fault_calls == 1);

  TwoLevelStageSpatialPreparer::TestBackend bad_time;
  bad_time.metadata.level = 0;
  auto bad_time_preparer =
      TwoLevelStageSpatialPreparer::create_for_test(bad_time);
  expect_throw([&] {
    bad_time_preparer.prepare_for_test(
        level_context(0),
        stage(StageKind::primary, step_clock_t(1, 2), 12.25), nullptr);
  });
  assert(bad_time.fault_calls == 1);
}

} // namespace

int main() {
  test_primary_level_zero_prepares_live_then_promotes();
  test_fractional_level_one_uses_exact_parent_theta_and_scratch();
  test_endpoint_probe_is_transaction_private_scratch();
  test_preflight_failures_do_not_prepare_or_promote_and_fault_once();
  test_fill_failure_never_promotes_validity_and_faults_once();
  test_metadata_and_stage_drift_fail_closed();
  std::cout << "subcycling stage spatial preparer unit tests passed\n";
}
