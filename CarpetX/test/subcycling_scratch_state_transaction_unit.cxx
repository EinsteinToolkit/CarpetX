#include "subcycling_dense_output.hxx"
#include "subcycling_dense_mfab_state.hxx"
#include "subcycling_dense_stencil.hxx"
#include "subcycling_scratch_state_transaction.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"
#include "transaction_level_step_session.hxx"

#include <AMReX_MultiFab.H>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using namespace CarpetX;

static_assert(!std::is_copy_constructible_v<ScratchStateToken>);
static_assert(!std::is_copy_assignable_v<ScratchStateToken>);
static_assert(std::is_nothrow_move_constructible_v<ScratchStateToken>);
static_assert(std::is_nothrow_move_assignable_v<ScratchStateToken>);
static_assert(!std::is_copy_constructible_v<ScratchStateTransaction>);

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

struct Source {
  int storage_token{0};
  ScratchGFKey key{};
  std::unique_ptr<amrex::MultiFab> multifab;
  std::vector<ScratchValidity> validity;
  bool grid_function_real{true};
};

Source make_source(const int group, const double base,
                   const int box_identity = 17,
                   const int distribution_identity = 19) {
  static amrex::Arena arena(23);
  static amrex::TestFabFactory factory(29);
  Source source;
  source.storage_token = group + 100;
  source.key = ScratchGFKey{7, 0, 0, group};
  source.multifab = std::make_unique<amrex::MultiFab>(
      amrex::BoxArray(box_identity),
      amrex::DistributionMapping(distribution_identity), 2,
      amrex::IntVect(1, 2, 1), &arena, factory);
  source.validity = {{true, false, true}, {false, true, false}};
  for (std::size_t region = 0; region < 3; ++region)
    for (int component = 0; component < 2; ++component)
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell)
        source.multifab->test_set(
            static_cast<amrex::TestRegion>(region), component, cell,
            base + 100.0 * static_cast<double>(region) +
                10.0 * static_cast<double>(component) +
                static_cast<double>(cell));
  return source;
}

struct Fixture {
  int requested_level{0};
  int hierarchy_level_count{1};
  int time_refinement_factor{1};
  std::int64_t observed_epoch{7};
  int level_identity{501};
  Source evolved{make_source(10, 10.0)};
  Source rhs{make_source(11, 2.0)};
  Source dependent{make_source(12, 50.0)};
  std::vector<Source *> entries{&evolved, &rhs, &dependent};

  explicit Fixture(const int level = 0, const int level_count = 1)
      : requested_level(level), hierarchy_level_count(level_count),
        time_refinement_factor(level == 0 ? 1 : 2) {
    for (auto *const source : entries)
      source->key.level = requested_level;
  }

  ScratchLevelFrame copy_live_frame() {
    UncertifiedScratchLevelManifest manifest{7, requested_level, {}};
    std::vector<ScratchGFView> views;
    for (const auto *source : entries) {
      manifest.entries.push_back(
          {source->key, &source->storage_token});
      views.push_back({source->key, 0, &source->storage_token,
                       source->multifab.get(), &source->validity});
    }
    return ScratchLevelFrame::copy_tl0(manifest, views);
  }

  ScratchStateTransactionFactoryMetadata metadata() {
    ScratchStateTransactionFactoryMetadata result{
        7,
        1,
        hierarchy_level_count,
        41,
        step_clock_t(1, 2),
        0.5,
        time_refinement_factor,
        false,
        {{10, 11}},
        {11, 12},
        [this] { return observed_epoch; },
        {}};
    for (auto *const source : entries) {
      result.live_entry_readers.push_back([this, source] {
        return ScratchLiveEntrySnapshot{
            source->key, &level_identity, source, &source->multifab,
            &source->storage_token, &source->validity,
            source->multifab.get(), source->validity,
            source->grid_function_real,
            [source](const amrex::MultiFab &state,
                     const std::vector<ScratchValidity> &validity) {
              amrex::MultiFab::Copy(*source->multifab, state, 0, 0,
                                    state.nComp(), state.nGrowVect());
              source->validity = validity;
            }};
      });
    }
    return result;
  }
};

ScratchScheduleExecutionReceipt successful_schedule(
    const ScratchSchedulePhase phase, const StepContext &,
    const ScratchStageCoordinates &) {
  return {phase, 1, 1};
}

std::unique_ptr<ScratchStateTransaction>
make_transaction(Fixture &fixture,
                 ScratchStateTransactionFactory::TestScheduleExecutor executor =
                     successful_schedule) {
  return ScratchStateTransactionFactory::create_for_test(
      fixture.copy_live_frame(), fixture.metadata(), std::move(executor));
}

void assert_source_equal(const amrex::MultiFab &actual,
                         const Source &expected) {
  for (std::size_t region = 0; region < 3; ++region)
    for (int component = 0; component < 2; ++component)
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell)
        assert(actual.test_value(static_cast<amrex::TestRegion>(region),
                                 component, cell) ==
               expected.multifab->test_value(
                   static_cast<amrex::TestRegion>(region), component, cell));
}

void fill_all(amrex::MultiFab &multifab, const double value) {
  for (std::size_t region = 0; region < 3; ++region)
    for (int component = 0; component < multifab.nComp(); ++component)
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell)
        multifab.test_set(static_cast<amrex::TestRegion>(region), component,
                          cell, value);
}

struct ExactSourceSnapshot {
  std::vector<std::uint64_t> value_bits;
  std::vector<ScratchValidity> validity;
};

std::uint64_t bits_of(const double value) {
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

ExactSourceSnapshot exact_snapshot(const Source &source) {
  ExactSourceSnapshot snapshot;
  for (std::size_t region = 0; region < 3; ++region)
    for (int component = 0; component < source.multifab->nComp(); ++component)
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell)
        snapshot.value_bits.push_back(bits_of(source.multifab->test_value(
            static_cast<amrex::TestRegion>(region), component, cell)));
  snapshot.validity = source.validity;
  return snapshot;
}

void assert_exact_snapshot(const Source &source,
                           const ExactSourceSnapshot &expected) {
  std::size_t value = 0;
  for (std::size_t region = 0; region < 3; ++region)
    for (int component = 0; component < source.multifab->nComp(); ++component)
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell)
        assert(bits_of(source.multifab->test_value(
                   static_cast<amrex::TestRegion>(region), component,
                   cell)) == expected.value_bits.at(value++));
  assert(value == expected.value_bits.size());
  assert(source.validity.size() == expected.validity.size());
  for (std::size_t component = 0; component < source.validity.size();
       ++component) {
    assert(source.validity[component].interior ==
           expected.validity[component].interior);
    assert(source.validity[component].outer ==
           expected.validity[component].outer);
    assert(source.validity[component].ghosts ==
           expected.validity[component].ghosts);
  }
}

void overwrite_live(Fixture &fixture, const double value) {
  for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry) {
    auto &source = *fixture.entries[entry];
    fill_all(*source.multifab, value + static_cast<double>(entry));
    for (auto &validity : source.validity)
      validity = {false, true, false};
  }
}

void test_mid_primary_failure_rolls_live_tl0_back_exactly_without_schedule() {
  Fixture fixture;
  int schedule_calls = 0;
  auto transaction = make_transaction(
      fixture,
      [&schedule_calls](const ScratchSchedulePhase phase,
                        const StepContext &,
                        const ScratchStageCoordinates &) {
        ++schedule_calls;
        return ScratchScheduleExecutionReceipt{phase, 1, 1};
      });
  auto rollback_state = transaction->capture_live_evolved();
  std::vector<ExactSourceSnapshot> expected;
  expected.reserve(fixture.entries.size());
  for (const auto *const source : fixture.entries)
    expected.push_back(exact_snapshot(*source));

  try {
    overwrite_live(fixture, -9000.0);
    throw std::runtime_error("injected primary-stage failure");
  } catch (const std::runtime_error &) {
    transaction->rollback_live_evolved(rollback_state);
  }

  for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
    assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
  assert(schedule_calls == 0);
  assert(transaction->discarded());
  assert(!transaction->faulted());
}

void test_successful_primary_endpoint_is_not_rolled_back() {
  Fixture fixture;
  int schedule_calls = 0;
  auto transaction = make_transaction(
      fixture,
      [&schedule_calls](const ScratchSchedulePhase phase,
                        const StepContext &,
                        const ScratchStageCoordinates &) {
        ++schedule_calls;
        return ScratchScheduleExecutionReceipt{phase, 1, 1};
      });
  auto rollback_state = transaction->capture_live_evolved();
  assert(rollback_state.valid());
  overwrite_live(fixture, 4200.0);
  std::vector<ExactSourceSnapshot> accepted_endpoint;
  accepted_endpoint.reserve(fixture.entries.size());
  for (const auto *const source : fixture.entries)
    accepted_endpoint.push_back(exact_snapshot(*source));

  transaction->discard();

  for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
    assert_exact_snapshot(*fixture.entries[entry], accepted_endpoint[entry]);
  assert(schedule_calls == 0);
  assert(transaction->discarded());
  assert(!transaction->faulted());
}

void test_live_rollback_fails_closed_for_wrong_kind_or_live_drift() {
  {
    Fixture fixture;
    auto transaction = make_transaction(fixture);
    auto rhs = transaction->capture_live_rhs();
    expect_throw<std::invalid_argument>(
        [&] { transaction->rollback_live_evolved(rhs); });
    assert(transaction->faulted());
    assert(transaction->discarded());
  }
  {
    Fixture fixture;
    bool drift = false;
    auto metadata = fixture.metadata();
    const auto original_reader = metadata.live_entry_readers.front();
    int replacement_identity = 999;
    metadata.live_entry_readers.front() =
        [original_reader, &replacement_identity, &drift] {
          auto snapshot = original_reader();
          if (drift)
            snapshot.storage_identity = &replacement_identity;
          return snapshot;
        };
    auto transaction = ScratchStateTransactionFactory::create_for_test(
        fixture.copy_live_frame(), std::move(metadata), successful_schedule);
    auto source = transaction->capture_live_evolved();
    drift = true;
    overwrite_live(fixture, -1300.0);
    const auto still_mutated = exact_snapshot(fixture.evolved);
    expect_throw<std::runtime_error>(
        [&] { transaction->rollback_live_evolved(source); });
    assert_exact_snapshot(fixture.evolved, still_mutated);
    assert(transaction->faulted());
    assert(transaction->discarded());
  }
}

void test_full_frame_round_trip_and_fixed_working_addresses() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  const auto state = transaction->capture_live_evolved();
  std::array<const void *, 3> addresses{};
  for (std::size_t entry = 0; entry < addresses.size(); ++entry)
    addresses[entry] =
        ScratchStateTransactionFactory::working_entry_address_for_test(
            *transaction, entry);

  for (int cycle = 0; cycle < 4; ++cycle) {
    for (std::size_t entry = 0; entry < 3; ++entry) {
      fill_all(ScratchStateTransactionFactory::mutable_working_multifab_for_test(
                   *transaction, entry),
               -1000.0 - cycle);
      auto &valid =
          ScratchStateTransactionFactory::mutable_working_validity_for_test(
              *transaction, entry);
      for (auto &component : valid)
        component = {false, false, false};
    }
    transaction->restore_state(state);
    assert_source_equal(
        ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                  0),
        fixture.evolved);
    assert_source_equal(
        ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                  1),
        fixture.rhs);
    assert_source_equal(
        ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                  2),
        fixture.dependent);
    for (std::size_t entry = 0; entry < 3; ++entry) {
      assert(ScratchStateTransactionFactory::working_entry_address_for_test(
                 *transaction, entry) == addresses[entry]);
      const auto &expected = fixture.entries[entry]->validity;
      const auto &actual =
          ScratchStateTransactionFactory::working_validity_for_test(
              *transaction, entry);
      assert(actual.size() == expected.size());
      for (std::size_t component = 0; component < actual.size(); ++component) {
        assert(actual[component].interior == expected[component].interior);
        assert(actual[component].outer == expected[component].outer);
        assert(actual[component].ghosts == expected[component].ghosts);
      }
    }
  }
}

void test_pair_order_rhs_canonicalization_and_mixed_linear_combination() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  assert(transaction->group_pairs().size() == 1);
  assert(transaction->group_pairs()[0].evolved_group == 10);
  assert(transaction->group_pairs()[0].rhs_group == 11);

  auto y = transaction->capture_live_evolved();
  auto f = transaction->capture_live_rhs();
  auto destination = transaction->clone_state(y);
  transaction->linear_combination(destination, 1.0, {{0.5, &f}});
  transaction->restore_state(destination);

  const auto &result =
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                0);
  for (int component = 0; component < 2; ++component)
    for (std::size_t cell = 0;
         cell < amrex::MultiFab::test_cells_per_region(); ++cell) {
      const double expected =
          fixture.evolved.multifab->test_value(amrex::TestRegion::valid,
                                               component, cell) +
          0.5 * fixture.rhs.multifab->test_value(amrex::TestRegion::valid,
                                                 component, cell);
      assert(result.test_value(amrex::TestRegion::valid, component, cell) ==
             expected);
      assert(result.test_value(amrex::TestRegion::outer, component, cell) ==
             fixture.evolved.multifab->test_value(amrex::TestRegion::outer,
                                                  component, cell));
      assert(result.test_value(amrex::TestRegion::ghosts, component, cell) ==
             fixture.evolved.multifab->test_value(amrex::TestRegion::ghosts,
                                                  component, cell));
    }

  auto accumulator = transaction->clone_state(f);
  transaction->linear_combination(accumulator, 1.0, {{2.0, &f}});
  auto from_accumulator = transaction->clone_state(y);
  transaction->linear_combination(from_accumulator, 0.0,
                                  {{1.0, &accumulator}});
  transaction->restore_state(from_accumulator);
  const auto &accumulated =
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                0);
  assert(accumulated.test_value(amrex::TestRegion::valid, 0, 0) ==
         3.0 * fixture.rhs.multifab->test_value(amrex::TestRegion::valid, 0,
                                                0));

  auto self_alias = transaction->clone_state(y);
  transaction->linear_combination(self_alias, 2.0,
                                  {{3.0, &self_alias}});
  transaction->restore_state(self_alias);
  assert(ScratchStateTransactionFactory::working_multifab_for_test(
             *transaction, 0)
             .test_value(amrex::TestRegion::valid, 0, 0) ==
         5.0 * fixture.evolved.multifab->test_value(
                   amrex::TestRegion::valid, 0, 0));

  auto validity_source = transaction->clone_state(y);
  transaction->set_state_valid(validity_source, ScratchStateRegion::interior,
                               true);
  auto validity_destination = transaction->clone_state(y);
  transaction->set_state_valid(validity_destination,
                               ScratchStateRegion::interior, false);
  transaction->linear_combination(validity_destination, 0.0,
                                  {{1.0, &validity_source}, {0.0, &f}});
  assert(transaction->state_valid(validity_destination,
                                  ScratchStateRegion::interior));
}

void test_restore_left_overlays_only_physical_rhs_interior_and_validity() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  auto evolved = transaction->capture_live_evolved();

  const auto original_rhs_outer =
      fixture.rhs.multifab->test_value(amrex::TestRegion::outer, 0, 0);
  const auto original_rhs_ghost =
      fixture.rhs.multifab->test_value(amrex::TestRegion::ghosts, 0, 0);
  const auto original_rhs_validity = fixture.rhs.validity;
  for (int component = 0; component < fixture.rhs.multifab->nComp();
       ++component) {
    for (std::size_t cell = 0;
         cell < amrex::MultiFab::test_cells_per_region(); ++cell) {
      fixture.rhs.multifab->test_set(amrex::TestRegion::valid, component,
                                     cell, 700.0 + component + cell);
      fixture.rhs.multifab->test_set(amrex::TestRegion::outer, component,
                                     cell, 800.0 + component + cell);
      fixture.rhs.multifab->test_set(amrex::TestRegion::ghosts, component,
                                     cell, 900.0 + component + cell);
    }
    fixture.rhs.validity[static_cast<std::size_t>(component)] =
        {true, true, true};
  }
  auto raw_rhs = transaction->capture_live_rhs();
  assert(transaction->state_valid(raw_rhs, ScratchStateRegion::interior));
  assert(!transaction->state_valid(raw_rhs, ScratchStateRegion::outer));
  assert(!transaction->state_valid(raw_rhs, ScratchStateRegion::ghosts));
  for (std::size_t entry = 0; entry < 3; ++entry) {
    fill_all(ScratchStateTransactionFactory::mutable_working_multifab_for_test(
                 *transaction, entry),
             -600.0);
    for (auto &validity :
         ScratchStateTransactionFactory::mutable_working_validity_for_test(
             *transaction, entry))
      validity = {false, false, false};
  }
  transaction->restore_left(evolved, raw_rhs);
  assert_source_equal(
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                0),
      fixture.evolved);
  const auto &restored_rhs =
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                1);
  assert(restored_rhs.test_value(amrex::TestRegion::valid, 0, 0) == 700.0);
  assert(restored_rhs.test_value(amrex::TestRegion::outer, 0, 0) ==
         original_rhs_outer);
  assert(restored_rhs.test_value(amrex::TestRegion::ghosts, 0, 0) ==
         original_rhs_ghost);
  assert_source_equal(
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                2),
      fixture.dependent);
  const auto &rhs_validity =
      ScratchStateTransactionFactory::working_validity_for_test(*transaction,
                                                                1);
  for (std::size_t component = 0; component < rhs_validity.size(); ++component) {
    assert(rhs_validity[component].interior);
    assert(rhs_validity[component].outer ==
           original_rhs_validity[component].outer);
    assert(rhs_validity[component].ghosts ==
           original_rhs_validity[component].ghosts);
  }
}

void test_hierarchy_epoch_drift_faults_before_state_mutation() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  auto state = transaction->capture_live_evolved();
  fixture.observed_epoch = 8;

  expect_throw<std::runtime_error>([&] { transaction->restore_state(state); });
  assert(transaction->faulted());
  assert(transaction->discarded());
  expect_throw<std::logic_error>([&] {
    transaction->linear_combination(state, 1.0, {});
  });
}

void test_clean_discard_remains_nonfaulting_when_dense_take_is_empty() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  auto state = transaction->capture_live_evolved();
  transaction->discard();
  expect_throw<std::logic_error>(
      [&] { transaction->rollback_live_evolved(state); });
  assert(transaction->take_committed_dense() == nullptr);
  assert(transaction->discarded());
  assert(!transaction->faulted());
}

void test_zero_scale_and_zero_terms_do_not_read_nan_sources() {
  Fixture fixture;
  fill_all(*fixture.evolved.multifab,
           std::numeric_limits<double>::quiet_NaN());
  auto transaction = make_transaction(fixture);
  auto destination = transaction->capture_live_evolved();
  fill_all(*fixture.evolved.multifab, 4.0);
  auto finite_source = transaction->capture_live_evolved();
  fill_all(*fixture.rhs.multifab,
           std::numeric_limits<double>::quiet_NaN());
  auto zero_term = transaction->capture_live_rhs();
  transaction->linear_combination(destination, 0.0,
                                  {{1.0, &finite_source},
                                   {0.0, &zero_term}});
  transaction->restore_state(destination);
  const auto &result =
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                0);
  for (int component = 0; component < 2; ++component)
    for (std::size_t cell = 0;
         cell < amrex::MultiFab::test_cells_per_region(); ++cell)
      assert(result.test_value(amrex::TestRegion::valid, component, cell) ==
             4.0);
}

void test_owner_epoch_kind_and_schema_rejections() {
  Fixture first_fixture;
  Fixture second_fixture;
  auto first = make_transaction(first_fixture);
  auto second = make_transaction(second_fixture);
  auto first_state = first->capture_live_evolved();
  auto second_state = second->capture_live_evolved();
  expect_throw<std::invalid_argument>(
      [&] { first->restore_state(second_state); });

  auto rhs = first->capture_live_rhs();
  expect_throw<std::invalid_argument>([&] { first->restore_state(rhs); });
  auto stale =
      ScratchStateTransactionFactory::stale_epoch_token_for_test(*first,
                                                                  first_state);
  expect_throw<std::invalid_argument>([&] { first->restore_state(stale); });
  auto incompatible =
      ScratchStateTransactionFactory::incompatible_schema_token_for_test(
          *first, first_state);
  expect_throw<std::invalid_argument>(
      [&] { first->restore_state(incompatible); });
}

void test_post_step_sidecars_and_schedule_coordinates() {
  Fixture fixture;
  ScratchStateTransaction *transaction_ptr = nullptr;
  std::vector<ScratchStageCoordinates> coordinates;
  auto transaction = make_transaction(
      fixture,
      [&](const ScratchSchedulePhase phase, const StepContext &,
          const ScratchStageCoordinates &stage) {
        coordinates.push_back(stage);
        if (phase == ScratchSchedulePhase::post_step) {
          const auto &evolved_validity =
              ScratchStateTransactionFactory::working_validity_for_test(
                  *transaction_ptr, 0);
          assert(evolved_validity[0].interior ==
                 fixture.evolved.validity[0].interior);
          assert(evolved_validity[1].interior ==
                 fixture.evolved.validity[1].interior);
          for (const auto entry : {1U, 2U})
            for (const auto &validity :
                 ScratchStateTransactionFactory::working_validity_for_test(
                     *transaction_ptr, entry))
              assert(!validity.interior && !validity.outer &&
                     !validity.ghosts);
          auto &dependent = ScratchStateTransactionFactory::
              mutable_working_multifab_for_test(*transaction_ptr, 2);
          dependent.test_set(amrex::TestRegion::valid, 0, 0, 777.0);
          ScratchStateTransactionFactory::mutable_working_validity_for_test(
              *transaction_ptr, 2)[0] = {true, true, true};
        }
        return ScratchScheduleExecutionReceipt{phase, 1, 1};
      });
  transaction_ptr = transaction.get();
  const StepContext context{0, step_clock_t(3, 2), step_clock_t(2, 1), 1.5,
                            2.0, SubcyclingODEMethod::rk4};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    const StagePoint stage{StageKind::fractional, 2, 4,
                           step_clock_t(1, 2), 1.75};
    transaction->post_step_after_update(context, stage);
    transaction->evaluate_rhs(context, stage);
  }
  assert(coordinates.size() == 2);
  for (const auto &stage : coordinates) {
    assert(stage.level_iteration == 41);
    assert(stage.stage_time == 1.75);
    assert(stage.base_delta_clock == step_clock_t(1, 2));
    assert(stage.base_delta_time == 0.5);
    assert(stage.time_refinement_factor == 1);
  }
  assert(transaction->rhs_evaluation_count() == 1);

  const auto post_state = transaction->capture_scratch_evolved();
  const auto base_state = transaction->capture_live_evolved();
  for (int repeat = 0; repeat < 3; ++repeat) {
    transaction->restore_state(base_state);
    assert(ScratchStateTransactionFactory::working_multifab_for_test(
               *transaction, 2)
               .test_value(amrex::TestRegion::valid, 0, 0) == 50.0);
    transaction->restore_state(post_state);
    assert(ScratchStateTransactionFactory::working_multifab_for_test(
               *transaction, 2)
               .test_value(amrex::TestRegion::valid, 0, 0) == 777.0);
    assert(ScratchStateTransactionFactory::working_validity_for_test(
               *transaction, 2)[0]
               .ghosts);
  }
}

void test_level_one_timefac_two_uses_exact_stage_clock() {
  Fixture fixture(1, 2);
  std::vector<ScratchStageCoordinates> coordinates;
  auto transaction = make_transaction(
      fixture,
      [&](const ScratchSchedulePhase phase, const StepContext &context,
          const ScratchStageCoordinates &stage) {
        assert(context.level == 1);
        coordinates.push_back(stage);
        return ScratchScheduleExecutionReceipt{phase, 1, 1};
      });
  const StepContext context{1, step_clock_t(1, 2), step_clock_t(3, 4), 1.0,
                            1.25, SubcyclingODEMethod::rk4};
  const StagePoint stage{StageKind::fractional, 1, 3, step_clock_t(1, 3),
                         1.0 + 0.25 / 3.0};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    transaction->post_step_after_update(context, stage);
    transaction->evaluate_rhs(context, stage);
  }
  assert(coordinates.size() == 2);
  for (const auto &observed : coordinates) {
    assert(observed.level_iteration == 41);
    assert(observed.stage_time == stage.stage_time);
    assert(observed.base_delta_clock == step_clock_t(1, 2));
    assert(observed.base_delta_time == 0.5);
    assert(observed.time_refinement_factor == 2);
  }
  assert(transaction->rhs_evaluation_count() == 1);
}

void test_level_zero_remains_available_in_a_two_level_hierarchy() {
  Fixture fixture(0, 2);
  auto transaction = make_transaction(fixture);
  const StepContext context{0, step_clock_t(1), step_clock_t(3, 2), 2.0,
                            2.5, SubcyclingODEMethod::rk4};
  const StagePoint stage{StageKind::endpoint_probe, 1, 1, step_clock_t(1),
                         2.5};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    transaction->evaluate_rhs(context, stage);
  }
  assert(transaction->rhs_evaluation_count() == 1);
  assert(!transaction->faulted());
  assert(!transaction->discarded());
}

void test_level_one_factory_and_context_mismatches_fail_closed() {
  {
    Fixture fixture(1, 1);
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(make_transaction(fixture));
    });
  }
  {
    Fixture fixture(1, 3);
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(make_transaction(fixture));
    });
  }
  {
    Fixture fixture(1, 2);
    auto metadata = fixture.metadata();
    metadata.time_refinement_factor = 1;
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
  {
    Fixture fixture(1, 2);
    for (auto *const source : fixture.entries)
      source->key.patch = 1;
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(make_transaction(fixture));
    });
  }
  {
    Fixture fixture(1, 2);
    auto metadata = fixture.metadata();
    const auto reader = metadata.live_entry_readers.front();
    metadata.live_entry_readers.front() = [reader] {
      auto snapshot = reader();
      snapshot.key.level = 0;
      return snapshot;
    };
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
  {
    Fixture fixture(1, 2);
    auto transaction = make_transaction(fixture);
    const StepContext wrong_context{0, step_clock_t(1, 2),
                                    step_clock_t(3, 4), 1.0, 1.25,
                                    SubcyclingODEMethod::rk4};
    const StagePoint stage{StageKind::primary, 1, 4, step_clock_t(0), 1.0};
    struct NoopPreparer final : StagePreparer {
      void prepare_stage(const StepContext &, const StagePoint &) override {}
    } preparer;
    ScopedStepContext scope(wrong_context, preparer, transaction.get());
    expect_throw<std::invalid_argument>(
        [&] { transaction->evaluate_rhs(wrong_context, stage); });
    assert(transaction->faulted());
    assert(transaction->discarded());
  }
  {
    Fixture fixture(1, 2);
    auto transaction = make_transaction(fixture);
    const StepContext wrong_clock{1, step_clock_t(0), step_clock_t(1, 2),
                                  0.0, 0.5,
                                  SubcyclingODEMethod::rk4};
    const StagePoint stage{StageKind::fractional, 1, 2,
                           step_clock_t(1, 2), 0.25};
    struct NoopPreparer final : StagePreparer {
      void prepare_stage(const StepContext &, const StagePoint &) override {}
    } preparer;
    ScopedStepContext scope(wrong_clock, preparer, transaction.get());
    expect_throw<std::invalid_argument>(
        [&] { transaction->evaluate_rhs(wrong_clock, stage); });
    assert(transaction->faulted());
    assert(transaction->discarded());
  }
}

DenseCapability capability(const SubcyclingODEMethod method) {
  TableauFingerprint fingerprint{};
  fingerprint[0] = static_cast<std::uint8_t>(method) + 1;
  const int controls = method == SubcyclingODEMethod::rk4
                           ? 5
                           : method == SubcyclingODEMethod::rkf78_order7 ? 8
                                                                        : 9;
  const int stages = method == SubcyclingODEMethod::rk4
                         ? 4
                         : method == SubcyclingODEMethod::rkf78_order7 ? 11
                                                                      : 13;
  const int order = controls - 1;
  return {method, fingerprint, order, order, stages, 0, controls, true, true};
}

void publish_dense_for_session(ScratchStateTransaction &transaction,
                               const StepContext &context) {
  const auto cap = capability(context.method);
  DenseOutputProvider provider(cap);
  const DenseIntervalId interval{
      context.level, context.begin_clock, context.end_clock,
      context.begin_time, context.end_time, context.method,
      cap.tableau_fingerprint};
  const auto &constraints =
      reference_dense_stencil(context.method).specification().constraints;
  std::vector<ScratchStateToken> states;
  std::vector<ScratchDenseSampleRef> samples;
  states.reserve(constraints.size());
  samples.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    const bool value = constraint.kind == DenseSampleKind::value;
    states.push_back(value ? transaction.capture_live_evolved()
                           : transaction.capture_live_rhs());
    samples.push_back({constraint.theta,
                       value ? ScratchDenseSampleKind::value
                             : ScratchDenseSampleKind::raw_derivative,
                       &states.back()});
  }
  transaction.commit_dense(context, provider, interval, samples);
}

struct NoopSessionPreparer final : StagePreparer {
  void prepare_stage(const StepContext &, const StagePoint &) override {}
};

void test_transaction_level_session_accepts_dense_and_nondense_endpoints() {
  {
    Fixture fixture;
    for (auto *const source : {&fixture.evolved, &fixture.rhs})
      for (auto &validity : source->validity)
        validity = {true, true, true};
    const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0,
                              0.5, SubcyclingODEMethod::rk4};
    int evolution_calls = 0;
    int commit_calls = 0;
    auto transaction = make_transaction(fixture);
    auto *const transaction_address = transaction.get();
    TransactionLevelStepSession session(
        std::move(transaction), true,
        [&](ScratchStateTransaction &active) {
          ++evolution_calls;
          publish_dense_for_session(active, context);
          overwrite_live(fixture, 6100.0);
        },
        [&] { ++commit_calls; });
    assert(session.transaction() == transaction_address);

    LevelAdvanceResult result;
    NoopSessionPreparer preparer;
    {
      ScopedStepContext scope(context, preparer, session.transaction());
      result = session.advance();
    }
    assert(result.dense_interval != nullptr);
    assert(result.dense_interval->id().begin_clock == context.begin_clock);
    assert(result.dense_interval->id().end_clock == context.end_clock);
    assert(evolution_calls == 1);
    assert(commit_calls == 0);
    expect_throw<std::logic_error>([&] { session.advance(); });

    std::vector<ExactSourceSnapshot> accepted;
    for (const auto *const source : fixture.entries)
      accepted.push_back(exact_snapshot(*source));
    session.commit();
    assert(commit_calls == 1);
    assert(session.transaction() == nullptr);
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], accepted[entry]);
    expect_throw<std::logic_error>([&] { session.commit(); });
    expect_throw<std::logic_error>([&] { session.advance(); });
  }

  {
    Fixture fixture;
    const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0,
                              0.5, SubcyclingODEMethod::rk4};
    int evolution_calls = 0;
    int commit_calls = 0;
    TransactionLevelStepSession session(
        make_transaction(fixture), false,
        [&](ScratchStateTransaction &) {
          ++evolution_calls;
          overwrite_live(fixture, 7100.0);
        },
        [&] { ++commit_calls; });
    expect_throw<std::logic_error>([&] { session.commit(); });
    NoopSessionPreparer preparer;
    LevelAdvanceResult result;
    {
      ScopedStepContext scope(context, preparer, session.transaction());
      result = session.advance();
    }
    assert(result.dense_interval == nullptr);
    session.commit();
    assert(evolution_calls == 1);
    assert(commit_calls == 1);
  }
}

void test_transaction_level_session_rolls_back_evolution_or_commit_failure() {
  {
    Fixture fixture;
    std::vector<ExactSourceSnapshot> expected;
    for (const auto *const source : fixture.entries)
      expected.push_back(exact_snapshot(*source));
    int evolution_calls = 0;
    int commit_calls = 0;
    TransactionLevelStepSession session(
        make_transaction(fixture), false,
        [&](ScratchStateTransaction &) {
          ++evolution_calls;
          overwrite_live(fixture, -8100.0);
          throw std::runtime_error("injected level evolution failure");
        },
        [&] { ++commit_calls; });
    expect_throw<std::runtime_error>([&] { session.advance(); });
    assert(evolution_calls == 1);
    assert(commit_calls == 0);
    assert(session.transaction() == nullptr);
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
    expect_throw<std::logic_error>([&] { session.advance(); });
    expect_throw<std::logic_error>([&] { session.commit(); });
  }

  {
    Fixture fixture;
    std::vector<ExactSourceSnapshot> expected;
    for (const auto *const source : fixture.entries)
      expected.push_back(exact_snapshot(*source));
    TransactionLevelStepSession session(
        make_transaction(fixture), false,
        [&](ScratchStateTransaction &) { overwrite_live(fixture, -9100.0); },
        [] { throw std::runtime_error("injected accepted-metadata failure"); });
    const auto result = session.advance();
    assert(result.dense_interval == nullptr);
    expect_throw<std::runtime_error>([&] { session.commit(); });
    assert(session.transaction() == nullptr);
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
  }
}

void test_transaction_level_session_rolls_back_terminal_transaction_failures() {
  const auto exercise_schedule_failure =
      [](const ScratchSchedulePhase failing_phase, const double poison) {
        Fixture fixture;
        std::vector<ExactSourceSnapshot> expected;
        for (const auto *const source : fixture.entries)
          expected.push_back(exact_snapshot(*source));
        int schedule_calls = 0;
        const StepContext context{0, step_clock_t(3, 2), step_clock_t(2, 1),
                                  1.5, 2.0, SubcyclingODEMethod::rk4};
        const StagePoint stage{StageKind::fractional, 2, 4,
                               step_clock_t(1, 2), 1.75};
        TransactionLevelStepSession session(
            make_transaction(
                fixture,
                [&](const ScratchSchedulePhase phase, const StepContext &,
                    const ScratchStageCoordinates &)
                    -> ScratchScheduleExecutionReceipt {
                  ++schedule_calls;
                  if (phase == failing_phase)
                    throw std::runtime_error("injected scratch schedule failure");
                  return {phase, 1, 1};
                }),
            false,
            [&](ScratchStateTransaction &active) {
              overwrite_live(fixture, poison);
              if (failing_phase == ScratchSchedulePhase::post_step)
                active.post_step_after_update(context, stage);
              else
                active.evaluate_rhs(context, stage);
            });
        NoopSessionPreparer preparer;
        {
          ScopedStepContext scope(context, preparer, session.transaction());
          expect_throw<std::runtime_error>([&] { session.advance(); });
        }
        assert(schedule_calls == 1);
        assert(session.transaction() == nullptr);
        for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
          assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
      };

  exercise_schedule_failure(ScratchSchedulePhase::post_step, -11100.0);
  exercise_schedule_failure(ScratchSchedulePhase::rhs, -11200.0);

  Fixture fixture;
  std::vector<ExactSourceSnapshot> expected;
  for (const auto *const source : fixture.entries)
    expected.push_back(exact_snapshot(*source));
  TransactionLevelStepSession session(
      make_transaction(fixture), false,
      [&](ScratchStateTransaction &active) {
        overwrite_live(fixture, -11300.0);
        active.discard();
        throw std::runtime_error("injected solver discard");
      });
  expect_throw<std::runtime_error>([&] { session.advance(); });
  assert(session.transaction() == nullptr);
  for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
    assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
}

void test_armed_live_rollback_disarms_and_fails_closed_on_stale_epoch() {
  {
    Fixture fixture;
    auto transaction = make_transaction(fixture);
    transaction->arm_live_evolved_rollback();
    overwrite_live(fixture, 11400.0);
    std::vector<ExactSourceSnapshot> accepted;
    for (const auto *const source : fixture.entries)
      accepted.push_back(exact_snapshot(*source));

    transaction->disarm_live_evolved_rollback();
    transaction->discard();
    expect_throw<std::logic_error>(
        [&] { transaction->rollback_live_evolved(); });
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], accepted[entry]);
    assert(transaction->discarded());
    assert(!transaction->faulted());
  }

  {
    Fixture fixture;
    auto transaction = make_transaction(fixture);
    transaction->arm_live_evolved_rollback();
    overwrite_live(fixture, -11500.0);
    std::vector<ExactSourceSnapshot> still_mutated;
    for (const auto *const source : fixture.entries)
      still_mutated.push_back(exact_snapshot(*source));
    fixture.observed_epoch = 8;

    expect_throw<std::runtime_error>(
        [&] { transaction->rollback_live_evolved(); });
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], still_mutated[entry]);
    assert(transaction->faulted());
    assert(transaction->discarded());

    fixture.observed_epoch = 7;
    expect_throw<std::logic_error>(
        [&] { transaction->rollback_live_evolved(); });
    for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
      assert_exact_snapshot(*fixture.entries[entry], still_mutated[entry]);
  }
}

void test_uncommitted_transaction_level_session_destructor_rolls_back() {
  Fixture fixture;
  std::vector<ExactSourceSnapshot> expected;
  for (const auto *const source : fixture.entries)
    expected.push_back(exact_snapshot(*source));
  int commit_calls = 0;
  {
    TransactionLevelStepSession session(
        make_transaction(fixture), false,
        [&](ScratchStateTransaction &) { overwrite_live(fixture, -10100.0); },
        [&] { ++commit_calls; });
    const auto result = session.advance();
    assert(result.dense_interval == nullptr);
    assert(commit_calls == 0);
  }
  assert(commit_calls == 0);
  for (std::size_t entry = 0; entry < fixture.entries.size(); ++entry)
    assert_exact_snapshot(*fixture.entries[entry], expected[entry]);
}

void test_dense_commit_is_atomic_and_uses_phase7_builder() {
  Fixture fixture(1, 2);
  for (auto *source : {&fixture.evolved, &fixture.rhs})
    for (auto &validity : source->validity)
      validity.interior = true;
  auto transaction = ScratchStateTransactionFactory::create_for_test(
      fixture.copy_live_frame(), fixture.metadata(), successful_schedule);
  const auto method = SubcyclingODEMethod::rk4;
  const auto cap = capability(method);
  DenseOutputProvider provider(cap);
  const StepContext context{1, step_clock_t(0), step_clock_t(1, 4), 0.0,
                            0.25, method};
  DenseIntervalId id{1, step_clock_t(0), step_clock_t(1, 4), 0.0, 0.25,
                     method, cap.tableau_fingerprint};

  const auto &constraints =
      reference_dense_stencil(method).specification().constraints;
  auto y0 = transaction->capture_live_evolved();
  auto raw_f = transaction->capture_live_rhs();
  std::vector<ScratchStateToken> owned_samples;
  owned_samples.reserve(constraints.size());
  std::vector<ScratchDenseSampleRef> refs;
  refs.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    const bool value = constraint.kind == DenseSampleKind::value;
    if (value) {
      auto sample = transaction->clone_state(y0);
      transaction->linear_combination(
          sample, 1.0, {{constraint.theta * 0.25, &raw_f}});
      owned_samples.push_back(std::move(sample));
    } else {
      owned_samples.push_back(transaction->clone_state(raw_f));
    }
    refs.push_back({constraint.theta,
                    value ? ScratchDenseSampleKind::value
                          : ScratchDenseSampleKind::raw_derivative,
                    &owned_samples.back()});
  }
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    transaction->commit_dense(context, provider, id, refs);
  }
  const auto interval = transaction->take_committed_dense();
  assert(interval != nullptr);
  assert(interval->id().level == 1);
  assert(interval->control_count() == 5);
  assert(transaction->take_committed_dense() == nullptr);

  auto destination = OwnedMultiFabDenseState::copy_of(
      {{{7, 0, 1, 10}, fixture.evolved.multifab.get()}});
  constexpr double theta = 0.37;
  interval->evaluate(theta, *destination);
  const double expected =
      fixture.evolved.multifab->test_value(amrex::TestRegion::valid, 0, 0) +
      theta * 0.25 *
          fixture.rhs.multifab->test_value(amrex::TestRegion::valid, 0, 0);
  assert(std::abs(destination->multifab(0).test_value(
                      amrex::TestRegion::valid, 0, 0) -
                  expected) < 1.0e-11);
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    expect_throw<std::logic_error>(
        [&] { transaction->commit_dense(context, provider, id, refs); });
  }

  Fixture stale_fixture(1, 2);
  for (auto *source : {&stale_fixture.evolved, &stale_fixture.rhs})
    for (auto &validity : source->validity)
      validity.interior = true;
  auto stale_transaction = ScratchStateTransactionFactory::create_for_test(
      stale_fixture.copy_live_frame(), stale_fixture.metadata(),
      successful_schedule);
  auto stale_y = stale_transaction->capture_live_evolved();
  auto stale_f = stale_transaction->capture_live_rhs();
  std::vector<ScratchStateToken> stale_owned;
  stale_owned.reserve(constraints.size());
  std::vector<ScratchDenseSampleRef> stale_refs;
  stale_refs.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    const bool value = constraint.kind == DenseSampleKind::value;
    stale_owned.push_back(stale_transaction->clone_state(value ? stale_y
                                                               : stale_f));
    stale_refs.push_back(
        {constraint.theta,
         value ? ScratchDenseSampleKind::value
               : ScratchDenseSampleKind::raw_derivative,
         &stale_owned.back()});
  }
  {
    ScopedStepContext scope(context, preparer, stale_transaction.get());
    stale_transaction->commit_dense(context, provider, id, stale_refs);
  }
  stale_fixture.observed_epoch = 8;
  assert(stale_transaction->take_committed_dense() == nullptr);
  assert(stale_transaction->faulted());
  assert(stale_transaction->discarded());

  Fixture invalid_fixture(1, 2);
  for (auto *source : {&invalid_fixture.evolved, &invalid_fixture.rhs})
    for (auto &validity : source->validity)
      validity.interior = true;
  auto invalid_transaction =
      ScratchStateTransactionFactory::create_for_test(
          invalid_fixture.copy_live_frame(), invalid_fixture.metadata(),
          successful_schedule);
  std::vector<ScratchStateToken> invalid_owned;
  invalid_owned.reserve(constraints.size());
  std::vector<ScratchDenseSampleRef> invalid_refs;
  invalid_refs.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    const bool value = constraint.kind == DenseSampleKind::value;
    invalid_owned.push_back(value
                                ? invalid_transaction->capture_live_evolved()
                                : invalid_transaction->capture_live_rhs());
    invalid_refs.push_back(
        {constraint.theta,
         value ? ScratchDenseSampleKind::value
               : ScratchDenseSampleKind::raw_derivative,
         &invalid_owned.back()});
  }
  invalid_refs.front().kind = ScratchDenseSampleKind::raw_derivative;
  {
    ScopedStepContext scope(context, preparer, invalid_transaction.get());
    expect_throw<std::invalid_argument>([&] {
      invalid_transaction->commit_dense(context, provider, id, invalid_refs);
    });
  }
  assert(invalid_transaction->take_committed_dense() == nullptr);
  assert(invalid_transaction->faulted());
}

void test_dense_commit_rejects_unrelated_interval_identity() {
  const auto method = SubcyclingODEMethod::rk4;
  const auto cap = capability(method);
  DenseOutputProvider provider(cap);
  const auto &constraints =
      reference_dense_stencil(method).specification().constraints;
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;

  for (int variant = 0; variant < 6; ++variant) {
    Fixture fixture;
    for (auto *source : {&fixture.evolved, &fixture.rhs})
      for (auto &validity : source->validity)
        validity.interior = true;
    auto transaction = make_transaction(fixture);
    StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                        method};
    DenseIntervalId id{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                       method, cap.tableau_fingerprint};
    if (variant == 0) {
      id.begin_clock = step_clock_t(1);
      id.end_clock = step_clock_t(3, 2);
    } else if (variant == 1) {
      id.begin_time = 1.0;
      id.end_time = 1.5;
    } else if (variant == 2) {
      context.level = 1;
      id.level = 1;
    } else if (variant == 3) {
      id.method = SubcyclingODEMethod::dp87_order8;
    } else if (variant == 4) {
      id.tableau_fingerprint[0] ^= 0xff;
    }

    std::vector<ScratchStateToken> owned;
    owned.reserve(constraints.size());
    std::vector<ScratchDenseSampleRef> refs;
    refs.reserve(constraints.size());
    for (const auto &constraint : constraints) {
      const bool value = constraint.kind == DenseSampleKind::value;
      owned.push_back(value ? transaction->capture_live_evolved()
                            : transaction->capture_live_rhs());
      refs.push_back({constraint.theta,
                      value ? ScratchDenseSampleKind::value
                            : ScratchDenseSampleKind::raw_derivative,
                      &owned.back()});
    }
    if (variant == 5) {
      expect_throw<std::logic_error>(
          [&] { transaction->commit_dense(context, provider, id, refs); });
    } else {
      ScopedStepContext scope(context, preparer, transaction.get());
      expect_throw<std::exception>(
          [&] { transaction->commit_dense(context, provider, id, refs); });
    }
    assert(transaction->take_committed_dense() == nullptr);
    assert(transaction->faulted());
  }
}

void test_dense_commit_rejects_invalid_mapped_interiors() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  const auto method = SubcyclingODEMethod::rk4;
  const auto cap = capability(method);
  DenseOutputProvider provider(cap);
  const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                            method};
  DenseIntervalId id{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                     method, cap.tableau_fingerprint};
  const auto &constraints =
      reference_dense_stencil(method).specification().constraints;
  std::vector<ScratchStateToken> owned;
  owned.reserve(constraints.size());
  std::vector<ScratchDenseSampleRef> refs;
  refs.reserve(constraints.size());
  for (const auto &constraint : constraints) {
    const bool value = constraint.kind == DenseSampleKind::value;
    owned.push_back(value ? transaction->capture_live_evolved()
                          : transaction->capture_live_rhs());
    refs.push_back({constraint.theta,
                    value ? ScratchDenseSampleKind::value
                          : ScratchDenseSampleKind::raw_derivative,
                    &owned.back()});
  }
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    expect_throw<std::invalid_argument>(
        [&] { transaction->commit_dense(context, provider, id, refs); });
  }
  assert(transaction->take_committed_dense() == nullptr);
  assert(transaction->faulted());
}

void test_fault_discards_transaction_without_partial_publication() {
  Fixture fixture;
  auto transaction = make_transaction(
      fixture,
      [](const ScratchSchedulePhase phase, const StepContext &,
         const ScratchStageCoordinates &) -> ScratchScheduleExecutionReceipt {
        if (phase == ScratchSchedulePhase::rhs)
          throw std::runtime_error("injected RHS failure");
        return {phase, 1, 1};
      });
  auto state = transaction->capture_live_evolved();
  const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                            SubcyclingODEMethod::rk4};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    expect_throw<std::runtime_error>(
        [&] {
          transaction->evaluate_rhs(
              context, {StageKind::fractional, 2, 4,
                        step_clock_t(1, 2), 0.25});
        });
  }
  assert(transaction->faulted());
  assert(transaction->discarded());
  assert(transaction->rhs_evaluation_count() == 0);
  expect_throw<std::logic_error>(
      [&] { transaction->rollback_live_evolved(state); });
  expect_throw<std::logic_error>([&] { transaction->restore_state(state); });
}

void test_live_capture_drift_and_factory_envelope_rejection() {
  {
    Fixture fixture;
    auto transaction = make_transaction(fixture);
    fixture.observed_epoch = 8;
    expect_throw<std::runtime_error>(
        [&] { static_cast<void>(transaction->capture_live_evolved()); });
  }
  {
    Fixture fixture;
    auto metadata = fixture.metadata();
    const auto original_reader = metadata.live_entry_readers[0];
    int replacement_identity = 999;
    bool drift = false;
    metadata.live_entry_readers[0] = [original_reader,
                                      &replacement_identity, &drift] {
      auto snapshot = original_reader();
      if (drift)
        snapshot.storage_identity = &replacement_identity;
      return snapshot;
    };
    auto transaction = ScratchStateTransactionFactory::create_for_test(
        fixture.copy_live_frame(), std::move(metadata), successful_schedule);
    drift = true;
    expect_throw<std::runtime_error>(
        [&] { static_cast<void>(transaction->capture_live_rhs()); });
  }
  {
    Fixture fixture;
    auto metadata = fixture.metadata();
    metadata.patch_count = 2;
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
  {
    Fixture fixture;
    auto metadata = fixture.metadata();
    metadata.time_refinement_factor = 3;
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
  {
    Fixture fixture;
    fixture.rhs.multifab = make_source(11, 2.0, 31).multifab;
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(make_transaction(fixture));
    });
  }
  {
    Fixture fixture;
    auto second_evolved = make_source(13, 70.0);
    auto second_rhs = make_source(14, 7.0);
    fixture.entries.push_back(&second_evolved);
    fixture.entries.push_back(&second_rhs);
    auto metadata = fixture.metadata();
    metadata.group_pairs = {{13, 14}, {10, 11}};
    metadata.dependent_groups.push_back(14);
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
  {
    Fixture fixture;
    auto metadata = fixture.metadata();
    metadata.dependent_groups.push_back(10);
    expect_throw<std::invalid_argument>([&] {
      static_cast<void>(ScratchStateTransactionFactory::create_for_test(
          fixture.copy_live_frame(), std::move(metadata),
          successful_schedule));
    });
  }
}

void test_scoped_context_borrows_and_restores_transaction_pointer() {
  Fixture fixture;
  Fixture other_fixture;
  auto transaction = make_transaction(fixture);
  auto other_transaction = make_transaction(other_fixture);
  const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                            SubcyclingODEMethod::rk4};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  assert(current_scratch_state_transaction() == nullptr);
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    assert(current_scratch_state_transaction() == transaction.get());
    expect_throw<std::logic_error>(
        [&] {
          other_transaction->evaluate_rhs(
              context, {StageKind::fractional, 2, 4,
                        step_clock_t(1, 2), 0.25});
        });
  }
  assert(current_scratch_state_transaction() == nullptr);
  {
    ScopedStepContext legacy_scope(context, preparer);
    assert(current_scratch_state_transaction() == nullptr);
  }
}

void test_exact_stage_fraction_rejects_a_matching_range_but_wrong_time() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  const StepContext context{0, step_clock_t(0), step_clock_t(1, 2), 0.0, 0.5,
                            SubcyclingODEMethod::rk4};
  struct NoopPreparer final : StagePreparer {
    void prepare_stage(const StepContext &, const StagePoint &) override {}
  } preparer;
  {
    ScopedStepContext scope(context, preparer, transaction.get());
    expect_throw<std::invalid_argument>([&] {
      transaction->evaluate_rhs(
          context, {StageKind::fractional, 2, 4,
                    step_clock_t(1, 4), 0.25});
    });
  }
  assert(transaction->faulted());
  assert(transaction->discarded());
  assert(transaction->rhs_evaluation_count() == 0);
}

void test_dependent_groups_are_exposed_in_canonical_order() {
  Fixture fixture;
  auto metadata = fixture.metadata();
  metadata.dependent_groups = {12, 11};
  auto transaction = ScratchStateTransactionFactory::create_for_test(
      fixture.copy_live_frame(), std::move(metadata), successful_schedule);
  assert((transaction->dependent_groups() == std::vector<int>{11, 12}));
}

} // namespace

int main() {
  test_mid_primary_failure_rolls_live_tl0_back_exactly_without_schedule();
  test_successful_primary_endpoint_is_not_rolled_back();
  test_live_rollback_fails_closed_for_wrong_kind_or_live_drift();
  test_full_frame_round_trip_and_fixed_working_addresses();
  test_pair_order_rhs_canonicalization_and_mixed_linear_combination();
  test_restore_left_overlays_only_physical_rhs_interior_and_validity();
  test_hierarchy_epoch_drift_faults_before_state_mutation();
  test_clean_discard_remains_nonfaulting_when_dense_take_is_empty();
  test_zero_scale_and_zero_terms_do_not_read_nan_sources();
  test_owner_epoch_kind_and_schema_rejections();
  test_post_step_sidecars_and_schedule_coordinates();
  test_level_one_timefac_two_uses_exact_stage_clock();
  test_level_zero_remains_available_in_a_two_level_hierarchy();
  test_level_one_factory_and_context_mismatches_fail_closed();
  test_transaction_level_session_accepts_dense_and_nondense_endpoints();
  test_transaction_level_session_rolls_back_evolution_or_commit_failure();
  test_transaction_level_session_rolls_back_terminal_transaction_failures();
  test_armed_live_rollback_disarms_and_fails_closed_on_stale_epoch();
  test_uncommitted_transaction_level_session_destructor_rolls_back();
  test_dense_commit_is_atomic_and_uses_phase7_builder();
  test_dense_commit_rejects_unrelated_interval_identity();
  test_dense_commit_rejects_invalid_mapped_interiors();
  test_fault_discards_transaction_without_partial_publication();
  test_live_capture_drift_and_factory_envelope_rejection();
  test_scoped_context_borrows_and_restores_transaction_pointer();
  test_exact_stage_fraction_rejects_a_matching_range_but_wrong_time();
  test_dependent_groups_are_exposed_in_canonical_order();
  std::cout << "Scratch-state transaction tests passed\n";
  return 0;
}
