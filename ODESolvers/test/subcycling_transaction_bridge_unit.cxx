#include "subcycling_transaction_bridge.hxx"

#include "subcycling_dense_output.hxx"
#include "subcycling_dense_stencil.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"

#include <AMReX_MultiFab.H>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using namespace CarpetX;
using ODESolvers::TransactionPrimaryObserver;
using ODESolvers::TransactionStateBackend;

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
};

Source make_source(const int group, const double value) {
  static amrex::Arena arena(31);
  static amrex::TestFabFactory factory(37);
  Source source;
  source.storage_token = group + 100;
  source.key = ScratchGFKey{9, 0, 0, group};
  source.multifab = std::make_unique<amrex::MultiFab>(
      amrex::BoxArray(41), amrex::DistributionMapping(43), 1,
      amrex::IntVect(1, 1, 1), &arena, factory);
  source.validity = {{true, true, true}};
  for (std::size_t region = 0; region < 3; ++region)
    for (std::size_t cell = 0;
         cell < amrex::MultiFab::test_cells_per_region(); ++cell)
      source.multifab->test_set(static_cast<amrex::TestRegion>(region), 0,
                                cell, value);
  return source;
}

struct Fixture {
  std::int64_t epoch{9};
  int level_identity{101};
  Source evolved{make_source(20, 10.0)};
  Source rhs{make_source(21, 2.0)};
  Source dependent{make_source(22, 50.0)};
  std::vector<Source *> entries{&evolved, &rhs, &dependent};

  ScratchLevelFrame frame() {
    UncertifiedScratchLevelManifest manifest{9, 0, {}};
    std::vector<ScratchGFView> views;
    for (auto *source : entries) {
      manifest.entries.push_back({source->key, &source->storage_token});
      views.push_back({source->key, 0, &source->storage_token,
                       source->multifab.get(), &source->validity});
    }
    return ScratchLevelFrame::copy_tl0(manifest, views);
  }

  ScratchStateTransactionFactoryMetadata metadata() {
    ScratchStateTransactionFactoryMetadata result{
        9,
        1,
        1,
        3,
        step_clock_t(1, 2),
        0.5,
        1,
        false,
        {{20, 21}},
        {21, 22},
        [this] { return epoch; },
        {}};
    for (auto *source : entries) {
      result.live_entry_readers.push_back([this, source] {
        return ScratchLiveEntrySnapshot{
            source->key, &level_identity, source, &source->multifab,
            &source->storage_token, &source->validity, source->multifab.get(),
            source->validity, true,
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

std::unique_ptr<ScratchStateTransaction> make_transaction(Fixture &fixture) {
  return ScratchStateTransactionFactory::create_for_test(
      fixture.frame(), fixture.metadata(),
      [](const ScratchSchedulePhase phase, const StepContext &,
         const ScratchStageCoordinates &) {
        return ScratchScheduleExecutionReceipt{phase, 1, 1};
      });
}

void test_backend_clone_arithmetic_and_owner_rejection() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  auto y = TransactionStateBackend::from_token(
      *transaction, transaction->capture_live_evolved());
  auto f = TransactionStateBackend::from_token(
      *transaction, transaction->capture_live_rhs());

  auto result = y.clone();
  result.linear_combination(0.0, {{1.0, &y}, {0.5, &f}});
  assert(result.kind() == ScratchStateKind::evolved);
  transaction->restore_state(result.token());
  const auto &working =
      ScratchStateTransactionFactory::working_multifab_for_test(*transaction,
                                                                 0);
  assert(working.test_value(static_cast<amrex::TestRegion>(0), 0, 0) ==
         11.0);

  Fixture foreign_fixture;
  auto foreign_transaction = make_transaction(foreign_fixture);
  auto foreign = TransactionStateBackend::from_token(
      *foreign_transaction, foreign_transaction->capture_live_evolved());
  expect_throw<std::invalid_argument>(
      [&] { result.linear_combination(0.0, {{1.0, &foreign}}); });
}

void test_primary_observer_captures_exact_transaction_states() {
  Fixture fixture;
  auto transaction = make_transaction(fixture);
  TransactionPrimaryObserver observer(*transaction);
  const int ignored = 0;
  observer.initial_state(1.5, ignored);
  observer.initial_rhs(1.5, ignored);
  observer.accepted_endpoint(2.0, ignored);

  auto capture = observer.take_complete();
  assert(capture.left_state.kind() == ScratchStateKind::evolved);
  assert(capture.left_rhs.kind() == ScratchStateKind::raw_rhs);
  assert(capture.accepted_endpoint.kind() == ScratchStateKind::evolved);
  expect_throw<std::logic_error>([&] { observer.take_complete(); });
}

} // namespace

int main() {
  test_backend_clone_arithmetic_and_owner_rejection();
  test_primary_observer_captures_exact_transaction_states();
}
