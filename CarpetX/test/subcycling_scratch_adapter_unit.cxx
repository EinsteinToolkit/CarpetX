#define private public
#include "subcycling_scratch_adapter.hxx"
#undef private
#include "subcycling_scratch_adapter_internal.hxx"

#include <AMReX_MultiFab.H>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using namespace CarpetX;
using namespace CarpetX::subcycling_detail;

struct ScriptedCollective final : CollectiveOps {
  struct Override {
    ScratchCollectivePhase phase;
    std::size_t ordinal;
    char operation;
    std::int64_t value;
  };
  std::vector<std::string> events;
  std::vector<Override> overrides;

  std::int64_t reduce_min(const ScratchCollectivePhase phase,
                          const std::size_t ordinal,
                          const std::int64_t local) override {
    events.push_back("min:" + std::to_string(static_cast<int>(phase)) + ":" +
                     std::to_string(ordinal));
    return replacement(phase, ordinal, 'm', local);
  }
  std::int64_t reduce_max(const ScratchCollectivePhase phase,
                          const std::size_t ordinal,
                          const std::int64_t local) override {
    events.push_back("max:" + std::to_string(static_cast<int>(phase)) + ":" +
                     std::to_string(ordinal));
    return replacement(phase, ordinal, 'M', local);
  }
  bool reduce_and(const ScratchCollectivePhase phase,
                  const std::size_t ordinal, const bool local) override {
    events.push_back("and:" + std::to_string(static_cast<int>(phase)) + ":" +
                     std::to_string(ordinal));
    return replacement(phase, ordinal, 'a', local ? 1 : 0) != 0;
  }

private:
  std::int64_t replacement(const ScratchCollectivePhase phase,
                           const std::size_t ordinal, const char operation,
                           const std::int64_t local) const {
    for (const auto &entry : overrides)
      if (entry.phase == phase && entry.ordinal == ordinal &&
          entry.operation == operation)
        return entry.value;
    return local;
  }
};

struct Source {
  amrex::MultiFab multifab;
  std::vector<ScratchValidity> validity;
  int token;
  Source(const int components, const bool source_ok, const int id)
      : multifab(components, true, source_ok, id), token(id) {
    for (int component = 0; component < components; ++component)
      validity.push_back(
          ScratchValidity{component % 2 == 0, component % 3 == 0,
                          component % 2 != 0});
  }
};

LocalScratchEntry make_entry(const int group, Source &source,
                             std::vector<std::int64_t> schema = {}) {
  if (schema.empty())
    schema = {77, 2, 0, group, source.multifab.nComp(),
              static_cast<std::int64_t>(source.validity.size()), group + 10};
  const auto schema_length = static_cast<std::int64_t>(schema.size());
  return LocalScratchEntry{ScratchGFKey{77, 2, 0, group}, &source.token,
                           &source.multifab, source.validity,
                           std::move(schema), schema_length};
}

LocalScratchBatch make_batch(std::vector<LocalScratchEntry> entries) {
  const auto entry_count = static_cast<std::int64_t>(entries.size());
  return LocalScratchBatch{ScratchAdapterFailure::none, 77, 0, 2,
                           std::move(entries), entry_count};
}

ScratchLevelFrame run(LocalScratchBatch &batch, ScriptedCollective &collective,
                      const std::int64_t current_epoch = 77,
                      const bool fail_copy = false) {
  return run_certification_transaction(
      batch, collective, [current_epoch] { return current_epoch; },
      [fail_copy](const UncertifiedScratchLevelManifest &manifest,
                  const std::vector<ScratchGFView> &views) {
        if (fail_copy)
          throw std::runtime_error("requested copy failure");
        return ScratchLevelFrame::copy_tl0(manifest, views);
      });
}

template <typename F> void expect_collective_rejection(F &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const ScratchAdapterCollectiveError &) {
    rejected = true;
  }
  assert(rejected);
}

void test_multiple_groups_and_exact_validity() {
  Source first(2, true, 11), second(3, true, 13);
  auto batch = make_batch(
      {make_entry(2, first), make_entry(7, second)});
  ScriptedCollective collective;
  auto frame = run(batch, collective);
  assert(frame.entry_count() == 2);
  assert(frame.key(0).group_index == 2);
  assert(frame.key(1).group_index == 7);
  assert(frame.validity(0).size() == first.validity.size());
  assert(frame.validity(1)[2].interior == second.validity[2].interior);
  assert(first.multifab.ok_call_count() == 1);
  assert(second.multifab.ok_call_count() == 1);
}

void test_preflight_failure_still_reaches_collective() {
  Source source(1, true, 17);
  auto batch = make_batch({make_entry(1, source)});
  batch.preflight_status = ScratchAdapterFailure::invalid_tl0;
  ScriptedCollective collective;
  expect_collective_rejection([&] { static_cast<void>(run(batch, collective)); });
  assert(collective.events.size() == 4);
}

void test_prevalidated_size_failures_reach_fixed_collective_block() {
  Source source(1, true, 18);
  {
    auto batch = make_batch({make_entry(1, source)});
    batch.prevalidated_entry_count = -1;
    ScriptedCollective collective;
    expect_collective_rejection(
        [&] { static_cast<void>(run(batch, collective)); });
    assert(collective.events.size() == 4);
  }
  {
    auto batch = make_batch({make_entry(1, source)});
    batch.entries[0].prevalidated_schema_length = -1;
    ScriptedCollective collective;
    expect_collective_rejection(
        [&] { static_cast<void>(run(batch, collective)); });
    assert(collective.events.size() == 4);
  }
  {
    auto batch = make_batch({make_entry(1, source)});
    batch.entries[0].prevalidated_schema_length += 1;
    ScriptedCollective collective;
    expect_collective_rejection(
        [&] { static_cast<void>(run(batch, collective)); });
    assert(collective.events.size() == 4);
  }
}

void test_count_and_schema_mismatch() {
  Source source(1, true, 19);
  auto batch = make_batch({make_entry(1, source)});
  ScriptedCollective count;
  count.overrides.push_back(
      {ScratchCollectivePhase::entry_count, 0, 'M', 2});
  expect_collective_rejection([&] { static_cast<void>(run(batch, count)); });

  ScriptedCollective schema;
  schema.overrides.push_back(
      {ScratchCollectivePhase::schema_fields, 3, 'M', 999});
  expect_collective_rejection([&] { static_cast<void>(run(batch, schema)); });
}

void test_box_processor_and_validity_mismatch() {
  Source source(2, true, 23);
  auto batch = make_batch({make_entry(1, source,
      {77, 0, 2, 1, 2, 2, 1, 0, 0, 7, 7, 7, 9, 9, 9,
       2, 0, 1, 2, 8, 8, 8, 0, 1, 1, 0, 1})});
  for (const std::size_t ordinal : {9U, 20U, 25U}) {
    ScriptedCollective collective;
    collective.overrides.push_back(
        {ScratchCollectivePhase::schema_fields, ordinal, 'M', 12345});
    expect_collective_rejection(
        [&] { static_cast<void>(run(batch, collective)); });
  }
}

void test_duplicate_omitted_and_extra_groups() {
  Source first(1, true, 29), second(1, true, 31);
  for (const auto failure : {ScratchAdapterFailure::duplicate_group,
                             ScratchAdapterFailure::omitted_group,
                             ScratchAdapterFailure::extra_group}) {
    auto batch = make_batch({make_entry(1, first), make_entry(2, second)});
    batch.preflight_status = failure;
    ScriptedCollective collective;
    expect_collective_rejection(
        [&] { static_cast<void>(run(batch, collective)); });
  }
}

void test_all_sources_checked_after_false() {
  Source first(1, false, 37), second(1, true, 41), third(1, false, 43);
  auto batch = make_batch(
      {make_entry(1, first), make_entry(2, second), make_entry(3, third)});
  ScriptedCollective collective;
  expect_collective_rejection([&] { static_cast<void>(run(batch, collective)); });
  assert(first.multifab.ok_call_count() == 1);
  assert(second.multifab.ok_call_count() == 1);
  assert(third.multifab.ok_call_count() == 1);
}

void test_epoch_change() {
  Source source(1, true, 47);
  auto batch = make_batch({make_entry(1, source)});
  ScriptedCollective collective;
  expect_collective_rejection(
      [&] { static_cast<void>(run(batch, collective, 78)); });
}

void test_copy_failure_discards_local_success() {
  Source source(1, true, 53);
  auto batch = make_batch({make_entry(1, source)});
  ScriptedCollective local_failure;
  expect_collective_rejection(
      [&] { static_cast<void>(run(batch, local_failure, 77, true)); });

  ScriptedCollective peer_failure;
  peer_failure.overrides.push_back(
      {ScratchCollectivePhase::copy_success, 0, 'a', 0});
  const int before = amrex::MultiFab::owned_live_count();
  expect_collective_rejection(
      [&] { static_cast<void>(run(batch, peer_failure)); });
  assert(amrex::MultiFab::owned_live_count() == before);
}

void test_stable_record_view_lifetimes() {
  Source first(1, true, 59), second(1, true, 61);
  auto batch = make_batch({make_entry(1, first), make_entry(2, second)});
  const auto *first_validity = batch.entries[0].validity.data();
  ScriptedCollective collective;
  auto frame = run(batch, collective);
  assert(batch.entries[0].validity.data() == first_validity);
  assert(frame.validity(0)[0].interior == first.validity[0].interior);
}

void test_move_contract() {
  static_assert(!std::is_copy_constructible_v<CertifiedScratchLevelFrame>);
  static_assert(!std::is_copy_assignable_v<CertifiedScratchLevelFrame>);
  static_assert(
      std::is_nothrow_move_constructible_v<CertifiedScratchLevelFrame>);
  static_assert(std::is_nothrow_move_assignable_v<CertifiedScratchLevelFrame>);
  Source source(1, true, 67);
  auto batch = make_batch({make_entry(1, source)});
  ScriptedCollective collective;
  auto certified = CertifiedScratchLevelFrame(0, run(batch, collective));
  auto moved = std::move(certified);
  assert(certified.patch() == -1);
  assert(certified.frame().hierarchy_epoch() == -1);
  assert(moved.patch() == 0);
  Source replacement_source(1, true, 69);
  auto replacement_batch = make_batch({make_entry(2, replacement_source)});
  ScriptedCollective replacement_collective;
  auto assigned = CertifiedScratchLevelFrame(
      0, run(replacement_batch, replacement_collective));
  assigned = std::move(moved);
  assert(moved.patch() == -1);
  auto *self = &assigned;
  assigned = std::move(*self);
  assert(assigned.patch() == 0);
}

void test_pointer_identities_never_enter_schema() {
  // Contract marker: pointer identities never enter schema.
  Source source(1, true, 71);
  auto entry = make_entry(1, source);
  for (const auto field : entry.schema_fields)
    assert(field != reinterpret_cast<std::intptr_t>(entry.source_storage_identity));
}

} // namespace

int main() {
  test_multiple_groups_and_exact_validity();
  test_preflight_failure_still_reaches_collective();
  test_prevalidated_size_failures_reach_fixed_collective_block();
  test_count_and_schema_mismatch();
  test_box_processor_and_validity_mismatch();
  test_duplicate_omitted_and_extra_groups();
  test_all_sources_checked_after_false();
  test_epoch_change();
  test_copy_failure_discards_local_success();
  test_stable_record_view_lifetimes();
  test_move_contract();
  test_pointer_identities_never_enter_schema();
  std::cout << "Phase 8B1 certified TL0 transaction tests passed\n";
}
