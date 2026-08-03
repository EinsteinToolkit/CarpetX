#include "subcycling_scratch_state.hxx"

#include <AMReX_MultiFab.H>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using CarpetX::ScratchGFKey;
using CarpetX::ScratchGFManifestEntry;
using CarpetX::ScratchGFView;
using CarpetX::ScratchLevelFrame;
using CarpetX::ScratchValidity;
using CarpetX::UncertifiedScratchLevelManifest;

static_assert(!std::is_copy_constructible_v<ScratchLevelFrame>);
static_assert(!std::is_copy_assignable_v<ScratchLevelFrame>);
static_assert(std::is_nothrow_move_constructible_v<ScratchLevelFrame>);
static_assert(std::is_nothrow_move_assignable_v<ScratchLevelFrame>);

template <typename Function>
void expect_invalid_argument(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  assert(threw);
}

template <typename Function> void expect_out_of_range(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const std::out_of_range &) {
    threw = true;
  }
  assert(threw);
}

template <typename Function> void expect_runtime_error(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const std::runtime_error &) {
    threw = true;
  }
  assert(threw);
}

ScratchGFKey key(const std::int64_t epoch = 41, const int level = 2,
                 const int patch = 0, const int group = 3) {
  return ScratchGFKey{epoch, level, patch, group};
}

bool same_key(const ScratchGFKey &left, const ScratchGFKey &right) {
  return left.hierarchy_epoch == right.hierarchy_epoch &&
         left.level == right.level && left.patch == right.patch &&
         left.group_index == right.group_index;
}

struct Source {
  std::unique_ptr<amrex::MultiFab> multifab;
  std::vector<ScratchValidity> validity;
};

Source make_source(amrex::Arena *arena, const int factory_setting,
                   const int box_array, const int distribution,
                   const int components, const amrex::IntVect &grow,
                   const bool defined = true) {
  amrex::TestFabFactory factory(factory_setting);
  Source result;
  result.multifab = std::make_unique<amrex::MultiFab>(
      amrex::BoxArray(box_array), amrex::DistributionMapping(distribution),
      components, grow, arena, factory, defined);
  if (components > 0) {
    result.validity.reserve(static_cast<std::size_t>(components));
    for (int component = 0; component < components; ++component) {
      result.validity.push_back(ScratchValidity{
          component % 2 == 0, component % 3 == 0, component % 2 != 0});
    }
  }
  return result;
}

void populate(amrex::MultiFab &multifab, const double base) {
  const std::array<amrex::TestRegion, 3> regions{
      amrex::TestRegion::valid, amrex::TestRegion::outer,
      amrex::TestRegion::ghosts};
  for (std::size_t region = 0; region < regions.size(); ++region) {
    for (int component = 0; component < multifab.nComp(); ++component) {
      for (std::size_t cell = 0;
           cell < amrex::MultiFab::test_cells_per_region(); ++cell) {
        multifab.test_set(regions[region], component, cell,
                          base + 100.0 * static_cast<double>(region) +
                              10.0 * static_cast<double>(component) +
                              static_cast<double>(cell));
      }
    }
  }
}

ScratchGFView view(const ScratchGFKey &entry_key, const void *token,
                   const Source &source, const int time_level = 0) {
  return ScratchGFView{entry_key, time_level, token, source.multifab.get(),
                       &source.validity};
}

UncertifiedScratchLevelManifest manifest(
    const std::int64_t epoch, const int level,
    const std::vector<ScratchGFKey> &keys,
    const std::vector<const void *> &tokens) {
  assert(keys.size() == tokens.size());
  UncertifiedScratchLevelManifest result{epoch, level, {}};
  result.entries.reserve(keys.size());
  for (std::size_t index = 0; index < keys.size(); ++index)
    result.entries.push_back(ScratchGFManifestEntry{keys[index], tokens[index]});
  return result;
}

void expect_prevalidation_rejection(
    const UncertifiedScratchLevelManifest &invalid_manifest,
    const std::vector<ScratchGFView> &invalid_views) {
  amrex::MultiFab::reset_deep_copy_log();
  const int live_before = amrex::MultiFab::owned_live_count();
  expect_invalid_argument([&] {
    static_cast<void>(
        ScratchLevelFrame::copy_tl0(invalid_manifest, invalid_views));
  });
  assert(amrex::MultiFab::deep_copy_attempt_count() == 0);
  assert(amrex::MultiFab::owned_live_count() == live_before);
}

void test_ordered_multi_patch_multi_group_copy_and_exact_metadata() {
  amrex::Arena first_arena(11);
  amrex::Arena second_arena(13);
  auto first = make_source(&first_arena, 101, 17, 19, 2,
                           amrex::IntVect(1, 2, 3));
  auto second = make_source(&first_arena, 103, 23, 29, 1,
                            amrex::IntVect(0, 4, 2));
  auto third = make_source(&second_arena, 107, 31, 37, 3,
                           amrex::IntVect(3, 1, 5));
  populate(*first.multifab, 1000.0);
  populate(*second.multifab, 2000.0);
  populate(*third.multifab, 3000.0);

  int first_token = 1;
  int second_token = 2;
  int third_token = 3;
  const std::vector<ScratchGFKey> keys{key(41, 2, 0, 3), key(41, 2, 0, 8),
                                       key(41, 2, 1, 2)};
  const std::vector<const void *> tokens{&first_token, &second_token,
                                         &third_token};
  const auto supplied_manifest = manifest(41, 2, keys, tokens);
  const std::vector<ScratchGFView> views{
      view(keys[0], tokens[0], first), view(keys[1], tokens[1], second),
      view(keys[2], tokens[2], third)};

  amrex::MultiFab::reset_deep_copy_log();
  const int live_before = amrex::MultiFab::owned_live_count();
  {
    auto frame = ScratchLevelFrame::copy_tl0(supplied_manifest, views);
    assert(frame.hierarchy_epoch() == 41);
    assert(frame.level() == 2);
    assert(frame.entry_count() == 3);
    assert(amrex::MultiFab::deep_copy_attempt_count() == 3);
    assert(amrex::MultiFab::owned_live_count() == live_before + 3);

    const std::array<const Source *, 3> sources{&first, &second, &third};
    const std::array<amrex::TestRegion, 3> regions{
        amrex::TestRegion::valid, amrex::TestRegion::outer,
        amrex::TestRegion::ghosts};
    for (std::size_t index = 0; index < sources.size(); ++index) {
      const auto &source = *sources[index]->multifab;
      const auto &owned = frame.multifab(index);
      assert(same_key(frame.key(index), keys[index]));
      assert(&owned != &source);
      assert(owned.test_storage_address() != source.test_storage_address());
      assert(owned.boxArray() == source.boxArray());
      assert(owned.DistributionMap() == source.DistributionMap());
      assert(owned.nComp() == source.nComp());
      assert(owned.nGrowVect() == source.nGrowVect());
      assert(owned.arena() == source.arena());
      const auto *source_factory =
          dynamic_cast<const amrex::TestFabFactory *>(&source.Factory());
      const auto *owned_factory =
          dynamic_cast<const amrex::TestFabFactory *>(&owned.Factory());
      assert(source_factory != nullptr);
      assert(owned_factory != nullptr);
      assert(owned_factory != source_factory);
      assert(owned_factory->setting() == source_factory->setting());
      assert(frame.validity(index).size() == sources[index]->validity.size());
      for (std::size_t component = 0;
           component < sources[index]->validity.size(); ++component) {
        assert(frame.validity(index)[component].interior ==
               sources[index]->validity[component].interior);
        assert(frame.validity(index)[component].outer ==
               sources[index]->validity[component].outer);
        assert(frame.validity(index)[component].ghosts ==
               sources[index]->validity[component].ghosts);
      }
      for (const auto region : regions) {
        for (int component = 0; component < source.nComp(); ++component) {
          for (std::size_t cell = 0;
               cell < amrex::MultiFab::test_cells_per_region(); ++cell) {
            assert(owned.test_value(region, component, cell) ==
                   source.test_value(region, component, cell));
          }
        }
      }
    }

    frame.mutable_multifab(0).test_set(amrex::TestRegion::valid, 0, 0, -10.0);
    frame.mutable_multifab(0).test_set(amrex::TestRegion::outer, 0, 0, -20.0);
    frame.mutable_multifab(0).test_set(amrex::TestRegion::ghosts, 0, 0, -30.0);
    assert(first.multifab->test_value(amrex::TestRegion::valid, 0, 0) ==
           1000.0);
    assert(first.multifab->test_value(amrex::TestRegion::outer, 0, 0) ==
           1100.0);
    assert(first.multifab->test_value(amrex::TestRegion::ghosts, 0, 0) ==
           1200.0);

    frame.mutable_validity(0)[0] = ScratchValidity{false, false, true};
    assert(first.validity[0].interior);
    assert(first.validity[0].outer);
    assert(!first.validity[0].ghosts);
  }
  assert(amrex::MultiFab::owned_live_count() == live_before);
}

void test_manifest_errors_fail_before_first_deep_copy() {
  amrex::Arena arena(17);
  auto first = make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1, 2, 3));
  auto second = make_source(&arena, 2, 5, 7, 1, amrex::IntVect(2, 1, 4));
  int token_a = 1;
  int token_b = 2;
  const auto first_key = key(41, 2, 0, 3);
  const auto second_key = key(41, 2, 1, 4);
  const std::vector<ScratchGFView> valid_views{
      view(first_key, &token_a, first), view(second_key, &token_b, second)};

  expect_prevalidation_rejection(manifest(-1, 2, {first_key}, {&token_a}),
                                 {valid_views[0]});
  expect_prevalidation_rejection(manifest(41, -1, {first_key}, {&token_a}),
                                 {valid_views[0]});
  expect_prevalidation_rejection(manifest(41, 2, {}, {}), {});

  for (const auto &invalid_key :
       std::vector<ScratchGFKey>{key(-1, 2, 0, 3), key(41, -1, 0, 3),
                                 key(41, 2, -1, 3), key(41, 2, 0, -1),
                                 key(42, 2, 0, 3), key(41, 3, 0, 3)}) {
    expect_prevalidation_rejection(
        manifest(41, 2, {invalid_key}, {&token_a}),
        {ScratchGFView{invalid_key, 0, &token_a, first.multifab.get(),
                       &first.validity}});
  }

  const auto later_group = key(41, 2, 0, 8);
  const auto earlier_group = key(41, 2, 0, 3);
  expect_prevalidation_rejection(
      manifest(41, 2, {later_group, earlier_group}, {&token_a, &token_b}),
      {view(later_group, &token_a, first),
       view(earlier_group, &token_b, second)});
  expect_prevalidation_rejection(
      manifest(41, 2, {first_key, first_key}, {&token_a, &token_b}),
      {view(first_key, &token_a, first), view(first_key, &token_b, second)});
  expect_prevalidation_rejection(manifest(41, 2, {first_key}, {nullptr}),
                                 {view(first_key, nullptr, first)});
  expect_prevalidation_rejection(
      manifest(41, 2, {first_key, second_key}, {&token_a, &token_a}),
      {view(first_key, &token_a, first),
       view(second_key, &token_a, second)});
}

void test_view_errors_fail_before_first_deep_copy() {
  amrex::Arena arena(19);
  auto first = make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1, 2, 3));
  auto second = make_source(&arena, 2, 5, 7, 2, amrex::IntVect(2, 1, 4));
  auto third = make_source(&arena, 3, 11, 13, 1, amrex::IntVect(3, 2, 1));
  int token_a = 1;
  int token_b = 2;
  int token_c = 3;
  int wrong_token = 4;
  const auto first_key = key(41, 2, 0, 3);
  const auto second_key = key(41, 2, 0, 8);
  const auto third_key = key(41, 2, 1, 2);
  const auto two_manifest =
      manifest(41, 2, {first_key, second_key}, {&token_a, &token_b});
  const std::vector<ScratchGFView> valid_views{
      view(first_key, &token_a, first), view(second_key, &token_b, second)};

  expect_prevalidation_rejection(two_manifest, {valid_views[0]});
  expect_prevalidation_rejection(
      two_manifest,
      {valid_views[0], valid_views[1], view(third_key, &token_c, third)});

  auto invalid_views = valid_views;
  invalid_views[1].key = key(41, 2, 1, 8);
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].source_storage_identity = &wrong_token;
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].time_level = 1;
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].multifab = nullptr;
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].multifab = first.multifab.get();
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].validity = nullptr;
  expect_prevalidation_rejection(two_manifest, invalid_views);
  invalid_views = valid_views;
  invalid_views[1].validity = &first.validity;
  expect_prevalidation_rejection(two_manifest, invalid_views);

  auto undefined =
      make_source(&arena, 5, 17, 19, 1, amrex::IntVect(1), false);
  expect_prevalidation_rejection(manifest(41, 2, {first_key}, {&token_a}),
                                 {view(first_key, &token_a, undefined)});
  auto zero_components =
      make_source(&arena, 7, 23, 29, 0, amrex::IntVect(1));
  expect_prevalidation_rejection(manifest(41, 2, {first_key}, {&token_a}),
                                 {view(first_key, &token_a, zero_components)});
  auto bad_validity =
      make_source(&arena, 11, 31, 37, 2, amrex::IntVect(1));
  bad_validity.validity.pop_back();
  expect_prevalidation_rejection(manifest(41, 2, {first_key}, {&token_a}),
                                 {view(first_key, &token_a, bad_validity)});
}

void test_distinct_alias_wrappers_with_one_storage_token_are_rejected() {
  amrex::Arena arena(23);
  auto first_wrapper =
      make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1, 2, 3));
  auto second_wrapper =
      make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1, 2, 3));
  assert(first_wrapper.multifab.get() != second_wrapper.multifab.get());
  int canonical_storage_token = 1;
  const auto first_key = key(41, 2, 0, 3);
  const auto second_key = key(41, 2, 0, 8);
  expect_prevalidation_rejection(
      manifest(41, 2, {first_key, second_key},
               {&canonical_storage_token, &canonical_storage_token}),
      {view(first_key, &canonical_storage_token, first_wrapper),
       view(second_key, &canonical_storage_token, second_wrapper)});
}

void test_late_invalid_last_view_prevalidates_the_complete_batch() {
  amrex::Arena arena(29);
  auto first = make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1));
  auto second = make_source(&arena, 2, 5, 7, 1, amrex::IntVect(2));
  auto third = make_source(&arena, 3, 11, 13, 2, amrex::IntVect(3));
  third.validity.pop_back();
  int token_a = 1;
  int token_b = 2;
  int token_c = 3;
  const auto first_key = key(41, 2, 0, 3);
  const auto second_key = key(41, 2, 0, 8);
  const auto third_key = key(41, 2, 1, 2);
  expect_prevalidation_rejection(
      manifest(41, 2, {first_key, second_key, third_key},
               {&token_a, &token_b, &token_c}),
      {view(first_key, &token_a, first), view(second_key, &token_b, second),
       view(third_key, &token_c, third)});
}

void test_first_and_later_deep_copy_exceptions_roll_back_owned_entries() {
  amrex::Arena arena(31);
  auto first = make_source(&arena, 1, 2, 3, 1, amrex::IntVect(1));
  auto second = make_source(&arena, 2, 5, 7, 1, amrex::IntVect(2));
  auto third = make_source(&arena, 3, 11, 13, 1, amrex::IntVect(3));
  int token_a = 1;
  int token_b = 2;
  int token_c = 3;
  const auto first_key = key(41, 2, 0, 3);
  const auto second_key = key(41, 2, 0, 8);
  const auto third_key = key(41, 2, 1, 2);
  const auto supplied_manifest =
      manifest(41, 2, {first_key, second_key, third_key},
               {&token_a, &token_b, &token_c});
  const std::vector<ScratchGFView> views{
      view(first_key, &token_a, first), view(second_key, &token_b, second),
      view(third_key, &token_c, third)};
  const int live_before = amrex::MultiFab::owned_live_count();

  amrex::MultiFab::reset_deep_copy_log();
  amrex::MultiFab::fail_deep_copy_after(0);
  expect_runtime_error([&] {
    static_cast<void>(ScratchLevelFrame::copy_tl0(supplied_manifest, views));
  });
  assert(amrex::MultiFab::deep_copy_attempt_count() == 1);
  assert(amrex::MultiFab::owned_live_count() == live_before);

  amrex::MultiFab::reset_deep_copy_log();
  amrex::MultiFab::fail_deep_copy_after(1);
  expect_runtime_error([&] {
    static_cast<void>(ScratchLevelFrame::copy_tl0(supplied_manifest, views));
  });
  assert(amrex::MultiFab::deep_copy_attempt_count() == 2);
  assert(amrex::MultiFab::owned_live_count() == live_before);
  amrex::MultiFab::reset_deep_copy_log();
}

ScratchLevelFrame make_single_entry_frame(const std::int64_t epoch,
                                          const int level,
                                          const double value) {
  amrex::Arena arena(37);
  auto source = make_source(&arena, 41, 43, 47, 1, amrex::IntVect(1, 2, 3));
  populate(*source.multifab, value);
  int token = 1;
  const auto entry_key = key(epoch, level, 0, 3);
  return ScratchLevelFrame::copy_tl0(
      manifest(epoch, level, {entry_key}, {&token}),
      {view(entry_key, &token, source)});
}

void assert_moved_from(const ScratchLevelFrame &frame) {
  assert(frame.hierarchy_epoch() == -1);
  assert(frame.level() == -1);
  assert(frame.entry_count() == 0);
  expect_out_of_range([&] { static_cast<void>(frame.key(0)); });
  expect_out_of_range([&] { static_cast<void>(frame.multifab(0)); });
  expect_out_of_range([&] { static_cast<void>(frame.validity(0)); });
}

void test_move_construction_assignment_self_move_and_bounds() {
  const int live_before = amrex::MultiFab::owned_live_count();
  {
    auto source = make_single_entry_frame(51, 4, 7000.0);
    auto moved = ScratchLevelFrame(std::move(source));
    assert_moved_from(source);
    assert(moved.hierarchy_epoch() == 51);
    assert(moved.level() == 4);
    assert(moved.entry_count() == 1);
    assert(moved.multifab(0).test_value(amrex::TestRegion::valid, 0, 0) ==
           7000.0);

    const void *storage_before_self_move =
        moved.multifab(0).test_storage_address();
    auto *self_alias = &moved;
    moved = std::move(*self_alias);
    assert(moved.hierarchy_epoch() == 51);
    assert(moved.level() == 4);
    assert(moved.entry_count() == 1);
    assert(moved.multifab(0).test_storage_address() == storage_before_self_move);

    auto assignment_source = make_single_entry_frame(61, 5, 8000.0);
    moved = std::move(assignment_source);
    assert_moved_from(assignment_source);
    assert(moved.hierarchy_epoch() == 61);
    assert(moved.level() == 5);
    assert(moved.entry_count() == 1);
    assert(moved.multifab(0).test_value(amrex::TestRegion::ghosts, 0, 1) ==
           8201.0);

    expect_out_of_range([&] { static_cast<void>(moved.key(1)); });
    expect_out_of_range([&] { static_cast<void>(moved.multifab(1)); });
    expect_out_of_range(
        [&] { static_cast<void>(moved.mutable_multifab(1)); });
    expect_out_of_range([&] { static_cast<void>(moved.validity(1)); });
    expect_out_of_range(
        [&] { static_cast<void>(moved.mutable_validity(1)); });
  }
  assert(amrex::MultiFab::owned_live_count() == live_before);
}

} // namespace

int main() {
  test_ordered_multi_patch_multi_group_copy_and_exact_metadata();
  test_manifest_errors_fail_before_first_deep_copy();
  test_view_errors_fail_before_first_deep_copy();
  test_distinct_alias_wrappers_with_one_storage_token_are_rejected();
  test_late_invalid_last_view_prevalidates_the_complete_batch();
  test_first_and_later_deep_copy_exceptions_roll_back_owned_entries();
  test_move_construction_assignment_self_move_and_bounds();
  std::cout << "Scratch level-state ownership tests passed\n";
  return 0;
}
