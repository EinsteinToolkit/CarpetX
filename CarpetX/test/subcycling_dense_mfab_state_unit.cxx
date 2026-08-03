#include "subcycling_dense_mfab_state.hxx"

#include <AMReX_MultiFab.H>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using CarpetX::DenseMFabKey;
using CarpetX::DenseMFabView;
using CarpetX::DenseStateVector;
using CarpetX::OwnedMultiFabDenseState;

static_assert(!std::is_copy_constructible_v<OwnedMultiFabDenseState>);
static_assert(!std::is_copy_assignable_v<OwnedMultiFabDenseState>);
static_assert(!std::is_move_constructible_v<OwnedMultiFabDenseState>);
static_assert(!std::is_move_assignable_v<OwnedMultiFabDenseState>);
static_assert(std::is_final_v<OwnedMultiFabDenseState>);

template <typename Function> void expect_invalid_argument(Function &&function) {
  bool threw = false;
  try {
    function();
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  assert(threw);
}

class OtherDenseState final : public DenseStateVector {
public:
  bool compatible(const DenseStateVector &) const noexcept override {
    return false;
  }

  void copy_from(const DenseStateVector &) override {}

  void linear_combination(
      const std::vector<double> &,
      const std::vector<const DenseStateVector *> &) override {}
};

DenseMFabKey key(const std::int64_t epoch = 7, const int patch = 1,
                 const int level = 2, const int group = 3) {
  return DenseMFabKey{epoch, patch, level, group};
}

std::unique_ptr<amrex::MultiFab>
borrowed_mfab(const int box_array = 11, const int distribution = 13,
              const int components = 2, const int grow = 1,
              const double value = 1.0) {
  auto result = std::make_unique<amrex::MultiFab>(
      amrex::BoxArray(box_array), amrex::DistributionMapping(distribution),
      components, amrex::IntVect(grow));
  result->setVal(value, 0, components, 0);
  for (int component = 0; component < components; ++component) {
    result->test_set_valid(component, 1,
                           value + 0.25 * static_cast<double>(component + 1));
    result->test_set_ghost(component, value + 1000.0 + component);
  }
  return result;
}

std::unique_ptr<OwnedMultiFabDenseState>
state(const DenseMFabKey &entry_key = key(), const int box_array = 11,
      const int distribution = 13, const int components = 2,
      const int grow = 1, const double value = 1.0) {
  auto borrowed =
      borrowed_mfab(box_array, distribution, components, grow, value);
  return OwnedMultiFabDenseState::copy_of(
      std::vector<DenseMFabView>{{entry_key, borrowed.get()}});
}

double first_value(const OwnedMultiFabDenseState &state) {
  return state.multifab(0).test_valid(0, 0);
}

void assert_first_value(const OwnedMultiFabDenseState &state,
                        const double expected) {
  assert(std::abs(first_value(state) - expected) < 1.0e-12);
}

void test_copy_of_owns_values_and_copies_only_valid_interiors() {
  amrex::MultiFab::reset_operation_log();
  auto borrowed = borrowed_mfab(11, 13, 2, 2, 4.0);
  const int live_before_copy = amrex::MultiFab::live_count();

  auto owned = OwnedMultiFabDenseState::copy_of(
      std::vector<DenseMFabView>{{key(), borrowed.get()}});

  assert(amrex::MultiFab::live_count() == live_before_copy + 1);
  assert(owned->entry_count() == 1);
  assert(owned->key(0).hierarchy_epoch == 7);
  assert(owned->key(0).patch == 1);
  assert(owned->key(0).level == 2);
  assert(owned->key(0).group_index == 3);
  assert(owned->multifab(0).test_valid(0, 0) == 4.0);
  assert(owned->multifab(0).test_valid(1, 1) == 4.5);
  assert(owned->multifab(0).test_ghost(0) ==
         amrex::MultiFab::uninitialized_value());
  assert(amrex::MultiFab::all_operations_interior_only());

  borrowed.reset();
  assert(amrex::MultiFab::live_count() == live_before_copy);
  assert(owned->multifab(0).test_valid(0, 0) == 4.0);
  assert(owned->multifab(0).test_valid(1, 1) == 4.5);

  const auto empty = OwnedMultiFabDenseState::empty_like(*owned);
  assert(empty->compatible(*owned));
  assert(owned->compatible(*empty));
}

void test_allocate_like_preserves_keys_and_layout_without_copying_values() {
  auto first = borrowed_mfab(11, 13, 2, 1, 4.0);
  auto second = borrowed_mfab(17, 19, 3, 2, 9.0);
  const std::vector<DenseMFabView> views{
      {key(7, 1, 2, 3), first.get()}, {key(7, 1, 2, 4), second.get()}};
  const auto copied = OwnedMultiFabDenseState::copy_of(views);

  amrex::MultiFab::reset_construction_count();
  const auto allocated = OwnedMultiFabDenseState::allocate_like(views);

  assert(amrex::MultiFab::construction_count() == 2);
  assert(allocated->entry_count() == 2);
  assert(allocated->key(0).hierarchy_epoch == 7);
  assert(allocated->key(0).patch == 1);
  assert(allocated->key(0).level == 2);
  assert(allocated->key(0).group_index == 3);
  assert(allocated->key(1).group_index == 4);
  assert(allocated->multifab(0).boxArray() == first->boxArray());
  assert(allocated->multifab(0).DistributionMap() ==
         first->DistributionMap());
  assert(allocated->multifab(0).nComp() == first->nComp());
  assert(allocated->multifab(0).nGrowVect() == first->nGrowVect());
  assert(allocated->multifab(1).boxArray() == second->boxArray());
  assert(allocated->multifab(1).DistributionMap() ==
         second->DistributionMap());
  assert(allocated->multifab(1).nComp() == second->nComp());
  assert(allocated->multifab(1).nGrowVect() == second->nGrowVect());
  assert(allocated->compatible(*copied));
  assert(copied->compatible(*allocated));
  assert(allocated->multifab(0).test_valid(0, 0) ==
         amrex::MultiFab::uninitialized_value());
  assert(allocated->multifab(0).test_ghost(0) ==
         amrex::MultiFab::uninitialized_value());
  assert(allocated->multifab(1).test_valid(0, 0) ==
         amrex::MultiFab::uninitialized_value());

  amrex::MultiFab::reset_construction_count();
  expect_invalid_argument(
      [] { OwnedMultiFabDenseState::allocate_like({}); });
  assert(amrex::MultiFab::construction_count() == 0);
}

void test_copy_of_rejects_invalid_views_before_allocation() {
  auto first = borrowed_mfab();
  auto second = borrowed_mfab();

  amrex::MultiFab::reset_construction_count();
  expect_invalid_argument(
      [] { OwnedMultiFabDenseState::copy_of({}); });
  assert(amrex::MultiFab::construction_count() == 0);

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(
        std::vector<DenseMFabView>{{key(), nullptr}});
  });
  assert(amrex::MultiFab::construction_count() == 0);

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(std::vector<DenseMFabView>{
        {key(7, 1, 2, 3), first.get()}, {key(7, 1, 2, 4), nullptr}});
  });
  assert(amrex::MultiFab::construction_count() == 0);

  for (const auto &invalid_key :
       std::vector<DenseMFabKey>{key(-1, 1, 2, 3), key(7, -1, 2, 3),
                                 key(7, 1, -1, 3), key(7, 1, 2, -1)}) {
    expect_invalid_argument([&] {
      OwnedMultiFabDenseState::copy_of(
          std::vector<DenseMFabView>{{invalid_key, first.get()}});
    });
    assert(amrex::MultiFab::construction_count() == 0);
  }

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(std::vector<DenseMFabView>{
        {key(7, 1, 2, 3), first.get()}, {key(7, 1, 2, -1), second.get()}});
  });
  assert(amrex::MultiFab::construction_count() == 0);

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(std::vector<DenseMFabView>{
        {key(7, 1, 2, 3), first.get()}, {key(7, 1, 2, 3), second.get()}});
  });
  assert(amrex::MultiFab::construction_count() == 0);

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(std::vector<DenseMFabView>{
        {key(7, 1, 2, 4), first.get()}, {key(7, 1, 2, 3), second.get()}});
  });
  assert(amrex::MultiFab::construction_count() == 0);

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::copy_of(std::vector<DenseMFabView>{
        {key(7, 1, 2, 3), first.get()}, {key(8, 1, 2, 3), second.get()}});
  });
  assert(amrex::MultiFab::construction_count() == 0);
}

void test_exact_key_and_layout_compatibility() {
  const auto reference = state();
  const auto identical = state();
  assert(reference->compatible(*identical));
  assert(identical->compatible(*reference));

  const std::vector<std::unique_ptr<OwnedMultiFabDenseState>> mismatches = [&] {
    std::vector<std::unique_ptr<OwnedMultiFabDenseState>> result;
    result.push_back(state(key(8, 1, 2, 3)));
    result.push_back(state(key(7, 2, 2, 3)));
    result.push_back(state(key(7, 1, 3, 3)));
    result.push_back(state(key(7, 1, 2, 4)));
    result.push_back(state(key(), 17, 13, 2, 1));
    result.push_back(state(key(), 11, 19, 2, 1));
    result.push_back(state(key(), 11, 13, 3, 1));
    result.push_back(state(key(), 11, 13, 2, 2));
    return result;
  }();

  for (const auto &mismatch : mismatches) {
    assert(!reference->compatible(*mismatch));
    assert(!mismatch->compatible(*reference));
  }

  OtherDenseState other_type;
  assert(!reference->compatible(other_type));
}

void test_copy_through_dense_state_interface() {
  auto destination = state(key(), 11, 13, 2, 1, -3.0);
  const auto source = state(key(), 11, 13, 2, 1, 9.0);
  DenseStateVector &destination_base = *destination;
  const DenseStateVector &source_base = *source;

  amrex::MultiFab::reset_operation_log();
  destination_base.copy_from(source_base);

  assert_first_value(*destination, 9.0);
  assert(destination->multifab(0).test_valid(1, 1) == 9.5);
  assert(amrex::MultiFab::all_operations_interior_only());
}

void exercise_base_combination(const std::size_t count,
                               const double expected) {
  std::vector<std::unique_ptr<OwnedMultiFabDenseState>> owned_sources;
  std::vector<const DenseStateVector *> base_sources;
  std::vector<double> weights;
  owned_sources.reserve(count);
  base_sources.reserve(count);
  weights.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    owned_sources.push_back(state(key(), 11, 13, 2, 1,
                                  static_cast<double>(i + 1)));
    base_sources.push_back(owned_sources.back().get());
    weights.push_back(static_cast<double>(i + 1));
  }

  auto destination = state(key(), 11, 13, 2, 1, -99.0);
  DenseStateVector &destination_base = *destination;
  amrex::MultiFab::reset_operation_log();
  destination_base.linear_combination(weights, base_sources);
  assert_first_value(*destination, expected);
  assert(amrex::MultiFab::all_operations_interior_only());
}

void test_five_eight_and_nine_source_combinations_through_base_interface() {
  exercise_base_combination(5, 55.0);
  exercise_base_combination(8, 204.0);
  exercise_base_combination(9, 285.0);

  std::vector<std::unique_ptr<OwnedMultiFabDenseState>> owned_sources;
  std::vector<const OwnedMultiFabDenseState *> concrete_sources;
  std::vector<double> weights{1.0, 2.0, 3.0, 4.0, 5.0};
  for (std::size_t i = 0; i < weights.size(); ++i) {
    owned_sources.push_back(state(key(), 11, 13, 2, 1,
                                  static_cast<double>(i + 1)));
    concrete_sources.push_back(owned_sources.back().get());
  }
  const auto result = OwnedMultiFabDenseState::linear_combination_of(
      weights, concrete_sources);
  assert_first_value(*result, 55.0);
}

template <typename Function>
void expect_failure_before_mutation(OwnedMultiFabDenseState &destination,
                                    Function &&function) {
  const double before = first_value(destination);
  expect_invalid_argument(std::forward<Function>(function));
  assert(first_value(destination) == before);
}

void test_copy_and_combination_failures_prevalidate_before_mutation() {
  auto destination = state(key(), 11, 13, 2, 1, 314.0);
  const auto compatible_source = state(key(), 11, 13, 2, 1, 2.0);
  const auto incompatible_source = state(key(7, 1, 2, 4), 11, 13, 2, 1, 3.0);
  OtherDenseState wrong_type;

  expect_failure_before_mutation(*destination,
                                 [&] { destination->copy_from(*destination); });
  expect_failure_before_mutation(
      *destination, [&] { destination->copy_from(wrong_type); });
  expect_failure_before_mutation(
      *destination, [&] { destination->copy_from(*incompatible_source); });

  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination({}, {});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0}, {compatible_source.get(), compatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination({1.0}, {nullptr});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0, 1.0}, {compatible_source.get(), nullptr});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination({1.0}, {&wrong_type});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0, 1.0}, {compatible_source.get(), &wrong_type});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination({1.0}, {incompatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0, 1.0}, {compatible_source.get(), incompatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {std::numeric_limits<double>::quiet_NaN()},
        {compatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0, std::numeric_limits<double>::quiet_NaN()},
        {compatible_source.get(), compatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {std::numeric_limits<double>::infinity()},
        {compatible_source.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination({1.0}, {destination.get()});
  });
  expect_failure_before_mutation(*destination, [&] {
    destination->linear_combination(
        {1.0, 1.0}, {compatible_source.get(), destination.get()});
  });

  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::linear_combination_of(
        {}, std::vector<const OwnedMultiFabDenseState *>{});
  });
  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::linear_combination_of(
        {1.0}, std::vector<const OwnedMultiFabDenseState *>{nullptr});
  });
  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::linear_combination_of(
        {std::numeric_limits<double>::quiet_NaN()},
        std::vector<const OwnedMultiFabDenseState *>{compatible_source.get()});
  });
  expect_invalid_argument([&] {
    OwnedMultiFabDenseState::linear_combination_of(
        {1.0},
        std::vector<const OwnedMultiFabDenseState *>{compatible_source.get(),
                                                      compatible_source.get()});
  });
}

} // namespace

int main() {
  test_copy_of_owns_values_and_copies_only_valid_interiors();
  test_allocate_like_preserves_keys_and_layout_without_copying_values();
  test_copy_of_rejects_invalid_views_before_allocation();
  test_exact_key_and_layout_compatibility();
  test_copy_through_dense_state_interface();
  test_five_eight_and_nine_source_combinations_through_base_interface();
  test_copy_and_combination_failures_prevalidate_before_mutation();
  std::cout << "Owned MultiFab dense-state tests passed\n";
  return 0;
}
