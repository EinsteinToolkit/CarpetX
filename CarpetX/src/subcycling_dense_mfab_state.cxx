#include "subcycling_dense_mfab_state.hxx"

#include <cmath>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace CarpetX {
namespace {

bool non_negative(const DenseMFabKey &key) noexcept {
  return key.hierarchy_epoch >= 0 && key.patch >= 0 && key.level >= 0 &&
         key.group_index >= 0;
}

auto ordered_key(const DenseMFabKey &key) noexcept {
  return std::tie(key.hierarchy_epoch, key.patch, key.level, key.group_index);
}

bool same_key(const DenseMFabKey &left, const DenseMFabKey &right) noexcept {
  return ordered_key(left) == ordered_key(right);
}

bool same_layout(const amrex::MultiFab &left,
                 const amrex::MultiFab &right) noexcept {
  return left.boxArray() == right.boxArray() &&
         left.DistributionMap() == right.DistributionMap() &&
         left.nComp() == right.nComp() &&
         left.nGrowVect() == right.nGrowVect();
}

void validate_views(const std::vector<DenseMFabView> &views) {
  if (views.empty())
    throw std::invalid_argument("dense MultiFab state requires at least one view");

  const std::int64_t hierarchy_epoch = views.front().key.hierarchy_epoch;
  for (std::size_t index = 0; index < views.size(); ++index) {
    const auto &view = views[index];
    if (view.multifab == nullptr)
      throw std::invalid_argument("dense MultiFab view must not be null");
    if (!non_negative(view.key))
      throw std::invalid_argument("dense MultiFab keys must be non-negative");
    if (view.key.hierarchy_epoch != hierarchy_epoch)
      throw std::invalid_argument("dense MultiFab views must share one epoch");
    if (index != 0 &&
        !(ordered_key(views[index - 1].key) < ordered_key(view.key)))
      throw std::invalid_argument(
          "dense MultiFab keys must be strictly increasing");
  }
}

void validate_concrete_combination(
    const std::vector<double> &weights,
    const std::vector<const OwnedMultiFabDenseState *> &sources) {
  if (sources.empty())
    throw std::invalid_argument(
        "dense MultiFab combination requires at least one source");
  if (weights.size() != sources.size())
    throw std::invalid_argument(
        "dense MultiFab weight and source counts differ");
  for (const double weight : weights) {
    if (!std::isfinite(weight))
      throw std::invalid_argument(
          "dense MultiFab combination weights must be finite");
  }
  if (sources.front() == nullptr)
    throw std::invalid_argument("dense MultiFab source must not be null");
  for (const auto *source : sources) {
    if (source == nullptr)
      throw std::invalid_argument("dense MultiFab source must not be null");
    if (!sources.front()->compatible(*source))
      throw std::invalid_argument(
          "dense MultiFab combination sources are incompatible");
  }
}

} // namespace

OwnedMultiFabDenseState::OwnedMultiFabDenseState(std::vector<Entry> entries)
    : entries_(std::move(entries)) {}

std::unique_ptr<OwnedMultiFabDenseState>
OwnedMultiFabDenseState::copy_of(const std::vector<DenseMFabView> &views) {
  validate_views(views);

  std::vector<Entry> entries;
  entries.reserve(views.size());
  for (const auto &view : views) {
    auto owned = std::make_unique<amrex::MultiFab>(
        view.multifab->boxArray(), view.multifab->DistributionMap(),
        view.multifab->nComp(), view.multifab->nGrowVect());
    amrex::MultiFab::Copy(*owned, *view.multifab, 0, 0,
                          view.multifab->nComp(), 0);
    entries.push_back(Entry{view.key, std::move(owned)});
  }
  return std::unique_ptr<OwnedMultiFabDenseState>(
      new OwnedMultiFabDenseState(std::move(entries)));
}

std::unique_ptr<OwnedMultiFabDenseState>
OwnedMultiFabDenseState::allocate_like(
    const std::vector<DenseMFabView> &views) {
  validate_views(views);

  std::vector<Entry> entries;
  entries.reserve(views.size());
  for (const auto &view : views) {
    auto owned = std::make_unique<amrex::MultiFab>(
        view.multifab->boxArray(), view.multifab->DistributionMap(),
        view.multifab->nComp(), view.multifab->nGrowVect());
    entries.push_back(Entry{view.key, std::move(owned)});
  }
  return std::unique_ptr<OwnedMultiFabDenseState>(
      new OwnedMultiFabDenseState(std::move(entries)));
}

std::unique_ptr<OwnedMultiFabDenseState>
OwnedMultiFabDenseState::empty_like(const OwnedMultiFabDenseState &source) {
  std::vector<Entry> entries;
  entries.reserve(source.entries_.size());
  for (const auto &source_entry : source.entries_) {
    const auto &source_mfab = *source_entry.multifab;
    auto owned = std::make_unique<amrex::MultiFab>(
        source_mfab.boxArray(), source_mfab.DistributionMap(),
        source_mfab.nComp(), source_mfab.nGrowVect());
    entries.push_back(Entry{source_entry.key, std::move(owned)});
  }
  return std::unique_ptr<OwnedMultiFabDenseState>(
      new OwnedMultiFabDenseState(std::move(entries)));
}

std::unique_ptr<OwnedMultiFabDenseState>
OwnedMultiFabDenseState::linear_combination_of(
    const std::vector<double> &weights,
    const std::vector<const OwnedMultiFabDenseState *> &sources) {
  validate_concrete_combination(weights, sources);

  auto result = empty_like(*sources.front());
  std::vector<const DenseStateVector *> base_sources;
  base_sources.reserve(sources.size());
  for (const auto *source : sources)
    base_sources.push_back(source);
  result->linear_combination(weights, base_sources);
  return result;
}

std::size_t OwnedMultiFabDenseState::entry_count() const noexcept {
  return entries_.size();
}

const DenseMFabKey &
OwnedMultiFabDenseState::key(const std::size_t index) const {
  return entries_.at(index).key;
}

const amrex::MultiFab &
OwnedMultiFabDenseState::multifab(const std::size_t index) const {
  return *entries_.at(index).multifab;
}

bool OwnedMultiFabDenseState::compatible(
    const DenseStateVector &other) const noexcept {
  const auto *typed = dynamic_cast<const OwnedMultiFabDenseState *>(&other);
  if (typed == nullptr || entries_.size() != typed->entries_.size())
    return false;
  for (std::size_t index = 0; index < entries_.size(); ++index) {
    if (!same_key(entries_[index].key, typed->entries_[index].key) ||
        !same_layout(*entries_[index].multifab,
                     *typed->entries_[index].multifab))
      return false;
  }
  return true;
}

void OwnedMultiFabDenseState::copy_from(const DenseStateVector &other) {
  const auto *typed = dynamic_cast<const OwnedMultiFabDenseState *>(&other);
  if (typed == nullptr)
    throw std::invalid_argument("dense MultiFab copy source has wrong type");
  if (typed == this)
    throw std::invalid_argument("dense MultiFab copy rejects aliasing");
  if (!compatible(*typed))
    throw std::invalid_argument("dense MultiFab copy source is incompatible");

  for (std::size_t index = 0; index < entries_.size(); ++index) {
    auto &destination = *entries_[index].multifab;
    const auto &source = *typed->entries_[index].multifab;
    amrex::MultiFab::Copy(destination, source, 0, 0,
                          destination.nComp(), 0);
  }
}

void OwnedMultiFabDenseState::linear_combination(
    const std::vector<double> &weights,
    const std::vector<const DenseStateVector *> &sources) {
  if (sources.empty())
    throw std::invalid_argument(
        "dense MultiFab combination requires at least one source");
  if (weights.size() != sources.size())
    throw std::invalid_argument(
        "dense MultiFab weight and source counts differ");
  for (const double weight : weights) {
    if (!std::isfinite(weight))
      throw std::invalid_argument(
          "dense MultiFab combination weights must be finite");
  }

  std::vector<const OwnedMultiFabDenseState *> typed_sources;
  typed_sources.reserve(sources.size());
  for (const auto *source : sources) {
    if (source == nullptr)
      throw std::invalid_argument("dense MultiFab source must not be null");
    if (source == this)
      throw std::invalid_argument(
          "dense MultiFab combination rejects destination aliasing");
    const auto *typed =
        dynamic_cast<const OwnedMultiFabDenseState *>(source);
    if (typed == nullptr)
      throw std::invalid_argument("dense MultiFab source has wrong type");
    if (!compatible(*typed))
      throw std::invalid_argument("dense MultiFab source is incompatible");
    typed_sources.push_back(typed);
  }

  for (std::size_t entry = 0; entry < entries_.size(); ++entry) {
    auto &destination = *entries_[entry].multifab;
    destination.setVal(0.0, 0, destination.nComp(), 0);
    for (std::size_t source = 0; source < typed_sources.size(); ++source) {
      amrex::MultiFab::Saxpy(
          destination, weights[source],
          *typed_sources[source]->entries_[entry].multifab, 0, 0,
          destination.nComp(), 0);
    }
  }
}

} // namespace CarpetX
