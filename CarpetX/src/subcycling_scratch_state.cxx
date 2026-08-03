#include "subcycling_scratch_state.hxx"

#include <AMReX_MultiFab.H>

#include <memory>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <utility>

namespace CarpetX {
namespace {

bool non_negative(const ScratchGFKey &key) noexcept {
  return key.hierarchy_epoch >= 0 && key.level >= 0 && key.patch >= 0 &&
         key.group_index >= 0;
}

auto patch_group(const ScratchGFKey &key) noexcept {
  return std::tie(key.patch, key.group_index);
}

bool same_key(const ScratchGFKey &left, const ScratchGFKey &right) noexcept {
  return left.hierarchy_epoch == right.hierarchy_epoch &&
         left.level == right.level && left.patch == right.patch &&
         left.group_index == right.group_index;
}

#ifndef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
bool same_transaction_layout(const amrex::MultiFab &left,
                             const amrex::MultiFab &right) noexcept {
  return left.boxArray() == right.boxArray() &&
         left.DistributionMap() == right.DistributionMap() &&
         left.nComp() == right.nComp() &&
         left.nGrowVect() == right.nGrowVect();
}
#endif

void validate_manifest(const UncertifiedScratchLevelManifest &manifest) {
  if (manifest.hierarchy_epoch < 0)
    throw std::invalid_argument("scratch manifest epoch must be non-negative");
  if (manifest.level < 0)
    throw std::invalid_argument("scratch manifest level must be non-negative");
  if (manifest.entries.empty())
    throw std::invalid_argument("scratch manifest must not be empty");

  std::unordered_set<const void *> unique_storage;
  unique_storage.reserve(manifest.entries.size());
  for (std::size_t index = 0; index < manifest.entries.size(); ++index) {
    const auto &entry = manifest.entries[index];
    if (!non_negative(entry.key))
      throw std::invalid_argument("scratch manifest keys must be non-negative");
    if (entry.key.hierarchy_epoch != manifest.hierarchy_epoch ||
        entry.key.level != manifest.level)
      throw std::invalid_argument(
          "scratch manifest key epoch and level must match the manifest");
    if (index != 0 &&
        !(patch_group(manifest.entries[index - 1].key) <
          patch_group(entry.key)))
      throw std::invalid_argument(
          "scratch manifest keys must be strictly ordered by patch and group");
    if (entry.source_storage_identity == nullptr)
      throw std::invalid_argument(
          "scratch manifest storage identity must not be null");
    if (!unique_storage.insert(entry.source_storage_identity).second)
      throw std::invalid_argument(
          "scratch manifest storage identities must be unique");
  }
}

void validate_views(const UncertifiedScratchLevelManifest &manifest,
                    const std::vector<ScratchGFView> &views) {
  if (views.size() != manifest.entries.size())
    throw std::invalid_argument(
        "scratch view count must exactly match the manifest");

  std::unordered_set<const amrex::MultiFab *> unique_multifabs;
  std::unordered_set<const std::vector<ScratchValidity> *> unique_validity;
  unique_multifabs.reserve(views.size());
  unique_validity.reserve(views.size());
  for (std::size_t index = 0; index < views.size(); ++index) {
    const auto &manifest_entry = manifest.entries[index];
    const auto &view = views[index];
    if (!same_key(view.key, manifest_entry.key))
      throw std::invalid_argument(
          "scratch view key must match the ordered manifest key");
    if (view.source_storage_identity !=
        manifest_entry.source_storage_identity)
      throw std::invalid_argument(
          "scratch view storage identity must match the ordered manifest");
    if (view.time_level != 0)
      throw std::invalid_argument("scratch view must select time level zero");
    if (view.multifab == nullptr)
      throw std::invalid_argument("scratch view MultiFab must not be null");
    if (!unique_multifabs.insert(view.multifab).second)
      throw std::invalid_argument("scratch view MultiFabs must be unique");
    if (view.validity == nullptr)
      throw std::invalid_argument("scratch view validity must not be null");
    if (!unique_validity.insert(view.validity).second)
      throw std::invalid_argument("scratch view validity must be unique");
    if (!view.multifab->isDefined())
      throw std::invalid_argument("scratch source MultiFab must be defined");
    if (view.multifab->nComp() <= 0)
      throw std::invalid_argument(
          "scratch source MultiFab must have positive components");
    if (view.validity->size() !=
        static_cast<std::size_t>(view.multifab->nComp()))
      throw std::invalid_argument(
          "scratch validity size must equal the MultiFab component count");
  }
}

} // namespace

struct ScratchLevelFrame::Entry {
  ScratchGFKey key;
  std::unique_ptr<amrex::MultiFab> multifab;
  std::vector<ScratchValidity> validity;
};

ScratchLevelFrame ScratchLevelFrame::copy_tl0(
    const UncertifiedScratchLevelManifest &manifest,
    const std::vector<ScratchGFView> &views) {
  validate_manifest(manifest);
  validate_views(manifest, views);

  std::vector<Entry> entries;
  entries.reserve(views.size());
  for (const auto &view : views) {
    auto owned =
        std::make_unique<amrex::MultiFab>(view.multifab->deepCopy());
    entries.push_back(Entry{view.key, std::move(owned), *view.validity});
  }
  return ScratchLevelFrame(manifest.hierarchy_epoch, manifest.level,
                           std::move(entries));
}

ScratchLevelFrame::ScratchLevelFrame(const std::int64_t hierarchy_epoch,
                                     const int level,
                                     std::vector<Entry> entries)
    : hierarchy_epoch_(hierarchy_epoch), level_(level),
      entries_(std::move(entries)) {}

ScratchLevelFrame::~ScratchLevelFrame() = default;

ScratchLevelFrame::ScratchLevelFrame(ScratchLevelFrame &&other) noexcept
    : hierarchy_epoch_(other.hierarchy_epoch_), level_(other.level_),
      entries_(std::move(other.entries_)) {
  other.hierarchy_epoch_ = -1;
  other.level_ = -1;
  other.entries_.clear();
}

ScratchLevelFrame &
ScratchLevelFrame::operator=(ScratchLevelFrame &&other) noexcept {
  if (this == &other)
    return *this;
  hierarchy_epoch_ = other.hierarchy_epoch_;
  level_ = other.level_;
  entries_ = std::move(other.entries_);
  other.hierarchy_epoch_ = -1;
  other.level_ = -1;
  other.entries_.clear();
  return *this;
}

std::int64_t ScratchLevelFrame::hierarchy_epoch() const noexcept {
  return hierarchy_epoch_;
}

int ScratchLevelFrame::level() const noexcept { return level_; }

std::size_t ScratchLevelFrame::entry_count() const noexcept {
  return entries_.size();
}

const ScratchGFKey &ScratchLevelFrame::key(const std::size_t index) const {
  return entries_.at(index).key;
}

const amrex::MultiFab &
ScratchLevelFrame::multifab(const std::size_t index) const {
  return *entries_.at(index).multifab;
}

amrex::MultiFab &
ScratchLevelFrame::mutable_multifab(const std::size_t index) {
  return *entries_.at(index).multifab;
}

const std::vector<ScratchValidity> &
ScratchLevelFrame::validity(const std::size_t index) const {
  return entries_.at(index).validity;
}

std::vector<ScratchValidity> &
ScratchLevelFrame::mutable_validity(const std::size_t index) {
  return entries_.at(index).validity;
}

#ifndef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
ScratchLevelFrame ScratchLevelFrame::clone_for_transaction() const {
  std::vector<Entry> cloned;
  cloned.reserve(entries_.size());
  for (const auto &entry : entries_) {
    auto multifab = std::make_unique<amrex::MultiFab>(entry.multifab->deepCopy());
    cloned.push_back(Entry{entry.key, std::move(multifab), entry.validity});
  }
  return ScratchLevelFrame(hierarchy_epoch_, level_, std::move(cloned));
}

void ScratchLevelFrame::restore_from_transaction(
    const ScratchLevelFrame &source) {
  if (hierarchy_epoch_ != source.hierarchy_epoch_ || level_ != source.level_ ||
      entries_.size() != source.entries_.size())
    throw std::invalid_argument("scratch transaction frame metadata differs");

  for (std::size_t index = 0; index < entries_.size(); ++index) {
    auto &destination = entries_[index];
    const auto &source_entry = source.entries_[index];
    if (!same_key(destination.key, source_entry.key) ||
        !same_transaction_layout(*destination.multifab,
                                 *source_entry.multifab) ||
        destination.validity.size() != source_entry.validity.size())
      throw std::invalid_argument("scratch transaction frame schema differs");
  }

  for (std::size_t index = 0; index < entries_.size(); ++index) {
    auto &destination = entries_[index];
    const auto &source_entry = source.entries_[index];
    amrex::MultiFab::Copy(*destination.multifab, *source_entry.multifab, 0, 0,
                          destination.multifab->nComp(),
                          destination.multifab->nGrowVect());
    for (std::size_t component = 0; component < destination.validity.size();
         ++component)
      destination.validity[component] = source_entry.validity[component];
  }
}

void ScratchLevelFrame::copy_interior_for_transaction(
    const std::size_t destination_entry, const ScratchLevelFrame &source,
    const std::size_t source_entry) {
  auto &destination = *entries_.at(destination_entry).multifab;
  const auto &source_multifab = *source.entries_.at(source_entry).multifab;
  if (!same_transaction_layout(destination, source_multifab))
    throw std::invalid_argument(
        "scratch transaction interior copy layout differs");
  amrex::MultiFab::Copy(destination, source_multifab, 0, 0,
                        destination.nComp(), 0);
}

void ScratchLevelFrame::copy_validity_for_transaction(
    const std::size_t destination_entry, const ScratchLevelFrame &source,
    const std::size_t source_entry) {
  auto &destination = entries_.at(destination_entry).validity;
  const auto &source_validity = source.entries_.at(source_entry).validity;
  if (destination.size() != source_validity.size())
    throw std::invalid_argument(
        "scratch transaction validity component count differs");
  for (std::size_t component = 0; component < destination.size(); ++component)
    destination[component] = source_validity[component];
}

void ScratchLevelFrame::copy_interior_validity_for_transaction(
    const std::size_t destination_entry, const ScratchLevelFrame &source,
    const std::size_t source_entry, const bool invalidate_noninterior) {
  auto &destination = entries_.at(destination_entry).validity;
  const auto &source_validity = source.entries_.at(source_entry).validity;
  if (destination.size() != source_validity.size())
    throw std::invalid_argument(
        "scratch transaction validity component count differs");
  for (std::size_t component = 0; component < destination.size(); ++component) {
    destination[component].interior = source_validity[component].interior;
    if (invalidate_noninterior) {
      destination[component].outer = false;
      destination[component].ghosts = false;
    }
  }
}

void ScratchLevelFrame::linear_combination_interior_for_transaction(
    const std::size_t destination_entry, const double destination_scale,
    const std::vector<double> &weights,
    const std::vector<const ScratchLevelFrame *> &sources,
    const std::vector<std::size_t> &source_entries) {
  if (weights.size() != sources.size() ||
      sources.size() != source_entries.size())
    throw std::invalid_argument(
        "scratch transaction combination metadata differs");
  auto &destination = *entries_.at(destination_entry).multifab;
  for (std::size_t source = 0; source < sources.size(); ++source) {
    if (sources[source] == nullptr)
      throw std::invalid_argument("scratch transaction source is null");
    const auto &source_multifab =
        *sources[source]->entries_.at(source_entries[source]).multifab;
    if (sources[source] == this &&
        source_entries[source] == destination_entry)
      throw std::invalid_argument(
          "scratch transaction combination rejects destination aliasing");
    if (!same_transaction_layout(destination, source_multifab))
      throw std::invalid_argument(
          "scratch transaction combination layout differs");
  }

  if (destination_scale == 0.0)
    destination.setVal(0.0, 0, destination.nComp(), 0);
  else
    destination.mult(destination_scale, 0, destination.nComp(), 0);
  for (std::size_t source = 0; source < sources.size(); ++source) {
    const auto &source_multifab =
        *sources[source]->entries_.at(source_entries[source]).multifab;
    amrex::MultiFab::Saxpy(destination, weights[source], source_multifab, 0, 0,
                          destination.nComp(), 0);
  }
}

const void *ScratchLevelFrame::entry_address_for_transaction(
    const std::size_t entry) const {
  return entries_.at(entry).multifab.get();
}
#endif

} // namespace CarpetX
