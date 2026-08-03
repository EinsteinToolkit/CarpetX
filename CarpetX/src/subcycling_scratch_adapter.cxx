#include "subcycling_scratch_adapter.hxx"

#include <AMReX_MultiFab.H>

#ifndef CARPETX_SUBCYCLING_ADAPTER_STANDALONE
#include "driver.hxx"

#include <AMReX_ParallelContext.H>
#endif

// driver.hxx provides the native Cactus feature macros (including CCTK_MPI).
// Include the internal declarations only after those macros are visible so the
// MPI collective declaration and its implementation use the same guard.
#include "subcycling_scratch_adapter_internal.hxx"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(CARPETX_SUBCYCLING_ADAPTER_MPI_UNIT) || defined(CCTK_MPI)
#include <mpi.h>
#endif

namespace CarpetX {

CertifiedScratchLevelFrame::CertifiedScratchLevelFrame(
    const int patch, ScratchLevelFrame frame)
    : patch_(patch), frame_(std::move(frame)) {}

CertifiedScratchLevelFrame::CertifiedScratchLevelFrame(
    CertifiedScratchLevelFrame &&other) noexcept
    : patch_(other.patch_), frame_(std::move(other.frame_)) {
  other.patch_ = -1;
}

CertifiedScratchLevelFrame &CertifiedScratchLevelFrame::operator=(
    CertifiedScratchLevelFrame &&other) noexcept {
  if (this == &other)
    return *this;
  patch_ = other.patch_;
  frame_ = std::move(other.frame_);
  other.patch_ = -1;
  return *this;
}

int CertifiedScratchLevelFrame::patch() const noexcept { return patch_; }

const ScratchLevelFrame &CertifiedScratchLevelFrame::frame() const noexcept {
  return frame_;
}

namespace subcycling_detail {
namespace {

[[noreturn]] void reject(const char *message) {
  throw ScratchAdapterCollectiveError(message);
}

bool size_to_i64(const std::size_t value, std::int64_t &converted) noexcept {
  if constexpr (sizeof(std::size_t) > sizeof(std::int64_t)) {
    if (value > static_cast<std::size_t>(
                    std::numeric_limits<std::int64_t>::max()))
      return false;
  }
  converted = static_cast<std::int64_t>(value);
  return true;
}

bool matches_prevalidated_size(const std::size_t actual,
                               const std::int64_t prevalidated) noexcept {
  std::int64_t converted = 0;
  return prevalidated >= 0 && size_to_i64(actual, converted) &&
         converted == prevalidated;
}

} // namespace

ScratchLevelFrame run_certification_transaction(
    LocalScratchBatch &batch, CollectiveOps &collectives,
    const EpochReader &epoch_reader, const CopyTL0 &copy_tl0) {
  auto effective_status = batch.preflight_status;
  if (!matches_prevalidated_size(batch.entries.size(),
                                 batch.prevalidated_entry_count))
    effective_status = ScratchAdapterFailure::invalid_schema;
  for (const auto &entry : batch.entries)
    if (!matches_prevalidated_size(entry.schema_fields.size(),
                                   entry.prevalidated_schema_length))
      effective_status = ScratchAdapterFailure::invalid_schema;

  const auto local_status = static_cast<std::int64_t>(effective_status);
  const auto local_count = batch.prevalidated_entry_count >= 0
                               ? batch.prevalidated_entry_count
                               : std::int64_t{0};

  const auto status_min = collectives.reduce_min(
      ScratchCollectivePhase::preflight_status, 0, local_status);
  const auto status_max = collectives.reduce_max(
      ScratchCollectivePhase::preflight_status, 0, local_status);
  const auto count_min = collectives.reduce_min(
      ScratchCollectivePhase::entry_count, 0, local_count);
  const auto count_max = collectives.reduce_max(
      ScratchCollectivePhase::entry_count, 0, local_count);
  if (status_min != 0 || status_max != 0)
    reject("canonical TL0 local preflight rejected collectively");
  if (count_min != count_max || count_min <= 0)
    reject("canonical TL0 entry count disagrees across ranks");

  bool lengths_agree = true;
  for (std::size_t entry_index = 0; entry_index < batch.entries.size();
       ++entry_index) {
    const auto local_length =
        batch.entries[entry_index].prevalidated_schema_length;
    const auto length_min = collectives.reduce_min(
        ScratchCollectivePhase::schema_lengths, entry_index, local_length);
    const auto length_max = collectives.reduce_max(
        ScratchCollectivePhase::schema_lengths, entry_index, local_length);
    lengths_agree = lengths_agree && length_min == length_max &&
                    length_min > 0;
  }
  if (!lengths_agree)
    reject("canonical TL0 schema descriptor lengths disagree across ranks");

  bool schema_agrees = true;
  std::size_t field_ordinal = 0;
  for (const auto &entry : batch.entries) {
    for (const auto local_field : entry.schema_fields) {
      const auto field_min = collectives.reduce_min(
          ScratchCollectivePhase::schema_fields, field_ordinal, local_field);
      const auto field_max = collectives.reduce_max(
          ScratchCollectivePhase::schema_fields, field_ordinal, local_field);
      schema_agrees =
          schema_agrees && field_min == field_max && field_min == local_field;
      ++field_ordinal;
    }
  }
  if (!schema_agrees)
    reject("canonical TL0 ordered schema disagrees across ranks");

  bool all_sources_ok = true;
  for (const auto &entry : batch.entries) {
    const bool entry_ok = entry.multifab->ok();
    all_sources_ok = all_sources_ok && entry_ok;
  }
  const bool sources_ok = collectives.reduce_and(
      ScratchCollectivePhase::source_ok, 0, all_sources_ok);
  if (!sources_ok)
    reject("canonical TL0 source allocation check rejected collectively");

  const std::int64_t observed_epoch = epoch_reader();
  const auto epoch_min = collectives.reduce_min(
      ScratchCollectivePhase::epoch_reader, 0, observed_epoch);
  const auto epoch_max = collectives.reduce_max(
      ScratchCollectivePhase::epoch_reader, 0, observed_epoch);
  if (epoch_min != epoch_max || epoch_min != batch.captured_epoch)
    reject("hierarchy epoch changed during canonical TL0 certification");

  UncertifiedScratchLevelManifest manifest{batch.captured_epoch, batch.level,
                                            {}};
  std::vector<ScratchGFView> views;
  manifest.entries.reserve(batch.entries.size());
  views.reserve(batch.entries.size());
  for (auto &entry : batch.entries) {
    manifest.entries.push_back(
        ScratchGFManifestEntry{entry.key, entry.source_storage_identity});
    views.push_back(ScratchGFView{entry.key, 0, entry.source_storage_identity,
                                  entry.multifab, &entry.validity});
  }

  std::optional<ScratchLevelFrame> local_frame;
  std::string local_diagnostic;
  try {
    local_frame.emplace(copy_tl0(manifest, views));
  } catch (const std::bad_alloc &) {
    throw;
  } catch (const std::exception &error) {
    local_diagnostic = error.what();
  }
  const bool all_copies_succeeded = collectives.reduce_and(
      ScratchCollectivePhase::copy_success, 0, local_frame.has_value());
  if (!all_copies_succeeded) {
    local_frame.reset();
    if (local_diagnostic.empty())
      reject("canonical TL0 owned copy failed on a peer rank");
    throw ScratchAdapterCollectiveError(
        "canonical TL0 owned copy failed collectively: " +
        local_diagnostic);
  }
  if (!local_frame.has_value())
    reject("canonical TL0 owned copy protocol lost its local result");
  return std::move(*local_frame);
}

#if defined(CARPETX_SUBCYCLING_ADAPTER_MPI_UNIT) || defined(CCTK_MPI)
MpiCollectiveOps::MpiCollectiveOps(void *const communicator_storage) noexcept
    : communicator_storage_(communicator_storage) {}

std::int64_t MpiCollectiveOps::reduce_min(
    const ScratchCollectivePhase, const std::size_t,
    const std::int64_t local) {
  std::int64_t global = 0;
  const MPI_Comm communicator =
      *static_cast<const MPI_Comm *>(communicator_storage_);
  if (MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_MIN, communicator) !=
      MPI_SUCCESS)
    throw std::runtime_error("MPI minimum reduction failed");
  return global;
}

std::int64_t MpiCollectiveOps::reduce_max(
    const ScratchCollectivePhase, const std::size_t,
    const std::int64_t local) {
  std::int64_t global = 0;
  const MPI_Comm communicator =
      *static_cast<const MPI_Comm *>(communicator_storage_);
  if (MPI_Allreduce(&local, &global, 1, MPI_INT64_T, MPI_MAX, communicator) !=
      MPI_SUCCESS)
    throw std::runtime_error("MPI maximum reduction failed");
  return global;
}

bool MpiCollectiveOps::reduce_and(const ScratchCollectivePhase,
                                  const std::size_t, const bool local) {
  const int local_integer = local ? 1 : 0;
  int global_integer = 0;
  const MPI_Comm communicator =
      *static_cast<const MPI_Comm *>(communicator_storage_);
  if (MPI_Allreduce(&local_integer, &global_integer, 1, MPI_INT, MPI_LAND,
                    communicator) != MPI_SUCCESS)
    throw std::runtime_error("MPI logical-and reduction failed");
  return global_integer != 0;
}
#endif

} // namespace subcycling_detail

#ifndef CARPETX_SUBCYCLING_ADAPTER_STANDALONE
namespace {

using subcycling_detail::LocalScratchBatch;
using subcycling_detail::LocalScratchEntry;
using subcycling_detail::ScratchAdapterFailure;
using subcycling_detail::ScratchCollectivePhase;

class SerialCollectiveOps final : public subcycling_detail::CollectiveOps {
public:
  std::int64_t reduce_min(ScratchCollectivePhase, std::size_t,
                          const std::int64_t local) override {
    return local;
  }
  std::int64_t reduce_max(ScratchCollectivePhase, std::size_t,
                          const std::int64_t local) override {
    return local;
  }
  bool reduce_and(ScratchCollectivePhase, std::size_t,
                  const bool local) override {
    return local;
  }
};

bool native_size_to_i64(const std::size_t value,
                        std::int64_t &converted) noexcept {
  if constexpr (sizeof(std::size_t) > sizeof(std::int64_t)) {
    if (value > static_cast<std::size_t>(
                    std::numeric_limits<std::int64_t>::max()))
      return false;
  }
  converted = static_cast<std::int64_t>(value);
  return true;
}

template <typename Integer>
bool integer_to_i64(const Integer value, std::int64_t &converted) noexcept {
  using Raw = std::remove_cv_t<Integer>;
  static_assert(std::is_integral_v<Raw>);
  if constexpr (std::is_signed_v<Raw>) {
    if constexpr (sizeof(Raw) > sizeof(std::int64_t)) {
      if (value < static_cast<Raw>(std::numeric_limits<std::int64_t>::min()) ||
          value > static_cast<Raw>(std::numeric_limits<std::int64_t>::max()))
        return false;
    }
  } else {
    if constexpr (sizeof(Raw) >= sizeof(std::int64_t)) {
      if (value > static_cast<Raw>(
                      std::numeric_limits<std::int64_t>::max()))
        return false;
    }
  }
  converted = static_cast<std::int64_t>(value);
  return true;
}

template <typename Integer>
bool append_integer(std::vector<std::int64_t> &fields, const Integer value) {
  std::int64_t converted = 0;
  if (!integer_to_i64(value, converted))
    return false;
  fields.push_back(converted);
  return true;
}

bool checked_add(const std::size_t left, const std::size_t right,
                 std::size_t &sum) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left)
    return false;
  sum = left + right;
  return true;
}

bool checked_multiply(const std::size_t left, const std::size_t right,
                      std::size_t &product) noexcept {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  product = left * right;
  return true;
}

bool descriptor_length(const std::size_t validity_count,
                       const std::size_t box_count,
                       const std::size_t processor_count,
                       std::size_t &length) noexcept {
  std::size_t validity_fields = 0;
  std::size_t box_fields = 0;
  std::size_t total = 18;
  return checked_multiply(validity_count, 3, validity_fields) &&
         checked_multiply(box_count, 6, box_fields) &&
         checked_add(total, validity_fields, total) &&
         checked_add(total, box_fields, total) &&
         checked_add(total, processor_count, length);
}

LocalScratchBatch inventory_canonical_tl0(const GHExt &ghext,
                                          const int patch,
                                          const int level) {
  LocalScratchBatch batch{ScratchAdapterFailure::none,
                          static_cast<std::int64_t>(CarpetX_GetEpoch()), patch,
                          level, {}, 0};
  if (patch < 0 || static_cast<std::size_t>(patch) >= ghext.patchdata.size()) {
    batch.preflight_status = ScratchAdapterFailure::invalid_patch;
    return batch;
  }
  const auto &patchdata = ghext.patchdata[static_cast<std::size_t>(patch)];
  if (level < 0 ||
      static_cast<std::size_t>(level) >= patchdata.leveldata.size()) {
    batch.preflight_status = ScratchAdapterFailure::invalid_level;
    return batch;
  }
  if (batch.captured_epoch < 0) {
    batch.preflight_status = ScratchAdapterFailure::invalid_epoch;
    return batch;
  }

  const auto &leveldata = patchdata.leveldata[static_cast<std::size_t>(level)];
  std::size_t complete_count = 0;
  for (const auto &group : leveldata.groupdata) {
    if (group != nullptr) {
      if (complete_count == std::numeric_limits<std::size_t>::max()) {
        batch.preflight_status = ScratchAdapterFailure::invalid_schema;
        return batch;
      }
      ++complete_count;
    }
  }
  if (!native_size_to_i64(complete_count,
                          batch.prevalidated_entry_count) ||
      complete_count > batch.entries.max_size()) {
    batch.preflight_status = ScratchAdapterFailure::invalid_schema;
    return batch;
  }
  batch.entries.reserve(complete_count);

  for (std::size_t group_index = 0;
       group_index < leveldata.groupdata.size(); ++group_index) {
    const auto &group_owner = leveldata.groupdata[group_index];
    if (group_owner == nullptr)
      continue;
    const auto &groupdata = *group_owner;
    if (group_index >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        groupdata.groupindex != static_cast<int>(group_index) ||
        groupdata.mfab.empty() || groupdata.valid.empty() ||
        groupdata.mfab.size() != groupdata.valid.size()) {
      batch.preflight_status = ScratchAdapterFailure::invalid_tl0;
      continue;
    }
    std::int64_t ignored_size = 0;
    if (!native_size_to_i64(groupdata.mfab.size(), ignored_size) ||
        !native_size_to_i64(groupdata.valid.size(), ignored_size)) {
      batch.preflight_status = ScratchAdapterFailure::invalid_schema;
      continue;
    }
    const auto *mfab = groupdata.mfab.at(0).get();
    const auto &why_valid = groupdata.valid.at(0);
    if (mfab == nullptr || !mfab->isDefined() || groupdata.numvars <= 0 ||
        mfab->nComp() <= 0 || mfab->nComp() != groupdata.numvars ||
        why_valid.size() != static_cast<std::size_t>(groupdata.numvars)) {
      batch.preflight_status = ScratchAdapterFailure::invalid_components;
      continue;
    }

    const auto &box_array = mfab->boxArray();
    const auto box_count = box_array.size();
    const auto &processor_map = mfab->DistributionMap().ProcessorMap();
    if (!box_array.ok() || box_count < 0 ||
        box_count > std::numeric_limits<int>::max() ||
        static_cast<std::size_t>(box_count) != processor_map.size() ||
        !native_size_to_i64(why_valid.size(), ignored_size) ||
        !native_size_to_i64(processor_map.size(), ignored_size)) {
      batch.preflight_status = ScratchAdapterFailure::invalid_schema;
      continue;
    }

    std::vector<ScratchValidity> validity;
    if (why_valid.size() > validity.max_size()) {
      batch.preflight_status = ScratchAdapterFailure::invalid_schema;
      continue;
    }
    validity.reserve(why_valid.size());
    for (const auto &why : why_valid) {
      const auto &bits = why.get();
      validity.push_back(ScratchValidity{bits.valid_int, bits.valid_outer,
                                         bits.valid_ghosts});
    }

    std::vector<std::int64_t> schema_fields;
    std::size_t expected_schema_length = 0;
    std::int64_t prevalidated_schema_length = 0;
    if (!descriptor_length(validity.size(),
                           static_cast<std::size_t>(box_count),
                           processor_map.size(), expected_schema_length) ||
        !native_size_to_i64(expected_schema_length,
                            prevalidated_schema_length) ||
        expected_schema_length > schema_fields.max_size()) {
      batch.preflight_status = ScratchAdapterFailure::invalid_schema;
      continue;
    }
    schema_fields.reserve(expected_schema_length);
    bool fields_representable = true;
    const auto append = [&](const auto value) {
      if (fields_representable)
        fields_representable = append_integer(schema_fields, value);
    };
    append(batch.captured_epoch);
    append(patch);
    append(level);
    append(group_index);
    append(groupdata.numvars);
    append(mfab->nComp());
    append(groupdata.mfab.size());
    append(mfab->size());
    append(groupdata.valid.size());
    append(why_valid.size());
    for (const auto &bits : validity) {
      append(bits.interior);
      append(bits.outer);
      append(bits.ghosts);
    }
    append(box_count);
    const auto index_type = box_array.ixType().toIntVect();
    for (int direction = 0; direction < 3; ++direction)
      append(index_type[direction]);
    for (int box_index = 0; box_index < static_cast<int>(box_count);
         ++box_index) {
      const auto box = box_array[box_index];
      for (int direction = 0; direction < 3; ++direction) {
        append(box.smallEnd(direction));
        append(box.bigEnd(direction));
      }
    }
    const auto grow = mfab->nGrowVect();
    for (int direction = 0; direction < 3; ++direction)
      append(grow[direction]);
    append(processor_map.size());
    for (const auto owner_rank : processor_map)
      append(owner_rank);
    if (!fields_representable ||
        schema_fields.size() != expected_schema_length) {
      batch.preflight_status = ScratchAdapterFailure::invalid_schema;
      continue;
    }

    batch.entries.push_back(LocalScratchEntry{
        ScratchGFKey{batch.captured_epoch, level, patch,
                     static_cast<int>(group_index)},
        static_cast<const void *>(mfab), mfab, std::move(validity),
        std::move(schema_fields), prevalidated_schema_length});
  }
  if (batch.entries.size() != complete_count &&
      batch.preflight_status == ScratchAdapterFailure::none)
    batch.preflight_status = ScratchAdapterFailure::omitted_group;
  return batch;
}

} // namespace

CertifiedScratchLevelFrame
copy_canonical_tl0_collective(const GHExt &ghext, const int patch,
                              const int level) {
  auto batch = inventory_canonical_tl0(ghext, patch, level);
#ifdef CCTK_MPI
  MPI_Comm communicator = amrex::ParallelContext::CommunicatorSub();
  subcycling_detail::MpiCollectiveOps collectives(&communicator);
#else
  SerialCollectiveOps collectives;
#endif
  auto frame = subcycling_detail::run_certification_transaction(
      batch, collectives,
      [] { return static_cast<std::int64_t>(CarpetX_GetEpoch()); },
      [](const UncertifiedScratchLevelManifest &manifest,
         const std::vector<ScratchGFView> &views) {
        return ScratchLevelFrame::copy_tl0(manifest, views);
      });
  return CertifiedScratchLevelFrame(patch, std::move(frame));
}
#endif

} // namespace CarpetX
