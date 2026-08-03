#ifndef CARPETX_SUBCYCLING_SCRATCH_ADAPTER_INTERNAL_HXX
#define CARPETX_SUBCYCLING_SCRATCH_ADAPTER_INTERNAL_HXX

#include "subcycling_scratch_state.hxx"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <vector>

namespace CarpetX::subcycling_detail {

enum class ScratchAdapterFailure : std::int64_t {
  none = 0,
  invalid_patch = 1,
  invalid_level = 2,
  invalid_epoch = 3,
  invalid_tl0 = 4,
  invalid_components = 5,
  invalid_validity = 6,
  invalid_schema = 7,
  duplicate_group = 8,
  omitted_group = 9,
  extra_group = 10,
};

enum class ScratchCollectivePhase : std::int64_t {
  preflight_status = 0,
  entry_count = 1,
  schema_lengths = 2,
  schema_fields = 3,
  source_ok = 4,
  epoch_reader = 5,
  copy_success = 6,
};

struct ScratchSchemaRow {
  std::vector<std::int64_t> schema_fields;
};

struct LocalScratchEntry {
  ScratchGFKey key;
  const void *source_storage_identity;
  const amrex::MultiFab *multifab;
  std::vector<ScratchValidity> validity;
  std::vector<std::int64_t> schema_fields;
  std::int64_t prevalidated_schema_length;
};

struct LocalScratchBatch {
  ScratchAdapterFailure preflight_status;
  std::int64_t captured_epoch;
  int patch;
  int level;
  std::vector<LocalScratchEntry> entries;
  std::int64_t prevalidated_entry_count;
};

class CollectiveOps {
public:
  virtual ~CollectiveOps() = default;
  virtual std::int64_t reduce_min(ScratchCollectivePhase phase,
                                  std::size_t ordinal,
                                  std::int64_t local) = 0;
  virtual std::int64_t reduce_max(ScratchCollectivePhase phase,
                                  std::size_t ordinal,
                                  std::int64_t local) = 0;
  virtual bool reduce_and(ScratchCollectivePhase phase, std::size_t ordinal,
                          bool local) = 0;
};

class ScratchAdapterCollectiveError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

using EpochReader = std::function<std::int64_t()>;
using CopyTL0 = std::function<ScratchLevelFrame(
    const UncertifiedScratchLevelManifest &,
    const std::vector<ScratchGFView> &)>;

[[nodiscard]] ScratchLevelFrame run_certification_transaction(
    LocalScratchBatch &batch, CollectiveOps &collectives,
    const EpochReader &epoch_reader, const CopyTL0 &copy_tl0);

#if defined(CARPETX_SUBCYCLING_ADAPTER_MPI_UNIT) || defined(CCTK_MPI)
class MpiCollectiveOps final : public CollectiveOps {
public:
  explicit MpiCollectiveOps(void *communicator_storage) noexcept;
  std::int64_t reduce_min(ScratchCollectivePhase phase, std::size_t ordinal,
                          std::int64_t local) override;
  std::int64_t reduce_max(ScratchCollectivePhase phase, std::size_t ordinal,
                          std::int64_t local) override;
  bool reduce_and(ScratchCollectivePhase phase, std::size_t ordinal,
                  bool local) override;

private:
  void *communicator_storage_;
};
#endif

} // namespace CarpetX::subcycling_detail

#endif // CARPETX_SUBCYCLING_SCRATCH_ADAPTER_INTERNAL_HXX
