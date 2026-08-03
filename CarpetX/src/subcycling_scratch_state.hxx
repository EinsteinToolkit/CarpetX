#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_HXX
#define CARPETX_SUBCYCLING_SCRATCH_STATE_HXX

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amrex {
class MultiFab;
}

namespace CarpetX {

class ScratchStateTransaction;
class ScratchStateTransactionFactory;
class ScratchStateTransactionCore;

struct ScratchGFKey {
  std::int64_t hierarchy_epoch;
  int level;
  int patch;
  int group_index;
};

struct ScratchValidity {
  bool interior;
  bool outer;
  bool ghosts;
};

struct ScratchGFView {
  ScratchGFKey key;
  int time_level;
  const void *source_storage_identity;
  const amrex::MultiFab *multifab;
  const std::vector<ScratchValidity> *validity;
};

struct ScratchGFManifestEntry {
  ScratchGFKey key;
  const void *source_storage_identity;
};

// This manifest is agreement supplied by a caller, not a certified inventory
// of the runtime GridFunction set. The identity token is opaque and is only
// compared; it is never dereferenced by this primitive.
struct UncertifiedScratchLevelManifest {
  std::int64_t hierarchy_epoch;
  int level;
  std::vector<ScratchGFManifestEntry> entries;
};

class ScratchLevelFrame {
public:
  // A defined source whose storage was created with allocation disabled is a
  // caller-contract violation in Phase 8A. Collective source certification is
  // deliberately deferred to Phase 8B.
  static ScratchLevelFrame
  copy_tl0(const UncertifiedScratchLevelManifest &manifest,
           const std::vector<ScratchGFView> &views);

  ~ScratchLevelFrame();
  ScratchLevelFrame(const ScratchLevelFrame &) = delete;
  ScratchLevelFrame &operator=(const ScratchLevelFrame &) = delete;
  ScratchLevelFrame(ScratchLevelFrame &&) noexcept;
  ScratchLevelFrame &operator=(ScratchLevelFrame &&) noexcept;

  std::int64_t hierarchy_epoch() const noexcept;
  int level() const noexcept;
  std::size_t entry_count() const noexcept;
  const ScratchGFKey &key(std::size_t index) const;
  const amrex::MultiFab &multifab(std::size_t index) const;
  amrex::MultiFab &mutable_multifab(std::size_t index);
  const std::vector<ScratchValidity> &validity(std::size_t index) const;
  std::vector<ScratchValidity> &mutable_validity(std::size_t index);

private:
  struct Entry;
  ScratchLevelFrame(std::int64_t hierarchy_epoch, int level,
                    std::vector<Entry> entries);

  ScratchLevelFrame clone_for_transaction() const;
  void restore_from_transaction(const ScratchLevelFrame &source);
  void copy_interior_for_transaction(std::size_t destination_entry,
                                     const ScratchLevelFrame &source,
                                     std::size_t source_entry);
  void copy_validity_for_transaction(std::size_t destination_entry,
                                     const ScratchLevelFrame &source,
                                     std::size_t source_entry);
  void copy_interior_validity_for_transaction(
      std::size_t destination_entry, const ScratchLevelFrame &source,
      std::size_t source_entry, bool invalidate_noninterior);
  void linear_combination_interior_for_transaction(
      std::size_t destination_entry, double destination_scale,
      const std::vector<double> &weights,
      const std::vector<const ScratchLevelFrame *> &sources,
      const std::vector<std::size_t> &source_entries);
  const void *entry_address_for_transaction(std::size_t entry) const;

  friend class ScratchStateTransaction;
  friend class ScratchStateTransactionFactory;
  friend class ScratchStateTransactionCore;

  std::int64_t hierarchy_epoch_;
  int level_;
  std::vector<Entry> entries_;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_STATE_HXX
