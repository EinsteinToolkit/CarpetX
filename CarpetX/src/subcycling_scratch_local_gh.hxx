#ifndef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_HXX
#define CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_HXX

#include "subcycling_scratch_adapter.hxx"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#if !defined(CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_STANDALONE) &&              \
    !defined(CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE)
#include <cctk.h>
#else
struct cGH;
#endif

namespace CarpetX {

class GHExt;
class CertifiedScratchStageExecutor;
class ScratchStateTransaction;
class ScratchStateTransactionFactory;
struct ScratchStageCoordinates;
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
class ScratchLocalLevelBindingTestAccess;
#endif

// Owns a certified scratch frame and non-executing local cGH views into it.
// The object cannot move because every cGH data pointer is tied to frame_.
class ScratchLocalLevelBinding final {
public:
  [[nodiscard]] static std::unique_ptr<ScratchLocalLevelBinding>
  bind(const GHExt &ghext, CertifiedScratchLevelFrame &&frame);

  ~ScratchLocalLevelBinding();
  ScratchLocalLevelBinding(const ScratchLocalLevelBinding &) = delete;
  ScratchLocalLevelBinding &
  operator=(const ScratchLocalLevelBinding &) = delete;
  ScratchLocalLevelBinding(ScratchLocalLevelBinding &&) = delete;
  ScratchLocalLevelBinding &operator=(ScratchLocalLevelBinding &&) = delete;

  int patch() const noexcept;
  int level() const noexcept;
  std::int64_t hierarchy_epoch() const noexcept;
  std::size_t local_context_count() const noexcept;

private:
  struct Storage;
  explicit ScratchLocalLevelBinding(CertifiedScratchLevelFrame &&frame);
  void build_contexts(const GHExt &ghext);
  cGH *context_for_executor(std::size_t ordinal) noexcept;
  ScratchLevelFrame &frame_for_executor() noexcept;
  ScratchLevelFrame &frame_for_transaction() noexcept { return frame_.frame_; }
  int level_for_executor() const noexcept;
  void claim_for_executor();
  void release_for_executor() noexcept;
  void fault_for_executor() noexcept;
  void install_stage_coordinates_for_executor(
      const ScratchStageCoordinates &coordinates) noexcept;

  friend class CertifiedScratchStageExecutor;
  friend class ScratchStateTransaction;
  friend class ScratchStateTransactionFactory;

#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
  friend class ScratchLocalLevelBindingTestAccess;
#endif

  // Contexts are destroyed first, then the frame their pointers address.
  CertifiedScratchLevelFrame frame_;
  std::vector<std::unique_ptr<Storage>> contexts_;
  std::atomic_bool executor_busy_{false};
  std::atomic_bool executor_faulted_{false};
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_HXX
