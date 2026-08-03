#ifndef CARPETX_SUBCYCLING_SCRATCH_ADAPTER_HXX
#define CARPETX_SUBCYCLING_SCRATCH_ADAPTER_HXX

#include "subcycling_scratch_state.hxx"

#include <utility>

namespace CarpetX {

class GHExt;
class ScratchLocalLevelBinding;
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
class ScratchLocalLevelBindingTestAccess;
#endif

class CertifiedScratchLevelFrame {
public:
  CertifiedScratchLevelFrame(const CertifiedScratchLevelFrame &) = delete;
  CertifiedScratchLevelFrame &
  operator=(const CertifiedScratchLevelFrame &) = delete;
  CertifiedScratchLevelFrame(CertifiedScratchLevelFrame &&) noexcept;
  CertifiedScratchLevelFrame &
  operator=(CertifiedScratchLevelFrame &&) noexcept;

  int patch() const noexcept;
  const ScratchLevelFrame &frame() const noexcept;

private:
  friend class ScratchLocalLevelBinding;
#ifdef CARPETX_SUBCYCLING_SCRATCH_LOCAL_GH_UNIT
  friend class ScratchLocalLevelBindingTestAccess;
#endif
  friend CertifiedScratchLevelFrame
  copy_canonical_tl0_collective(const GHExt &, int, int);
  CertifiedScratchLevelFrame(int patch, ScratchLevelFrame frame);

  int patch_;
  ScratchLevelFrame frame_;
};

// Collective on amrex::ParallelContext::CommunicatorSub(). All participating
// ranks must call at one quiescent hierarchy point with the same patch/level.
// There is deliberately no caller-selected communicator, manifest, or subset.
[[nodiscard]] CertifiedScratchLevelFrame
copy_canonical_tl0_collective(const GHExt &ghext, int patch, int level);

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_ADAPTER_HXX
