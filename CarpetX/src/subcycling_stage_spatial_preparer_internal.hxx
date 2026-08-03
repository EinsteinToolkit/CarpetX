#ifndef CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_INTERNAL_HXX
#define CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_INTERNAL_HXX

#include "subcycling_stage_spatial_preparer.hxx"

#include <cstddef>
#include <cstdint>

namespace CarpetX {

struct StageSpatialMetadata {
  int patch_count;
  int level_count;
  int spatial_refinement_ratio;
  int level;
  std::int64_t transaction_epoch;
  std::int64_t observed_epoch;
  bool global_sync_ownership_conflict;
  std::size_t evolved_group_count;
};

namespace detail {

class StageSpatialPreparerBackend {
public:
  virtual ~StageSpatialPreparerBackend() = default;
  virtual StageSpatialMetadata
  inspect(ScratchStateTransaction *transaction) = 0;
  virtual void prepare_level_zero(ScratchStateTransaction *transaction,
                                  StageSpatialTarget target,
                                  step_clock_t stage_clock) = 0;
  virtual void prepare_level_one(ScratchStateTransaction *transaction,
                                 StageSpatialTarget target,
                                 step_clock_t stage_clock,
                                 step_clock_t parent_theta,
                                 const DenseInterval &parent_dense) = 0;
  virtual void promote(ScratchStateTransaction *transaction,
                       StageSpatialTarget target) = 0;
  virtual void fault(ScratchStateTransaction *transaction) noexcept = 0;
};

} // namespace detail
} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_INTERNAL_HXX
