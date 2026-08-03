#ifndef CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_HXX
#define CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_HXX

#include "subcycling_dense_output.hxx"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

#ifdef CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_UNIT
#include <string>
#include <vector>
#endif

namespace CarpetX {

class GHExt;
class ScratchStateTransaction;
namespace detail {
class StageSpatialPreparerBackend;
}

enum class StageSpatialScheduleOwnership : std::uint8_t {
  certified_local_no_global_sync,
  conflicting_global_sync,
};

enum class StageSpatialTarget : std::uint8_t {
  primary_live_tl0,
  transaction_scratch,
};

struct StageSpatialPreparationReceipt {
  StageSpatialTarget target;
  int patch;
  int level;
  step_clock_t stage_clock;
  std::optional<step_clock_t> parent_theta;
  std::size_t evolved_group_count;
};

class TwoLevelStageSpatialPreparer final {
public:
  TwoLevelStageSpatialPreparer(
      GHExt &ghext, StageSpatialScheduleOwnership schedule_ownership);
  ~TwoLevelStageSpatialPreparer();

  TwoLevelStageSpatialPreparer(const TwoLevelStageSpatialPreparer &) = delete;
  TwoLevelStageSpatialPreparer &
  operator=(const TwoLevelStageSpatialPreparer &) = delete;
  TwoLevelStageSpatialPreparer(TwoLevelStageSpatialPreparer &&) noexcept;
  TwoLevelStageSpatialPreparer &
  operator=(TwoLevelStageSpatialPreparer &&) noexcept;

  StageSpatialPreparationReceipt
  prepare(ScratchStateTransaction &transaction, const StepContext &context,
          const StagePoint &stage_point, const DenseInterval *parent_dense);

#ifdef CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_UNIT
  struct TestMetadata {
    int patch_count{1};
    int level_count{1};
    int spatial_refinement_ratio{2};
    int level{0};
    std::int64_t transaction_epoch{7};
    std::int64_t observed_epoch{7};
    bool global_sync_ownership_conflict{false};
    std::size_t evolved_group_count{1};
  };

  struct TestBackend {
    TestMetadata metadata;
    bool throw_during_prepare{false};
    bool throw_during_promote{false};
    int prepare_calls{0};
    int promote_calls{0};
    int fault_calls{0};
    std::optional<step_clock_t> observed_parent_theta;
    std::vector<std::string> events;
  };

  static TwoLevelStageSpatialPreparer create_for_test(TestBackend &backend);
  StageSpatialPreparationReceipt
  prepare_for_test(const StepContext &context, const StagePoint &stage_point,
                   const DenseInterval *parent_dense);
#endif

private:
  explicit TwoLevelStageSpatialPreparer(
      std::unique_ptr<detail::StageSpatialPreparerBackend> backend);

  StageSpatialPreparationReceipt
  prepare_impl(ScratchStateTransaction *transaction,
               const StepContext &context, const StagePoint &stage_point,
               const DenseInterval *parent_dense);

  std::unique_ptr<detail::StageSpatialPreparerBackend> backend_;
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_STAGE_SPATIAL_PREPARER_HXX
