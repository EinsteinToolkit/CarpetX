#ifndef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_FACTORY_HXX
#define CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_FACTORY_HXX

#include "subcycling_scratch_adapter.hxx"
#include "subcycling_scratch_schedule_executor.hxx"
#include "subcycling_scratch_state_transaction.hxx"

#include <AMReX_MultiFab.H>

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace CarpetX {

class CertifiedLocalScheduleRegistry;
class GHExt;

using ScratchLiveEntryRestorer = std::function<void(
    const amrex::MultiFab &, const std::vector<ScratchValidity> &)>;

struct ScratchLiveEntrySnapshot {
  ScratchGFKey key;
  const void *level_identity;
  const void *group_identity;
  const void *tl0_identity;
  const void *storage_identity;
  const void *validity_identity;
  const amrex::MultiFab *multifab;
  std::vector<ScratchValidity> validity;
  bool grid_function_real;
  ScratchLiveEntryRestorer restore;
};

using ScratchLiveEntryReader =
    std::function<ScratchLiveEntrySnapshot()>;

struct ScratchStateTransactionFactoryMetadata {
  std::int64_t hierarchy_epoch;
  int patch_count;
  int level_count;
  int level_iteration;
  step_clock_t base_delta_clock;
  double base_delta_time;
  int time_refinement_factor;
  bool poison_undefined_values;
  std::vector<ScratchGroupPair> group_pairs;
  std::vector<int> dependent_groups;
  std::function<std::int64_t()> epoch_reader;
  std::vector<ScratchLiveEntryReader> live_entry_readers;
};

class ScratchStateTransactionFactory final {
public:
  static std::unique_ptr<ScratchStateTransaction> create_native(
      const GHExt &ghext, CertifiedLocalScheduleRegistry &registry,
      CertifiedScratchLevelFrame &&frame,
      ScratchStateTransactionFactoryMetadata metadata);

#ifdef CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_UNIT
  using TestScheduleExecutor = std::function<ScratchScheduleExecutionReceipt(
      ScratchSchedulePhase, const StepContext &,
      const ScratchStageCoordinates &)>;

  static std::unique_ptr<ScratchStateTransaction> create_for_test(
      ScratchLevelFrame working_frame,
      ScratchStateTransactionFactoryMetadata metadata,
      TestScheduleExecutor executor);

  static const void *working_entry_address_for_test(
      const ScratchStateTransaction &transaction, std::size_t entry);
  static const amrex::MultiFab &working_multifab_for_test(
      const ScratchStateTransaction &transaction, std::size_t entry);
  static amrex::MultiFab &mutable_working_multifab_for_test(
      ScratchStateTransaction &transaction, std::size_t entry);
  static const std::vector<ScratchValidity> &working_validity_for_test(
      const ScratchStateTransaction &transaction, std::size_t entry);
  static std::vector<ScratchValidity> &mutable_working_validity_for_test(
      ScratchStateTransaction &transaction, std::size_t entry);
  static ScratchStateToken stale_epoch_token_for_test(
      const ScratchStateTransaction &transaction,
      const ScratchStateToken &source);
  static ScratchStateToken incompatible_schema_token_for_test(
      const ScratchStateTransaction &transaction,
      const ScratchStateToken &source);
#endif
};

} // namespace CarpetX

#endif // CARPETX_SUBCYCLING_SCRATCH_STATE_TRANSACTION_FACTORY_HXX
