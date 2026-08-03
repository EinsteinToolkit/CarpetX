#ifndef CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_HXX
#define CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_HXX

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#if !defined(CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_STANDALONE) &&         \
    !defined(CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_STANDALONE)
#include <cctk.h>
#else
struct cGH;
#endif

namespace CarpetX {

class CertifiedScratchStageExecutor;
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_UNIT
class CertifiedLocalScheduleRegistryExecutorTestAccess;
#endif

enum class SubcyclingScheduleTarget : std::uint8_t { rhs, post_step };
enum class SpecialScheduleGroupState : std::uint8_t { missing };
enum class CanonicalScheduleMode : std::uint8_t { default_local, local };

struct SchedulePatchProvenance {
  std::string phase;
  std::string sha256;
  friend bool operator==(const SchedulePatchProvenance &,
                         const SchedulePatchProvenance &) noexcept;
};

struct ScheduleBuildProvenance {
  std::uint32_t schema_version{};
  std::uint32_t inventory_abi_version{};
  std::string cactus_flesh_parent;
  std::string cactus_inventory_patch_sha256;
  std::string carpetx_parent;
  std::vector<SchedulePatchProvenance> ordered_carpetx_patches;
  std::string compiler_id;
  std::string compiler_version;
  std::string target_triple;
  std::string cxx_abi;
  std::string byte_order;
  std::vector<std::pair<int, int>> cctk_type_sizes;
  std::string executable_sha256;
  friend bool operator==(const ScheduleBuildProvenance &,
                         const ScheduleBuildProvenance &) noexcept;
};

struct CanonicalTag {
  std::string key;
  int type_code{};
  std::uint64_t element_count{};
  std::vector<std::uint8_t> value_bytes;
  friend bool operator==(const CanonicalTag &, const CanonicalTag &) noexcept;
};

struct CanonicalAccess {
  std::string variable_name;
  int timelevel{};
  int read_mask{};
  int write_mask{};
  int invalidate_mask{};
  friend bool operator==(const CanonicalAccess &,
                         const CanonicalAccess &) noexcept;
};

struct CanonicalGroupTimelevel {
  std::string group;
  int timelevel{};
  friend bool operator==(const CanonicalGroupTimelevel &,
                         const CanonicalGroupTimelevel &) noexcept;
};

struct CanonicalExecutionFlags {
  int meta{};
  int meta_early{};
  int meta_late{};
  int global{};
  int global_early{};
  int global_late{};
  int level{};
  int singlemap{};
  int local{};
  int loop_meta{};
  int loop_global{};
  int loop_level{};
  int loop_singlemap{};
  int loop_local{};
  friend bool operator==(const CanonicalExecutionFlags &,
                         const CanonicalExecutionFlags &) noexcept;
};

struct CanonicalScheduleItem {
  std::string description;
  std::string implementation;
  std::string where;
  std::string thorn;
  std::string routine;
  int language{};
  int function_type{};
  CanonicalExecutionFlags execution_flags;
  CanonicalScheduleMode effective_mode{CanonicalScheduleMode::default_local};
  std::vector<CanonicalTag> tags;
  std::vector<std::string> reads_clauses;
  std::vector<std::string> writes_clauses;
  std::vector<std::string> invalidates_clauses;
  std::vector<CanonicalAccess> accesses;
  std::vector<CanonicalGroupTimelevel> storage_groups;
  std::vector<std::string> communication_groups;
  std::vector<std::string> sync_groups;
  std::vector<std::string> trigger_groups;
  std::vector<std::string> ifs;
  std::vector<std::string> whiles;
  friend bool operator==(const CanonicalScheduleItem &,
                         const CanonicalScheduleItem &) noexcept;
};

struct CanonicalScheduleScope {
  std::uint64_t item_ordinal{};
  std::uint64_t parent_local_ordinal{};
  CanonicalScheduleItem item;
  friend bool operator==(const CanonicalScheduleScope &,
                         const CanonicalScheduleScope &) noexcept;
};

struct CanonicalScheduleRecord {
  std::uint64_t traversal_ordinal{};
  std::uint64_t item_ordinal{};
  std::uint64_t parent_local_ordinal{};
  CanonicalScheduleItem item;
  std::vector<CanonicalScheduleScope> ancestry;
  friend bool operator==(const CanonicalScheduleRecord &,
                         const CanonicalScheduleRecord &) noexcept;
};

struct CanonicalTargetSchedule {
  SubcyclingScheduleTarget target{SubcyclingScheduleTarget::rhs};
  std::string body_root;
  SpecialScheduleGroupState entry_state{SpecialScheduleGroupState::missing};
  SpecialScheduleGroupState exit_state{SpecialScheduleGroupState::missing};
  bool entry_exit_included{};
  std::vector<CanonicalScheduleRecord> records;
  friend bool operator==(const CanonicalTargetSchedule &,
                         const CanonicalTargetSchedule &) noexcept;
};

struct CanonicalScheduleBundle final {
  std::array<CanonicalTargetSchedule, 2> targets;
  bool operator==(const CanonicalScheduleBundle &) const noexcept;
};

struct ScheduleCertificationExpectation {
  CanonicalScheduleBundle manifest;
  ScheduleBuildProvenance provenance;
};

enum class ScheduleCertificationErrorCode : std::uint8_t {
  inventory_error,
  special_group_present,
  abi_mismatch,
  malformed_record,
  callback_internal_error,
  unsupported_metadata,
  name_resolution_error,
  tag_error,
  manifest_mismatch,
  provenance_mismatch,
};

struct ScheduleCertificationFailure {
  ScheduleCertificationErrorCode code;
  std::optional<SubcyclingScheduleTarget> target;
  std::optional<std::uint64_t> traversal_ordinal;
  std::string field;
  std::string detail;
};

class CertifiedLocalScheduleRegistry;
struct ScheduleCertificationResult;

// Rank/process-local, read-only certification. The caller is responsible for
// supplying an independently reviewed expectation and trusted build
// provenance; this API intentionally offers no observed-to-expected approval
// path. It inventories metadata only and never executes a scheduled routine.
[[nodiscard]] ScheduleCertificationResult certify_local_subcycling_schedules(
    const ScheduleCertificationExpectation &expected,
    const ScheduleBuildProvenance &observed_provenance);

class CertifiedLocalScheduleRegistry final {
public:
  ~CertifiedLocalScheduleRegistry();
  CertifiedLocalScheduleRegistry(const CertifiedLocalScheduleRegistry &) =
      delete;
  CertifiedLocalScheduleRegistry &
  operator=(const CertifiedLocalScheduleRegistry &) = delete;
  CertifiedLocalScheduleRegistry(CertifiedLocalScheduleRegistry &&) = delete;
  CertifiedLocalScheduleRegistry &
  operator=(CertifiedLocalScheduleRegistry &&) = delete;

  const CanonicalScheduleBundle &manifest() const noexcept;
  const ScheduleBuildProvenance &provenance() const noexcept;
  std::size_t size(SubcyclingScheduleTarget target) const noexcept;

private:
  struct Storage;
  explicit CertifiedLocalScheduleRegistry(std::unique_ptr<Storage> storage);
  const CanonicalScheduleRecord &
  record_for_executor(SubcyclingScheduleTarget target,
                      std::size_t ordinal) const;
  int invoke_for_executor(SubcyclingScheduleTarget target,
                          std::size_t ordinal, cGH *scratch_gh) const;
  std::unique_ptr<Storage> storage_;
  friend class CertifiedScratchStageExecutor;
#ifdef CARPETX_SUBCYCLING_SCRATCH_SCHEDULE_EXECUTOR_UNIT
  friend class CertifiedLocalScheduleRegistryExecutorTestAccess;
#endif
  friend ScheduleCertificationResult certify_local_subcycling_schedules(
      const ScheduleCertificationExpectation &expected,
      const ScheduleBuildProvenance &observed_provenance);
};

struct ScheduleCertificationResult {
  std::unique_ptr<CertifiedLocalScheduleRegistry> registry;
  std::optional<ScheduleCertificationFailure> failure;
  explicit operator bool() const noexcept { return registry != nullptr; }
};

} // namespace CarpetX

#endif
