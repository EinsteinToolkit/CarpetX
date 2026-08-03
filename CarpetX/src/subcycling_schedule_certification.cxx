#include "subcycling_schedule_certification.hxx"
#include "subcycling_schedule_certification_native_gate.hxx"

#ifndef CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_STANDALONE
#include <cctk_Constants.h>
#include <cctk_Groups.h>
#include <cctk_Schedule.h>
#include <util_Table.h>
#endif

#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

namespace CarpetX {

bool operator==(const SchedulePatchProvenance &a,
                const SchedulePatchProvenance &b) noexcept {
  return std::tie(a.phase, a.sha256) == std::tie(b.phase, b.sha256);
}

bool operator==(const ScheduleBuildProvenance &a,
                const ScheduleBuildProvenance &b) noexcept {
  return std::tie(a.schema_version, a.inventory_abi_version,
                  a.cactus_flesh_parent, a.cactus_inventory_patch_sha256,
                  a.carpetx_parent, a.ordered_carpetx_patches, a.compiler_id,
                  a.compiler_version, a.target_triple, a.cxx_abi, a.byte_order,
                  a.cctk_type_sizes, a.executable_sha256) ==
         std::tie(b.schema_version, b.inventory_abi_version,
                  b.cactus_flesh_parent, b.cactus_inventory_patch_sha256,
                  b.carpetx_parent, b.ordered_carpetx_patches, b.compiler_id,
                  b.compiler_version, b.target_triple, b.cxx_abi, b.byte_order,
                  b.cctk_type_sizes, b.executable_sha256);
}

bool operator==(const CanonicalTag &a, const CanonicalTag &b) noexcept {
  return std::tie(a.key, a.type_code, a.element_count, a.value_bytes) ==
         std::tie(b.key, b.type_code, b.element_count, b.value_bytes);
}

bool operator==(const CanonicalAccess &a, const CanonicalAccess &b) noexcept {
  return std::tie(a.variable_name, a.timelevel, a.read_mask, a.write_mask,
                  a.invalidate_mask) ==
         std::tie(b.variable_name, b.timelevel, b.read_mask, b.write_mask,
                  b.invalidate_mask);
}

bool operator==(const CanonicalGroupTimelevel &a,
                const CanonicalGroupTimelevel &b) noexcept {
  return std::tie(a.group, a.timelevel) == std::tie(b.group, b.timelevel);
}

bool operator==(const CanonicalExecutionFlags &a,
                const CanonicalExecutionFlags &b) noexcept {
  return std::tie(a.meta, a.meta_early, a.meta_late, a.global, a.global_early,
                  a.global_late, a.level, a.singlemap, a.local, a.loop_meta,
                  a.loop_global, a.loop_level, a.loop_singlemap,
                  a.loop_local) ==
         std::tie(b.meta, b.meta_early, b.meta_late, b.global, b.global_early,
                  b.global_late, b.level, b.singlemap, b.local, b.loop_meta,
                  b.loop_global, b.loop_level, b.loop_singlemap,
                  b.loop_local);
}

bool operator==(const CanonicalScheduleItem &a,
                const CanonicalScheduleItem &b) noexcept {
  return std::tie(a.description, a.implementation, a.where, a.thorn, a.routine,
                  a.language, a.function_type, a.execution_flags,
                  a.effective_mode, a.tags, a.reads_clauses, a.writes_clauses,
                  a.invalidates_clauses, a.accesses, a.storage_groups,
                  a.communication_groups, a.sync_groups, a.trigger_groups,
                  a.ifs, a.whiles) ==
         std::tie(b.description, b.implementation, b.where, b.thorn, b.routine,
                  b.language, b.function_type, b.execution_flags,
                  b.effective_mode, b.tags, b.reads_clauses, b.writes_clauses,
                  b.invalidates_clauses, b.accesses, b.storage_groups,
                  b.communication_groups, b.sync_groups, b.trigger_groups,
                  b.ifs, b.whiles);
}

bool operator==(const CanonicalScheduleScope &a,
                const CanonicalScheduleScope &b) noexcept {
  return std::tie(a.item_ordinal, a.parent_local_ordinal, a.item) ==
         std::tie(b.item_ordinal, b.parent_local_ordinal, b.item);
}

bool operator==(const CanonicalScheduleRecord &a,
                const CanonicalScheduleRecord &b) noexcept {
  return std::tie(a.traversal_ordinal, a.item_ordinal, a.parent_local_ordinal,
                  a.item, a.ancestry) ==
         std::tie(b.traversal_ordinal, b.item_ordinal, b.parent_local_ordinal,
                  b.item, b.ancestry);
}

bool operator==(const CanonicalTargetSchedule &a,
                const CanonicalTargetSchedule &b) noexcept {
  return std::tie(a.target, a.body_root, a.entry_state, a.exit_state,
                  a.entry_exit_included, a.records) ==
         std::tie(b.target, b.body_root, b.entry_state, b.exit_state,
                  b.entry_exit_included, b.records);
}

bool CanonicalScheduleBundle::operator==(
    const CanonicalScheduleBundle &other) const noexcept {
  return targets == other.targets;
}

namespace {

constexpr std::size_t max_ancestry_depth = 256;
constexpr std::uint32_t provenance_schema_version = 1;

struct TableHandle final {
  int handle{-1};
  TableHandle() = default;
  explicit TableHandle(const int value) : handle(value) {}
  ~TableHandle() {
    if (handle >= 0)
      (void)Util_TableDestroy(handle);
  }
  TableHandle(const TableHandle &) = delete;
  TableHandle &operator=(const TableHandle &) = delete;
  TableHandle(TableHandle &&other) noexcept : handle(other.handle) {
    other.handle = -1;
  }
  TableHandle &operator=(TableHandle &&other) noexcept {
    if (this != &other) {
      if (handle >= 0)
        (void)Util_TableDestroy(handle);
      handle = other.handle;
      other.handle = -1;
    }
    return *this;
  }
};

struct TableIterator final {
  int handle{-1};
  explicit TableIterator(const int value) : handle(value) {}
  ~TableIterator() {
    if (handle >= 0)
      (void)Util_TableItDestroy(handle);
  }
  TableIterator(const TableIterator &) = delete;
  TableIterator &operator=(const TableIterator &) = delete;
};

[[nodiscard]] ScheduleCertificationFailure
make_failure(const ScheduleCertificationErrorCode code,
             const std::optional<SubcyclingScheduleTarget> target,
             const std::optional<std::uint64_t> ordinal,
             std::string field, std::string detail) {
  return {code, target, ordinal, std::move(field), std::move(detail)};
}

[[nodiscard]] std::string copy_string(const char *const value) {
  return value == nullptr ? std::string{} : std::string(value);
}

[[nodiscard]] std::vector<char> copy_buffer(const char *const value) {
  if (value == nullptr)
    return {'\0'};
  const auto length = std::strlen(value);
  std::vector<char> result(length + 1);
  std::memcpy(result.data(), value, length + 1);
  return result;
}

[[nodiscard]] std::string ascii_lower(std::string value) {
  for (char &c : value) {
    if (c >= 'A' && c <= 'Z')
      c = static_cast<char>(c - 'A' + 'a');
  }
  return value;
}

[[nodiscard]] bool is_hex_string(const std::string &value,
                                 const std::size_t expected_size) {
  if (value.size() != expected_size)
    return false;
  return std::all_of(value.begin(), value.end(), [](const char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') ||
           (c >= 'A' && c <= 'F');
  });
}

[[nodiscard]] bool normalize_provenance(
    const ScheduleBuildProvenance &input, ScheduleBuildProvenance &output,
    std::string &field, std::string &detail) {
  output = input;
  if (output.schema_version != provenance_schema_version) {
    field = "provenance.schema_version";
    detail = "unsupported provenance schema";
    return false;
  }
  if (output.inventory_abi_version != CCTK_SCHEDULE_INVENTORY_ABI_VERSION) {
    field = "provenance.inventory_abi_version";
    detail = "schedule inventory ABI does not match this build";
    return false;
  }
  const auto normalize_git = [&](std::string &value,
                                 const char *const name) -> bool {
    if (!is_hex_string(value, 40)) {
      field = name;
      detail = "expected exactly 40 hexadecimal characters";
      return false;
    }
    value = ascii_lower(std::move(value));
    return true;
  };
  const auto normalize_sha = [&](std::string &value,
                                 const char *const name) -> bool {
    if (!is_hex_string(value, 64)) {
      field = name;
      detail = "expected exactly 64 hexadecimal characters";
      return false;
    }
    value = ascii_lower(std::move(value));
    return true;
  };
  if (!normalize_git(output.cactus_flesh_parent,
                     "provenance.cactus_flesh_parent") ||
      !normalize_sha(output.cactus_inventory_patch_sha256,
                     "provenance.cactus_inventory_patch_sha256") ||
      !normalize_git(output.carpetx_parent, "provenance.carpetx_parent"))
    return false;
  if (output.ordered_carpetx_patches.empty()) {
    field = "provenance.ordered_carpetx_patches";
    detail = "ordered overlay list must not be empty";
    return false;
  }
  for (std::size_t i = 0; i < output.ordered_carpetx_patches.size(); ++i) {
    auto &patch = output.ordered_carpetx_patches[i];
    if (patch.phase.empty()) {
      field = "provenance.ordered_carpetx_patches[" + std::to_string(i) +
              "].phase";
      detail = "phase must not be empty";
      return false;
    }
    if (!normalize_sha(patch.sha256,
                       ("provenance.ordered_carpetx_patches[" +
                        std::to_string(i) + "].sha256")
                           .c_str()))
      return false;
  }
  if (output.compiler_id.empty() || output.compiler_version.empty() ||
      output.target_triple.empty() || output.cxx_abi.empty() ||
      output.byte_order.empty()) {
    field = "provenance.compiler_abi";
    detail = "compiler, target, ABI, and byte order fields must be nonempty";
    return false;
  }
  if (output.byte_order != "little" && output.byte_order != "big") {
    field = "provenance.byte_order";
    detail = "byte order must be 'little' or 'big'";
    return false;
  }
  if (output.cctk_type_sizes.empty()) {
    field = "provenance.cctk_type_sizes";
    detail = "CCTK type-size tuple must not be empty";
    return false;
  }
  std::map<int, int> type_sizes;
  for (std::size_t i = 0; i < output.cctk_type_sizes.size(); ++i) {
    const auto &entry = output.cctk_type_sizes[i];
    if (entry.second <= 0 || !type_sizes.emplace(entry).second) {
      field = "provenance.cctk_type_sizes";
      detail = "type codes must be unique and sizes positive";
      return false;
    }
    const int live_size = CCTK_VarTypeSize(entry.first);
    if (live_size <= 0 || entry.second != live_size) {
      field = "provenance.cctk_type_sizes[" + std::to_string(i) + "].size";
      detail = "type size does not match CCTK_VarTypeSize in this build";
      return false;
    }
  }
  if (!normalize_sha(output.executable_sha256,
                     "provenance.executable_sha256"))
    return false;
  return true;
}

[[nodiscard]] bool valid_count_and_pointer(const int count, const void *pointer) {
  return count >= 0 && (count == 0 || pointer != nullptr);
}

[[nodiscard]] CanonicalExecutionFlags
execution_flags(const cFunctionData &data) {
  return {data.meta,
          data.meta_early,
          data.meta_late,
          data.global,
          data.global_early,
          data.global_late,
          data.level,
          data.singlemap,
          data.local,
          data.loop_meta,
          data.loop_global,
          data.loop_level,
          data.loop_singlemap,
          data.loop_local};
}

[[nodiscard]] bool validate_local_scope(
    const cFunctionData &data, CanonicalScheduleMode &mode,
    ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &prefix) {
  const auto flags = execution_flags(data);
  if (flags.meta != 0 || flags.meta_early != 0 || flags.meta_late != 0 ||
      flags.global != 0 || flags.global_early != 0 ||
      flags.global_late != 0 || flags.level != 0 || flags.singlemap != 0 ||
      flags.local < 0 || flags.local > 1 || flags.loop_meta != 0 ||
      flags.loop_global != 0 || flags.loop_level != 0 ||
      flags.loop_singlemap != 0 || flags.loop_local != 0) {
    failure = make_failure(ScheduleCertificationErrorCode::unsupported_metadata,
                           target, ordinal, prefix + ".execution_mode",
                           "only default-local or explicit local is supported");
    return false;
  }
  mode = data.local == 0 ? CanonicalScheduleMode::default_local
                         : CanonicalScheduleMode::local;
  return true;
}

[[nodiscard]] bool canonicalize_table(
    const int source_handle, const bool clone_leaf,
    std::vector<CanonicalTag> &tags, TableHandle &owned_clone,
    ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &prefix) {
  const int flags = Util_TableQueryFlags(source_handle);
  if (flags != UTIL_TABLE_FLAGS_CASE_INSENSITIVE) {
    failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                           ordinal, prefix + ".tags.flags",
                           "tag table must be exactly case-insensitive");
    return false;
  }
  const int key_count = Util_TableQueryNKeys(source_handle);
  const int max_key_length = Util_TableQueryMaxKeyLength(source_handle);
  if (key_count < 0 || max_key_length < 0) {
    failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                           ordinal, prefix + ".tags.shape",
                           "failed to query tag table shape");
    return false;
  }
  TableIterator iterator(Util_TableItCreate(source_handle));
  if (iterator.handle < 0) {
    failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                           ordinal, prefix + ".tags.iterator",
                           "failed to create tag-table iterator");
    return false;
  }
  tags.clear();
  tags.reserve(static_cast<std::size_t>(key_count));
  std::map<std::string, bool> normalized_keys;
  for (int i = 0; i < key_count; ++i) {
    std::vector<char> key(static_cast<std::size_t>(max_key_length) + 1U, '\0');
    CCTK_INT type_code = 0;
    CCTK_INT element_count = 0;
    const int key_length = Util_TableItQueryKeyValueInfo(
        iterator.handle, static_cast<int>(key.size()), key.data(), &type_code,
        &element_count);
    if (key_length < 0) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.iterator",
                             "failed to query tag-table iterator entry");
      return false;
    }
    if (key_length > max_key_length || element_count <= 0) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.entry",
                             "malformed tag-table entry");
      return false;
    }
    std::string normalized = ascii_lower(std::string(key.data()));
    if (!normalized_keys.emplace(normalized, true).second) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.key",
                             "duplicate normalized tag key");
      return false;
    }
    const int type = static_cast<int>(type_code);
    if (type == CCTK_VARIABLE_STRING || type == CCTK_VARIABLE_POINTER ||
        type == CCTK_VARIABLE_POINTER_TO_CONST ||
        type == CCTK_VARIABLE_FPOINTER) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.type",
                             "pointer-valued tags are unsupported");
      return false;
    }
    const int type_size = CCTK_VarTypeSize(type);
    if (type_size <= 0) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.type",
                             "unknown or zero-size tag type");
      return false;
    }
    const auto count = static_cast<std::uint64_t>(element_count);
    if (count > std::numeric_limits<std::size_t>::max() /
                    static_cast<std::size_t>(type_size)) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.size",
                             "tag byte size overflows size_t");
      return false;
    }
    std::vector<std::uint8_t> bytes(
        static_cast<std::size_t>(count) * static_cast<std::size_t>(type_size));
    const int got = Util_TableGetGenericArray(source_handle, type,
                                              static_cast<int>(element_count),
                                              bytes.data(), key.data());
    if (got != static_cast<int>(element_count)) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.value",
                             "tag element count changed while reading");
      return false;
    }
    tags.push_back({std::move(normalized), type, count, std::move(bytes)});
    const int advance = Util_TableItAdvance(iterator.handle);
    const int expected_advance = i + 1 < key_count ? 1 : 0;
    if (advance != expected_advance) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.iterator",
                             "iterator cardinality differs from key count");
      return false;
    }
  }
  if (Util_TableItQueryIsNull(iterator.handle) != 1) {
    failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                           ordinal, prefix + ".tags.iterator",
                           "iterator did not end after exact key count");
    return false;
  }
  std::sort(tags.begin(), tags.end(), [](const CanonicalTag &a,
                                         const CanonicalTag &b) {
    return a.key < b.key;
  });
  if (clone_leaf) {
    const int clone = Util_TableClone(source_handle);
    if (clone < 0) {
      failure = make_failure(ScheduleCertificationErrorCode::tag_error, target,
                             ordinal, prefix + ".tags.clone",
                             "failed to clone leaf tag table");
      return false;
    }
    owned_clone = TableHandle(clone);
  }
  return true;
}

[[nodiscard]] bool canonical_group_name(
    const int index, std::string &name, ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &field) {
  const char *const raw = CCTK_FullGroupName(index);
  if (raw == nullptr) {
    failure = make_failure(ScheduleCertificationErrorCode::name_resolution_error,
                           target, ordinal, field,
                           "group index has no fully qualified name");
    return false;
  }
  name = raw;
  if (CCTK_GroupIndex(name.c_str()) != index) {
    failure = make_failure(ScheduleCertificationErrorCode::name_resolution_error,
                           target, ordinal, field,
                           "group full-name roundtrip failed");
    return false;
  }
  return true;
}

[[nodiscard]] bool canonical_variable_name(
    const int index, std::string &name, ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &field) {
  const char *const raw = CCTK_FullVarName(index);
  if (raw == nullptr) {
    failure = make_failure(ScheduleCertificationErrorCode::name_resolution_error,
                           target, ordinal, field,
                           "variable index has no fully qualified name");
    return false;
  }
  name = raw;
  if (CCTK_VarIndex(name.c_str()) != index) {
    failure = make_failure(ScheduleCertificationErrorCode::name_resolution_error,
                           target, ordinal, field,
                           "variable full-name roundtrip failed");
    return false;
  }
  return true;
}

[[nodiscard]] bool canonicalize_accesses(
    const cFunctionData &data, std::vector<CanonicalAccess> &accesses,
    ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &prefix) {
  if (!valid_count_and_pointer(data.n_RDWR, data.RDWR)) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".rdwr",
                           "invalid RDWR count/pointer pair");
    return false;
  }
  using AccessKey = std::pair<std::string, int>;
  std::map<AccessKey, CanonicalAccess> merged;
  constexpr int valid_regions = CCTK_VALID_GHOSTS | CCTK_VALID_BOUNDARY |
                                CCTK_VALID_INTERIOR;
  for (int i = 0; i < data.n_RDWR; ++i) {
    const auto &entry = data.RDWR[i];
    if (entry.timelevel != 0 || entry.where_rd < 0 || entry.where_wr < 0 ||
        entry.where_inv < 0 ||
        ((entry.where_rd | entry.where_wr | entry.where_inv) &
         ~valid_regions) != 0) {
      failure = make_failure(ScheduleCertificationErrorCode::unsupported_metadata,
                             target, ordinal, prefix + ".accesses",
                             "only timelevel-zero standard region masks are supported");
      return false;
    }
    std::string variable_name;
    if (!canonical_variable_name(entry.varindex, variable_name, failure, target,
                                 ordinal, prefix + ".accesses.variable_name"))
      return false;
    const int group_index = CCTK_GroupIndexFromVarI(entry.varindex);
    std::string group_name;
    if (group_index < 0 ||
        !canonical_group_name(group_index, group_name, failure, target, ordinal,
                              prefix + ".accesses.group_name"))
      return false;
    cGroup group{};
    if (CCTK_GroupData(group_index, &group) != 0 ||
        group.grouptype != CCTK_GF || group.vartype != CCTK_VARIABLE_REAL ||
        CCTK_VarTypeI(entry.varindex) != CCTK_VARIABLE_REAL ||
        CCTK_DeclaredTimeLevelsVI(entry.varindex) <= 0) {
      failure = make_failure(ScheduleCertificationErrorCode::unsupported_metadata,
                             target, ordinal, prefix + ".accesses.variable",
                             "only GF REAL variables with declared TL0 are supported");
      return false;
    }
    const AccessKey key{variable_name, entry.timelevel};
    auto found = merged.find(key);
    if (found == merged.end()) {
      found = merged
                  .emplace(key, CanonicalAccess{variable_name, entry.timelevel,
                                                0, 0, 0})
                  .first;
    }
    found->second.read_mask |= entry.where_rd;
    found->second.write_mask |= entry.where_wr;
    found->second.invalidate_mask |= entry.where_inv;
  }
  accesses.clear();
  for (auto &entry : merged)
    accesses.push_back(std::move(entry.second));
  const bool normalized_reads = std::any_of(
      accesses.begin(), accesses.end(),
      [](const CanonicalAccess &a) { return a.read_mask != 0; });
  const bool normalized_writes = std::any_of(
      accesses.begin(), accesses.end(),
      [](const CanonicalAccess &a) { return a.write_mask != 0; });
  const bool normalized_invalidates = std::any_of(
      accesses.begin(), accesses.end(),
      [](const CanonicalAccess &a) { return a.invalidate_mask != 0; });
  if ((data.n_ReadsClauses > 0) != normalized_reads ||
      (data.n_WritesClauses > 0) != normalized_writes ||
      (data.n_InvalidatesClauses > 0) != normalized_invalidates) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".accesses.raw_normalized",
                           "raw clauses contradict normalized access masks");
    return false;
  }
  return true;
}

[[nodiscard]] bool copy_string_array(
    const int count, const char *const *const values,
    std::vector<std::string> &output, ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &field) {
  if (!valid_count_and_pointer(count, values)) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, field,
                           "invalid string-array count/pointer pair");
    return false;
  }
  output.clear();
  output.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    if (values[i] == nullptr) {
      failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                             target, ordinal, field,
                             "string-array element is null");
      return false;
    }
    output.emplace_back(values[i]);
  }
  return true;
}

[[nodiscard]] bool canonicalize_group_array(
    const int count, const int *const indices, std::vector<std::string> &output,
    ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &field) {
  if (!valid_count_and_pointer(count, indices)) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, field,
                           "invalid group-array count/pointer pair");
    return false;
  }
  output.clear();
  output.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    std::string name;
    if (!canonical_group_name(indices[i], name, failure, target, ordinal,
                              field + ".name"))
      return false;
    output.push_back(std::move(name));
  }
  return true;
}

[[nodiscard]] bool canonicalize_item(
    const char *const description, const char *const implementation,
    const cFunctionData *const function_data, const int n_storage_groups,
    const int *const storage_groups, const int *const storage_timelevels,
    const int n_communication_groups, const int *const communication_groups,
    const int n_ifs, const char *const *const ifs, const int n_whiles,
    const char *const *const whiles, const bool leaf,
    CanonicalScheduleItem &item, TableHandle &leaf_tag_clone,
    ScheduleCertificationFailure &failure,
    const SubcyclingScheduleTarget target, const std::uint64_t ordinal,
    const std::string &prefix) {
  if (function_data == nullptr) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".function_data",
                           "function data is null");
    return false;
  }
  const cFunctionData &data = *function_data;
  if (description == nullptr || implementation == nullptr ||
      data.where == nullptr || data.thorn == nullptr ||
      data.routine == nullptr) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".identity",
                           "schedule identity strings must not be null");
    return false;
  }
  item.description = copy_string(description);
  item.implementation = copy_string(implementation);
  item.where = copy_string(data.where);
  item.thorn = copy_string(data.thorn);
  item.routine = copy_string(data.routine);
  item.language = static_cast<int>(data.language);
  item.function_type = static_cast<int>(data.type);
  item.execution_flags = execution_flags(data);

  if (!canonicalize_table(data.tags, leaf, item.tags, leaf_tag_clone, failure,
                          target, ordinal, prefix))
    return false;

  if (!validate_local_scope(data, item.effective_mode, failure, target,
                            ordinal, prefix))
    return false;
  if (leaf &&
      (data.language != LangC || data.type != FunctionStandard)) {
    failure = make_failure(ScheduleCertificationErrorCode::unsupported_metadata,
                           target, ordinal, prefix + ".call_convention",
                           "leaf must be LangC and FunctionStandard");
    return false;
  }

  if (!valid_count_and_pointer(data.n_SyncGroups, data.SyncGroups) ||
      !valid_count_and_pointer(data.n_TriggerGroups, data.TriggerGroups)) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".hazard_groups",
                           "invalid sync/trigger count/pointer pair");
    return false;
  }
  if (!canonicalize_group_array(data.n_SyncGroups, data.SyncGroups,
                                item.sync_groups, failure, target, ordinal,
                                prefix + ".sync_groups") ||
      !canonicalize_group_array(data.n_TriggerGroups, data.TriggerGroups,
                                item.trigger_groups, failure, target, ordinal,
                                prefix + ".trigger_groups"))
    return false;

  if (!copy_string_array(data.n_ReadsClauses, data.ReadsClauses,
                         item.reads_clauses, failure, target, ordinal,
                         prefix + ".reads_clauses") ||
      !copy_string_array(data.n_WritesClauses, data.WritesClauses,
                         item.writes_clauses, failure, target, ordinal,
                         prefix + ".writes_clauses") ||
      !copy_string_array(data.n_InvalidatesClauses,
                         data.InvalidatesClauses,
                         item.invalidates_clauses, failure, target, ordinal,
                         prefix + ".invalidates_clauses") ||
      !canonicalize_accesses(data, item.accesses, failure, target, ordinal,
                             prefix))
    return false;

  if (!valid_count_and_pointer(n_storage_groups, storage_groups) ||
      !valid_count_and_pointer(n_storage_groups, storage_timelevels)) {
    failure = make_failure(ScheduleCertificationErrorCode::malformed_record,
                           target, ordinal, prefix + ".storage_groups",
                           "invalid storage count/pointer pair");
    return false;
  }
  item.storage_groups.clear();
  item.storage_groups.reserve(static_cast<std::size_t>(n_storage_groups));
  for (int i = 0; i < n_storage_groups; ++i) {
    std::string name;
    if (!canonical_group_name(storage_groups[i], name, failure, target,
                              ordinal, prefix + ".storage_groups.name"))
      return false;
    item.storage_groups.push_back({std::move(name), storage_timelevels[i]});
  }
  if (!canonicalize_group_array(n_communication_groups, communication_groups,
                                item.communication_groups, failure, target,
                                ordinal, prefix + ".communication_groups") ||
      !copy_string_array(n_ifs, ifs, item.ifs, failure, target, ordinal,
                         prefix + ".ifs") ||
      !copy_string_array(n_whiles, whiles, item.whiles, failure, target,
                         ordinal, prefix + ".whiles"))
    return false;

  if (!item.storage_groups.empty() || !item.communication_groups.empty() ||
      !item.sync_groups.empty() || !item.trigger_groups.empty() ||
      !item.ifs.empty() || !item.whiles.empty()) {
    failure = make_failure(ScheduleCertificationErrorCode::unsupported_metadata,
                           target, ordinal, prefix + ".hazards",
                           "storage, communication, sync, trigger, IF, and WHILE are unsupported");
    return false;
  }
  return true;
}

struct OwnedFunctionData final {
  cFunctionData data{};
  void *opaque_call_handle{};
  TableHandle tags;
  std::vector<int> sync_groups;
  std::vector<int> trigger_groups;
  std::vector<RDWR_entry> rdwr;
  std::vector<std::vector<char>> reads_buffers;
  std::vector<std::vector<char>> writes_buffers;
  std::vector<std::vector<char>> invalidates_buffers;
  std::vector<const char *> reads_pointers;
  std::vector<const char *> writes_pointers;
  std::vector<const char *> invalidates_pointers;
  std::vector<char> where;
  std::vector<char> routine;
  std::vector<char> thorn;
  std::vector<int> storage_group_indices;
  std::vector<int> storage_timelevels;
  std::vector<int> communication_group_indices;
  std::vector<std::string> if_strings;
  std::vector<std::string> while_strings;
  CanonicalScheduleRecord canonical;

  OwnedFunctionData() = default;
  OwnedFunctionData(const OwnedFunctionData &) = delete;
  OwnedFunctionData &operator=(const OwnedFunctionData &) = delete;
  OwnedFunctionData(OwnedFunctionData &&) = delete;
  OwnedFunctionData &operator=(OwnedFunctionData &&) = delete;

  void capture(const cFunctionData &source, TableHandle source_tags,
               void *const handle, const cScheduleInventoryRecord &record) {
    data = source;
    // The call-convention bridge is metadata only in this phase. Preserve its
    // exact value with the rest of cFunctionData, but never invoke it.
    data.FortranCaller = source.FortranCaller;
    opaque_call_handle = handle;
    tags = std::move(source_tags);
    if (source.n_SyncGroups > 0)
      sync_groups.assign(source.SyncGroups,
                         source.SyncGroups + source.n_SyncGroups);
    if (source.n_TriggerGroups > 0)
      trigger_groups.assign(source.TriggerGroups,
                            source.TriggerGroups + source.n_TriggerGroups);
    if (source.n_RDWR > 0)
      rdwr.assign(source.RDWR, source.RDWR + source.n_RDWR);
    copy_clause_buffers(source.n_ReadsClauses, source.ReadsClauses,
                        reads_buffers, reads_pointers);
    copy_clause_buffers(source.n_WritesClauses, source.WritesClauses,
                        writes_buffers, writes_pointers);
    copy_clause_buffers(source.n_InvalidatesClauses,
                        source.InvalidatesClauses, invalidates_buffers,
                        invalidates_pointers);
    where = copy_buffer(source.where);
    routine = copy_buffer(source.routine);
    thorn = copy_buffer(source.thorn);
    if (record.n_storage_groups > 0) {
      storage_group_indices.assign(
          record.storage_groups,
          record.storage_groups + record.n_storage_groups);
      storage_timelevels.assign(
          record.storage_timelevels,
          record.storage_timelevels + record.n_storage_groups);
    }
    if (record.n_communication_groups > 0)
      communication_group_indices.assign(
          record.communication_groups,
          record.communication_groups + record.n_communication_groups);
    for (int i = 0; i < record.n_ifs; ++i)
      if_strings.emplace_back(record.ifs[i]);
    for (int i = 0; i < record.n_whiles; ++i)
      while_strings.emplace_back(record.whiles[i]);
    rebind();
  }

private:
  static void copy_clause_buffers(
      const int count, const char *const *const source,
      std::vector<std::vector<char>> &buffers,
      std::vector<const char *> &pointers) {
    buffers.clear();
    pointers.clear();
    buffers.reserve(static_cast<std::size_t>(count));
    pointers.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i)
      buffers.push_back(copy_buffer(source[i]));
    for (const auto &buffer : buffers)
      pointers.push_back(buffer.data());
  }

  void rebind() noexcept {
    data.SyncGroups = sync_groups.empty() ? nullptr : sync_groups.data();
    data.TriggerGroups =
        trigger_groups.empty() ? nullptr : trigger_groups.data();
    data.RDWR = rdwr.empty() ? nullptr : rdwr.data();
    data.ReadsClauses =
        reads_pointers.empty() ? nullptr : reads_pointers.data();
    data.WritesClauses =
        writes_pointers.empty() ? nullptr : writes_pointers.data();
    data.InvalidatesClauses =
        invalidates_pointers.empty() ? nullptr : invalidates_pointers.data();
    data.where = where.data();
    data.routine = routine.data();
    data.thorn = thorn.data();
    data.tags = tags.handle;
  }
};

struct CandidateTarget final {
  CanonicalTargetSchedule manifest;
  std::vector<std::unique_ptr<OwnedFunctionData>> functions;
};

using ObservedParentPath = std::vector<std::uint64_t>;

struct ObservedScheduleNode final {
  CanonicalScheduleItem item;
  ObservedParentPath parent_path;
  std::uint64_t parent_local_ordinal{};
  bool is_scope{};
};

struct CollectContext final {
  SubcyclingScheduleTarget target;
  std::string root;
  CandidateTarget *candidate{};
  std::optional<ScheduleCertificationFailure> failure;
  std::uint64_t next_traversal_ordinal{};
  std::optional<std::uint64_t> last_item_ordinal;
  // This flag is the only state touched from a catch block in the noexcept C
  // callback. Constructing the diagnostic is deliberately deferred until the
  // inventory call has returned to the ordinary C++ boundary.
  bool callback_exception{};
  std::optional<std::uint64_t> last_first_seen_item_ordinal;
  std::map<std::uint64_t, ObservedScheduleNode> nodes_by_item_ordinal;
  std::map<ObservedParentPath, std::map<std::uint64_t, std::uint64_t>>
      children_by_parent_path;
};

[[nodiscard]] bool validate_observed_node(
    CollectContext &context, const std::uint64_t item_ordinal,
    const std::uint64_t parent_local_ordinal,
    const CanonicalScheduleItem &item, const ObservedParentPath &parent_path,
    const bool is_scope, const std::uint64_t traversal,
    const std::string &prefix) {
  const auto known_node = context.nodes_by_item_ordinal.find(item_ordinal);
  if (known_node != context.nodes_by_item_ordinal.end()) {
    if (known_node->second.is_scope != is_scope) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          traversal, prefix + ".node_kind",
          "repeated item ordinal changed between scope and leaf");
      return false;
    }
    if (!(known_node->second.item == item)) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          traversal, prefix + ".item_identity",
          "repeated item ordinal changed canonical item metadata");
      return false;
    }
    if (known_node->second.parent_path != parent_path) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          traversal, prefix + ".parent_path",
          "repeated item ordinal changed its observed parent path");
      return false;
    }
    if (known_node->second.parent_local_ordinal != parent_local_ordinal) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          traversal, prefix + ".parent_local_ordinal",
          "repeated item ordinal changed its sibling-local ordinal");
      return false;
    }
  } else if (context.last_first_seen_item_ordinal.has_value() &&
             item_ordinal <= *context.last_first_seen_item_ordinal) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::malformed_record, context.target,
        traversal, prefix + ".item_ordinal",
        "first-seen item ordinals must be strictly increasing");
    return false;
  }

  auto &children = context.children_by_parent_path[parent_path];
  const auto known_slot = children.find(parent_local_ordinal);
  if (known_slot != children.end() && known_slot->second != item_ordinal) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::malformed_record, context.target,
        traversal, prefix + ".parent_local_ordinal",
        "parent path and sibling-local ordinal identify multiple items");
    return false;
  }
  if (known_slot == children.end()) {
    const auto expected_local_ordinal =
        static_cast<std::uint64_t>(children.size());
    if (parent_local_ordinal != expected_local_ordinal) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          traversal, prefix + ".parent_local_ordinal",
          "new observed siblings must have contiguous local ordinals");
      return false;
    }
    children.emplace(parent_local_ordinal, item_ordinal);
  }
  if (known_node == context.nodes_by_item_ordinal.end())
    context.nodes_by_item_ordinal.emplace(
        item_ordinal,
        ObservedScheduleNode{item, parent_path, parent_local_ordinal,
                             is_scope});
  if (known_node == context.nodes_by_item_ordinal.end())
    context.last_first_seen_item_ordinal = item_ordinal;
  return true;
}

[[nodiscard]] bool validate_scope_shape(
    const cScheduleInventoryScope &scope, const std::uint64_t traversal,
    const SubcyclingScheduleTarget target, const std::size_t index,
    ScheduleCertificationFailure &failure) {
  if (scope.abi_version != CCTK_SCHEDULE_INVENTORY_ABI_VERSION ||
      scope.struct_size != sizeof(cScheduleInventoryScope)) {
    failure = make_failure(ScheduleCertificationErrorCode::abi_mismatch, target,
                           traversal,
                           "record.ancestry[" + std::to_string(index) +
                               "].abi",
                           "scope ABI version or size mismatch");
    return false;
  }
  return true;
}

[[nodiscard]] bool collect_one_record(CollectContext &context,
                                      const cScheduleInventoryRecord &record,
                                      void *const opaque_call_handle) {
  const auto ordinal = context.next_traversal_ordinal;
  if (record.abi_version != CCTK_SCHEDULE_INVENTORY_ABI_VERSION ||
      record.struct_size != sizeof(cScheduleInventoryRecord)) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::abi_mismatch, context.target, ordinal,
        "record.abi", "record ABI version or size mismatch");
    return false;
  }
  if (record.root_group == nullptr || context.root != record.root_group ||
      record.entry_exit_included != 0) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::malformed_record, context.target,
        ordinal, "record.root",
        "record root mismatch or ENTRY/EXIT was included");
    return false;
  }
  if (record.traversal_ordinal != ordinal ||
      (context.last_item_ordinal.has_value() &&
       record.item_ordinal <= *context.last_item_ordinal)) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::malformed_record, context.target,
        ordinal, "record.ordinals",
        "traversal ordinals must be contiguous and item ordinals increasing");
    return false;
  }
  if (record.n_enclosing_groups > max_ancestry_depth ||
      (record.n_enclosing_groups > 0 && record.enclosing_groups == nullptr)) {
    context.failure = make_failure(
        ScheduleCertificationErrorCode::malformed_record, context.target,
        ordinal, "record.ancestry", "invalid ancestry depth/pointer pair");
    return false;
  }

  auto owned = std::make_unique<OwnedFunctionData>();
  owned->canonical.traversal_ordinal = record.traversal_ordinal;
  owned->canonical.item_ordinal = record.item_ordinal;
  owned->canonical.parent_local_ordinal = record.parent_local_ordinal;
  TableHandle leaf_tag_clone;
  ScheduleCertificationFailure failure{};
  if (!canonicalize_item(
          record.description, record.implementation, record.function_data,
          record.n_storage_groups, record.storage_groups,
          record.storage_timelevels, record.n_communication_groups,
          record.communication_groups, record.n_ifs, record.ifs,
          record.n_whiles, record.whiles, true, owned->canonical.item,
          leaf_tag_clone, failure, context.target, ordinal, "record.item")) {
    context.failure = std::move(failure);
    return false;
  }
  std::optional<std::uint64_t> previous_scope_item;
  ObservedParentPath parent_path;
  owned->canonical.ancestry.reserve(record.n_enclosing_groups);
  for (std::size_t i = 0; i < record.n_enclosing_groups; ++i) {
    const auto &scope = record.enclosing_groups[i];
    if (!validate_scope_shape(scope, ordinal, context.target, i, failure)) {
      context.failure = std::move(failure);
      return false;
    }
    if ((previous_scope_item.has_value() &&
         scope.item_ordinal <= *previous_scope_item) ||
        scope.item_ordinal >= record.item_ordinal) {
      context.failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context.target,
          ordinal,
          "record.ancestry[" + std::to_string(i) + "].item_ordinal",
          "ancestry item ordinals are not outer-to-inner increasing");
      return false;
    }
    CanonicalScheduleScope canonical_scope;
    canonical_scope.item_ordinal = scope.item_ordinal;
    canonical_scope.parent_local_ordinal = scope.parent_local_ordinal;
    TableHandle no_clone;
    if (!canonicalize_item(
            scope.description, scope.implementation, scope.function_data,
            scope.n_storage_groups, scope.storage_groups,
            scope.storage_timelevels, scope.n_communication_groups,
            scope.communication_groups, scope.n_ifs, scope.ifs,
            scope.n_whiles, scope.whiles, false, canonical_scope.item,
            no_clone, failure, context.target, ordinal,
            "record.ancestry[" + std::to_string(i) + "].item")) {
      context.failure = std::move(failure);
      return false;
    }
    const std::string scope_prefix =
        "record.ancestry[" + std::to_string(i) + "]";
    if (!validate_observed_node(context, canonical_scope.item_ordinal,
                                canonical_scope.parent_local_ordinal,
                                canonical_scope.item, parent_path, true, ordinal,
                                scope_prefix))
      return false;
    previous_scope_item = scope.item_ordinal;
    parent_path.push_back(scope.item_ordinal);
    owned->canonical.ancestry.push_back(std::move(canonical_scope));
  }
  if (!validate_observed_node(
          context, owned->canonical.item_ordinal,
          owned->canonical.parent_local_ordinal, owned->canonical.item,
          parent_path, false, ordinal, "record"))
    return false;
  owned->capture(*record.function_data, std::move(leaf_tag_clone),
                 opaque_call_handle, record);
  context.last_item_ordinal = record.item_ordinal;
  ++context.next_traversal_ordinal;
  context.candidate->manifest.records.push_back(owned->canonical);
  context.candidate->functions.push_back(std::move(owned));
  return true;
}

extern "C" int certification_inventory_callback(
    const cScheduleInventoryRecord *const record, void *const opaque_call_handle,
    void *const user_data) noexcept {
  auto *const context = static_cast<CollectContext *>(user_data);
  if (context == nullptr)
    return 1;
  if (context->failure.has_value())
    return 1;
  try {
    if (record == nullptr) {
      context->failure = make_failure(
          ScheduleCertificationErrorCode::malformed_record, context->target,
          context->next_traversal_ordinal, "record", "record pointer is null");
      return 1;
    }
    return collect_one_record(*context, *record, opaque_call_handle) ? 0 : 1;
  } catch (...) {
    context->callback_exception = true;
    return 1;
  }
}

extern "C" int reject_special_callback(
    const cScheduleInventoryRecord *, void *, void *) noexcept {
  return 0;
}

[[nodiscard]] std::string first_string_vector_difference(
    const std::string &prefix, const std::vector<std::string> &expected,
    const std::vector<std::string> &observed) {
  if (expected.size() != observed.size())
    return prefix + ".size";
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (expected[i] != observed[i])
      return prefix + "[" + std::to_string(i) + "]";
  }
  return {};
}

[[nodiscard]] std::string first_flags_difference(
    const std::string &prefix, const CanonicalExecutionFlags &a,
    const CanonicalExecutionFlags &b) {
#define CARPETX_COMPARE_FLAG(field_name)                                       \
  if (a.field_name != b.field_name)                                            \
    return prefix + "." #field_name
  CARPETX_COMPARE_FLAG(meta);
  CARPETX_COMPARE_FLAG(meta_early);
  CARPETX_COMPARE_FLAG(meta_late);
  CARPETX_COMPARE_FLAG(global);
  CARPETX_COMPARE_FLAG(global_early);
  CARPETX_COMPARE_FLAG(global_late);
  CARPETX_COMPARE_FLAG(level);
  CARPETX_COMPARE_FLAG(singlemap);
  CARPETX_COMPARE_FLAG(local);
  CARPETX_COMPARE_FLAG(loop_meta);
  CARPETX_COMPARE_FLAG(loop_global);
  CARPETX_COMPARE_FLAG(loop_level);
  CARPETX_COMPARE_FLAG(loop_singlemap);
  CARPETX_COMPARE_FLAG(loop_local);
#undef CARPETX_COMPARE_FLAG
  return {};
}

[[nodiscard]] std::string first_tags_difference(
    const std::string &prefix, const std::vector<CanonicalTag> &a,
    const std::vector<CanonicalTag> &b) {
  if (a.size() != b.size())
    return prefix + ".size";
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto item = prefix + "[" + std::to_string(i) + "]";
    if (a[i].key != b[i].key)
      return item + ".key";
    if (a[i].type_code != b[i].type_code)
      return item + ".type_code";
    if (a[i].element_count != b[i].element_count)
      return item + ".element_count";
    if (a[i].value_bytes.size() != b[i].value_bytes.size())
      return item + ".value_bytes.size";
    for (std::size_t j = 0; j < a[i].value_bytes.size(); ++j) {
      if (a[i].value_bytes[j] != b[i].value_bytes[j])
        return item + ".value_bytes[" + std::to_string(j) + "]";
    }
  }
  return {};
}

[[nodiscard]] std::string first_accesses_difference(
    const std::string &prefix, const std::vector<CanonicalAccess> &a,
    const std::vector<CanonicalAccess> &b) {
  if (a.size() != b.size())
    return prefix + ".size";
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto item = prefix + "[" + std::to_string(i) + "]";
    if (a[i].variable_name != b[i].variable_name)
      return item + ".variable_name";
    if (a[i].timelevel != b[i].timelevel)
      return item + ".timelevel";
    if (a[i].read_mask != b[i].read_mask)
      return item + ".read_mask";
    if (a[i].write_mask != b[i].write_mask)
      return item + ".write_mask";
    if (a[i].invalidate_mask != b[i].invalidate_mask)
      return item + ".invalidate_mask";
  }
  return {};
}

[[nodiscard]] std::string first_storage_difference(
    const std::string &prefix,
    const std::vector<CanonicalGroupTimelevel> &a,
    const std::vector<CanonicalGroupTimelevel> &b) {
  if (a.size() != b.size())
    return prefix + ".size";
  for (std::size_t i = 0; i < a.size(); ++i) {
    const auto item = prefix + "[" + std::to_string(i) + "]";
    if (a[i].group != b[i].group)
      return item + ".group";
    if (a[i].timelevel != b[i].timelevel)
      return item + ".timelevel";
  }
  return {};
}

[[nodiscard]] std::string first_item_difference(
    const std::string &prefix, const CanonicalScheduleItem &a,
    const CanonicalScheduleItem &b) {
#define CARPETX_COMPARE_ITEM_FIELD(field_name)                                 \
  if (a.field_name != b.field_name)                                            \
    return prefix + "." #field_name
  CARPETX_COMPARE_ITEM_FIELD(description);
  CARPETX_COMPARE_ITEM_FIELD(implementation);
  CARPETX_COMPARE_ITEM_FIELD(where);
  CARPETX_COMPARE_ITEM_FIELD(thorn);
  CARPETX_COMPARE_ITEM_FIELD(routine);
  CARPETX_COMPARE_ITEM_FIELD(language);
  CARPETX_COMPARE_ITEM_FIELD(function_type);
#undef CARPETX_COMPARE_ITEM_FIELD
  if (const auto difference =
          first_flags_difference(prefix + ".execution_flags",
                                 a.execution_flags, b.execution_flags);
      !difference.empty())
    return difference;
  if (a.effective_mode != b.effective_mode)
    return prefix + ".effective_mode";
  if (const auto difference =
          first_tags_difference(prefix + ".tags", a.tags, b.tags);
      !difference.empty())
    return difference;
  if (const auto difference = first_string_vector_difference(
          prefix + ".reads_clauses", a.reads_clauses, b.reads_clauses);
      !difference.empty())
    return difference;
  if (const auto difference = first_string_vector_difference(
          prefix + ".writes_clauses", a.writes_clauses, b.writes_clauses);
      !difference.empty())
    return difference;
  if (const auto difference = first_string_vector_difference(
          prefix + ".invalidates_clauses", a.invalidates_clauses,
          b.invalidates_clauses);
      !difference.empty())
    return difference;
  if (const auto difference = first_accesses_difference(
          prefix + ".accesses", a.accesses, b.accesses);
      !difference.empty())
    return difference;
  if (const auto difference = first_storage_difference(
          prefix + ".storage_groups", a.storage_groups, b.storage_groups);
      !difference.empty())
    return difference;
  for (const auto &entry :
       {std::pair<const char *, const std::vector<std::string> *>{
            "communication_groups", &a.communication_groups},
        {"sync_groups", &a.sync_groups},
        {"trigger_groups", &a.trigger_groups}, {"ifs", &a.ifs},
        {"whiles", &a.whiles}}) {
    const std::vector<std::string> *other = nullptr;
    const std::string name = entry.first;
    if (name == "communication_groups")
      other = &b.communication_groups;
    else if (name == "sync_groups")
      other = &b.sync_groups;
    else if (name == "trigger_groups")
      other = &b.trigger_groups;
    else if (name == "ifs")
      other = &b.ifs;
    else
      other = &b.whiles;
    if (const auto difference = first_string_vector_difference(
            prefix + "." + name, *entry.second, *other);
        !difference.empty())
      return difference;
  }
  return {};
}

[[nodiscard]] std::string first_manifest_difference(
    const CanonicalScheduleBundle &expected,
    const CanonicalScheduleBundle &observed) {
  for (std::size_t ti = 0; ti < expected.targets.size(); ++ti) {
    const auto prefix = "manifest.targets[" + std::to_string(ti) + "]";
    const auto &a = expected.targets[ti];
    const auto &b = observed.targets[ti];
    if (a.target != b.target)
      return prefix + ".target";
    if (a.body_root != b.body_root)
      return prefix + ".body_root";
    if (a.entry_state != b.entry_state)
      return prefix + ".entry_state";
    if (a.exit_state != b.exit_state)
      return prefix + ".exit_state";
    if (a.entry_exit_included != b.entry_exit_included)
      return prefix + ".entry_exit_included";
    if (a.records.size() != b.records.size())
      return prefix + ".records.size";
    for (std::size_t ri = 0; ri < a.records.size(); ++ri) {
      const auto record_prefix =
          prefix + ".records[" + std::to_string(ri) + "]";
      const auto &ar = a.records[ri];
      const auto &br = b.records[ri];
      if (ar.traversal_ordinal != br.traversal_ordinal)
        return record_prefix + ".traversal_ordinal";
      if (ar.item_ordinal != br.item_ordinal)
        return record_prefix + ".item_ordinal";
      if (ar.parent_local_ordinal != br.parent_local_ordinal)
        return record_prefix + ".parent_local_ordinal";
      if (const auto difference = first_item_difference(
              record_prefix + ".item", ar.item, br.item);
          !difference.empty())
        return difference;
      if (ar.ancestry.size() != br.ancestry.size())
        return record_prefix + ".ancestry.size";
      for (std::size_t si = 0; si < ar.ancestry.size(); ++si) {
        const auto scope_prefix =
            record_prefix + ".ancestry[" + std::to_string(si) + "]";
        if (ar.ancestry[si].item_ordinal != br.ancestry[si].item_ordinal)
          return scope_prefix + ".item_ordinal";
        if (ar.ancestry[si].parent_local_ordinal !=
            br.ancestry[si].parent_local_ordinal)
          return scope_prefix + ".parent_local_ordinal";
        if (const auto difference = first_item_difference(
                scope_prefix + ".item", ar.ancestry[si].item,
                br.ancestry[si].item);
            !difference.empty())
          return difference;
      }
    }
  }
  return {};
}

[[nodiscard]] bool compare_manifest(const CanonicalScheduleBundle &expected,
                                    const CanonicalScheduleBundle &observed,
                                    std::string &field) {
  if (expected == observed) {
    field.clear();
    return true;
  }
  field = first_manifest_difference(expected, observed);
  if (field.empty())
    field = "manifest.unknown";
  return false;
}

[[nodiscard]] std::string first_provenance_difference(
    const ScheduleBuildProvenance &a, const ScheduleBuildProvenance &b) {
#define CARPETX_COMPARE_PROVENANCE_FIELD(field_name)                           \
  if (a.field_name != b.field_name)                                            \
    return "provenance." #field_name
  CARPETX_COMPARE_PROVENANCE_FIELD(schema_version);
  CARPETX_COMPARE_PROVENANCE_FIELD(inventory_abi_version);
  CARPETX_COMPARE_PROVENANCE_FIELD(cactus_flesh_parent);
  CARPETX_COMPARE_PROVENANCE_FIELD(cactus_inventory_patch_sha256);
  CARPETX_COMPARE_PROVENANCE_FIELD(carpetx_parent);
  if (a.ordered_carpetx_patches.size() != b.ordered_carpetx_patches.size())
    return "provenance.ordered_carpetx_patches.size";
  for (std::size_t i = 0; i < a.ordered_carpetx_patches.size(); ++i) {
    const auto prefix = "provenance.ordered_carpetx_patches[" +
                        std::to_string(i) + "]";
    if (a.ordered_carpetx_patches[i].phase !=
        b.ordered_carpetx_patches[i].phase)
      return prefix + ".phase";
    if (a.ordered_carpetx_patches[i].sha256 !=
        b.ordered_carpetx_patches[i].sha256)
      return prefix + ".sha256";
  }
  CARPETX_COMPARE_PROVENANCE_FIELD(compiler_id);
  CARPETX_COMPARE_PROVENANCE_FIELD(compiler_version);
  CARPETX_COMPARE_PROVENANCE_FIELD(target_triple);
  CARPETX_COMPARE_PROVENANCE_FIELD(cxx_abi);
  CARPETX_COMPARE_PROVENANCE_FIELD(byte_order);
  if (a.cctk_type_sizes.size() != b.cctk_type_sizes.size())
    return "provenance.cctk_type_sizes.size";
  for (std::size_t i = 0; i < a.cctk_type_sizes.size(); ++i) {
    if (a.cctk_type_sizes[i].first != b.cctk_type_sizes[i].first)
      return "provenance.cctk_type_sizes[" + std::to_string(i) + "].type";
    if (a.cctk_type_sizes[i].second != b.cctk_type_sizes[i].second)
      return "provenance.cctk_type_sizes[" + std::to_string(i) + "].size";
  }
  CARPETX_COMPARE_PROVENANCE_FIELD(executable_sha256);
#undef CARPETX_COMPARE_PROVENANCE_FIELD
  return "provenance.unknown";
}

[[nodiscard]] bool validate_and_compare_provenance(
    const ScheduleBuildProvenance &expected,
    const ScheduleBuildProvenance &observed,
    ScheduleBuildProvenance &normalized_observed, std::string &field,
    std::string &detail) {
  ScheduleBuildProvenance normalized_expected;
  if (!normalize_provenance(expected, normalized_expected, field, detail) ||
      !normalize_provenance(observed, normalized_observed, field, detail))
    return false;
  if (!(normalized_expected == normalized_observed)) {
    field = first_provenance_difference(normalized_expected,
                                        normalized_observed);
    detail = "expected and observed build provenance differ";
    return false;
  }
  return true;
}

[[nodiscard]] bool validate_provenance_tag_type_coverage(
    const ScheduleBuildProvenance &provenance,
    const CanonicalScheduleBundle &manifest, std::string &field,
    std::string &detail) {
  std::map<int, bool> described_types;
  for (const auto &entry : provenance.cctk_type_sizes)
    described_types.emplace(entry.first, true);
  std::map<int, bool> used_types;
  const auto observe_item = [&](const CanonicalScheduleItem &item) {
    for (const auto &tag : item.tags)
      used_types.emplace(tag.type_code, true);
  };
  for (const auto &target : manifest.targets) {
    for (const auto &record : target.records) {
      observe_item(record.item);
      for (const auto &scope : record.ancestry)
        observe_item(scope.item);
    }
  }
  for (const auto &entry : used_types) {
    if (described_types.find(entry.first) == described_types.end()) {
      field = "provenance.cctk_type_sizes.missing[" +
              std::to_string(entry.first) + "]";
      detail = "canonical tag type is absent from build provenance";
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::optional<ScheduleCertificationFailure>
collect_native_gate_observation(
    const ScheduleBuildProvenance &normalized_provenance,
    CanonicalScheduleBundle &observed_manifest) {
  constexpr std::array<SubcyclingScheduleTarget, 2> target_ids{
      SubcyclingScheduleTarget::rhs, SubcyclingScheduleTarget::post_step};
  constexpr std::array<const char *, 2> roots{"ODESolvers_RHS",
                                               "ODESolvers_PostStep"};
  std::array<CandidateTarget, 2> candidates;
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    const auto target = target_ids[i];
    const std::string root = roots[i];
    candidates[i].manifest.target = target;
    candidates[i].manifest.body_root = root;
    candidates[i].manifest.entry_state = SpecialScheduleGroupState::missing;
    candidates[i].manifest.exit_state = SpecialScheduleGroupState::missing;
    candidates[i].manifest.entry_exit_included = false;

    const std::string entry = root + "$ENTRY";
    const int entry_result = CCTK_ScheduleInventoryDry(
        entry.c_str(), reject_special_callback, nullptr);
    if (entry_result == 0)
      return make_failure(
          ScheduleCertificationErrorCode::special_group_present, target,
          std::nullopt, "entry_state",
          "ENTRY group exists; missing named group was required");
    if (entry_result != CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND)
      return make_failure(
          ScheduleCertificationErrorCode::inventory_error, target,
          std::nullopt, "entry_inventory",
          "ENTRY inventory did not return named-missing");

    CollectContext context{target, root, &candidates[i], std::nullopt, 0,
                           std::nullopt, false, std::nullopt, {}, {}};
    const int body_result = CCTK_ScheduleInventoryDry(
        root.c_str(), certification_inventory_callback, &context);
    if (context.callback_exception)
      return make_failure(
          ScheduleCertificationErrorCode::callback_internal_error, target,
          context.next_traversal_ordinal, "record.exception",
          "exception contained while copying borrowed inventory metadata");
    if (context.failure.has_value())
      return std::move(context.failure);
    if (body_result != 0)
      return make_failure(
          ScheduleCertificationErrorCode::inventory_error, target,
          std::nullopt, "body_inventory",
          "body inventory failed or body group is missing");

    const std::string exit = root + "$EXIT";
    const int exit_result = CCTK_ScheduleInventoryDry(
        exit.c_str(), reject_special_callback, nullptr);
    if (exit_result == 0)
      return make_failure(
          ScheduleCertificationErrorCode::special_group_present, target,
          std::nullopt, "exit_state",
          "EXIT group exists; missing named group was required");
    if (exit_result != CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND)
      return make_failure(
          ScheduleCertificationErrorCode::inventory_error, target,
          std::nullopt, "exit_inventory",
          "EXIT inventory did not return named-missing");
  }

  for (std::size_t i = 0; i < candidates.size(); ++i)
    observed_manifest.targets[i] = candidates[i].manifest;
  std::string field;
  std::string detail;
  if (!validate_provenance_tag_type_coverage(
          normalized_provenance, observed_manifest, field, detail))
    return make_failure(ScheduleCertificationErrorCode::provenance_mismatch,
                        std::nullopt, std::nullopt, std::move(field),
                        std::move(detail));
  return std::nullopt;
}

} // namespace

struct CertifiedLocalScheduleRegistry::Storage final {
  CanonicalScheduleBundle manifest;
  ScheduleBuildProvenance provenance;
  std::array<std::vector<std::unique_ptr<OwnedFunctionData>>, 2> functions;
};

CertifiedLocalScheduleRegistry::CertifiedLocalScheduleRegistry(
    std::unique_ptr<Storage> storage)
    : storage_(std::move(storage)) {}

CertifiedLocalScheduleRegistry::~CertifiedLocalScheduleRegistry() = default;

const CanonicalScheduleBundle &
CertifiedLocalScheduleRegistry::manifest() const noexcept {
  return storage_->manifest;
}

const ScheduleBuildProvenance &
CertifiedLocalScheduleRegistry::provenance() const noexcept {
  return storage_->provenance;
}

std::size_t CertifiedLocalScheduleRegistry::size(
    const SubcyclingScheduleTarget target) const noexcept {
  const auto index = target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  return storage_->functions[index].size();
}

const CanonicalScheduleRecord &
CertifiedLocalScheduleRegistry::record_for_executor(
    const SubcyclingScheduleTarget target, const std::size_t ordinal) const {
  const auto index = target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  return storage_->functions.at(index).at(ordinal)->canonical;
}

int CertifiedLocalScheduleRegistry::invoke_for_executor(
    const SubcyclingScheduleTarget target, const std::size_t ordinal,
    cGH *const scratch_gh) const {
#ifdef CARPETX_SUBCYCLING_SCHEDULE_CERTIFICATION_STANDALONE
  (void)target;
  (void)ordinal;
  (void)scratch_gh;
  throw std::logic_error(
      "schedule invocation is unavailable in certification-only tests");
#else
  const auto index = target == SubcyclingScheduleTarget::rhs ? 0U : 1U;
  auto &function = *storage_->functions.at(index).at(ordinal);
  return CCTK_CallFunction(function.opaque_call_handle, &function.data,
                           scratch_gh);
#endif
}

NativeGateScheduleObservationResult
observe_local_subcycling_schedules_for_native_gate(
    const ScheduleBuildProvenance &observed_provenance) {
  ScheduleBuildProvenance normalized;
  std::string field;
  std::string detail;
  if (!normalize_provenance(observed_provenance, normalized, field, detail))
    return {std::nullopt,
            make_failure(ScheduleCertificationErrorCode::provenance_mismatch,
                         std::nullopt, std::nullopt, std::move(field),
                         std::move(detail))};

  CanonicalScheduleBundle manifest;
  if (auto failure =
          collect_native_gate_observation(normalized, manifest))
    return {std::nullopt, std::move(failure)};
  return {ScheduleCertificationExpectation{std::move(manifest),
                                            std::move(normalized)},
          std::nullopt};
}

ScheduleCertificationResult certify_local_subcycling_schedules(
    const ScheduleCertificationExpectation &expected,
    const ScheduleBuildProvenance &observed_provenance) {
  ScheduleBuildProvenance normalized_observed;
  std::string provenance_field;
  std::string provenance_detail;
  if (!validate_and_compare_provenance(
          expected.provenance, observed_provenance, normalized_observed,
          provenance_field, provenance_detail)) {
    return {nullptr,
            make_failure(ScheduleCertificationErrorCode::provenance_mismatch,
                         std::nullopt, std::nullopt,
                         std::move(provenance_field),
                         std::move(provenance_detail))};
  }

  constexpr std::array<SubcyclingScheduleTarget, 2> target_ids{
      SubcyclingScheduleTarget::rhs, SubcyclingScheduleTarget::post_step};
  constexpr std::array<const char *, 2> roots{"ODESolvers_RHS",
                                               "ODESolvers_PostStep"};
  std::array<CandidateTarget, 2> candidates;
  for (std::size_t i = 0; i < candidates.size(); ++i) {
    const auto target = target_ids[i];
    const std::string root = roots[i];
    candidates[i].manifest.target = target;
    candidates[i].manifest.body_root = root;
    candidates[i].manifest.entry_state = SpecialScheduleGroupState::missing;
    candidates[i].manifest.exit_state = SpecialScheduleGroupState::missing;
    candidates[i].manifest.entry_exit_included = false;

    const std::string entry = root + "$ENTRY";
    const int entry_result = CCTK_ScheduleInventoryDry(
        entry.c_str(), reject_special_callback, nullptr);
    if (entry_result == 0) {
      return {nullptr,
              make_failure(
                  ScheduleCertificationErrorCode::special_group_present,
                  target, std::nullopt, "entry_state",
                  "ENTRY group exists; missing named group was required")};
    }
    if (entry_result != CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND) {
      return {nullptr,
              make_failure(ScheduleCertificationErrorCode::inventory_error,
                           target, std::nullopt, "entry_inventory",
                           "ENTRY inventory did not return named-missing")};
    }

    CollectContext context{target, root, &candidates[i], std::nullopt, 0,
                           std::nullopt, false, std::nullopt, {}, {}};
    const int body_result = CCTK_ScheduleInventoryDry(
        root.c_str(), certification_inventory_callback, &context);
    if (context.callback_exception) {
      return {nullptr,
              make_failure(
                  ScheduleCertificationErrorCode::callback_internal_error,
                  target,
                  context.next_traversal_ordinal, "record.exception",
                  "exception contained while copying borrowed inventory metadata")};
    }
    if (context.failure.has_value())
      return {nullptr, std::move(context.failure)};
    if (body_result != 0) {
      return {nullptr,
              make_failure(ScheduleCertificationErrorCode::inventory_error,
                           target, std::nullopt, "body_inventory",
                           "body inventory failed or body group is missing")};
    }

    const std::string exit = root + "$EXIT";
    const int exit_result = CCTK_ScheduleInventoryDry(
        exit.c_str(), reject_special_callback, nullptr);
    if (exit_result == 0) {
      return {nullptr,
              make_failure(
                  ScheduleCertificationErrorCode::special_group_present,
                  target, std::nullopt, "exit_state",
                  "EXIT group exists; missing named group was required")};
    }
    if (exit_result != CCTK_SCHEDULE_INVENTORY_GROUP_NOT_FOUND) {
      return {nullptr,
              make_failure(ScheduleCertificationErrorCode::inventory_error,
                           target, std::nullopt, "exit_inventory",
                           "EXIT inventory did not return named-missing")};
    }
  }

  CanonicalScheduleBundle observed_manifest;
  for (std::size_t i = 0; i < candidates.size(); ++i)
    observed_manifest.targets[i] = candidates[i].manifest;
  if (!validate_provenance_tag_type_coverage(
          normalized_observed, observed_manifest, provenance_field,
          provenance_detail)) {
    return {nullptr,
            make_failure(ScheduleCertificationErrorCode::provenance_mismatch,
                         std::nullopt, std::nullopt,
                         std::move(provenance_field),
                         std::move(provenance_detail))};
  }
  std::string manifest_field;
  if (!compare_manifest(expected.manifest, observed_manifest, manifest_field)) {
    return {nullptr,
            make_failure(ScheduleCertificationErrorCode::manifest_mismatch,
                         std::nullopt, std::nullopt,
                         std::move(manifest_field),
                         "expected and observed typed manifests differ")};
  }

  std::unique_ptr<CertifiedLocalScheduleRegistry::Storage> storage =
      std::make_unique<CertifiedLocalScheduleRegistry::Storage>();
  storage->manifest = std::move(observed_manifest);
  storage->provenance = std::move(normalized_observed);
  for (std::size_t i = 0; i < candidates.size(); ++i)
    storage->functions[i] = std::move(candidates[i].functions);
  return {std::unique_ptr<CertifiedLocalScheduleRegistry>(
              new CertifiedLocalScheduleRegistry(std::move(storage))),
          std::nullopt};
}

} // namespace CarpetX
