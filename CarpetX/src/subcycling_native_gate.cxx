#include "subcycling_native_gate.hxx"

#include <cmath>
#include <stdexcept>

#ifndef CARPETX_SUBCYCLING_NATIVE_GATE_CONTRACT_UNIT
#include "driver.hxx"
#include "schedule.hxx"
#include "subcycling_dense_output.hxx"
#include "subcycling_schedule_certification.hxx"
#include "subcycling_schedule_certification_native_gate.hxx"
#include "subcycling_scratch_adapter.hxx"
#include "subcycling_scratch_state_transaction_factory.hxx"

#include <cctk.h>
#include <cctk_Groups.h>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#endif

namespace CarpetX {

NativeGateMethodContract
native_gate_method_contract(const int evolution_iteration) {
  switch (evolution_iteration) {
  case 1:
    return {SubcyclingODEMethod::rk4, "RK4", 4, 5};
  case 2:
    return {SubcyclingODEMethod::rkf78_order7, "RKF78", 23, 8};
  case 3:
    return {SubcyclingODEMethod::dp87_order8, "DP87", 39, 9};
  default:
    throw std::invalid_argument(
        "native gate requires evolution iteration 1, 2, or 3");
  }
}

namespace native_gate_detail {

StepContext make_patch_step_context(const PatchContextInput &input,
                                    const SubcyclingODEMethod method) {
  if (input.active_min_level != 0 || input.active_max_level != 1 ||
      input.active_min_patch != 0 || input.active_max_patch != 1)
    throw std::invalid_argument(
        "native gate active hierarchy is not patch zero, level zero");
  if (input.iteration < 0 || input.time_refinement_factor != 1 ||
      input.level_iteration != step_clock_t(input.iteration) ||
      !(step_clock_t(0) < input.delta_iteration) ||
      !std::isfinite(input.time) || !std::isfinite(input.delta_time) ||
      !(input.delta_time > 0.0))
    throw std::invalid_argument(
        "native gate patch clock envelope is unsupported");
  return {0, input.level_iteration - input.delta_iteration,
          input.level_iteration, input.time - input.delta_time, input.time,
          method, true,
          static_cast<std::uint64_t>(input.iteration)};
}

} // namespace native_gate_detail

#ifndef CARPETX_SUBCYCLING_NATIVE_GATE_CONTRACT_UNIT
namespace {

constexpr std::uint32_t expectation_schema = 1;

YAML::Node encode_strings(const std::vector<std::string> &values) {
  YAML::Node node(YAML::NodeType::Sequence);
  for (const auto &value : values)
    node.push_back(value);
  return node;
}

std::vector<std::string> decode_strings(const YAML::Node &node) {
  if (!node || !node.IsSequence())
    throw std::invalid_argument("native gate string list is malformed");
  std::vector<std::string> result;
  result.reserve(node.size());
  for (const auto &value : node)
    result.push_back(value.as<std::string>());
  return result;
}

YAML::Node encode_flags(const CanonicalExecutionFlags &flags) {
  YAML::Node node;
#define CARPETX_NATIVE_GATE_FLAG(name) node[#name] = flags.name
  CARPETX_NATIVE_GATE_FLAG(meta);
  CARPETX_NATIVE_GATE_FLAG(meta_early);
  CARPETX_NATIVE_GATE_FLAG(meta_late);
  CARPETX_NATIVE_GATE_FLAG(global);
  CARPETX_NATIVE_GATE_FLAG(global_early);
  CARPETX_NATIVE_GATE_FLAG(global_late);
  CARPETX_NATIVE_GATE_FLAG(level);
  CARPETX_NATIVE_GATE_FLAG(singlemap);
  CARPETX_NATIVE_GATE_FLAG(local);
  CARPETX_NATIVE_GATE_FLAG(loop_meta);
  CARPETX_NATIVE_GATE_FLAG(loop_global);
  CARPETX_NATIVE_GATE_FLAG(loop_level);
  CARPETX_NATIVE_GATE_FLAG(loop_singlemap);
  CARPETX_NATIVE_GATE_FLAG(loop_local);
#undef CARPETX_NATIVE_GATE_FLAG
  return node;
}

CanonicalExecutionFlags decode_flags(const YAML::Node &node) {
  if (!node || !node.IsMap())
    throw std::invalid_argument("native gate execution flags are malformed");
  CanonicalExecutionFlags flags;
#define CARPETX_NATIVE_GATE_FLAG(name) flags.name = node[#name].as<int>()
  CARPETX_NATIVE_GATE_FLAG(meta);
  CARPETX_NATIVE_GATE_FLAG(meta_early);
  CARPETX_NATIVE_GATE_FLAG(meta_late);
  CARPETX_NATIVE_GATE_FLAG(global);
  CARPETX_NATIVE_GATE_FLAG(global_early);
  CARPETX_NATIVE_GATE_FLAG(global_late);
  CARPETX_NATIVE_GATE_FLAG(level);
  CARPETX_NATIVE_GATE_FLAG(singlemap);
  CARPETX_NATIVE_GATE_FLAG(local);
  CARPETX_NATIVE_GATE_FLAG(loop_meta);
  CARPETX_NATIVE_GATE_FLAG(loop_global);
  CARPETX_NATIVE_GATE_FLAG(loop_level);
  CARPETX_NATIVE_GATE_FLAG(loop_singlemap);
  CARPETX_NATIVE_GATE_FLAG(loop_local);
#undef CARPETX_NATIVE_GATE_FLAG
  return flags;
}

YAML::Node encode_tag(const CanonicalTag &tag) {
  YAML::Node node;
  node["key"] = tag.key;
  node["type_code"] = tag.type_code;
  node["element_count"] = tag.element_count;
  YAML::Node bytes(YAML::NodeType::Sequence);
  for (const auto byte : tag.value_bytes)
    bytes.push_back(static_cast<unsigned int>(byte));
  node["value_bytes"] = std::move(bytes);
  return node;
}

CanonicalTag decode_tag(const YAML::Node &node) {
  CanonicalTag tag;
  tag.key = node["key"].as<std::string>();
  tag.type_code = node["type_code"].as<int>();
  tag.element_count = node["element_count"].as<std::uint64_t>();
  const auto bytes = node["value_bytes"];
  if (!bytes || !bytes.IsSequence())
    throw std::invalid_argument("native gate tag bytes are malformed");
  tag.value_bytes.reserve(bytes.size());
  for (const auto &entry : bytes) {
    const auto value = entry.as<unsigned int>();
    if (value > 255U)
      throw std::invalid_argument("native gate tag byte is out of range");
    tag.value_bytes.push_back(static_cast<std::uint8_t>(value));
  }
  return tag;
}

YAML::Node encode_access(const CanonicalAccess &access) {
  YAML::Node node;
  node["variable_name"] = access.variable_name;
  node["timelevel"] = access.timelevel;
  node["read_mask"] = access.read_mask;
  node["write_mask"] = access.write_mask;
  node["invalidate_mask"] = access.invalidate_mask;
  return node;
}

CanonicalAccess decode_access(const YAML::Node &node) {
  return {node["variable_name"].as<std::string>(),
          node["timelevel"].as<int>(), node["read_mask"].as<int>(),
          node["write_mask"].as<int>(),
          node["invalidate_mask"].as<int>()};
}

YAML::Node encode_group_timelevel(const CanonicalGroupTimelevel &entry) {
  YAML::Node node;
  node["group"] = entry.group;
  node["timelevel"] = entry.timelevel;
  return node;
}

CanonicalGroupTimelevel decode_group_timelevel(const YAML::Node &node) {
  return {node["group"].as<std::string>(), node["timelevel"].as<int>()};
}

YAML::Node encode_item(const CanonicalScheduleItem &item) {
  YAML::Node node;
  node["description"] = item.description;
  node["implementation"] = item.implementation;
  node["where"] = item.where;
  node["thorn"] = item.thorn;
  node["routine"] = item.routine;
  node["language"] = item.language;
  node["function_type"] = item.function_type;
  node["execution_flags"] = encode_flags(item.execution_flags);
  node["effective_mode"] = static_cast<int>(item.effective_mode);
  YAML::Node tags(YAML::NodeType::Sequence);
  for (const auto &tag : item.tags)
    tags.push_back(encode_tag(tag));
  node["tags"] = std::move(tags);
  node["reads_clauses"] = encode_strings(item.reads_clauses);
  node["writes_clauses"] = encode_strings(item.writes_clauses);
  node["invalidates_clauses"] = encode_strings(item.invalidates_clauses);
  YAML::Node accesses(YAML::NodeType::Sequence);
  for (const auto &access : item.accesses)
    accesses.push_back(encode_access(access));
  node["accesses"] = std::move(accesses);
  YAML::Node storage(YAML::NodeType::Sequence);
  for (const auto &entry : item.storage_groups)
    storage.push_back(encode_group_timelevel(entry));
  node["storage_groups"] = std::move(storage);
  node["communication_groups"] = encode_strings(item.communication_groups);
  node["sync_groups"] = encode_strings(item.sync_groups);
  node["trigger_groups"] = encode_strings(item.trigger_groups);
  node["ifs"] = encode_strings(item.ifs);
  node["whiles"] = encode_strings(item.whiles);
  return node;
}

CanonicalScheduleItem decode_item(const YAML::Node &node) {
  CanonicalScheduleItem item;
  item.description = node["description"].as<std::string>();
  item.implementation = node["implementation"].as<std::string>();
  item.where = node["where"].as<std::string>();
  item.thorn = node["thorn"].as<std::string>();
  item.routine = node["routine"].as<std::string>();
  item.language = node["language"].as<int>();
  item.function_type = node["function_type"].as<int>();
  item.execution_flags = decode_flags(node["execution_flags"]);
  item.effective_mode = static_cast<CanonicalScheduleMode>(
      node["effective_mode"].as<int>());
  for (const auto &tag : node["tags"])
    item.tags.push_back(decode_tag(tag));
  item.reads_clauses = decode_strings(node["reads_clauses"]);
  item.writes_clauses = decode_strings(node["writes_clauses"]);
  item.invalidates_clauses = decode_strings(node["invalidates_clauses"]);
  for (const auto &access : node["accesses"])
    item.accesses.push_back(decode_access(access));
  for (const auto &entry : node["storage_groups"])
    item.storage_groups.push_back(decode_group_timelevel(entry));
  item.communication_groups = decode_strings(node["communication_groups"]);
  item.sync_groups = decode_strings(node["sync_groups"]);
  item.trigger_groups = decode_strings(node["trigger_groups"]);
  item.ifs = decode_strings(node["ifs"]);
  item.whiles = decode_strings(node["whiles"]);
  return item;
}

YAML::Node encode_scope(const CanonicalScheduleScope &scope) {
  YAML::Node node;
  node["item_ordinal"] = scope.item_ordinal;
  node["parent_local_ordinal"] = scope.parent_local_ordinal;
  node["item"] = encode_item(scope.item);
  return node;
}

CanonicalScheduleScope decode_scope(const YAML::Node &node) {
  return {node["item_ordinal"].as<std::uint64_t>(),
          node["parent_local_ordinal"].as<std::uint64_t>(),
          decode_item(node["item"])};
}

YAML::Node encode_record(const CanonicalScheduleRecord &record) {
  YAML::Node node;
  node["traversal_ordinal"] = record.traversal_ordinal;
  node["item_ordinal"] = record.item_ordinal;
  node["parent_local_ordinal"] = record.parent_local_ordinal;
  node["item"] = encode_item(record.item);
  YAML::Node ancestry(YAML::NodeType::Sequence);
  for (const auto &scope : record.ancestry)
    ancestry.push_back(encode_scope(scope));
  node["ancestry"] = std::move(ancestry);
  return node;
}

CanonicalScheduleRecord decode_record(const YAML::Node &node) {
  CanonicalScheduleRecord record;
  record.traversal_ordinal = node["traversal_ordinal"].as<std::uint64_t>();
  record.item_ordinal = node["item_ordinal"].as<std::uint64_t>();
  record.parent_local_ordinal =
      node["parent_local_ordinal"].as<std::uint64_t>();
  record.item = decode_item(node["item"]);
  for (const auto &scope : node["ancestry"])
    record.ancestry.push_back(decode_scope(scope));
  return record;
}

YAML::Node encode_target(const CanonicalTargetSchedule &target) {
  YAML::Node node;
  node["target"] = static_cast<int>(target.target);
  node["body_root"] = target.body_root;
  node["entry_state"] = static_cast<int>(target.entry_state);
  node["exit_state"] = static_cast<int>(target.exit_state);
  node["entry_exit_included"] = target.entry_exit_included;
  YAML::Node records(YAML::NodeType::Sequence);
  for (const auto &record : target.records)
    records.push_back(encode_record(record));
  node["records"] = std::move(records);
  return node;
}

CanonicalTargetSchedule decode_target(const YAML::Node &node) {
  CanonicalTargetSchedule target;
  target.target =
      static_cast<SubcyclingScheduleTarget>(node["target"].as<int>());
  target.body_root = node["body_root"].as<std::string>();
  target.entry_state =
      static_cast<SpecialScheduleGroupState>(node["entry_state"].as<int>());
  target.exit_state =
      static_cast<SpecialScheduleGroupState>(node["exit_state"].as<int>());
  target.entry_exit_included = node["entry_exit_included"].as<bool>();
  for (const auto &record : node["records"])
    target.records.push_back(decode_record(record));
  return target;
}

YAML::Node encode_manifest(const CanonicalScheduleBundle &manifest) {
  YAML::Node targets(YAML::NodeType::Sequence);
  for (const auto &target : manifest.targets)
    targets.push_back(encode_target(target));
  return targets;
}

CanonicalScheduleBundle decode_manifest(const YAML::Node &node) {
  if (!node || !node.IsSequence() || node.size() != 2)
    throw std::invalid_argument("native gate manifest must have two targets");
  CanonicalScheduleBundle manifest;
  for (std::size_t i = 0; i < manifest.targets.size(); ++i)
    manifest.targets[i] = decode_target(node[i]);
  return manifest;
}

YAML::Node encode_provenance(const ScheduleBuildProvenance &provenance) {
  YAML::Node node;
  node["schema_version"] = provenance.schema_version;
  node["inventory_abi_version"] = provenance.inventory_abi_version;
  node["cactus_flesh_parent"] = provenance.cactus_flesh_parent;
  node["cactus_inventory_patch_sha256"] =
      provenance.cactus_inventory_patch_sha256;
  node["carpetx_parent"] = provenance.carpetx_parent;
  YAML::Node patches(YAML::NodeType::Sequence);
  for (const auto &patch : provenance.ordered_carpetx_patches) {
    YAML::Node entry;
    entry["phase"] = patch.phase;
    entry["sha256"] = patch.sha256;
    patches.push_back(std::move(entry));
  }
  node["ordered_carpetx_patches"] = std::move(patches);
  node["compiler_id"] = provenance.compiler_id;
  node["compiler_version"] = provenance.compiler_version;
  node["target_triple"] = provenance.target_triple;
  node["cxx_abi"] = provenance.cxx_abi;
  node["byte_order"] = provenance.byte_order;
  YAML::Node sizes(YAML::NodeType::Sequence);
  for (const auto &entry : provenance.cctk_type_sizes) {
    YAML::Node size;
    size["type"] = entry.first;
    size["size"] = entry.second;
    sizes.push_back(std::move(size));
  }
  node["cctk_type_sizes"] = std::move(sizes);
  node["executable_sha256"] = provenance.executable_sha256;
  return node;
}

ScheduleBuildProvenance decode_provenance(const YAML::Node &node) {
  if (!node || !node.IsMap())
    throw std::invalid_argument("native gate provenance is missing");
  ScheduleBuildProvenance provenance;
  provenance.schema_version = node["schema_version"].as<std::uint32_t>();
  provenance.inventory_abi_version =
      node["inventory_abi_version"].as<std::uint32_t>();
  provenance.cactus_flesh_parent =
      node["cactus_flesh_parent"].as<std::string>();
  provenance.cactus_inventory_patch_sha256 =
      node["cactus_inventory_patch_sha256"].as<std::string>();
  provenance.carpetx_parent = node["carpetx_parent"].as<std::string>();
  for (const auto &entry : node["ordered_carpetx_patches"])
    provenance.ordered_carpetx_patches.push_back(
        {entry["phase"].as<std::string>(),
         entry["sha256"].as<std::string>()});
  provenance.compiler_id = node["compiler_id"].as<std::string>();
  provenance.compiler_version = node["compiler_version"].as<std::string>();
  provenance.target_triple = node["target_triple"].as<std::string>();
  provenance.cxx_abi = node["cxx_abi"].as<std::string>();
  provenance.byte_order = node["byte_order"].as<std::string>();
  for (const auto &entry : node["cctk_type_sizes"])
    provenance.cctk_type_sizes.emplace_back(entry["type"].as<int>(),
                                            entry["size"].as<int>());
  provenance.executable_sha256 =
      node["executable_sha256"].as<std::string>();
  return provenance;
}

std::string expectation_path() {
  const char *const path = std::getenv("CARPETX_NATIVE_GATE_EXPECTATION");
  if (path == nullptr || *path == '\0')
    throw std::runtime_error(
        "CARPETX_NATIVE_GATE_EXPECTATION is not configured");
  return path;
}

YAML::Node load_document(const std::string &path) {
  const auto root = YAML::LoadFile(path);
  if (!root || !root.IsMap())
    throw std::invalid_argument("native gate expectation is malformed");
  if (root["schema"].as<std::uint32_t>() != expectation_schema)
    throw std::invalid_argument("native gate expectation schema differs");
  return root;
}

ScheduleCertificationExpectation load_expectation(const std::string &path) {
  const auto root = load_document(path);
  return {decode_manifest(root["manifest"]),
          decode_provenance(root["provenance"])};
}

void write_expectation(const std::string &path,
                       const ScheduleCertificationExpectation &expectation) {
  YAML::Node root;
  root["schema"] = expectation_schema;
  root["provenance"] = encode_provenance(expectation.provenance);
  root["manifest"] = encode_manifest(expectation.manifest);
  YAML::Emitter emitter;
  emitter << root;
  if (!emitter.good())
    throw std::runtime_error("native gate expectation serialization failed");
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output)
    throw std::runtime_error("native gate expectation cannot be written");
  output.write(emitter.c_str(), static_cast<std::streamsize>(emitter.size()));
  output.put('\n');
  if (!output)
    throw std::runtime_error("native gate expectation write failed");
}

constexpr std::uint32_t rotate_right(const std::uint32_t value,
                                     const int count) noexcept {
  return (value >> count) | (value << (32 - count));
}

void append_u64(std::vector<std::uint8_t> &bytes,
                const std::uint64_t value) {
  for (int shift = 56; shift >= 0; shift -= 8)
    bytes.push_back(static_cast<std::uint8_t>(value >> shift));
}

std::array<std::uint8_t, 32>
sha256(const std::vector<std::uint8_t> &message) {
  static constexpr std::array<std::uint32_t, 64> constants{{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U,
      0x3956c25bU, 0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U,
      0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U,
      0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
      0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U,
      0xc6e00bf3U, 0xd5a79147U, 0x06ca6351U, 0x14292967U,
      0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
      0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U,
      0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U,
      0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU, 0x682e6ff3U,
      0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U}};
  std::array<std::uint32_t, 8> hash{{
      0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
      0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}};
  std::vector<std::uint8_t> padded(message);
  const auto bits = static_cast<std::uint64_t>(message.size()) * 8U;
  padded.push_back(0x80U);
  while (padded.size() % 64U != 56U)
    padded.push_back(0U);
  append_u64(padded, bits);
  for (std::size_t block = 0; block < padded.size(); block += 64U) {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16U; ++i) {
      const auto offset = block + 4U * i;
      words[i] = (static_cast<std::uint32_t>(padded[offset]) << 24) |
                 (static_cast<std::uint32_t>(padded[offset + 1]) << 16) |
                 (static_cast<std::uint32_t>(padded[offset + 2]) << 8) |
                 static_cast<std::uint32_t>(padded[offset + 3]);
    }
    for (std::size_t i = 16; i < words.size(); ++i) {
      const auto s0 = rotate_right(words[i - 15], 7) ^
                      rotate_right(words[i - 15], 18) ^
                      (words[i - 15] >> 3);
      const auto s1 = rotate_right(words[i - 2], 17) ^
                      rotate_right(words[i - 2], 19) ^
                      (words[i - 2] >> 10);
      words[i] = words[i - 16] + s0 + words[i - 7] + s1;
    }
    auto a = hash[0];
    auto b = hash[1];
    auto c = hash[2];
    auto d = hash[3];
    auto e = hash[4];
    auto f = hash[5];
    auto g = hash[6];
    auto h = hash[7];
    for (std::size_t i = 0; i < words.size(); ++i) {
      const auto s1 = rotate_right(e, 6) ^ rotate_right(e, 11) ^
                      rotate_right(e, 25);
      const auto choice = (e & f) ^ ((~e) & g);
      const auto t1 = h + s1 + choice + constants[i] + words[i];
      const auto s0 = rotate_right(a, 2) ^ rotate_right(a, 13) ^
                      rotate_right(a, 22);
      const auto majority = (a & b) ^ (a & c) ^ (b & c);
      const auto t2 = s0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }
    hash[0] += a;
    hash[1] += b;
    hash[2] += c;
    hash[3] += d;
    hash[4] += e;
    hash[5] += f;
    hash[6] += g;
    hash[7] += h;
  }
  std::array<std::uint8_t, 32> result{};
  for (std::size_t i = 0; i < hash.size(); ++i)
    for (std::size_t byte = 0; byte < 4; ++byte)
      result[4 * i + byte] = static_cast<std::uint8_t>(
          hash[i] >> static_cast<int>(24 - 8 * byte));
  return result;
}

std::string current_executable_sha256() {
  std::ifstream input("/proc/self/exe", std::ios::binary);
  if (!input)
    throw std::runtime_error("native gate cannot read /proc/self/exe");
  std::vector<std::uint8_t> bytes;
  std::array<char, 1 << 20> buffer{};
  while (input.read(buffer.data(), static_cast<std::streamsize>(buffer.size())) ||
         input.gcount() > 0) {
    const auto count = static_cast<std::size_t>(input.gcount());
    bytes.insert(bytes.end(), buffer.begin(), buffer.begin() + count);
  }
  if (!input.eof())
    throw std::runtime_error("native gate executable read failed");
  const auto digest = sha256(bytes);
  std::ostringstream output;
  output << std::hex << std::setfill('0');
  for (const auto byte : digest)
    output << std::setw(2) << static_cast<unsigned int>(byte);
  return output.str();
}

[[noreturn]] void certification_failure(
    const ScheduleCertificationFailure &failure) {
  throw std::runtime_error("native gate schedule certification failed at " +
                           failure.field + ": " + failure.detail);
}

class NativeGateStagePreparer final : public StagePreparer {
public:
  explicit NativeGateStagePreparer(const StepContext &context)
      : context_(context) {}

  void prepare_stage(const StepContext &context,
                     const StagePoint &stage_point) override {
    const auto scale = std::max(
        {1.0, std::abs(context_.begin_time), std::abs(context_.end_time)});
    const auto tolerance =
        32.0 * std::numeric_limits<double>::epsilon() * scale;
    if (context.level != context_.level ||
        context.begin_clock != context_.begin_clock ||
        context.end_clock != context_.end_clock ||
        context.method != context_.method ||
        !std::isfinite(stage_point.stage_time) ||
        stage_point.stage_fraction < step_clock_t(0) ||
        stage_point.stage_fraction > step_clock_t(1) ||
        stage_point.stage_time < context_.begin_time - tolerance ||
        stage_point.stage_time > context_.end_time + tolerance)
      throw std::invalid_argument(
          "native gate received an invalid one-level stage request");
  }

private:
  StepContext context_;
};

int required_group(const char *const name) {
  const int group = CCTK_GroupIndex(name);
  if (group < 0)
    throw std::runtime_error(std::string("native gate group is missing: ") +
                             name);
  return group;
}

} // namespace

struct NativeGateSession::Storage {
  NativeGateMethodContract contract;
  StepContext context;
  std::unique_ptr<CertifiedLocalScheduleRegistry> registry;
  std::unique_ptr<ScratchStateTransaction> transaction;
  NativeGateStagePreparer preparer;
  std::unique_ptr<ScopedStepContext> scope;

  Storage(NativeGateMethodContract contract_, StepContext context_,
          std::unique_ptr<CertifiedLocalScheduleRegistry> registry_,
          std::unique_ptr<ScratchStateTransaction> transaction_)
      : contract(contract_), context(std::move(context_)),
        registry(std::move(registry_)), transaction(std::move(transaction_)),
        preparer(context) {
    scope = std::make_unique<ScopedStepContext>(context, preparer,
                                                transaction.get());
  }
};

NativeGateSession::NativeGateSession() noexcept = default;
NativeGateSession::~NativeGateSession() = default;
NativeGateSession::NativeGateSession(NativeGateSession &&) noexcept = default;
NativeGateSession &
NativeGateSession::operator=(NativeGateSession &&) noexcept = default;
NativeGateSession::NativeGateSession(std::unique_ptr<Storage> storage) noexcept
    : storage_(std::move(storage)) {}
bool NativeGateSession::active() const noexcept { return storage_ != nullptr; }

void write_native_gate_inventory(cGH *const cctkGH,
                                 const int cctk_itlast) {
  if (cctkGH == nullptr || cctk_itlast != 0)
    throw std::invalid_argument(
        "native gate inventory requires cctk_itlast=0");
  const auto path = expectation_path();
  auto root = load_document(path);
  auto provenance = decode_provenance(root["provenance"]);
  provenance.executable_sha256 = current_executable_sha256();
  auto observed =
      observe_local_subcycling_schedules_for_native_gate(provenance);
  if (!observed) {
    if (!observed.failure)
      throw std::runtime_error("native gate inventory failed without detail");
    certification_failure(*observed.failure);
  }
  write_expectation(path, *observed.expectation);
}

NativeGateSession begin_native_gate(cGH *const cctkGH,
                                    const int cctk_itlast) {
  if (cctkGH == nullptr || cctk_itlast != 3 || CCTK_nProcs(cctkGH) != 1)
    throw std::invalid_argument(
        "native gate requires one rank and cctk_itlast=3");
  if (!ghext || ghext->patchdata.size() != 1 ||
      ghext->patchdata.front().leveldata.size() != 1)
    throw std::invalid_argument(
        "native gate hierarchy is not exactly one patch and one level");
  if (!active_levels || active_levels->min_level != 0 ||
      active_levels->max_level != 1 || active_levels->min_patch != 0 ||
      active_levels->max_patch != 1)
    throw std::invalid_argument(
        "native gate active hierarchy is not patch zero, level zero");

  const auto &level = ghext->patchdata.front().leveldata.front();
  const auto *const patch_cctkGH = level.get_patch_cctkGH();
  if (patch_cctkGH == nullptr)
    throw std::invalid_argument("native gate patch GH is absent");

  auto expected = load_expectation(expectation_path());
  auto observed_provenance = expected.provenance;
  observed_provenance.executable_sha256 = current_executable_sha256();
  auto certification =
      certify_local_subcycling_schedules(expected, observed_provenance);
  if (!certification) {
    if (!certification.failure)
      throw std::runtime_error(
          "native gate certification failed without detail");
    certification_failure(*certification.failure);
  }

  const auto contract =
      native_gate_method_contract(patch_cctkGH->cctk_iteration);
  const StepContext context = native_gate_detail::make_patch_step_context(
      {active_levels->min_level,
       active_levels->max_level,
       active_levels->min_patch,
       active_levels->max_patch,
       patch_cctkGH->cctk_iteration,
       patch_cctkGH->cctk_timefac,
       step_clock_t(level.iteration.num, level.iteration.den),
       step_clock_t(level.delta_iteration.num, level.delta_iteration.den),
       patch_cctkGH->cctk_time,
       patch_cctkGH->cctk_delta_time},
      contract.method);

  const int state_group = required_group("TestODESolvers2::state");
  const int rhs_group = required_group("TestODESolvers2::rhs");
  const int dependent_group =
      required_group("TestODESolvers2::gate_dependent");
  ScratchStateTransactionFactoryMetadata metadata{
      0,
      0,
      0,
      0,
      step_clock_t(1),
      1.0,
      1,
      false,
      {{state_group, rhs_group}},
      {rhs_group, dependent_group},
      {},
      {}};
  auto frame = copy_canonical_tl0_collective(*ghext, 0, 0);
  auto transaction = ScratchStateTransactionFactory::create_native(
      *ghext, *certification.registry, std::move(frame),
      std::move(metadata));
  return NativeGateSession(std::make_unique<NativeGateSession::Storage>(
      contract, context, std::move(certification.registry),
      std::move(transaction)));
}

NativeGateReceipt end_native_gate(NativeGateSession &&session) {
  if (!session.storage_ || !session.storage_->transaction ||
      !session.storage_->scope)
    throw std::invalid_argument("native gate session is not active");
  auto &storage = *session.storage_;
  const auto rhs_count = storage.transaction->rhs_evaluation_count();
  auto interval = storage.transaction->take_committed_dense();
  if (!interval)
    throw std::runtime_error(
        "native gate ODESolvers did not publish a dense interval");
  const auto &id = interval->id();
  const auto &capability = interval->capability();
  if (id.level != storage.context.level ||
      id.begin_clock != storage.context.begin_clock ||
      id.end_clock != storage.context.end_clock ||
      id.begin_time != storage.context.begin_time ||
      id.end_time != storage.context.end_time ||
      id.method != storage.context.method ||
      capability.method != storage.context.method ||
      capability.extra_rhs_evaluations !=
          static_cast<int>(storage.contract.extra_rhs_evaluations) ||
      rhs_count != storage.contract.extra_rhs_evaluations ||
      interval->control_count() != storage.contract.control_count)
    throw std::runtime_error(
        "native gate dense interval metadata or count differs");
  const NativeGateReceipt receipt{storage.context.method, rhs_count,
                                  interval->control_count()};
  storage.scope.reset();
  session.storage_.reset();
  return receipt;
}
#endif

} // namespace CarpetX
