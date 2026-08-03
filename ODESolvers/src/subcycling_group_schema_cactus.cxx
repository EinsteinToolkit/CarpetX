#include "subcycling_group_schema_builder.hxx"

#include <cctk.h>
#include <cctk_Groups.h>
#include <cctk_Parameter.h>
#include <util_Table.h>

#include <cstring>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace ODESolvers {
namespace {

class CactusGroupSchemaCatalog final : public GroupSchemaCatalog {
public:
  int group_count() const override { return CCTK_NumGroups(); }

  GroupSchemaGroupType group_type(const int group) const override {
    switch (CCTK_GroupTypeI(group)) {
    case CCTK_GF:
      return GroupSchemaGroupType::grid_function;
    case CCTK_ARRAY:
      return GroupSchemaGroupType::array;
    case CCTK_SCALAR:
      return GroupSchemaGroupType::scalar;
    default:
      throw std::invalid_argument("CCTK group has an unknown group type");
    }
  }

  int group_numvars(const int group) const override {
    const int count = CCTK_NumVarsInGroupI(group);
    if (count < 0)
      throw std::invalid_argument(
          "could not read the CCTK group variable count");
    return count;
  }

  std::string full_group_name(const int group) const override {
    const char *const name = CCTK_FullGroupName(group);
    if (name == nullptr)
      throw std::invalid_argument("could not read the CCTK full group name");
    return name;
  }

  std::optional<std::string>
  group_tag(const int group, const std::string_view key) const override {
    const int tags = CCTK_GroupTagsTableI(group);
    if (tags < 0)
      throw std::invalid_argument("could not read the CCTK group tag table");
    const std::string key_string(key);
    std::vector<char> buffer(1000);
    const int status = Util_TableGetString(
        tags, static_cast<int>(buffer.size()), buffer.data(),
        key_string.c_str());
    if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY)
      return std::nullopt;
    if (status < 0)
      throw std::invalid_argument("could not read CCTK group tag \"" +
                                  key_string + "\"");
    return std::string(buffer.data());
  }

  int resolve_group(const int owner,
                    const std::string_view name) const override {
    std::string qualified(name);
    if (qualified.find(':') == std::string::npos) {
      const char *const thorn_or_impl = CCTK_GroupImplementationI(owner);
      if (thorn_or_impl == nullptr)
        throw std::invalid_argument(
            "could not read the evolved group's implementation");
      const char *const implementation =
          CCTK_ThornImplementation(thorn_or_impl);
      const char *const thorn = CCTK_ImplementationThorn(thorn_or_impl);
      if (implementation == nullptr && thorn == nullptr)
        throw std::invalid_argument(
            "could not resolve the evolved group's namespace");
      const char *prefix = nullptr;
      if (implementation == nullptr) {
        prefix = thorn;
      } else if (thorn == nullptr) {
        prefix = implementation;
      } else {
        if (std::strcmp(implementation, thorn) != 0)
          throw std::invalid_argument(
              "CCTK implementation and thorn namespace differ");
        prefix = implementation;
      }
      std::ostringstream stream;
      stream << prefix << "::" << qualified;
      qualified = stream.str();
    }
    return CCTK_GroupIndex(qualified.c_str());
  }
};

} // namespace

bool carpetx_subcycling_enabled() {
  const auto *const enabled = static_cast<const CCTK_INT *>(
      CCTK_ParameterGet("use_subcycling_wip", "CarpetX", nullptr));
  if (enabled == nullptr)
    throw std::logic_error(
        "CarpetX::use_subcycling_wip parameter is unavailable");
  return *enabled != 0;
}

GroupSchemaBuild
build_cactus_group_schema(const bool require_grid_functions) {
  const CactusGroupSchemaCatalog catalog;
  return build_group_schema(catalog, require_grid_functions);
}

std::optional<GroupSchemaBuild> maybe_build_cactus_subcycling_group_schema(
    const bool subcycling_enabled) {
  if (!subcycling_enabled)
    return std::nullopt;
  const CactusGroupSchemaCatalog catalog;
  return maybe_build_subcycling_group_schema(true, catalog);
}

} // namespace ODESolvers
