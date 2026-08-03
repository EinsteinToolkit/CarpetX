#include "subcycling_group_schema_builder.hxx"

#include <algorithm>
#include <cctype>
#include <limits>
#include <stdexcept>
#include <string>

namespace ODESolvers {
namespace {

std::vector<std::string> split_group_names(const std::string &names) {
  std::vector<std::string> result;
  std::size_t position = 0;
  for (;;) {
    while (position < names.size() &&
           std::isspace(static_cast<unsigned char>(names[position])))
      ++position;
    if (position == names.size())
      break;
    const auto begin = position;
    while (position < names.size() &&
           !std::isspace(static_cast<unsigned char>(names[position])))
      ++position;
    result.push_back(names.substr(begin, position - begin));
  }
  return result;
}

void require_valid_index(const int group, const int group_count,
                         const char *const role) {
  if (group < 0 || group >= group_count)
    throw std::invalid_argument(std::string(role) +
                                " group index is outside the CCTK catalog");
}

void require_grid_function(const GroupSchemaCatalog &catalog,
                           const int group, const char *const role) {
  if (catalog.group_type(group) != GroupSchemaGroupType::grid_function)
    throw std::invalid_argument(
        std::string(role) + " group \"" + catalog.full_group_name(group) +
        "\" is not a CCTK grid-function group");
}

int resolve_required_group(const GroupSchemaCatalog &catalog,
                           const int owner, const std::string_view name,
                           const int group_count, const char *const role) {
  const int resolved = catalog.resolve_group(owner, name);
  if (resolved < 0)
    throw std::invalid_argument(
        "variable group \"" + catalog.full_group_name(owner) +
        "\" declares a " + role + " group \"" + std::string(name) +
        "\" that does not exist");
  require_valid_index(resolved, group_count, role);
  return resolved;
}

void sort_unique_or_reject(std::vector<int> &groups,
                           const char *const role) {
  std::sort(groups.begin(), groups.end());
  const auto last = std::unique(groups.begin(), groups.end());
  if (last != groups.end())
    throw std::invalid_argument(std::string("duplicate ") + role +
                                " group in ODE schema");
}

} // namespace

GroupSchemaBuild build_group_schema(const GroupSchemaCatalog &catalog,
                                    const bool require_grid_functions) {
  const int group_count = catalog.group_count();
  if (group_count < 0)
    throw std::invalid_argument("CCTK group count is negative");

  GroupSchemaBuild build;
  for (int evolved_group = 0; evolved_group < group_count; ++evolved_group) {
    const auto type = catalog.group_type(evolved_group);
    if (type != GroupSchemaGroupType::grid_function &&
        !require_grid_functions)
      continue;

    const auto rhs_tag = catalog.group_tag(evolved_group, "rhs");
    if (!rhs_tag || rhs_tag->empty())
      continue;

    if (type != GroupSchemaGroupType::grid_function)
      require_grid_function(catalog, evolved_group, "evolved");

    const int rhs_group = resolve_required_group(
        catalog, evolved_group, *rhs_tag, group_count, "RHS");
    if (rhs_group == evolved_group)
      throw std::invalid_argument(
          "an evolved group cannot also be its own RHS group");
    require_grid_function(catalog, rhs_group, "RHS");
    if (catalog.group_numvars(evolved_group) !=
        catalog.group_numvars(rhs_group))
      throw std::invalid_argument(
          "evolved and RHS groups have different variable counts");

    if (build.evolved_variable_count >
        std::numeric_limits<int>::max() -
            catalog.group_numvars(evolved_group))
      throw std::overflow_error("evolved variable count overflow");
    build.evolved_variable_count += catalog.group_numvars(evolved_group);
    build.evolved_groups.push_back(evolved_group);
    build.rhs_groups.push_back(rhs_group);
    build.contract.ordered_group_pairs.push_back(
        {evolved_group, rhs_group});

    const auto dependents_tag = catalog.group_tag(evolved_group, "dependents");
    if (!dependents_tag)
      continue;
    for (const auto &name : split_group_names(*dependents_tag)) {
      const int dependent_group = resolve_required_group(
          catalog, evolved_group, name, group_count, "dependent");
      require_grid_function(catalog, dependent_group, "dependent");
      build.contract.dependent_groups.push_back(dependent_group);
    }
  }

  sort_unique_or_reject(build.evolved_groups, "evolved");
  sort_unique_or_reject(build.rhs_groups, "RHS");

  build.contract.dependent_groups.insert(
      build.contract.dependent_groups.end(), build.rhs_groups.begin(),
      build.rhs_groups.end());
  std::sort(build.contract.dependent_groups.begin(),
            build.contract.dependent_groups.end());
  build.contract.dependent_groups.erase(
      std::unique(build.contract.dependent_groups.begin(),
                  build.contract.dependent_groups.end()),
      build.contract.dependent_groups.end());

  for (const int evolved_group : build.evolved_groups)
    if (std::binary_search(build.contract.dependent_groups.begin(),
                           build.contract.dependent_groups.end(),
                           evolved_group))
      throw std::invalid_argument(
          "an evolved group is also in the dependent-group closure");
  for (const int rhs_group : build.rhs_groups)
    if (std::binary_search(build.evolved_groups.begin(),
                           build.evolved_groups.end(), rhs_group))
      throw std::invalid_argument(
          "an RHS group is also an evolved group");

  return build;
}

std::optional<GroupSchemaBuild>
maybe_build_subcycling_group_schema(const bool subcycling_enabled,
                                    const GroupSchemaCatalog &catalog) {
  if (!subcycling_enabled)
    return std::nullopt;
  return build_group_schema(catalog, true);
}

} // namespace ODESolvers
