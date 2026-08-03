#ifndef ODESOLVERS_SUBCYCLING_GROUP_SCHEMA_BUILDER_HXX
#define ODESOLVERS_SUBCYCLING_GROUP_SCHEMA_BUILDER_HXX

#include <subcycling_method_contract.hxx>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace ODESolvers {

enum class GroupSchemaGroupType { grid_function, array, scalar };

// Narrow metadata seam shared by the pure builder tests and the Cactus
// adapter. Implementations must expose groups in CCTK group-index order.
class GroupSchemaCatalog {
public:
  virtual ~GroupSchemaCatalog() = default;

  virtual int group_count() const = 0;
  virtual GroupSchemaGroupType group_type(int group) const = 0;
  virtual int group_numvars(int group) const = 0;
  virtual std::string full_group_name(int group) const = 0;
  virtual std::optional<std::string>
  group_tag(int group, std::string_view key) const = 0;
  virtual int resolve_group(int owner, std::string_view name) const = 0;
};

struct GroupSchemaBuild {
  CarpetX::SubcyclingGroupSchema contract;
  std::vector<int> evolved_groups;
  std::vector<int> rhs_groups;
  int evolved_variable_count{0};
};

// Build the single canonical evolved/RHS/dependent schema. When
// require_grid_functions is false, evolved arrays/scalars are ignored exactly
// as in the legacy solve path. When true, every evolved, RHS, and dependent
// group must be a CCTK grid-function group.
GroupSchemaBuild build_group_schema(const GroupSchemaCatalog &catalog,
                                    bool require_grid_functions);

// This guard is intentionally outside the catalog: disabled subcycling must
// return before any Cactus metadata access.
std::optional<GroupSchemaBuild>
maybe_build_subcycling_group_schema(bool subcycling_enabled,
                                    const GroupSchemaCatalog &catalog);

// Production adapter used by both ODESolvers_Solve and the PARAMCHECK
// method-contract handoff. The optional form is a complete no-op when disabled.
bool carpetx_subcycling_enabled();
GroupSchemaBuild
build_cactus_group_schema(bool require_grid_functions);
std::optional<GroupSchemaBuild>
maybe_build_cactus_subcycling_group_schema(bool subcycling_enabled);

} // namespace ODESolvers

#endif // ODESOLVERS_SUBCYCLING_GROUP_SCHEMA_BUILDER_HXX
