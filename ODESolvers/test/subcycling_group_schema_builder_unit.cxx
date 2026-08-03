#include "subcycling_group_schema_builder.hxx"

#include <cassert>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using ODESolvers::GroupSchemaCatalog;
using ODESolvers::GroupSchemaGroupType;

template <class Exception, class Function>
void require_throws(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Exception &) {
    threw = true;
  }
  assert(threw);
}

struct FakeGroup {
  GroupSchemaGroupType type{GroupSchemaGroupType::grid_function};
  int numvars{1};
  std::string full_name;
  std::optional<std::string> rhs;
  std::optional<std::string> dependents;
};

class FakeCatalog final : public GroupSchemaCatalog {
public:
  mutable int calls{0};
  std::vector<FakeGroup> groups;
  std::map<std::pair<int, std::string>, int> resolutions;

  int group_count() const override {
    ++calls;
    return static_cast<int>(groups.size());
  }

  GroupSchemaGroupType group_type(const int group) const override {
    ++calls;
    return groups.at(group).type;
  }

  int group_numvars(const int group) const override {
    ++calls;
    return groups.at(group).numvars;
  }

  std::string full_group_name(const int group) const override {
    ++calls;
    return groups.at(group).full_name;
  }

  std::optional<std::string>
  group_tag(const int group, const std::string_view key) const override {
    ++calls;
    if (key == "rhs")
      return groups.at(group).rhs;
    if (key == "dependents")
      return groups.at(group).dependents;
    throw std::invalid_argument("unexpected fake group tag");
  }

  int resolve_group(const int owner,
                    const std::string_view name) const override {
    ++calls;
    const auto found = resolutions.find({owner, std::string(name)});
    return found == resolutions.end() ? -1 : found->second;
  }
};

FakeCatalog valid_catalog() {
  FakeCatalog catalog;
  catalog.groups = {
      {GroupSchemaGroupType::grid_function, 2, "Impl::evolved_z",
       std::string("rhs_z"), std::string("aux_z aux_a aux_z")},
      {GroupSchemaGroupType::grid_function, 1, "Impl::unused", std::nullopt,
       std::nullopt},
      {GroupSchemaGroupType::grid_function, 1, "Other::evolved_a",
       std::string("Other::rhs_a"), std::string("Other::aux_a")},
      {GroupSchemaGroupType::grid_function, 2, "Impl::rhs_z", std::nullopt,
       std::nullopt},
      {GroupSchemaGroupType::grid_function, 1, "Other::rhs_a", std::nullopt,
       std::nullopt},
      {GroupSchemaGroupType::grid_function, 1, "Other::aux_a", std::nullopt,
       std::nullopt},
      {GroupSchemaGroupType::grid_function, 1, "Impl::aux_z", std::nullopt,
       std::nullopt},
      {GroupSchemaGroupType::grid_function, 1, "Impl::aux_a", std::nullopt,
       std::nullopt},
  };
  catalog.resolutions = {
      {{0, "rhs_z"}, 3},       {{0, "aux_z"}, 6},
      {{0, "aux_a"}, 7},       {{2, "Other::rhs_a"}, 4},
      {{2, "Other::aux_a"}, 5},
  };
  return catalog;
}

void test_order_and_exact_dependent_closure() {
  auto catalog = valid_catalog();
  const auto build = ODESolvers::build_group_schema(catalog, true);

  assert(build.contract.ordered_group_pairs.size() == 2);
  assert(build.contract.ordered_group_pairs[0].evolved_group == 0);
  assert(build.contract.ordered_group_pairs[0].rhs_group == 3);
  assert(build.contract.ordered_group_pairs[1].evolved_group == 2);
  assert(build.contract.ordered_group_pairs[1].rhs_group == 4);
  assert((build.evolved_groups == std::vector<int>{0, 2}));
  assert((build.rhs_groups == std::vector<int>{3, 4}));
  assert((build.contract.dependent_groups ==
          std::vector<int>{3, 4, 5, 6, 7}));
  assert(build.evolved_variable_count == 3);
}

void test_duplicate_rhs_and_numvar_mismatch_fail_closed() {
  {
    auto catalog = valid_catalog();
    catalog.groups[2].rhs = "Impl::rhs_z";
    catalog.resolutions[{2, "Impl::rhs_z"}] = 3;
    require_throws<std::invalid_argument>(
        [&] { (void)ODESolvers::build_group_schema(catalog, true); });
  }
  {
    auto catalog = valid_catalog();
    catalog.groups[4].numvars = 2;
    require_throws<std::invalid_argument>(
        [&] { (void)ODESolvers::build_group_schema(catalog, true); });
  }
}

void test_subcycling_requires_grid_functions() {
  {
    auto catalog = valid_catalog();
    catalog.groups[0].type = GroupSchemaGroupType::scalar;
    require_throws<std::invalid_argument>(
        [&] { (void)ODESolvers::build_group_schema(catalog, true); });
  }
  {
    auto catalog = valid_catalog();
    catalog.groups[3].type = GroupSchemaGroupType::array;
    require_throws<std::invalid_argument>(
        [&] { (void)ODESolvers::build_group_schema(catalog, true); });
  }
  {
    auto catalog = valid_catalog();
    catalog.groups[6].type = GroupSchemaGroupType::scalar;
    require_throws<std::invalid_argument>(
        [&] { (void)ODESolvers::build_group_schema(catalog, true); });
  }
}

void test_legacy_mode_ignores_non_gf_evolved_groups() {
  auto catalog = valid_catalog();
  catalog.groups.push_back({GroupSchemaGroupType::scalar, 1,
                            "Impl::global_evolved", std::string("rhs_z"),
                            std::string("aux_z")});
  catalog.resolutions[{8, "rhs_z"}] = 3;
  catalog.resolutions[{8, "aux_z"}] = 6;

  const auto build = ODESolvers::build_group_schema(catalog, false);
  assert(build.contract.ordered_group_pairs.size() == 2);
  assert((build.evolved_groups == std::vector<int>{0, 2}));
}

void test_disabled_subcycling_is_a_complete_noop() {
  FakeCatalog catalog;
  catalog.groups.push_back({GroupSchemaGroupType::scalar, 1,
                            "Impl::bad_global", std::string("missing"),
                            std::nullopt});
  const auto build =
      ODESolvers::maybe_build_subcycling_group_schema(false, catalog);
  assert(!build.has_value());
  assert(catalog.calls == 0);
}

} // namespace

int main() {
  test_order_and_exact_dependent_closure();
  test_duplicate_rhs_and_numvar_mismatch_fail_closed();
  test_subcycling_requires_grid_functions();
  test_legacy_mode_ignores_non_gf_evolved_groups();
  test_disabled_subcycling_is_a_complete_noop();
}
