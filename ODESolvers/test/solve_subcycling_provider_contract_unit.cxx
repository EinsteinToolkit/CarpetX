#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void check(const bool condition, const char *const message) {
  if (!condition)
    throw std::runtime_error(message);
}

std::string read_file(const char *const path) {
  std::ifstream input(path);
  if (!input)
    throw std::runtime_error("could not open ODESolvers solve source");
  std::ostringstream contents;
  contents << input.rdbuf();
  return contents.str();
}

std::size_t require_find(const std::string &source,
                         const std::string &needle) {
  const auto position = source.find(needle);
  if (position == std::string::npos)
    throw std::runtime_error("missing active-provider contract: " + needle);
  return position;
}

} // namespace

int main(const int argc, char **const argv) {
  const char *const source_path =
      argc == 2 ? argv[1] : "ODESolvers/src/solve.cxx";
  const auto source = read_file(source_path);

  require_find(source, "#include \"subcycling_ode_provider_registry.hxx\"");
  require_find(source, "#include \"subcycling_group_schema_builder.hxx\"");
  const auto active_begin =
      require_find(source, "if (CarpetX::step_context_active()) {");
  const auto active_end = require_find(source, "static bool did_output");
  check(active_begin < active_end, "active preflight is not before setup");
  const auto active_preflight =
      source.substr(active_begin, active_end - active_begin);
  require_find(active_preflight,
               "require_subcycling_ode_provider(std::string_view(method))");
  require_find(active_preflight,
               "active_context->method != candidate.dense.method");
  require_find(active_preflight,
               "candidate.dense.tableau_fingerprint !=");
  require_find(active_preflight,
               "explicit_rk_tableau_fingerprint(candidate.method)");
  require_find(active_preflight,
               "candidate.dense.endpoint_order != tableau.endpoint_order");
  require_find(active_preflight,
               "candidate.dense.stage_count !=");
  check(active_preflight.find("CCTK_EQUALS(method") == std::string::npos,
        "active StepContext still manually dispatches the method string");

  const auto primary_setup = require_find(source, "statecomp_t var, rhs;");
  check(active_begin < primary_setup && active_end < primary_setup,
        "provider validation can occur after primary-state setup");

  require_find(source,
               "ExplicitRKMethod explicit_method = ExplicitRKMethod::rk4;");
  require_find(source, "explicit_method = active_provider->method;");
  require_find(source,
               "CarpetX::DenseOutputProvider provider(active_provider->dense)");
  require_find(source, "active_provider->dense.tableau_fingerprint");
  check(source.find("make_reference_dense_provider(explicit_method)") ==
            std::string::npos,
        "solve still constructs dense providers outside the registry");
  check(source.find("explicit_rk_tableau_fingerprint(explicit_method)") ==
            std::string::npos,
        "solve still constructs interval fingerprints outside the registry");

  require_find(source,
               "build_cactus_group_schema(carpetx_subcycling_enabled())");
  require_find(source,
               "group_schema.contract.ordered_group_pairs");
  require_find(source,
               "group_schema.contract.dependent_groups");
  check(source.find("int groupindex(") == std::string::npos,
        "solve still owns group-name resolution");
  check(source.find("int get_group_rhs(") == std::string::npos,
        "solve still owns RHS tag parsing");
  check(source.find("get_group_dependents(") == std::string::npos,
        "solve still owns dependent tag parsing");
  check(source.find("Util_TableGetString") == std::string::npos,
        "solve still reads CCTK group tags directly");

  std::cout << "ODESolvers active subcycling provider contract tests passed\n";
  return 0;
}
