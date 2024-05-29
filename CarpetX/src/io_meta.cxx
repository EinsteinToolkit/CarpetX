#include "io_meta.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <fstream>
#include <iomanip>
#include <utility>

namespace CarpetX {

const std::array<int, 3> cactus_metadata_version{1, 0, 0};

YAML::Emitter &operator<<(YAML::Emitter &yaml, const group_type_t &group_type) {
  yaml << YAML::LocalTag("group_type-1.0.0");
  switch (group_type) {
  case group_type_t::none:
    yaml << "none";
    break;
  case group_type_t::gf:
    yaml << "gf";
    break;
  case group_type_t::array:
    yaml << "array";
    break;
  default:
    assert(0);
  }
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const reduction_t &reduction) {
  yaml << YAML::LocalTag("reduction-1.0.0");
  switch (reduction) {
  case reduction_t::volume:
    yaml << "volume";
    break;
  case reduction_t::sum:
    yaml << "sum";
    break;
  case reduction_t::sum_abs:
    yaml << "sum_abs";
    break;
  case reduction_t::sum_squared:
    yaml << "sum_squared";
    break;
  case reduction_t::average:
    yaml << "average";
    break;
  case reduction_t::standard_deviation:
    yaml << "standard_deviation";
    break;
  case reduction_t::minimum:
    yaml << "minimum";
    break;
  case reduction_t::maximum:
    yaml << "maximum";
    break;
  case reduction_t::maximum_abs:
    yaml << "maximum_abs";
    break;
  case reduction_t::norm1:
    yaml << "norm1";
    break;
  case reduction_t::norm2:
    yaml << "norm2";
    break;
  case reduction_t::norm_inf:
    yaml << "norm_inf";
    break;
  case reduction_t::minimum_location:
    yaml << "minimum_location";
    break;
  case reduction_t::maximum_location:
    yaml << "maximum_location";
    break;
  case reduction_t::sum_location:
    yaml << "sum_location";
    break;
  default:
    assert(0);
  }
  return yaml;
}

struct real_output_file_description_t {
  std::string filename;                //
  std::string description;             // human readable
  std::string writer_thorn;            //
  std::vector<std::string> variables;  //
  std::vector<int> iterations;         //
  group_type_t group_type;             // set automatically
  int variable_dimensions;             // set automatically
  std::vector<int> output_directions;  //
  std::vector<reduction_t> reductions; //
  std::string format_name;             //
  std::array<int, 3> format_version;   //

  real_output_file_description_t() = delete;

  real_output_file_description_t(const real_output_file_description_t &) =
      default;
  real_output_file_description_t(real_output_file_description_t &&) = default;
  real_output_file_description_t &
  operator=(const real_output_file_description_t &) = default;
  real_output_file_description_t &
  operator=(real_output_file_description_t &&) = default;

  real_output_file_description_t(
      output_file_description_t output_file_description);
};

real_output_file_description_t::real_output_file_description_t(
    output_file_description_t ofd) {
  if (iterations.empty())
    return; // CCTK_ERROR("Empty iterations");
  if (variables.empty())
    return; // CCTK_ERROR("Empty variables");

  filename = ofd.filename;
  // TODO: Check whether file exists
  if (filename == "")
    CCTK_ERROR("Empty filename");

  description = ofd.description;

  writer_thorn = ofd.writer_thorn;
  const int is_active = CCTK_IsThornActive(writer_thorn.c_str());
  if (is_active <= 0)
    CCTK_ERROR("Writer thorn not active");
  for (auto &ch : writer_thorn)
    ch = std::tolower(ch);

  variables = ofd.variables;
  for (auto &var : variables) {
    const int varindex = CCTK_VarIndex(var.c_str());
    if (varindex < 0)
      CCTK_VERROR("Variable \"%s\" does not exist", var.c_str());
    for (auto &ch : var)
      ch = std::tolower(ch);
  }

  iterations = ofd.iterations;
  for (const auto &iter : iterations)
    if (iter < 0)
      CCTK_VERROR("Bad iteration number %d", iter);

  if (variables.empty()) {
    group_type = group_type_t::none;
    variable_dimensions = -1;
  } else {
    const int groupindex0 = CCTK_GroupIndexFromVar(variables.front().c_str());
    assert(groupindex0 >= 0);
    const int grouptype0 = CCTK_GroupTypeI(groupindex0);
    switch (grouptype0) {
    case CCTK_GF:
      group_type = group_type_t::gf;
      break;
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      group_type = group_type_t::array;
      break;
    default:
      assert(0);
    }
    variable_dimensions = CCTK_GroupDimI(groupindex0);
    for (auto &var : variables) {
      const int varindex = CCTK_VarIndex(var.c_str());
      if (varindex < 0)
        CCTK_VERROR("Variable \"%s\" does not exist", var.c_str());
      for (auto &ch : var)
        ch = std::tolower(ch);
    }
  }

  iterations = ofd.iterations;
  for (const auto &iter : iterations)
    if (iter < 0)
      CCTK_VERROR("Bad iteration number %d", iter);

  if (!variables.empty()) {
    const int groupindex0 = CCTK_GroupIndexFromVar(variables.front().c_str());
    assert(groupindex0 >= 0);
    const int grouptype0 = CCTK_GroupTypeI(groupindex0);
    switch (grouptype0) {
    case CCTK_GF:
      group_type = group_type_t::gf;
      break;
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      group_type = group_type_t::array;
      break;
    default:
      assert(0);
    }
    variable_dimensions = CCTK_GroupDimI(groupindex0);
    for (auto &var : variables) {
      const int groupindex = CCTK_GroupIndexFromVar(var.c_str());
      assert(groupindex >= 0);
      const int grouptype = CCTK_GroupTypeI(groupindex);
      group_type_t group_type1;
      switch (grouptype) {
      case CCTK_GF:
        group_type1 = group_type_t::gf;
        break;
      case CCTK_SCALAR:
      case CCTK_ARRAY:
        group_type1 = group_type_t::array;
        break;
      default:
        assert(0);
      }
      if (group_type1 != group_type)
        CCTK_VERROR("Variable \"%s\" has wrong group type", var.c_str());
      const int group_dim1 = CCTK_GroupDimI(groupindex);
      if (group_dim1 != variable_dimensions)
        CCTK_VERROR("Variable \"%s\" has wrong dimension", var.c_str());
    }
  }

  // TODO: Check whether directions are unique
  output_directions = ofd.output_directions;
  for (const auto &dir : output_directions)
    if (dir < 0 || dir >= variable_dimensions)
      CCTK_VERROR("Bad output direction %d", dir);

  reductions = ofd.reductions;

  if (!output_directions.empty() && !reductions.empty())
    CCTK_ERROR("Cannot have both output directions and reductions");

  // TODO: Maintain a list of known output formats?
  format_name = ofd.format_name;
  if (format_name == "")
    CCTK_VERROR("Bad format name \"%s\"", format_name.c_str());

  // TODO: Maintain a list of known versions for each output format?
  format_version = ofd.format_version;
  for (const auto &ver : format_version)
    if (ver < 0)
      CCTK_VERROR("Bad format version %d.%d.%d", format_version[0],
                  format_version[1], format_version[2]);
}

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const real_output_file_description_t &rofd) {
  yaml << YAML::LocalTag("cactus-metadata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "filename" << YAML::Value << rofd.filename;
  yaml << YAML::Key << "description" << YAML::Value << rofd.description;
  yaml << YAML::Key << "writer_thorn" << YAML::Value << rofd.writer_thorn;
  yaml << YAML::Key << "variables" << YAML::Value << rofd.variables;
  if (!rofd.variables.empty()) {
    yaml << YAML::Key << "iterations" << YAML::Value << YAML::Flow
         << rofd.iterations;
    yaml << YAML::Key << "group_type" << YAML::Value << rofd.group_type;
    yaml << YAML::Key << "variable_dimensions" << YAML::Value
         << rofd.variable_dimensions;
    yaml << YAML::Key << "real_output_directions" << YAML::Value << YAML::Flow
         << rofd.output_directions;
    yaml << YAML::Key << "reductions" << YAML::Value << YAML::Flow
         << rofd.reductions;
  } // if !rofd.variables.empty
  yaml << YAML::Key << "format_name" << YAML::Value << rofd.format_name;
  yaml << YAML::Key << "format_version" << YAML::Value << YAML::Flow
       << YAML::BeginSeq;
  for (const auto &elt : rofd.format_version)
    yaml << elt;
  yaml << YAML::EndSeq;
  yaml << YAML::EndMap;
  return yaml;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<real_output_file_description_t> real_output_file_descriptions;

void OutputMeta_RegisterOutputFile(
    output_file_description_t output_file_description) {
  // Only the root process outputs metadata
  if (CCTK_MyProc(nullptr) != 0)
    return;

  real_output_file_descriptions.emplace_back(
      std::move(output_file_description));
}

void OutputMeta(const cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // Only the root process outputs metadata
  if (CCTK_MyProc(nullptr) != 0)
    return;

  // Skip metadata output if no files were written
  if (real_output_file_descriptions.empty())
    return;

  if (out_metadata) {
    YAML::Emitter yaml;
    yaml << YAML::Comment("Cactus Metadata");
    yaml << YAML::BeginDoc;
    yaml << real_output_file_descriptions;
    yaml << YAML::EndDoc;

    std::ostringstream buf;
    buf << out_dir << "/cactus-metadata"
        << ".it" << std::setw(8) << std::setfill('0') << cctkGH->cctk_iteration
        << ".yaml";
    const std::string filename = buf.str();
    std::ofstream file(filename);
    file << yaml.c_str() << "\n";
    file.close();
  }

  real_output_file_descriptions.clear();
}

} // namespace CarpetX
