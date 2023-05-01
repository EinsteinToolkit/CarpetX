#include "io_norm.hxx"

#include "driver.hxx"
#include "io_meta.hxx"
#include "loop.hxx"
#include "reduction.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <limits>
#include <mutex>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace CarpetX {

void OutputNorms(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_norm_vars[0] == '\0')
    return;

  static Timer timer("OutputNorms");
  Interval interval(timer);

  // Find output groups
  const std::vector<bool> group_enabled = [&] {
    std::vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          std::vector<bool> &enabled = *static_cast<std::vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(out_norm_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    if (verbose) {
      CCTK_VINFO("TSV output for groups:");
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
        if (group_enabled.at(gi))
          CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
    }
    return enabled;
  }();
  const auto num_out_groups =
      std::count(group_enabled.begin(), group_enabled.end(), true);
  if (num_out_groups == 0)
    return;

  static std::vector<bool> previous_group_enabled;
  const bool group_enabled_changed = group_enabled != previous_group_enabled;
  if (group_enabled_changed)
    previous_group_enabled = group_enabled;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  const std::string sep = "\t";
  const int tl = 0;
  const int numgroups = CCTK_NumGroups();
  std::vector<std::string> reductions = {"min",    "max",    "sum",
                                         "avg",    "stddev", "volume",
                                         "L1norm", "L2norm", "maxabs"};
  if (!out_norm_omit_unstable) {
    for (int d = 0; d < dim; ++d) {
      std::ostringstream buf;
      buf << "minloc[" << d << "]";
      reductions.push_back(buf.str());
    }
    for (int d = 0; d < dim; ++d) {
      std::ostringstream buf;
      buf << "maxloc[" << d << "]";
      reductions.push_back(buf.str());
    }
  }

  if (is_root) {
    static std::once_flag create_directory;
    std::call_once(create_directory, [&]() {
      const int mode = 0755;
      int ierr = CCTK_CreateDirectory(mode, out_dir);
      assert(ierr >= 0);
      std::ostringstream buf;
      buf << out_dir << "/norms";
      ierr = CCTK_CreateDirectory(mode, buf.str().c_str());
      assert(ierr >= 0);
    });
  }

  if (is_root) {
    if (group_enabled_changed) {
      for (int gi = 0; gi < numgroups; ++gi) {
        if (!group_enabled.at(gi))
          continue;
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;

        const int patch = 0;
        const GHExt::PatchData &restrict patchdata = ghext->patchdata.at(patch);
        const int level = 0;
        const GHExt::PatchData::LevelData &restrict leveldata =
            patchdata.leveldata.at(level);
        const GHExt::PatchData::LevelData::GroupData &restrict groupdata =
            *leveldata.groupdata.at(gi);

        bool group_has_valid_variables = false;
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          // Only output variables with valid data
          if (!groupdata.valid.at(tl).at(vi).get().valid_int)
            continue;
          group_has_valid_variables = true;
        }
        if (!group_has_valid_variables)
          continue;

        std::ostringstream buf;
        std::string groupname = CCTK_FullGroupName(gi);
        groupname = std::regex_replace(groupname, std::regex("::"), "-");
        for (auto &ch : groupname)
          ch = std::tolower(ch);
        buf << out_dir << "/norms/" << groupname << ".tsv";
        const std::string filename = buf.str();
        std::ofstream file;
        file.open(filename, std::ios_base::app);

        int col = 0;
        file << "# " << ++col << ":iteration";
        file << sep << ++col << ":time";

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          // Only output variables with valid data
          if (!groupdata.valid.at(tl).at(vi).get().valid_int)
            continue;

          std::string varname = CCTK_FullVarName(groupdata.firstvarindex + vi);
          for (auto &ch : varname)
            ch = std::tolower(ch);
          for (const auto &reduction : reductions)
            file << sep << ++col << ":" << varname << "." << reduction;
        }
        file << "\n";

        file.close();
      }
    }
  }

  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    const int patch = 0;
    const GHExt::PatchData &restrict patchdata = ghext->patchdata.at(patch);
    const int level = 0;
    const GHExt::PatchData::LevelData &restrict leveldata =
        patchdata.leveldata.at(level);
    const GHExt::PatchData::LevelData::GroupData &restrict groupdata =
        *leveldata.groupdata.at(gi);

    bool group_has_valid_variables = false;
    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      // Only output variables with valid data
      if (!groupdata.valid.at(tl).at(vi).get().valid_int)
        continue;
      group_has_valid_variables = true;
    }
    if (!group_has_valid_variables)
      continue;

    std::ofstream file;
    output_file_description_t ofd;
    if (is_root) {
      std::ostringstream buf;
      std::string groupname = CCTK_FullGroupName(gi);
      groupname = std::regex_replace(groupname, std::regex("::"), "-");
      for (auto &ch : groupname)
        ch = std::tolower(ch);
      buf << out_dir << "/norms/" << groupname << ".tsv";
      const std::string filename = buf.str();
      file.open(filename, std::ios_base::app);
      ofd.filename = filename;

      // get more precision for floats, could also use
      // https://stackoverflow.com/a/30968371
      file << setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1)
           << scientific;

      file << cctk_iteration << sep << cctk_time;
    }

    for (int vi = 0; vi < groupdata.numvars; ++vi) {
      // Only output variables with valid data
      if (!groupdata.valid.at(tl).at(vi).get().valid_int)
        continue;

      ofd.variables.push_back(CCTK_FullVarName(groupdata.firstvarindex + vi));

      const reduction<CCTK_REAL, dim> red = reduce(gi, vi, tl);

      if (is_root) {
        file << sep << red.min << sep << red.max << sep << red.sum << sep
             << red.avg() << sep << red.sdv() << sep << red.norm0() << sep
             << red.norm1() << sep << red.norm2() << sep << red.norm_inf();
        if (!out_norm_omit_unstable) {
          for (int d = 0; d < dim; ++d)
            file << sep << red.minloc[d];
          for (int d = 0; d < dim; ++d)
            file << sep << red.maxloc[d];
        }
      }
    }

    if (is_root) {
      file << "\n";
      file.close();

      ofd.description = "CarpetX TSV norms output";
      ofd.writer_thorn = CCTK_THORNSTRING;
      ofd.iterations = {cctk_iteration};
      ofd.reductions = {
          reduction_t::minimum,
          reduction_t::maximum,
          reduction_t::sum,
          reduction_t::average,
          reduction_t::standard_deviation,
          reduction_t::volume,
          reduction_t::norm1,
          reduction_t::norm2,
          reduction_t::norm_inf,
      };
      if (!out_norm_omit_unstable) {
        ofd.reductions.push_back(reduction_t::minimum_location);
        ofd.reductions.push_back(reduction_t::maximum_location);
      }
      ofd.format_name = "CarpetX/norms/TSV";
      ofd.format_version = {1, 0, 0};

      OutputMeta_RegisterOutputFile(std::move(ofd));
    }
  }
}

} // namespace CarpetX
