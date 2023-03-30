#include "driver.hxx"
#include "io.hxx"
#include "io_adios2.hxx"
#include "io_meta.hxx"
#include "io_norm.hxx"
#include "io_openpmd.hxx"
#include "io_silo.hxx"
#include "io_tsv.hxx"
#include "schedule.hxx"
#include "timer.hxx"

#include <CactusBase/IOUtil/src/ioGH.h>
#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <AMReX.H>
#include <AMReX_Orientation.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <regex>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

namespace {
vector<bool> find_groups(const char *const method, const char *const out_vars) {
  DECLARE_CCTK_PARAMETERS;

  vector<bool> enabled(CCTK_NumGroups(), false);
  const auto callback{
      [](const int index, const char *const optstring, void *const arg) {
        vector<bool> &enabled = *static_cast<vector<bool> *>(arg);
        enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
      }};
  CCTK_TraverseString(out_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
  if (verbose) {
    CCTK_VINFO("%s output for groups:", method);
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
      if (enabled.at(gi))
        CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
  }
  return enabled;
}

string get_parameter_filename() {
  vector<char> buf(10000);
  int ilen = CCTK_ParameterFilename(buf.size(), buf.data());
  assert(ilen < int(buf.size() - 1));
  string parfilename(buf.data());
  // Remove directory prefix, if any
  auto slash = parfilename.rfind('/');
  if (slash != string::npos)
    parfilename = parfilename.substr(slash + 1);
  // Remove suffix, if it is there
  auto suffix = parfilename.rfind('.');
  if (suffix != string::npos && parfilename.substr(suffix) == ".par")
    parfilename = parfilename.substr(0, suffix);
  return parfilename;
}

string get_simulation_name() {
  string name = get_parameter_filename();
  const size_t last_slash = name.rfind('/');
  if (last_slash != string::npos && last_slash < name.length())
    name = name.substr(last_slash + 1);
  const size_t last_dot = name.rfind('.');
  if (last_dot != string::npos && last_dot > 0)
    name = name.substr(0, last_dot);
  return name;
}
} // namespace

////////////////////////////////////////////////////////////////////////////////

// Recovering

int recover_iteration = -1;

extern "C" int CarpetX_RecoverParameters() {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("RecoverParameters");
  Interval interval(timer);

  if (CCTK_EQUALS(recover_method, "openpmd")) {

#ifdef HAVE_CAPABILITY_openPMD_api
    recover_iteration = InputOpenPMDParameters(recover_dir, recover_file);
#else
    CCTK_VERROR(
        "CarpetX::recover_method is set to \"openpmd\", but openPMD_api "
        "is not enabled");
#endif

  } else if (CCTK_EQUALS(recover_method, "silo")) {

#ifdef HAVE_CAPABILITY_Silo
    if (recover_iteration < 0)
      recover_iteration = InputSiloParameters(recover_dir, recover_file);
#else
    CCTK_VERROR("CarpetX::recover_method is set to \"silo\", but Silo "
                "is not enabled");
#endif

  } else {
    CCTK_ERROR("unknown value for paramater CarpetX::recover_method");
  }

  return recover_iteration >= 0;
}

////////////////////////////////////////////////////////////////////////////////

void RecoverGridStructure(cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(recover_method, "openpmd")) {

#ifdef HAVE_CAPABILITY_openPMD_api
    InputOpenPMDGridStructure(cctkGH, recover_dir, recover_file,
                              recover_iteration);
#else
    CCTK_VERROR(
        "CarpetX::recover_method is set to \"openpmd\", but openPMD_api "
        "is not enabled");
#endif

  } else if (CCTK_EQUALS(recover_method, "silo")) {

#ifdef HAVE_CAPABILITY_Silo
    InputSiloGridStructure(cctkGH, recover_dir, recover_file,
                           recover_iteration);
#else
    CCTK_VERROR("CarpetX::recover_method is set to \"silo\", but Silo "
                "is not enabled");
#endif

  } else {
    CCTK_ERROR("unknown value for paramater CarpetX::recover_method");
  }
}

////////////////////////////////////////////////////////////////////////////////

void RecoverGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("RecoverGH");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root) {
    const int runtime = CCTK_RunTime(); // seconds
    CCTK_VINFO("RecoverGH: iteration %d, time %f, run time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600));
  }

  // Find input groups
  const vector<bool> group_enabled = [&] {
    vector<bool> enabled(CCTK_NumGroups(), false);
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      const GHExt::GlobalData &globaldata = ghext->globaldata;
      if (globaldata.arraygroupdata.at(gi)) {
        // grid array
        enabled.at(gi) = globaldata.arraygroupdata.at(gi)->do_checkpoint;
      } else {
        // grid function
        const int patch = 0;
        const int level = 0;
        const GHExt::PatchData::LevelData &leveldata =
            ghext->patchdata.at(patch).leveldata.at(level);
        assert(leveldata.groupdata.at(gi));
        enabled.at(gi) = leveldata.groupdata.at(gi)->do_checkpoint;
      }
    }
    return enabled;
  }();

  if (CCTK_EQUALS(recover_method, "openpmd")) {

#ifdef HAVE_CAPABILITY_openPMD_api
    InputOpenPMD(cctkGH, group_enabled, recover_dir, recover_file);
#else
    CCTK_VERROR(
        "CarpetX::recover_method is set to \"openpmd\", but openPMD_api "
        "is not enabled");
#endif

  } else if (CCTK_EQUALS(recover_method, "silo")) {

#ifdef HAVE_CAPABILITY_Silo
    InputSilo(cctkGH, group_enabled, recover_dir, recover_file);
#else
    CCTK_VERROR("CarpetX::recover_method is set to \"silo\", but Silo "
                "is not enabled");
#endif

  } else {
    CCTK_ERROR("unknown value for paramater CarpetX::recover_method");
  }
}

void InputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("InputGH");
  Interval interval(timer);

  if (filereader_ID_files[0] == '\0')
    return;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root) {
    const int runtime = CCTK_RunTime(); // seconds
    CCTK_VINFO("InputGH: iteration %d, time %f, run time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600.0));
  }

  if (CCTK_EQUALS(filereader_method, "openpmd")) {

#ifdef HAVE_CAPABILITY_openPMD_api
    // TODO: handle multiple file names in `filereader_ID_files`
    const vector<bool> input_group = find_groups("openPMD", filereader_ID_vars);
    const string simulation_name = get_simulation_name();
    InputOpenPMD(cctkGH, input_group, filereader_ID_dir, filereader_ID_files);
#else
    // TODO: Check this at paramcheck
    CCTK_VERROR(
        "CarpetX::filereader_method is set to \"openpmd\", but openPMD_api "
        "is not enabled");
#endif

  } else if (CCTK_EQUALS(filereader_method, "silo")) {

    const vector<bool> input_group = find_groups("Silo", filereader_ID_vars);
#ifdef HAVE_CAPABILITY_Silo
    // TODO: Stop at paramcheck time when Silo input parameters are
    // set, but Silo is not available
    // TODO: handle multiple file names in `filereader_ID_files`
    const string simulation_name = get_simulation_name();
    InputSilo(cctkGH, input_group, filereader_ID_dir, filereader_ID_files);
#else
    // TODO: Check this at paramcheck
    CCTK_VERROR("CarpetX::filereader_method is set to \"silo\", but Silo is "
                "not enabled");
#endif

  } else {
    CCTK_ERROR("unknown value for paramater CarpetX::filereader_method");
  }
}

////////////////////////////////////////////////////////////////////////////////

void OutputPlotfile(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputPlotfile");
  Interval interval(timer);

  const int numgroups = CCTK_NumGroups();
  vector<bool> group_enabled(numgroups, false);
  auto enable_group{[](int index, const char *optstring, void *callback) {
    vector<bool> &group_enabled = *static_cast<vector<bool> *>(callback);
    // CCTK_TraverseString expands the given groups into their variables
    // so I condense it down again to groups only
    group_enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
  }};
  CCTK_TraverseString(out_plotfile_groups, enable_group, &group_enabled,
                      CCTK_GROUP_OR_VAR);

  // TODO: Set the number of output files (the number of files written
  // in parallel) to reduce the total number of files
  // VisMF::SetNOutFiles(saveNFiles);

  for (int gi = 0; gi < numgroups; ++gi) {
    if (!group_enabled.at(gi))
      continue;
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;

    auto &restrict groupdata0 =
        *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      string groupname = CCTK_FullGroupName(gi);
      groupname = regex_replace(groupname, regex("::"), "-");
      for (auto &c : groupname)
        c = tolower(c);

      for (const auto &restrict patchdata : ghext->patchdata) {

        const string basename = [&]() {
          ostringstream buf;
          buf << groupname << ".it" << setw(8) << setfill('0')
              << cctk_iteration;
          if (ghext->num_patches() > 1)
            buf << ".m" << setw(2) << setfill('0') << patchdata.patch;
          return buf.str();
        }();
        const string filename = string(out_dir) + "/" + basename;

        amrex::Vector<string> varnames(groupdata0.numvars);
        for (int vi = 0; vi < groupdata0.numvars; ++vi) {
          ostringstream buf;
          buf << CCTK_VarName(groupdata0.firstvarindex + vi);
          for (int i = 0; i < tl; ++i)
            buf << "_p";
          varnames.at(vi) = buf.str();
        }

        amrex::Vector<const amrex::MultiFab *> mfabs(
            patchdata.leveldata.size());
        amrex::Vector<amrex::Geometry> geoms(patchdata.leveldata.size());
        amrex::Vector<int> iters(patchdata.leveldata.size());
        amrex::Vector<amrex::IntVect> reffacts(patchdata.leveldata.size());
        for (const auto &restrict leveldata : patchdata.leveldata) {
          mfabs.at(leveldata.level) = &*leveldata.groupdata.at(gi)->mfab.at(tl);
          geoms.at(leveldata.level) = patchdata.amrcore->Geom(leveldata.level);
          iters.at(leveldata.level) = cctk_iteration;
          reffacts.at(leveldata.level) = amrex::IntVect{2, 2, 2};
        }

        // TODO: Output all groups into a single file
        WriteMultiLevelPlotfile(filename, mfabs.size(), mfabs, varnames, geoms,
                                cctk_time, iters, reffacts);

        const bool is_root = CCTK_MyProc(nullptr) == 0;
        if (is_root) {
          const string visitname = [&]() {
            ostringstream buf;
            buf << out_dir << "/" << groupname << ".visit";
            return buf.str();
          }();
          ofstream visit(visitname, ios::app);
          assert(visit.good());
          visit << basename << "/Header\n";
          visit.close();
        } // if is_root
      }   // for patchdata
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

struct parameters {};
YAML::Emitter &operator<<(YAML::Emitter &yaml, parameters) {
  // Collect all parameters and their values
  std::vector<std::pair<std::string, const cParamData *> > parameter_values;
  int first = 1;
  for (;;) {
    const cParamData *data;
    // This call is not thread safe
    const int iret = CCTK_ParameterWalk(first, nullptr, nullptr, &data);
    assert(iret >= 0);
    if (iret)
      break;
    first = 0;
    // Ignore inactive thorns
    if (!CCTK_IsThornActive(data->thorn))
      continue;
    const std::string fullname = std::string(data->thorn) + "::" + data->name;
    parameter_values.emplace_back(fullname, data);
  }

  // Sort parameters alphabetically
  std::sort(parameter_values.begin(), parameter_values.end());

  // Output parameters
  yaml << YAML::LocalTag("parameters-1.0.0");
  yaml << YAML::BeginMap;
  for (const auto &parameter_value : parameter_values) {
    const std::string &fullname = parameter_value.first;
    const cParamData *const data = parameter_value.second;
    yaml << YAML::Key << fullname << YAML::Value << YAML::Flow;
    int type;
    const void *const pvalue =
        CCTK_ParameterGet(data->name, data->thorn, &type);
    switch (type) {
    case PARAMETER_KEYWORD:
    case PARAMETER_STRING:
      yaml << *static_cast<const char *const *>(pvalue);
      break;
    case PARAMETER_INT:
      yaml << *static_cast<const CCTK_INT *>(pvalue);
      break;
    case PARAMETER_REAL:
      yaml << *static_cast<const CCTK_REAL *>(pvalue);
      break;
    case PARAMETER_BOOLEAN:
      yaml << static_cast<bool>(*static_cast<const CCTK_INT *>(pvalue));
      break;
    default:
      assert(0);
    }
  }
  yaml << YAML::EndMap;

  return yaml;
}

void OutputMetadata(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!out_metadata)
    return;

  static Timer timer("OutputMetadata");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (!is_root)
    return;

  {
    YAML::Emitter yaml;
    yaml << YAML::Comment("CarpetX Metadata");
    yaml << YAML::BeginDoc;
    yaml << YAML::LocalTag("carpetx-metadata-1.0.0");
    yaml << YAML::BeginMap;
    yaml << YAML::Key << "nghostzones";
    yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
    for (int d = 0; d < dim; ++d)
      yaml << cctk_nghostzones[d];
    yaml << YAML::EndSeq;
    // Note: origin_space and delta_space are nan at this point
    // yaml << YAML::Key << "origin_space";
    // yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
    // for (int d = 0; d < dim; ++d)
    //   yaml << cctk_origin_space[d];
    // yaml << YAML::EndSeq;
    // yaml << YAML::Key << "delta_space";
    // yaml << YAML::Value << YAML::Flow << YAML::BeginSeq;
    // for (int d = 0; d < dim; ++d)
    //   yaml << cctk_delta_space[d];
    // yaml << YAML::EndSeq;
    yaml << YAML::Key << "iteration";
    yaml << YAML::Value << cctk_iteration;
    yaml << YAML::Key << "time";
    yaml << YAML::Value << cctk_time;
    yaml << YAML::Key << "delta_time";
    yaml << YAML::Value << cctk_delta_time;
    yaml << YAML::Key << "ghext";
    yaml << YAML::Value << *ghext;
    // yaml << YAML::Key << "parameters";
    // yaml << YAML::Value << parameters();
    yaml << YAML::EndMap;
    yaml << YAML::EndDoc;

    ostringstream buf;
    buf << out_dir << "/carpetx-metadata"
        << ".it" << setw(8) << setfill('0') << cctk_iteration << ".yaml";
    const string filename = buf.str();
    ofstream file(filename.c_str(), std::ofstream::out);
    file << yaml.c_str();
  }

  {
    YAML::Emitter yaml;
    yaml << YAML::Comment("Parameters");
    yaml << YAML::BeginDoc;
    yaml << parameters();
    yaml << YAML::EndDoc;

    ostringstream buf;
    buf << out_dir << "/parameters"
        << ".it" << setw(8) << setfill('0') << cctk_iteration << ".yaml";
    const string filename = buf.str();
    ofstream file(filename.c_str(), std::ofstream::out);
    file << yaml.c_str();
  }
}

int OutputGH(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("OutputGH");
  Interval interval(timer);

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root)
    CCTK_VINFO("OutputGH: iteration %d, time %f, run time %d s", cctk_iteration,
               double(cctk_time), CCTK_RunTime());

  {
    const int every = out_metadata_every == -1 ? out_every : out_metadata_every;
    if (every > 0 && cctk_iteration % every == 0) {
      OutputMetadata(cctkGH);
    }
  }

  {
    const int every = out_norm_every == -1 ? out_every : out_norm_every;
    if (every > 0 && cctk_iteration % every == 0) {
      OutputNorms(cctkGH);
    }
  }

  {
    const int every = out_adios2_every == -1 ? out_every : out_adios2_every;
    if (every > 0 && cctk_iteration % every == 0) {
      const vector<bool> group_enabled = find_groups("ADIOS2", out_adios2_vars);
#ifdef HAVE_CAPABILITY_ADIOS2
      // TODO: Stop at paramcheck time when ADIOS2 output parameters
      // are set, but ADIOS2 is not available
      const string simulation_name = get_simulation_name();
      OutputADIOS2(cctkGH, group_enabled, out_dir, simulation_name);
#else
      if (strlen(out_adios2_vars) != 0)
        CCTK_VERROR("ADIOS2 is not enabled. The parameter "
                    "CarpetX::out_adios2_vars must be empty.");
#endif
    }

    {
      const int every = out_openpmd_every == -1 ? out_every : out_openpmd_every;
      if (every > 0 && cctk_iteration % every == 0) {
        const vector<bool> group_enabled =
            find_groups("openPMD", out_openpmd_vars);
#ifdef HAVE_CAPABILITY_openPMD_api
        // TODO: Stop at paramcheck time when openPMD output parameters
        // are set, but openPMD is not available
        const string simulation_name = get_simulation_name();
        OutputOpenPMD(cctkGH, group_enabled, out_dir, simulation_name);
#else
        if (strlen(out_openpmd_vars) != 0)
          CCTK_VERROR("openPMD is not enabled. The parameter "
                      "CarpetX::out_openpmd_vars must be empty.");
#endif
      }
    }

    {
      const int every =
          out_plotfile_every == -1 ? out_every : out_plotfile_every;
      if (every > 0 && cctk_iteration % every == 0) {
        OutputPlotfile(cctkGH);
      }
    }

    {
      const int every = out_silo_every == -1 ? out_every : out_silo_every;
      if (every > 0 && cctk_iteration % every == 0) {
        const vector<bool> group_enabled = find_groups("Silo", out_silo_vars);
#ifdef HAVE_CAPABILITY_Silo
        // TODO: Stop at paramcheck time when Silo output parameters are
        // set, but Silo is not available
        const string simulation_name = get_simulation_name();
        OutputSilo(cctkGH, group_enabled, out_dir, simulation_name);
#else
        if (strlen(out_silo_vars) != 0)
          CCTK_VERROR("Silo is not enabled. The parameter "
                      "CarpetX::out_silo_vars must be empty.");
#endif
      }
    }

    OutputTSVold(cctkGH);

    {
      const int every = out_tsv_every == -1 ? out_every : out_tsv_every;
      if (every > 0 && cctk_iteration % every == 0) {
        OutputTSV(cctkGH);
      }
    }

    // Describe all output files
    OutputMeta(cctkGH);

    if (is_root)
      CCTK_VINFO("OutputGH done.");
  }

  // TODO: This should be the number of variables output
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

void Checkpoint(const cGH *const restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static int last_checkpoint_iteration = -1;

  if (cctkGH->cctk_iteration == recover_iteration) {
    CCTK_VINFO(
        "Recovered from checkpoint at iteration %d; skipping checkpointing",
        cctkGH->cctk_iteration);
    return;
  }

  if (cctkGH->cctk_iteration <= last_checkpoint_iteration) {
    CCTK_VINFO(
        "Already wrote checkpoint at iteration %d; skipping checkpointing",
        cctkGH->cctk_iteration);
    return;
  }
  last_checkpoint_iteration = cctkGH->cctk_iteration;

  static Timer timer("Checkpoint");
  Interval interval(timer);

  if (CCTK_EQUALS(checkpoint_method, "openpmd")) {

#ifdef HAVE_CAPABILITY_openPMD_api
    const vector<bool> checkpoint_group = [&] {
      vector<bool> enabled(CCTK_NumGroups(), false);
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        const GHExt::GlobalData &globaldata = ghext->globaldata;
        if (globaldata.arraygroupdata.at(gi)) {
          // grid array
          enabled.at(gi) = globaldata.arraygroupdata.at(gi)->do_checkpoint;
        } else {
          // grid function
          const int patch = 0;
          const GHExt::PatchData &patchdata = ghext->patchdata.at(patch);
          const int level = 0;
          const GHExt::PatchData::LevelData &leveldata =
              patchdata.leveldata.at(level);
          assert(leveldata.groupdata.at(gi));
          enabled.at(gi) = leveldata.groupdata.at(gi)->do_checkpoint;
        }
      }
      return enabled;
    }();
    OutputOpenPMD(cctkGH, checkpoint_group, checkpoint_dir, checkpoint_file);
#else
    // TODO: Check this at paramcheck
    CCTK_VERROR(
        "CarpetX::checkpoint_method is set to \"openpmd\", but openPMD_api is "
        "not enabled");
#endif

  } else if (CCTK_EQUALS(checkpoint_method, "silo")) {

#ifdef HAVE_CAPABILITY_Silo
    const vector<bool> checkpoint_group = [&] {
      vector<bool> enabled(CCTK_NumGroups(), false);
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        const GHExt::GlobalData &globaldata = ghext->globaldata;
        if (globaldata.arraygroupdata.at(gi)) {
          // grid array
          enabled.at(gi) = globaldata.arraygroupdata.at(gi)->do_checkpoint;
        } else {
          // grid function
          const int patch = 0;
          const GHExt::PatchData &patchdata = ghext->patchdata.at(patch);
          const int level = 0;
          const GHExt::PatchData::LevelData &leveldata =
              patchdata.leveldata.at(level);
          assert(leveldata.groupdata.at(gi));
          enabled.at(gi) = leveldata.groupdata.at(gi)->do_checkpoint;
        }
      }
      return enabled;
    }();
    OutputSilo(cctkGH, checkpoint_group, checkpoint_dir, checkpoint_file);
#else
    // TODO: Check this at paramcheck
    CCTK_VERROR("CarpetX::checkpoint_method is set to \"silo\", but Silo is "
                "not enabled");
#endif

  } else {
    CCTK_ERROR("unknown value for paramater CarpetX::checkpoint_method");
  }
}

int last_checkpoint_runtime = -1; // seconds

extern "C" void CarpetX_CheckpointInitial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int runtime = CCTK_RunTime(); // seconds

  if (checkpoint_ID) {
    CCTK_VINFO("Checkpointing initial conditions at iteration %d, time %f, run "
               "time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600.0));
    Checkpoint(cctkGH);
    last_checkpoint_runtime = runtime;
  }
}

extern "C" void CarpetX_Checkpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int runtime = CCTK_RunTime(); // seconds

  const bool checkpoint_by_iteration =
      checkpoint_every > 0 && cctk_iteration % checkpoint_every == 0;
  const bool checkpoint_by_walltime =
      checkpoint_every_walltime_hours > 0 &&
      runtime >= last_checkpoint_runtime +
                     lrint(checkpoint_every_walltime_hours * 3600);

  if (checkpoint_by_iteration || checkpoint_by_walltime) {
    CCTK_VINFO("Checkpointing at iteration %d, time %f, run time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600.0));
    Checkpoint(cctkGH);
    last_checkpoint_runtime = runtime;
  }
}

extern "C" void CarpetX_CheckpointTerminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (checkpoint_on_terminate) {
    const int runtime = CCTK_RunTime(); // seconds
    CCTK_VINFO("Checkpointing before terminating at iteration %d, time %f, run "
               "time %.2f h",
               cctk_iteration, double(cctk_time), double(runtime / 3600.0));
    Checkpoint(cctkGH);
  }
}

} // namespace CarpetX
