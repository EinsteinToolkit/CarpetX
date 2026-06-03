#include "io_silo.hxx"

#include "driver.hxx"
#include "io_meta.hxx"
#include "mpi_types.hxx"
#include "timer.hxx"

#include <tuple.hxx>

#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef HAVE_CAPABILITY_Silo

#include <AMReX.H>
#include <AMReX_IntVect.H>

#include <mpi.h>

#include <silo.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;
#else
#include <filesystem>
namespace filesystem = std::filesystem;
#endif
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {

namespace {

constexpr bool io_verbose = true;

// // Compile-time if-then-else expression
// template <bool cond, typename T, typename F>
// constexpr std::enable_if_t<cond, T> ifelse(T &&t, F &&f) {
//   return std::forward<T>(t);
// }
// template <bool cond, typename T, typename F>
// constexpr std::enable_if_t<!cond, F> ifelse(T &&t, F &&f) {
//   return std::forward<F>(f);
// }

template <typename T> struct db_datatype;
template <> struct db_datatype<char> : std::integral_constant<int, DB_CHAR> {};
template <>
struct db_datatype<short> : std::integral_constant<int, DB_SHORT> {};
template <> struct db_datatype<int> : std::integral_constant<int, DB_INT> {};
template <> struct db_datatype<long> : std::integral_constant<int, DB_LONG> {};
template <>
struct db_datatype<long long> : std::integral_constant<int, DB_LONG_LONG> {};
template <>
struct db_datatype<float> : std::integral_constant<int, DB_FLOAT> {};
template <>
struct db_datatype<double> : std::integral_constant<int, DB_DOUBLE> {};
template <typename T> constexpr int db_datatype_v = db_datatype<T>::value;

// Slice membership + slab bounds for component `c`, shared by the data-write
// and metadata passes so they skip the same components in the same order.
// Returns the thickness-1 slab (or nullopt if `c` misses the plane), tested
// against the interior box so a box-boundary plane belongs to one component.
inline std::optional<std::pair<Arith::vect<int, 3>, Arith::vect<int, 3> > >
slice_slab(const amrex::MultiFab &mfab, const int c,
           const amrex::Real *const x0, const amrex::Real *const dx,
           const amrex::IndexType &indextype, const slice_t &slice) {
  const amrex::Box &fabbox = mfab.fabbox(c); // exterior
  const Arith::vect<int, 3> ext_lo{fabbox.smallEnd(0), fabbox.smallEnd(1),
                                   fabbox.smallEnd(2)};
  const Arith::vect<int, 3> ext_hi{fabbox.bigEnd(0) + 1, fabbox.bigEnd(1) + 1,
                                   fabbox.bigEnd(2) + 1};
  const amrex::Box &validbox = mfab.box(c); // interior (variable-typed)
  // enclosedCells maps a nodal or cell valid box to the same cell range, so
  // membership (which box owns the plane) is centering-independent.
  const amrex::Box cvalidbox = amrex::enclosedCells(validbox);
  const Arith::vect<int, 3> cval_lo{
      cvalidbox.smallEnd(0), cvalidbox.smallEnd(1), cvalidbox.smallEnd(2)};
  const Arith::vect<int, 3> cval_hi{cvalidbox.bigEnd(0) + 1,
                                    cvalidbox.bigEnd(1) + 1,
                                    cvalidbox.bigEnd(2) + 1};
  const Arith::vect<CCTK_REAL, 3> sx0{x0[0], x0[1], x0[2]};
  const Arith::vect<CCTK_REAL, 3> sdx{dx[0], dx[1], dx[2]};
  const Arith::vect<bool, 3> is_cell_centred{indextype.cellCentered(0),
                                             indextype.cellCentered(1),
                                             indextype.cellCentered(2)};
  return slice.restrict_component(ext_lo, ext_hi, cval_lo, cval_hi, sx0, sdx,
                                  is_cell_centred);
}

struct mesh_props_t {
  std::array<int, dim> nghosts;

  friend bool operator==(const mesh_props_t &p, const mesh_props_t &q) {
    return p.nghosts == q.nghosts;
  }
  friend bool operator<(const mesh_props_t &p, const mesh_props_t &q) {
    return p.nghosts < q.nghosts;
  }
};

std::string make_subdirname(const std::string &file_name, const int iteration) {
  std::ostringstream buf;
  buf << file_name                                     //
      << ".it" << setw(8) << setfill('0') << iteration //
      << ".silo.dir";
  return buf.str();
}

std::string make_filename(const std::string &file_name, const int iteration,
                          const int ioserver = -1) {
  std::ostringstream buf;
  buf << file_name //
      << ".it" << setw(8) << setfill('0') << iteration;
  if (ioserver >= 0)
    buf << ".p" << setw(6) << setfill('0') << ioserver;
  buf << ".silo";
  return buf.str();
}

int match_filename(const std::string &file_name) {
  const std::regex filename_regex("[.]it(\\d+)[.]silo$",
                                  regex_constants::ECMAScript);
  std::smatch sm;
  regex_search(file_name, sm, filename_regex);
  if (sm.empty())
    return -1;
  return std::stoi(sm.str(1));
}

// The optional `band` tag namespaces subcycling consumer-band meshes/vars
// apart from the regular tl=0 data (see subcycling checkpoint/recovery). An
// empty tag (the default) leaves regular-data names byte-identical.
std::string make_meshname(const std::array<int, dim> &nghosts,
                          const int patch = -1, const int reflevel = -1,
                          const int component = -1,
                          const std::string &band = std::string()) {
  assert((patch == -1) == (reflevel == -1));
  assert((patch == -1) == (component == -1));
  std::ostringstream buf;
  if (patch >= 0)
    buf << "box";
  else
    buf << "gh";
  buf << ".ghosts";
  for (int d = 0; d < dim; ++d) {
    if (d > 0)
      buf << "_";
    buf << setw(2) << setfill('0') << nghosts[d];
  }
  if (patch >= 0)
    buf << ".m" << setw(4) << setfill('0') << patch     //
        << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  if (!band.empty())
    buf << ".band." << band;
  return DB::legalize_name(buf.str());
}

std::string make_varname(const int gi, const int vi, const int patch = -1,
                         const int reflevel = -1, const int component = -1,
                         const std::string &band = std::string()) {
  assert((patch == -1) == (reflevel == -1));
  assert((patch == -1) == (component == -1));
  std::string varname;
  if (vi < 0) {
    assert(0);
    varname = CCTK_FullGroupName(gi);
  } else {
    const int v0 = CCTK_FirstVarIndexI(gi);
    varname = CCTK_FullVarName(v0 + vi);
    varname = regex_replace(varname, regex("::"), "-");
    for (auto &ch : varname)
      ch = tolower(ch);
  }
  std::ostringstream buf;
  buf << varname;
  if (reflevel >= 0)
    buf << ".m" << setw(4) << setfill('0') << patch     //
        << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  if (!band.empty())
    buf << ".band." << band;
  return DB::legalize_name(buf.str());
}

const std::string driver_name = "CarpetX";

std::string make_fabarraybasename(const int patch, const int reflevel) {
  std::ostringstream buf;
  buf << "FabArrayBase"                           //
      << ".m" << setw(4) << setfill('0') << patch //
      << ".rl" << setw(2) << setfill('0') << reflevel;
  return DB::legalize_name(buf.str());
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

int InputSiloParameters(const std::string &input_dir,
                        const std::string &input_file) {
  DECLARE_CCTK_PARAMETERS;

  assert(!input_dir.empty());
  assert(!input_file.empty());

  // Set up timers
  static Timer timer("InputSiloParameters");
  Interval interval(timer);

  static Timer timer_setup("InputSiloParameters.setup");
  auto interval_setup = std::make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(nullptr);

  const int metafile_ioproc = 0;
  const bool read_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);
  DBSetAllowEmptyObjects(1);

  // TODO: directories instead of carefully chosen names

  static Timer timer_parameters("InputSiloParameters.parameters");
  auto interval_parameters = std::make_unique<Interval>(timer_parameters);

  // Find checkpoint input iteration
  int input_iteration = -1;
  if (read_metafile) {

    // Find latest iteration (if any)
    try {
      for (const auto &direntry : filesystem::directory_iterator(input_dir)) {
        const auto &filename = direntry.path().filename().string();
        const int iter = match_filename(filename);
        input_iteration = max(input_iteration, iter);
      }
    } catch (const filesystem::filesystem_error &) {
      // do nothing if directory does not exist
    }
  }

  MPI_Bcast(&input_iteration, 1, MPI_INT, metafile_ioproc, mpi_comm);

  if (input_iteration < 0) {
    // Did not find a checkpoint file
    if (io_verbose) {
      CCTK_VINFO("Not recovering parameters:");
      CCTK_VINFO("  Could not find a Silo checkpoint file \"%s/%s.*.silo\"",
                 input_dir.c_str(), input_file.c_str());
    }
    return input_iteration;
  }

  // Read metadata
  std::string parameters;
  if (read_metafile) {
    CCTK_VINFO(
        "Recovering parameters from checkpoint file \"%s\" iteration %d",
        make_filename(input_dir + "/" + input_file, input_iteration).c_str(),
        input_iteration);

    const std::string metafilename =
        input_dir + "/" + make_filename(input_file, input_iteration);
    // We could use DB_UNKNOWN instead of DB_HDF5
    // assert(metafilename.size() < 256);
    const DB::ptr<DBfile> metafile =
        DB::make(DBOpen(metafilename.c_str(), DB_HDF5, DB_READ));
    assert(metafile);

    // Read parameters
    {
      // TODO: Put this into a "drivername" subdirectory?
      const int type = DBGetVarType(metafile.get(), "AllParameters");
      assert(type == DB_CHAR);

      int length;
      int ndims = DBGetVarDims(metafile.get(), "AllParameters", 1, &length);
      assert(ndims == 1);

      char *const data =
          static_cast<char *>(DBGetVar(metafile.get(), "AllParameters"));
      assert(data);
      parameters = std::string(data, length);
      std::free(data);
    }

    int length = parameters.size();
    MPI_Bcast(&length, 1, MPI_INT, metafile_ioproc, mpi_comm);
    MPI_Bcast(parameters.data(), length, MPI_CHAR, metafile_ioproc, mpi_comm);

  } else {

    int length;
    MPI_Bcast(&length, 1, MPI_INT, metafile_ioproc, mpi_comm);
    parameters = std::string(length, ' ');
    MPI_Bcast(parameters.data(), length, MPI_CHAR, metafile_ioproc, mpi_comm);
  }

  IOUtil_SetAllParameters(parameters.data());

  interval_parameters = nullptr;

  return input_iteration;
}

void InputSiloGridStructure(cGH *restrict const cctkGH,
                            const std::string &input_dir,
                            const std::string &input_file,
                            const int input_iteration) {
  DECLARE_CCTK_PARAMETERS;

  assert(!input_dir.empty());
  assert(!input_file.empty());
  assert(input_iteration >= 0);

  int ierr;

  // Set up timers
  static Timer timer("InputSiloGridStructure");
  Interval interval(timer);

  static Timer timer_setup("InputSiloGridStructure.setup");
  auto interval_setup = std::make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);

  const int metafile_ioproc = 0;
  const bool read_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_meta("InputSiloGridStructure.meta");
  auto interval_meta = std::make_unique<Interval>(timer_meta);

  // TODO: Recover grid structure for grid arrays; call SetupGlobals for them?

  // Read metadata

  DB::ptr<DBfile> metafile;
  if (read_metafile) {
    const std::string metafilename =
        input_dir + "/" + make_filename(input_file, input_iteration);
    // We could use DB_UNKNOWN instead of DB_HDF5
    metafile = DB::make(DBOpen(metafilename.c_str(), DB_HDF5, DB_READ));
    assert(metafile);
  }

  if (read_metafile) {
    const int type = DBGetVarType(metafile.get(), "cycle");
    assert(type == DB_INT);
    int dim;
    int ndims = DBGetVarDims(metafile.get(), "cycle", 1, &dim);
    assert(ndims == 1);
    assert(dim == 1);
    ierr = DBReadVar(metafile.get(), "cycle", &cctkGH->cctk_iteration);
    assert(!ierr);
  }
  MPI_Bcast(&cctkGH->cctk_iteration, 1, MPI_INT, metafile_ioproc, mpi_comm);

  double dtime;
  if (read_metafile) {
    const int type = DBGetVarType(metafile.get(), "dtime");
    assert(type == DB_DOUBLE);
    int dim;
    int ndims = DBGetVarDims(metafile.get(), "dtime", 1, &dim);
    assert(ndims == 1);
    assert(dim == 1);
    ierr = DBReadVar(metafile.get(), "dtime", &dtime);
    assert(!ierr);
  }
  MPI_Bcast(&dtime, 1, MPI_DOUBLE, metafile_ioproc, mpi_comm);
  cctkGH->cctk_time = dtime;

  // Read internal driver state
  const std::string dirname = DB::legalize_name(driver_name);

  int npatches;
  std::vector<int> nlevels;
  if (read_metafile) {
    // Read number of levels
    const std::string varname = dirname + "/" + DB::legalize_name("nlevels");
    const int vartype = DBGetVarType(metafile.get(), varname.c_str());
    assert(vartype == DB_INT);
    const int varlength = DBGetVarLength(metafile.get(), varname.c_str());
    assert(varlength >= 0);
    npatches = varlength;
    nlevels.resize(varlength);
    ierr = DBReadVar(metafile.get(), varname.c_str(), nlevels.data());
    assert(!ierr);
  }
  MPI_Bcast(&npatches, 1, MPI_INT, metafile_ioproc, mpi_comm);
  CCTK_VINFO("Found %d patches", npatches);
  if (npatches != ghext->num_patches())
    CCTK_VERROR(
        "Wrong number of patches: Expected %d, found %d in the checkpoint file",
        ghext->num_patches(), npatches);
  if (!read_metafile)
    nlevels.resize(npatches);
  MPI_Bcast(nlevels.data(), npatches, MPI_INT, metafile_ioproc, mpi_comm);

  ghext->recovered_level_iterations.resize(ghext->num_patches());

  // Read FabArrayBase (component positions and shapes)
  for (int patch = 0; patch < npatches; ++patch) {
    CCTK_VINFO("Reading patch %d...", patch);

    CCTK_VINFO("  Found %d levels on patch %d", nlevels.at(patch), patch);
    auto &patchdata = ghext->patchdata.at(patch);
    patchdata.amrcore->SetFinestLevel(nlevels.at(patch) - 1);
    ghext->recovered_level_iterations.at(patch).resize(nlevels.at(patch));

    for (int level = 0; level < nlevels.at(patch); ++level) {
      CCTK_VINFO("  Reading level %d...", level);

      const std::string varname =
          dirname + "/" + make_fabarraybasename(patch, level);

      int nfabs;
      if (read_metafile) {
        const int vartype = DBGetVarType(metafile.get(), varname.c_str());
        assert(vartype == DB_INT);
        int vardims[2];
        const int varndims =
            DBGetVarDims(metafile.get(), varname.c_str(), 2, vardims);
        assert(varndims >= 0);
        assert(varndims == 2);
        assert(vardims[1] == 2 * ndims);
        nfabs = vardims[0];
        assert(nfabs >= 0);
      }
      MPI_Bcast(&nfabs, 1, MPI_INT, metafile_ioproc, mpi_comm);

      std::vector<int> data(2 * ndims * nfabs);
      if (read_metafile) {
        ierr = DBReadVar(metafile.get(), varname.c_str(), data.data());
        assert(!ierr);
      }
      MPI_Bcast(data.data(), data.size(), MPI_INT, metafile_ioproc, mpi_comm);

      amrex::Vector<amrex::Box> levboxes(nfabs);
      for (int component = 0; component < nfabs; ++component) {
        const amrex::IntVect small(&data.at(2 * ndims * component));
        const amrex::IntVect big(&data.at(ndims + 2 * ndims * component));
        levboxes.at(component) = amrex::Box(small, big);
      }

      // Don't set coarse level domain; this is already set by the driver
      if (level > 0) {
        amrex::Geometry geom = patchdata.amrcore->Geom(level - 1);
        geom.refine({2, 2, 2});
        patchdata.amrcore->SetGeometry(level, geom);
      }

      amrex::BoxList boxlist(std::move(levboxes));
      amrex::BoxArray boxarray(std::move(boxlist));
      patchdata.amrcore->SetBoxArray(level, boxarray);

      amrex::DistributionMapping dm(boxarray);
      patchdata.amrcore->SetDistributionMap(level, dm);

      patchdata.amrcore->SetupLevel(level, boxarray, dm,
                                    []() { return "Recovering"; });

      // Read per-level iteration if present (new checkpoint format)
      {
        const std::string varname_num =
            dirname + "/" +
            DB::legalize_name("iteration_num.m" + std::to_string(patch) +
                              ".rl" + std::to_string(level));
        const std::string varname_den =
            dirname + "/" +
            DB::legalize_name("iteration_den.m" + std::to_string(patch) +
                              ".rl" + std::to_string(level));

        long long iter_num = 0, iter_den = 0;
        bool have_iteration = false;
        if (read_metafile) {
          if (DBInqVarExists(metafile.get(), varname_num.c_str()) &&
              DBInqVarExists(metafile.get(), varname_den.c_str())) {
            int ierr;
            ierr = DBReadVar(metafile.get(), varname_num.c_str(), &iter_num);
            assert(!ierr);
            ierr = DBReadVar(metafile.get(), varname_den.c_str(), &iter_den);
            assert(!ierr);
            have_iteration = true;
          }
        }
        MPI_Bcast(&have_iteration, 1, MPI_C_BOOL, metafile_ioproc, mpi_comm);
        if (have_iteration) {
          MPI_Bcast(&iter_num, 1, MPI_LONG_LONG, metafile_ioproc, mpi_comm);
          MPI_Bcast(&iter_den, 1, MPI_LONG_LONG, metafile_ioproc, mpi_comm);
          ghext->recovered_level_iterations.at(patch).at(level) =
              rat64(iter_num, iter_den);
        }
        // else: old checkpoint without per-level iteration — leave as nullopt
      }
    } // for level
  } // for patch

  interval_meta = nullptr;
}

void InputSilo(const cGH *restrict const cctkGH,
               const std::vector<bool> &input_group,
               const std::string &input_dir, const std::string &input_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  assert(!input_dir.empty());
  assert(!input_file.empty());

  // Set up timers
  static Timer timer("InputSilo");
  Interval interval(timer);

  if (std::count(input_group.begin(), input_group.end(), true) == 0)
    return;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root) {
    CCTK_VINFO("Silo input for groups:");
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
      if (input_group.at(gi))
        CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
  }

  static Timer timer_setup("InputSilo.setup");
  auto interval_setup = std::make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  const int ioproc_every = [&]() {
    if (CCTK_EQUALS(out_mode, "proc"))
      return 1;
    if (CCTK_EQUALS(out_mode, "np"))
      return out_proc_every;
    if (CCTK_EQUALS(out_mode, "onefile"))
      return nprocs;
    assert(0);
  }();
  assert(ioproc_every > 0);

  const int myioproc = myproc / ioproc_every * ioproc_every;
  assert(myioproc <= myproc);
  assert(myproc < myioproc + ioproc_every);
  const bool read_file = myproc % ioproc_every == 0;
  assert((myioproc == myproc) == read_file);

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_data("InputSilo.data");
  auto interval_data = std::make_unique<Interval>(timer_data);

  // Read data
  {
    DB::ptr<DBfile> file;
    if (read_file) {
      const std::string subdirname =
          make_subdirname(input_file, cctk_iteration);
      const std::string filename =
          input_dir + "/" + subdirname + "/" +
          make_filename(input_file, cctk_iteration, myproc / ioproc_every);
      // We could use DB_UNKNOWN instead of DB_HDF5
      file = DB::make(DBOpen(filename.c_str(), DB_HDF5, DB_READ));
      assert(file);
    }

    // Bands are written all-or-nothing (only at unsynchronized checkpoints), so
    // detect their presence once with a global flag. This keeps the band read
    // symmetric with the regular read; a per-variable file probe inside the
    // loop would desync the band owner's MPI_Recv.
    bool file_has_bands = false;
    if (ghext->use_subcycling) {
      int local_has_bands = 0;
      if (read_file) {
        for (const auto &patchdata : ghext->patchdata) {
          for (const auto &leveldata : patchdata.leveldata) {
            for (int gi = 0; gi < CCTK_NumGroups() && !local_has_bands; ++gi) {
              if (!input_group.at(gi) || CCTK_GroupTypeI(gi) != CCTK_GF)
                continue;
              const auto &groupdata = *leveldata.groupdata.at(gi);
              if (groupdata.mfab.empty())
                continue;
              const auto probe = [&](const amrex::MultiFab *const band,
                                     const band_kind kind, const int stage) {
                if (local_has_bands || !band || band->empty())
                  return;
                const std::string band_tag = subcycling_band_tag(kind, stage);
                const amrex::DistributionMapping &bdm = band->DistributionMap();
                for (int c = 0; c < bdm.size(); ++c) {
                  if (bdm[c] / ioproc_every * ioproc_every != myproc)
                    continue; // not read by this ioproc
                  const std::string varname = make_varname(
                      gi, 0, patchdata.patch, leveldata.level, c, band_tag);
                  if (DBInqVarExists(file.get(), varname.c_str())) {
                    local_has_bands = 1;
                    return;
                  }
                }
              };
              for (int s = 0; s < max_num_rk_stages; ++s)
                probe(groupdata.ks_consumer_band[s].get(),
                      band_kind::ks_consumer, s);
              probe(groupdata.old_consumer_band.get(), band_kind::old_consumer,
                    -1);
            }
          }
        }
      }
      int global_has_bands = 0;
      MPI_Allreduce(&local_has_bands, &global_has_bands, 1, MPI_INT, MPI_MAX,
                    mpi_comm);
      file_has_bands = global_has_bands != 0;
    }

    // Loop over patches and levels
    for (const auto &patchdata : ghext->patchdata) {
      for (const auto &leveldata : patchdata.leveldata) {
        if (io_verbose)
          CCTK_VINFO("Reading patch %d level %d", patchdata.patch,
                     leveldata.level);

        // Loop over groups
        // std::set<mesh_props_t> have_meshes;
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          if (!input_group.at(gi))
            continue;
          // TODO: read grid arrays
          if (CCTK_GroupTypeI(gi) != CCTK_GF)
            continue;
          if (io_verbose)
            CCTK_VINFO("  Reading group %s", CCTK_FullGroupName(gi));

          auto &groupdata = *leveldata.groupdata.at(gi);
          if (groupdata.mfab.empty())
            continue;
          // TODO: Check that group has the default number of ghost zones
          const int numvars = groupdata.numvars;
          const int tl = 0;
          amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::IndexType &indextype = mfab.ixType();
          const amrex::DistributionMapping &dm = mfab.DistributionMap();

          const std::array<int, dim> &nghosts = groupdata.nghostzones;

          // Loop over components (AMReX boxes)
          const int ncomponents = dm.size();
          for (int component = 0; component < ncomponents; ++component) {
            if (io_verbose)
              CCTK_VINFO("    Reading component %d", component);

            const int proc = dm[component];
            const int ioproc = proc / ioproc_every * ioproc_every;
            const bool recv_this_fab = proc == myproc;
            const bool read_this_fab = ioproc == myproc;
            if (!(recv_this_fab || read_this_fab))
              continue;

            const amrex::Box &fabbox = mfab.fabbox(component); // exterior

            std::array<int, ndims> dims;
            for (int d = 0; d < ndims; ++d)
              dims[d] = fabbox.length(d);
            std::ptrdiff_t zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              zonecount *= dims[d];
            assert(zonecount >= 0 && zonecount <= INT_MAX);

            // Communicate variable, part 1
            static Timer timer_mpi("InputSilo.mpi");
            auto interval_mpi = std::make_unique<Interval>(timer_mpi);
            const int mpi_tag = 22901; // randomly chosen
            std::vector<CCTK_REAL> buffer;
            MPI_Request mpi_req;
            CCTK_REAL *data = nullptr;
            if (recv_this_fab && read_this_fab) {
              amrex::FArrayBox &fab = mfab[component];
              data = fab.dataPtr();
            } else if (recv_this_fab) {
              amrex::FArrayBox &fab = mfab[component];
              assert(numvars * zonecount <= INT_MAX);
              MPI_Irecv(fab.dataPtr(), numvars * zonecount,
                        mpi_datatype_v<CCTK_REAL>, ioproc, mpi_tag, mpi_comm,
                        &mpi_req);
            } else {
              buffer.resize(numvars * zonecount);
              assert(numvars * zonecount <= INT_MAX);
              data = buffer.data();
            }
            interval_mpi = nullptr;

            // Read variable
            if (read_file) {
              static Timer timer_var("InputSilo.var");
              Interval interval_var(timer_var);

              const std::string meshname = make_meshname(
                  nghosts, patchdata.patch, leveldata.level, component);

              const int centering = [&]() {
                const int rank = indextype.cellCentered(0) +
                                 indextype.cellCentered(1) +
                                 indextype.cellCentered(2);
                switch (rank) {
                case 0:
                  return DB_NODECENT;
                case 1:
                  return DB_EDGECENT;
                case 2:
                  return DB_FACECENT;
                case 3:
                  return DB_ZONECENT;
                }
                assert(0);
              }();

              if (centering == DB_EDGECENT || centering == DB_FACECENT) {
                // Need to find the other 2 edge- or face-centered
                // variables, and output them as well. Maybe input all 3
                // when the x- or xy-centered value is input, for those
                // which should be input? Maybe add a "sibling" tag to
                // such grid functions to find these other components?
                assert(0);
              }

              for (int vi = 0; vi < numvars; ++vi) {
                const std::string varname = make_varname(
                    gi, vi, patchdata.patch, leveldata.level, component);
                if (io_verbose)
                  CCTK_VINFO("      Reading variable %s", varname.c_str());

                const DB::ptr<DBquadvar> quadvar =
                    DB::make(DBGetQuadvar(file.get(), varname.c_str()));
                assert(quadvar);

                assert(quadvar->ndims == ndims);
                assert(ndims <= 3);
                for (int d = 0; d < ndims; ++d)
                  assert(quadvar->dims[d] == dims[d]);
                assert(quadvar->datatype == db_datatype_v<CCTK_REAL>);
                assert(quadvar->centering == centering);
                assert(quadvar->nvals == 1);

                // TODO: check DBOPT_COORDSYS:
                // int cartesian = DB_CARTESIAN;

                const int column_major = 0;
                assert(quadvar->major_order == column_major);

                const void *const read_ptr = quadvar->vals[0];

                void *const data_ptr = data + vi * zonecount;
                memcpy(data_ptr, read_ptr, zonecount * sizeof(CCTK_REAL));
              } // for vi
            } // if read_file

            // Communicate variable, part 2
            static Timer timer_wait("InputSilo.wait");
            auto interval_wait = std::make_unique<Interval>(timer_wait);
            if (recv_this_fab && read_this_fab) {
              // do nothing
            } else if (recv_this_fab) {
              MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
            } else {
              assert(std::ptrdiff_t(buffer.size()) == numvars * zonecount);
              MPI_Send(buffer.data(), numvars * zonecount,
                       mpi_datatype_v<CCTK_REAL>, proc, mpi_tag, mpi_comm);
            }
            if (recv_this_fab)
              for (int vi = 0; vi < numvars; ++vi)
                groupdata.valid.at(tl).at(vi).set_all(
                    make_valid_all(), []() { return "read from Silo file"; });
            interval_wait = nullptr;

          } // for component

          // Subcycling consumer bands: zero-ghost MultiFabs with their own
          // DistributionMap, each read in its own component loop mirroring the
          // Phase 1 write (reversed MPI). When file_has_bands, every rebuilt
          // non-empty band is present in full.
          if (file_has_bands) {
            const int centering = [&]() {
              const int rank = indextype.cellCentered(0) +
                               indextype.cellCentered(1) +
                               indextype.cellCentered(2);
              switch (rank) {
              case 0:
                return DB_NODECENT;
              case 1:
                return DB_EDGECENT;
              case 2:
                return DB_FACECENT;
              case 3:
                return DB_ZONECENT;
              }
              assert(0);
            }();
            assert(centering != DB_EDGECENT && centering != DB_FACECENT);

            const auto read_band = [&](amrex::MultiFab *const band,
                                       const band_kind kind, const int stage) {
              if (!band || band->empty())
                return;
              assert(band->nGrowVect() == 0);
              const std::string band_tag = subcycling_band_tag(kind, stage);
              const amrex::DistributionMapping &bdm = band->DistributionMap();
              const int bncomponents = bdm.size();

              for (int component = 0; component < bncomponents; ++component) {
                const int proc = bdm[component];
                const int ioproc = proc / ioproc_every * ioproc_every;
                const bool recv_this_fab = proc == myproc;
                const bool read_this_fab = ioproc == myproc;
                if (!(recv_this_fab || read_this_fab))
                  continue;

                const amrex::Box &fabbox =
                    band->fabbox(component); // zero ghost

                std::array<int, ndims> dims;
                for (int d = 0; d < ndims; ++d)
                  dims[d] = fabbox.length(d);
                std::ptrdiff_t zonecount = 1;
                for (int d = 0; d < ndims; ++d)
                  zonecount *= dims[d];
                assert(zonecount >= 0 && zonecount <= INT_MAX);

                // Communicate the band fab from the I/O process, part 1.
                const int mpi_tag = 22903; // randomly chosen, band path
                std::vector<CCTK_REAL> buffer;
                MPI_Request mpi_req;
                CCTK_REAL *data = nullptr;
                if (recv_this_fab && read_this_fab) {
                  amrex::FArrayBox &fab = (*band)[component];
                  data = fab.dataPtr();
                } else if (recv_this_fab) {
                  amrex::FArrayBox &fab = (*band)[component];
                  assert(numvars * zonecount <= INT_MAX);
                  MPI_Irecv(fab.dataPtr(), numvars * zonecount,
                            mpi_datatype_v<CCTK_REAL>, ioproc, mpi_tag,
                            mpi_comm, &mpi_req);
                } else {
                  buffer.resize(numvars * zonecount);
                  assert(numvars * zonecount <= INT_MAX);
                  data = buffer.data();
                }

                // Read the band variables.
                if (read_file) {
                  for (int vi = 0; vi < numvars; ++vi) {
                    const std::string varname =
                        make_varname(gi, vi, patchdata.patch, leveldata.level,
                                     component, band_tag);
                    const DB::ptr<DBquadvar> quadvar =
                        DB::make(DBGetQuadvar(file.get(), varname.c_str()));
                    assert(quadvar);
                    assert(quadvar->ndims == ndims);
                    for (int d = 0; d < ndims; ++d)
                      assert(quadvar->dims[d] == dims[d]);
                    assert(quadvar->datatype == db_datatype_v<CCTK_REAL>);
                    assert(quadvar->centering == centering);
                    assert(quadvar->nvals == 1);
                    const int column_major = 0;
                    assert(quadvar->major_order == column_major);
                    const void *const read_ptr = quadvar->vals[0];
                    void *const data_ptr = data + vi * zonecount;
                    memcpy(data_ptr, read_ptr, zonecount * sizeof(CCTK_REAL));
                  } // for vi
                } // if read_file

                // Communicate the band fab, part 2.
                if (recv_this_fab && read_this_fab) {
                  // do nothing
                } else if (recv_this_fab) {
                  MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
                } else {
                  assert(std::ptrdiff_t(buffer.size()) == numvars * zonecount);
                  MPI_Send(buffer.data(), numvars * zonecount,
                           mpi_datatype_v<CCTK_REAL>, proc, mpi_tag, mpi_comm);
                }
              } // for band component
            };

            for (int s = 0; s < max_num_rk_stages; ++s)
              read_band(groupdata.ks_consumer_band[s].get(),
                        band_kind::ks_consumer, s);
            read_band(groupdata.old_consumer_band.get(),
                      band_kind::old_consumer, -1);
          } // if file_has_bands

        } // for gi
      } // for leveldata
    } // for patchdata
  } // read data

  interval_data = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void OutputSilo(const cGH *restrict const cctkGH,
                const std::vector<bool> &output_group,
                const std::string &output_dir, const std::string &output_file,
                std::optional<slice_t> slice) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("OutputSilo");
  Interval interval(timer);

  if (std::count(output_group.begin(), output_group.end(), true) == 0)
    return;

  if (io_verbose)
    CCTK_VINFO("OutputSilo...");

  static Timer timer_setup("OutputSilo.setup");
  auto interval_setup = std::make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);
  const int nprocs = CCTK_nProcs(cctkGH);

  const int ioproc_every = [&]() {
    if (CCTK_EQUALS(out_mode, "proc"))
      return 1;
    if (CCTK_EQUALS(out_mode, "np"))
      return out_proc_every;
    if (CCTK_EQUALS(out_mode, "onefile"))
      return nprocs;
    assert(0);
  }();
  assert(ioproc_every > 0);

  const int myioproc = myproc / ioproc_every * ioproc_every;
  assert(myioproc <= myproc);
  assert(myproc < myioproc + ioproc_every);
  const bool write_file = myproc % ioproc_every == 0;
  assert((myioproc == myproc) == write_file);
  // If process 1 exists and if it is not an I/O process, then output
  // metadata there. Else output metadata on process 0.
  const int metafile_ioproc = nprocs == 1 || ioproc_every == 1 ? 0 : 1;
  const bool write_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);
  // DBSetAllowEmptyObjects(1);
  DBSetCompression(out_silo_compression_options);
  DBSetEnableChecksums(1);

  // TODO: directories instead of carefully chosen names

  constexpr int ndims = dim;

  // At an unsynchronized subcycling checkpoint the fine consumer bands hold
  // mid-cycle state that exists nowhere else and must be serialized. Otherwise
  // the on-disk format is unchanged.
  const bool write_bands = !all_levels_synchronized();

  // 2D slice output maps a physical coord to an index via (coord - x0)/dx,
  // which is undefined on a non-Cartesian patch. Reject up front so every loop
  // below can assume slice ⟹ all-Cartesian. is_cartesian is identical on every
  // rank, so this scan is collective-safe.
  if (slice) {
    for (const auto &patchdata : ghext->patchdata)
      if (!patchdata.is_cartesian)
        CCTK_VERROR("2D Silo slice output is only supported for "
                    "Cartesian patches");
  }

  interval_setup = nullptr;

  static Timer timer_data("OutputSilo.data");
  auto interval_data = std::make_unique<Interval>(timer_data);

  // Global coordinate extents for each component
  std::vector<std::vector<std::vector<std::array<CCTK_REAL, dim> > > >
      coordinate_minima, coordinate_maxima;

  // Write data
  {
    static Timer timer_mkdir("OutputSilo.mkdir");
    auto interval_mkdir = std::make_unique<Interval>(timer_mkdir);
    const std::string subdirname = make_subdirname(output_file, cctk_iteration);
    const std::string pathname = std::string(output_dir) + "/" + subdirname;
    const int mode = 0755;
    static once_flag create_directory;
    call_once(create_directory, [&]() {
      const int ierr = CCTK_CreateDirectory(mode, output_dir.c_str());
      assert(ierr >= 0);
    });
    ierr = CCTK_CreateDirectory(mode, pathname.c_str());
    assert(ierr >= 0);
    interval_mkdir = nullptr;

    DB::ptr<DBfile> file;
    if (write_file) {
      const std::string filename =
          pathname + "/" +
          make_filename(output_file, cctk_iteration, myproc / ioproc_every);
      // assert(filename.size() < 256);
      // assert(output_file.size() < 256);
      file = DB::make(DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL,
                               output_file.c_str(), DB_HDF5));
      assert(file);
    }

    if (write_file) {
      {
        char *const data = IOUtil_GetAllParameters(cctkGH, 1 /*all*/);
        const std::string parameters(data);
        std::free(data);
        const int dims = parameters.length();
        ierr = DBWrite(file.get(), "AllParameters", parameters.data(), &dims, 1,
                       DB_CHAR);
        assert(!ierr);
      }

      { // Tell VisIt that the mesh structure may change over time
        const int dims = 1;
        const int value = 1;
        ierr = DBWrite(file.get(), "MetadataIsTimeVarying", &value, &dims, 1,
                       DB_INT);
        assert(!ierr);
      }
    }

    coordinate_minima.resize(ghext->patchdata.size());
    coordinate_maxima.resize(ghext->patchdata.size());

    // Loop over patches
    for (const auto &patchdata : ghext->patchdata) {

      if (!patchdata.is_cartesian) {
        coordinate_minima.at(patchdata.patch)
            .resize(patchdata.leveldata.size());
        coordinate_maxima.at(patchdata.patch)
            .resize(patchdata.leveldata.size());
      }

      // Loop over levels
      for (const auto &leveldata : patchdata.leveldata) {

        const amrex::DistributionMapping &dm = leveldata.fab->DistributionMap();
        const int ncomponents = dm.size();

        if (!patchdata.is_cartesian) {
          coordinate_minima.at(patchdata.patch)
              .at(leveldata.level)
              .resize(ncomponents);
          coordinate_maxima.at(patchdata.patch)
              .at(leveldata.level)
              .resize(ncomponents);
        }

        // Loop over groups
        std::set<mesh_props_t> have_meshes;
        for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
          if (!output_group.at(gi))
            continue;
          // TODO: Output grid arrays
          if (CCTK_GroupTypeI(gi) != CCTK_GF)
            continue;

          const auto &groupdata = *leveldata.groupdata.at(gi);
          if (groupdata.mfab.empty())
            continue;
          // TODO: Check that group has the default number of ghost zones
          const int numvars = groupdata.numvars;
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::IndexType &indextype = mfab.ixType();
          const amrex::DistributionMapping &dm = mfab.DistributionMap();

          const std::array<int, dim> &nghosts = groupdata.nghostzones;
          const mesh_props_t mesh_props{nghosts};
          const bool have_mesh = have_meshes.count(mesh_props);

          // Loop over components (AMReX boxes)
          for (int component = 0; component < ncomponents; ++component) {
            const int proc = dm[component];
            const int ioproc = proc / ioproc_every * ioproc_every;
            const bool send_this_fab = proc == myproc;
            const bool write_this_fab = ioproc == myproc;
            if (!(send_this_fab || write_this_fab))
              continue;

            // TODO: Check whether data are valid
            const amrex::Box &fabbox = mfab.fabbox(component); // exterior

            // Exterior (source) half-open bounds. The variable `data` array is
            // always shaped by this box, whether the local fab pointer or an
            // MPI receive buffer.
            const Arith::vect<int, 3> ext_lo{
                fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)};
            const Arith::vect<int, 3> ext_hi{fabbox.bigEnd(0) + 1,
                                             fabbox.bigEnd(1) + 1,
                                             fabbox.bigEnd(2) + 1};

            // Number of zones in the full exterior box. Used for MPI send/recv
            // sizing and for indexing into the (un-sliced) `data` array.
            std::ptrdiff_t src_zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              src_zonecount *= ext_hi[d] - ext_lo[d];
            assert(src_zonecount >= 0 && src_zonecount <= INT_MAX);

            // Geometry for this (patch, level): index 0 sits at the Cartesian
            // coordinate x0[d], with cell size dx[d].
            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const amrex::Real *const x0 = geom.ProbLo();
            const amrex::Real *const dx = geom.CellSize();

            // Restrict to a thickness-1 slab when slicing. sub_lo/sub_hi
            // default to the full exterior box, so the non-slice path is
            // unchanged.
            Arith::vect<int, 3> sub_lo = ext_lo, sub_hi = ext_hi;
            const int slice_n = slice ? slice->normal_dir : -1;
            if (slice) {
              const auto r =
                  slice_slab(mfab, component, x0, dx, indextype, *slice);
              if (!r)
                continue; // component does not intersect the plane
              sub_lo = r->first;
              sub_hi = r->second;
            }

            // Mesh dimensions: the slab shape (== exterior shape when not
            // slicing).
            std::array<int, ndims> dims;
            for (int d = 0; d < ndims; ++d)
              dims[d] = sub_hi[d] - sub_lo[d];
            std::ptrdiff_t zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              zonecount *= dims[d];
            assert(zonecount >= 0 && zonecount <= INT_MAX);

            if (write_file && !have_mesh) {
              static Timer timer_mesh("OutputSilo.mesh");
              Interval interval_mesh(timer_mesh);

              const std::string meshname = make_meshname(
                  nghosts, patchdata.patch, leveldata.level, component);

              std::array<int, ndims> dims_vc;
              for (int d = 0; d < ndims; ++d)
                dims_vc[d] = dims[d] + int(indextype.cellCentered(d));

              std::array<std::vector<CCTK_REAL>, ndims> coords;
              std::array<const void *, ndims> coord_ptrs;
              if (patchdata.is_cartesian) {
                for (int d = 0; d < ndims; ++d) {
                  coords[d].resize(dims_vc[d]);
                  // sub_lo[d] == fabbox.smallEnd(d) for in-plane axes; for the
                  // sliced normal axis it is the plane index, so the slab's
                  // single coordinate plane is placed correctly.
                  for (int i = 0; i < dims_vc[d]; ++i)
                    coords[d][i] = x0[d] + (sub_lo[d] + i) * dx[d];
                }
                for (int d = 0; d < ndims; ++d)
                  coord_ptrs[d] = coords[d].data();
              } else {
                // Find the vertex-centred coordinate grid functions
                // Silo always wants vertex-centred coordinates.
                const int coordinate_group =
                    CCTK_GroupIndex("CoordinatesX::vertex_coords");
                assert(coordinate_group >= 0);
                // TODO: We may need to send the coordinates across
                // the network. We don't do this yet, so assert that
                // we actually have the coordinates stored locally.
                assert(send_this_fab && write_this_fab);
                assert(leveldata.groupdata.at(coordinate_group)->nghostzones ==
                       groupdata.nghostzones);
                const auto &coord_mfab =
                    *leveldata.groupdata.at(coordinate_group)->mfab.at(0);
                const auto &coord_fabbox = coord_mfab[component];
                for (int d = 0; d < ndims; ++d) {
                  coord_ptrs[d] = coord_fabbox.dataPtr(d);
                  assert(coord_ptrs[d]);
                }

                std::array<CCTK_REAL, dim> xmin, xmax;
                for (int d = 0; d < dim; ++d) {
#ifdef AMREX_USE_GPU
                  constexpr auto run_on = amrex::RunOn::Device;
#else
                  constexpr auto run_on = amrex::RunOn::Host;
#endif
                  const auto xminmax = coord_fabbox.minmax<run_on>(d);
                  xmin[d] = xminmax.first;
                  xmax[d] = xminmax.second;
                }
                coordinate_minima.at(patchdata.patch)
                    .at(leveldata.level)
                    .at(component) = xmin;
                coordinate_maxima.at(patchdata.patch)
                    .at(leveldata.level)
                    .at(component) = xmax;
              }

              const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
              assert(optlist);

              if (patchdata.is_cartesian) {
                int cartesian = DB_CARTESIAN;
                ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
                assert(!ierr);
              }

              int cycle = cctk_iteration;
              ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
              assert(!ierr);

              std::array<int, ndims> min_index, max_index;
              for (int d = 0; d < ndims; ++d) {
                // The thickness-1 slab carries no ghosts along the normal.
                const int ng = (slice && d == slice_n) ? 0 : nghosts[d];
                min_index[d] = ng;
                max_index[d] = ng;
              }
              ierr =
                  DBAddOption(optlist.get(), DBOPT_LO_OFFSET, min_index.data());
              assert(!ierr);
              ierr =
                  DBAddOption(optlist.get(), DBOPT_HI_OFFSET, max_index.data());
              assert(!ierr);

              int column_major = 0;
              ierr =
                  DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
              assert(!ierr);

              // float time = cctk_time;
              // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
              // assert(!ierr);
              double dtime = cctk_time;
              ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
              assert(!ierr);

              int hide_from_gui = 1;
              ierr = DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
                                 &hide_from_gui);
              assert(!ierr);

              const int coordtype =
                  patchdata.is_cartesian ? DB_COLLINEAR : DB_NONCOLLINEAR;

              ierr = DBPutQuadmesh(file.get(), meshname.c_str(), nullptr,
                                   coord_ptrs.data(), dims_vc.data(), ndims,
                                   db_datatype_v<CCTK_REAL>, coordtype,
                                   optlist.get());
              assert(!ierr);

              have_meshes.insert(mesh_props);
            } // if write mesh

            // Communicate variable
            static Timer timer_mpi("OutputSilo.mpi");
            auto interval_mpi = std::make_unique<Interval>(timer_mpi);
            const int mpi_tag = 22900; // randomly chosen
            std::vector<CCTK_REAL> buffer;
            const CCTK_REAL *data = nullptr;
            if (send_this_fab && write_this_fab) {
              const amrex::FArrayBox &fab = mfab[component];
              data = fab.dataPtr();
            } else if (send_this_fab) {
              const amrex::FArrayBox &fab = mfab[component];
              assert(numvars * src_zonecount <= INT_MAX);
              MPI_Send(fab.dataPtr(), numvars * src_zonecount,
                       mpi_datatype_v<CCTK_REAL>, ioproc, mpi_tag, mpi_comm);
            } else {
              buffer.resize(numvars * src_zonecount);
              assert(numvars * src_zonecount <= INT_MAX);
              MPI_Recv(buffer.data(), numvars * src_zonecount,
                       mpi_datatype_v<CCTK_REAL>, proc, mpi_tag, mpi_comm,
                       MPI_STATUS_IGNORE);
              data = buffer.data();
            }
            interval_mpi = nullptr;

            // Write variable
            if (write_file) {
              static Timer timer_var("OutputSilo.var");
              Interval interval_var(timer_var);

              const std::string meshname = make_meshname(
                  nghosts, patchdata.patch, leveldata.level, component);

              const int centering = [&]() {
                const int rank = indextype.cellCentered(0) +
                                 indextype.cellCentered(1) +
                                 indextype.cellCentered(2);
                switch (rank) {
                case 0:
                  return DB_NODECENT;
                case 1:
                  return DB_EDGECENT;
                case 2:
                  return DB_FACECENT;
                case 3:
                  return DB_ZONECENT;
                }
                assert(0);
              }();

              const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
              assert(optlist);

              int cartesian = DB_CARTESIAN;
              ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
              assert(!ierr);

              int cycle = cctk_iteration;
              ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
              assert(!ierr);

              int column_major = 0;
              ierr =
                  DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
              assert(!ierr);

              // float time = cctk_time;
              // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
              // assert(!ierr);
              double dtime = cctk_time;
              ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
              assert(!ierr);

              int hide_from_gui = 1;
              ierr = DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
                                 &hide_from_gui);
              assert(!ierr);

              if (centering == DB_EDGECENT || centering == DB_FACECENT) {
                // Need to find the other 2 edge- or face-centered
                // variables, and output them as well. Maybe output all
                // 3 when the x- or xy-centered value is output? Maybe
                // add a "sibling" tag to such grid functions to find
                // these other components?
                assert(0);
              }

              // A thickness-1 slab is contiguous only for the xy plane (normal
              // = z, the slowest axis); xz/yz planes are strided. Extract each
              // variable's slab from the exterior-shaped `data` into a
              // contiguous buffer before handing it to Silo. `zonecount` is the
              // slab count (== src_zonecount when not slicing).
              std::vector<CCTK_REAL> slab;
              if (slice)
                slab.resize(zonecount);
              for (int vi = 0; vi < numvars; ++vi) {
                const std::string varname = make_varname(
                    gi, vi, patchdata.patch, leveldata.level, component);

                const void *data_ptr;
                if (slice) {
                  extract_subbox(slab.data(), data + vi * src_zonecount, ext_lo,
                                 ext_hi, sub_lo, sub_hi);
                  data_ptr = slab.data();
                } else {
                  data_ptr = data + vi * src_zonecount;
                }

                ierr = DBPutQuadvar1(
                    file.get(), varname.c_str(), meshname.c_str(), data_ptr,
                    dims.data(), ndims, nullptr, 0, db_datatype_v<CCTK_REAL>,
                    centering, optlist.get());
              } // for vi
            } // if write_file

          } // for component

          // Subcycling consumer bands: zero-ghost MultiFabs in this level's
          // index space with their own DistributionMap, so each gets its own
          // component loop. Geometry is rebuilt on recovery; we serialize only
          // the data, namespaced by a band tag. 2D slice output covers only the
          // main grid data, so bands are skipped entirely when slicing.
          if (write_bands && !slice) {
            const int centering = [&]() {
              const int rank = indextype.cellCentered(0) +
                               indextype.cellCentered(1) +
                               indextype.cellCentered(2);
              switch (rank) {
              case 0:
                return DB_NODECENT;
              case 1:
                return DB_EDGECENT;
              case 2:
                return DB_FACECENT;
              case 3:
                return DB_ZONECENT;
              }
              assert(0);
            }();
            assert(centering != DB_EDGECENT && centering != DB_FACECENT);

            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const amrex::Real *const x0 = geom.ProbLo();
            const amrex::Real *const dx = geom.CellSize();
            const std::array<int, dim> band_nghosts{0, 0, 0};

            const auto write_band = [&](const amrex::MultiFab *const band,
                                        const band_kind kind, const int stage) {
              if (!band || band->empty())
                return;
              assert(band->nGrowVect() == 0);
              const std::string band_tag = subcycling_band_tag(kind, stage);
              const amrex::DistributionMapping &bdm = band->DistributionMap();
              const int bncomponents = bdm.size();

              for (int component = 0; component < bncomponents; ++component) {
                const int proc = bdm[component];
                const int ioproc = proc / ioproc_every * ioproc_every;
                const bool send_this_fab = proc == myproc;
                const bool write_this_fab = ioproc == myproc;
                if (!(send_this_fab || write_this_fab))
                  continue;

                const amrex::Box &fabbox =
                    band->fabbox(component); // zero ghost

                std::array<int, ndims> dims;
                for (int d = 0; d < ndims; ++d)
                  dims[d] = fabbox.length(d);
                std::ptrdiff_t zonecount = 1;
                for (int d = 0; d < ndims; ++d)
                  zonecount *= dims[d];
                assert(zonecount >= 0 && zonecount <= INT_MAX);

                // Write the band mesh (zero ghost => no LO/HI offset).
                if (write_file) {
                  const std::string meshname =
                      make_meshname(band_nghosts, patchdata.patch,
                                    leveldata.level, component, band_tag);

                  std::array<int, ndims> dims_vc;
                  for (int d = 0; d < ndims; ++d)
                    dims_vc[d] = dims[d] + int(indextype.cellCentered(d));

                  std::array<std::vector<CCTK_REAL>, ndims> coords;
                  std::array<const void *, ndims> coord_ptrs;
                  for (int d = 0; d < ndims; ++d) {
                    coords[d].resize(dims_vc[d]);
                    for (int i = 0; i < dims_vc[d]; ++i)
                      coords[d][i] = x0[d] + (fabbox.smallEnd(d) + i) * dx[d];
                    coord_ptrs[d] = coords[d].data();
                  }

                  const DB::ptr<DBoptlist> optlist =
                      DB::make(DBMakeOptlist(10));
                  assert(optlist);
                  int cartesian = DB_CARTESIAN;
                  ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
                  assert(!ierr);
                  int cycle = cctk_iteration;
                  ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
                  assert(!ierr);
                  std::array<int, ndims> zero_offset;
                  for (int d = 0; d < ndims; ++d)
                    zero_offset[d] = 0;
                  ierr = DBAddOption(optlist.get(), DBOPT_LO_OFFSET,
                                     zero_offset.data());
                  assert(!ierr);
                  ierr = DBAddOption(optlist.get(), DBOPT_HI_OFFSET,
                                     zero_offset.data());
                  assert(!ierr);
                  int column_major = 0;
                  ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER,
                                     &column_major);
                  assert(!ierr);
                  double dtime = cctk_time;
                  ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
                  assert(!ierr);
                  int hide_from_gui = 1;
                  ierr = DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
                                     &hide_from_gui);
                  assert(!ierr);

                  ierr = DBPutQuadmesh(file.get(), meshname.c_str(), nullptr,
                                       coord_ptrs.data(), dims_vc.data(), ndims,
                                       db_datatype_v<CCTK_REAL>, DB_COLLINEAR,
                                       optlist.get());
                  assert(!ierr);
                }

                // Communicate the band fab to the I/O process.
                const int mpi_tag = 22902; // randomly chosen, band path
                std::vector<CCTK_REAL> buffer;
                const CCTK_REAL *data = nullptr;
                if (send_this_fab && write_this_fab) {
                  const amrex::FArrayBox &fab = (*band)[component];
                  data = fab.dataPtr();
                } else if (send_this_fab) {
                  const amrex::FArrayBox &fab = (*band)[component];
                  assert(numvars * zonecount <= INT_MAX);
                  MPI_Send(fab.dataPtr(), numvars * zonecount,
                           mpi_datatype_v<CCTK_REAL>, ioproc, mpi_tag,
                           mpi_comm);
                } else {
                  buffer.resize(numvars * zonecount);
                  assert(numvars * zonecount <= INT_MAX);
                  MPI_Recv(buffer.data(), numvars * zonecount,
                           mpi_datatype_v<CCTK_REAL>, proc, mpi_tag, mpi_comm,
                           MPI_STATUS_IGNORE);
                  data = buffer.data();
                }

                // Write the band variables.
                if (write_file) {
                  const std::string meshname =
                      make_meshname(band_nghosts, patchdata.patch,
                                    leveldata.level, component, band_tag);

                  const DB::ptr<DBoptlist> optlist =
                      DB::make(DBMakeOptlist(10));
                  assert(optlist);
                  int cartesian = DB_CARTESIAN;
                  ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
                  assert(!ierr);
                  int cycle = cctk_iteration;
                  ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
                  assert(!ierr);
                  int column_major = 0;
                  ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER,
                                     &column_major);
                  assert(!ierr);
                  double dtime = cctk_time;
                  ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
                  assert(!ierr);
                  int hide_from_gui = 1;
                  ierr = DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI,
                                     &hide_from_gui);
                  assert(!ierr);

                  for (int vi = 0; vi < numvars; ++vi) {
                    const std::string varname =
                        make_varname(gi, vi, patchdata.patch, leveldata.level,
                                     component, band_tag);
                    const void *const data_ptr = data + vi * zonecount;
                    ierr = DBPutQuadvar1(
                        file.get(), varname.c_str(), meshname.c_str(), data_ptr,
                        dims.data(), ndims, nullptr, 0,
                        db_datatype_v<CCTK_REAL>, centering, optlist.get());
                    assert(!ierr);
                  } // for vi
                } // if write_file
              } // for band component
            };

            for (int s = 0; s < max_num_rk_stages; ++s)
              write_band(groupdata.ks_consumer_band[s].get(),
                         band_kind::ks_consumer, s);
            write_band(groupdata.old_consumer_band.get(),
                       band_kind::old_consumer, -1);
          } // if write_bands

        } // for gi
      } // for leveldata
    } // for patchdata
  } // write data

  for (const auto &patchdata : ghext->patchdata) {
    if (!patchdata.is_cartesian) {
      for (const auto &leveldata : patchdata.leveldata) {
        const int patch = patchdata.patch;
        const int level = leveldata.level;
        auto &minima = coordinate_minima.at(patch).at(level);
        auto &maxima = coordinate_maxima.at(patch).at(level);
        if (!minima.empty()) {
          const bool isroot = myproc == metafile_ioproc;
          const void *const minima_sendptr =
              isroot ? MPI_IN_PLACE : minima.data();
          void *const minima_recvptr = isroot ? minima.data() : nullptr;
          MPI_Reduce(minima_sendptr, minima_recvptr, dim * minima.size(),
                     mpi_datatype_v<CCTK_REAL>, MPI_MIN, metafile_ioproc,
                     mpi_comm);
          const void *const maxima_sendptr =
              isroot ? MPI_IN_PLACE : maxima.data();
          void *const maxima_recvptr = isroot ? maxima.data() : nullptr;
          MPI_Reduce(maxima_sendptr, maxima_recvptr, dim * maxima.size(),
                     mpi_datatype_v<CCTK_REAL>, MPI_MAX, metafile_ioproc,
                     mpi_comm);
        }
      }
    }
  }

  interval_data = nullptr;

  static Timer timer_meta("OutputSilo.meta");
  auto interval_meta = std::make_unique<Interval>(timer_meta);

  // Write metadata
  if (write_metafile) {

    const std::string metafilename = std::string(output_dir) + "/" +
                                     make_filename(output_file, cctk_iteration);
    const DB::ptr<DBfile> metafile =
        DB::make(DBCreate(metafilename.c_str(), DB_CLOBBER, DB_LOCAL,
                          output_file.c_str(), DB_HDF5));
    assert(metafile);

    {
      char *const data = IOUtil_GetAllParameters(cctkGH, 1 /*all*/);
      const std::string parameters(data);
      std::free(data);
      const int dims = parameters.length();
      ierr = DBWrite(metafile.get(), "AllParameters", parameters.data(), &dims,
                     1, DB_CHAR);
      assert(!ierr);
    }

    {
      // Tell VisIt that the mesh structure may change over time
      const int dims = 1;
      const int value = 1;
      ierr = DBWrite(metafile.get(), "MetadataIsTimeVarying", &value, &dims, 1,
                     DB_INT);
      assert(!ierr);
    }

    // Loop over groups
    // Map each shared multimesh to its block count so the multivar block below
    // can check that every slice multivar references the same number of blocks.
    std::map<mesh_props_t, std::size_t> mesh_nblocks;
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      if (!output_group.at(gi))
        continue;
      // TODO: Output grid arrays
      if (CCTK_GroupTypeI(gi) != CCTK_GF)
        continue;

      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const int numvars = groupdata0.numvars;
      const int tl = 0;
      const amrex::MultiFab &mfab0 = *groupdata0.mfab[tl];
      const amrex::IndexType &indextype = mfab0.ixType();

      const std::array<int, dim> &nghosts = groupdata0.nghostzones;
      const mesh_props_t mesh_props{nghosts};
      const bool have_mesh = mesh_nblocks.count(mesh_props);

      if (!have_mesh) {
        const std::string multimeshname = make_meshname(nghosts);

        // The AMR MRG tree enumerates all components unconditionally; a slice
        // writes only intersecting ones, so the tree would dangle. Slices get
        // a flat multimesh with no MRG tree.
        if (!slice) {
          // Count components
          const int nlevels = ghext->num_levels();
          const int npatches = ghext->num_patches();
          std::vector<int> comp0_level(nlevels);
          std::vector<int> ncomps_level(nlevels);
          std::vector<std::vector<int> > comp0_level_patch(nlevels);
          std::vector<std::vector<int> > ncomps_level_patch(nlevels);
          for (int level = 0; level < nlevels; ++level) {
            comp0_level_patch.at(level).resize(npatches, 0);
            ncomps_level_patch.at(level).resize(npatches, 0);
          }
          int ncomps_total = 0;
          for (int level = 0; level < nlevels; ++level) {
            comp0_level.at(level) = ncomps_total;
            for (int patch = 0; patch < npatches; ++patch) {
              const auto &patchdata = ghext->patchdata.at(patch);
              if (level < int(patchdata.leveldata.size())) {
                const auto &leveldata = patchdata.leveldata.at(level);
                comp0_level_patch.at(level).at(patch) = ncomps_total;
                const amrex::DistributionMapping &dm =
                    leveldata.fab->DistributionMap();
                const int ncomponents = dm.size();
                ncomps_total += ncomponents;
                ncomps_level_patch.at(level).at(patch) =
                    ncomps_total - comp0_level_patch.at(level).at(patch);
              }
            }
            ncomps_level.at(level) = ncomps_total - comp0_level.at(level);
          }

          // Describe which components belong to which level
          // Question: Can this name be changed?
          const std::string levelmaps_name =
              multimeshname + "_wmrgtree_lvlMaps";
          {
            std::vector<int> segment_types;
            std::vector<std::vector<int> > segment_data;
            segment_types.reserve(nlevels);
            segment_data.reserve(nlevels);
            for (int l = 0; l < nlevels; ++l) {
              const int comp0 = comp0_level.at(l);
              const int ncomps = ncomps_level.at(l);
              std::vector<int> data;
              data.reserve(ncomps);
              for (int c = 0; c < ncomps; ++c)
                data.push_back(comp0 + c);
              segment_types.push_back(DB_BLOCKCENT);
              segment_data.push_back(std::move(data));
            }
            std::vector<int> segment_lengths;
            std::vector<const int *> segment_data_ptrs;
            segment_lengths.reserve(segment_data.size());
            segment_data_ptrs.reserve(segment_data.size());
            for (const auto &data : segment_data) {
              segment_lengths.push_back(data.size());
              segment_data_ptrs.push_back(data.data());
            }

            ierr = DBPutGroupelmap(
                metafile.get(), levelmaps_name.c_str(), nlevels,
                segment_types.data(), segment_lengths.data(), nullptr,
                segment_data_ptrs.data(), nullptr, 0, nullptr);
            assert(!ierr);
          }

          // Describe which components are children of which other components
          // Question: Can this name be changed?
          const std::string childmaps_name =
              multimeshname + "_wmrgtree_chldMaps";
          std::vector<int> num_children;
          {
            std::vector<int> segment_types;
            std::vector<std::vector<int> > segment_data;
            segment_types.reserve(ncomps_total);
            segment_data.reserve(ncomps_total);

            for (const auto &patchdata : ghext->patchdata) {
              const int patch = patchdata.patch;
              for (const auto &leveldata : patchdata.leveldata) {
                const int level = leveldata.level;
                const int fine_level = level + 1;
                if (fine_level < int(patchdata.leveldata.size())) {
                  const auto &groupdata = *leveldata.groupdata.at(gi);
                  const amrex::MultiFab &mfab = *groupdata.mfab[tl];
                  const auto &fine_leveldata =
                      patchdata.leveldata.at(fine_level);
                  const auto &fine_groupdata = *fine_leveldata.groupdata.at(gi);
                  const amrex::MultiFab &fine_mfab = *fine_groupdata.mfab[tl];

                  const int ncomps = ncomps_level_patch.at(level).at(patch);
                  const int fine_comp0 =
                      comp0_level_patch.at(fine_level).at(patch);
                  const amrex::BoxArray &fine_boxarray = fine_mfab.boxarray;

                  for (int component = 0; component < ncomps; ++component) {
                    const amrex::Box &box = mfab.box(component); // interior
                    amrex::Box refined_box(box);
                    refined_box.refine(2);
                    const std::vector<pair<int, amrex::Box> > child_boxes =
                        fine_boxarray.intersections(refined_box);
                    std::vector<int> children;
                    children.reserve(child_boxes.size());
                    for (const auto &ib : child_boxes) {
                      const int fine_component = ib.first;
                      children.push_back(fine_comp0 + fine_component);
                    }

                    segment_types.push_back(DB_BLOCKCENT);
                    segment_data.push_back(std::move(children));
                  }

                } else {
                  // no finer level, hence no children
                  const int ncomps = ncomps_level.at(level);
                  for (int component = 0; component < ncomps; ++component) {
                    segment_types.push_back(DB_BLOCKCENT);
                    segment_data.emplace_back();
                  }
                }
              } // loop over all levels
            } // loop over all patches

            std::vector<int> &segment_lengths = num_children;
            std::vector<const int *> segment_data_ptrs;
            segment_lengths.reserve(segment_data.size());
            segment_data_ptrs.reserve(segment_data.size());
            for (const auto &data : segment_data) {
              segment_lengths.push_back(data.size());
              segment_data_ptrs.push_back(data.data());
            }

            ierr = DBPutGroupelmap(
                metafile.get(), childmaps_name.c_str(), ncomps_total,
                segment_types.data(), segment_lengths.data(), nullptr,
                segment_data_ptrs.data(), nullptr, 0, nullptr);
            assert(!ierr);
          }

          // Create Mrgtree
          {
            const int max_mgrtree_children = 2;
            const DB::ptr<DBmrgtree> mrgtree = DB::make(
                DBMakeMrgtree(DB_MULTIMESH, 0, max_mgrtree_children, nullptr));
            assert(mrgtree);

            // Describe AMR configuration
            const int max_amr_decomp_children = 2;
            ierr = DBAddRegion(mrgtree.get(), "amr_decomp", 0,
                               max_amr_decomp_children, nullptr, 0, nullptr,
                               nullptr, nullptr, nullptr);
            assert(!ierr);
            ierr = DBSetCwr(mrgtree.get(), "amr_decomp");
            assert(ierr >= 0);

            // Describe AMR levels
            {
              ierr = DBAddRegion(mrgtree.get(), "levels", 0, nlevels, nullptr,
                                 0, nullptr, nullptr, nullptr, nullptr);
              assert(!ierr);

              ierr = DBSetCwr(mrgtree.get(), "levels");
              assert(ierr >= 0);

              const std::vector<std::string> region_names{"@level%d@n"};
              std::vector<const char *> region_name_ptrs;
              region_name_ptrs.reserve(region_names.size());
              for (const std::string &name : region_names)
                region_name_ptrs.push_back(name.c_str());

              std::vector<int> segment_ids;
              std::vector<int> segment_types;
              segment_ids.reserve(nlevels);
              segment_types.reserve(nlevels);
              for (int l = 0; l < nlevels; ++l) {
                segment_ids.push_back(l);
                segment_types.push_back(DB_BLOCKCENT);
              }

              ierr = DBAddRegionArray(
                  mrgtree.get(), nlevels, region_name_ptrs.data(), 0,
                  levelmaps_name.c_str(), 1, segment_ids.data(),
                  ncomps_level.data(), segment_types.data(), nullptr);
              assert(!ierr);

              ierr = DBSetCwr(mrgtree.get(), "..");
              assert(ierr >= 0);
            }

            // Describe AMR children
            {
              ierr =
                  DBAddRegion(mrgtree.get(), "patches", 0, ncomps_total,
                              nullptr, 0, nullptr, nullptr, nullptr, nullptr);
              assert(ierr >= 0);

              ierr = DBSetCwr(mrgtree.get(), "patches");
              assert(ierr >= 0);

              const std::vector<std::string> region_names{"@patch%d@n"};
              std::vector<const char *> region_name_ptrs;
              region_name_ptrs.reserve(region_names.size());
              for (const std::string &name : region_names)
                region_name_ptrs.push_back(name.c_str());

              std::vector<int> segment_types;
              std::vector<int> segment_ids;
              segment_types.reserve(ncomps_total);
              segment_ids.reserve(ncomps_total);
              for (int c = 0; c < ncomps_total; ++c) {
                segment_ids.push_back(c);
                segment_types.push_back(DB_BLOCKCENT);
              }

              ierr = DBAddRegionArray(
                  mrgtree.get(), ncomps_total, region_name_ptrs.data(), 0,
                  childmaps_name.c_str(), 1, segment_ids.data(),
                  num_children.data(), segment_types.data(), nullptr);

              ierr = DBSetCwr(mrgtree.get(), "..");
              assert(ierr >= 0);
            }

            {
              const std::vector<std::string> mrgv_onames{
                  multimeshname + "_wmrgtree_lvlRatios",
                  multimeshname + "_wmrgtree_ijkExts",
                  multimeshname + "_wmrgtree_xyzExts", "rank"};
              std::vector<const char *> mrgv_oname_ptrs;
              mrgv_oname_ptrs.reserve(mrgv_onames.size() + 1);
              for (const std::string &name : mrgv_onames)
                mrgv_oname_ptrs.push_back(name.c_str());
              mrgv_oname_ptrs.push_back(nullptr);

              const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
              assert(optlist);

              ierr = DBAddOption(optlist.get(), DBOPT_MRGV_ONAMES,
                                 mrgv_oname_ptrs.data());
              assert(!ierr);

              ierr = DBPutMrgtree(metafile.get(), "mrgTree", "amr_mesh",
                                  mrgtree.get(), optlist.get());
              assert(!ierr);
            }
          }

          // Write refinement ratios
          {
            const std::string levelrationame =
                multimeshname + "_wmrgtree_lvlRatios";

            const std::vector<std::string> compnames{"iRatio", "jRatio",
                                                     "kRatio"};
            std::vector<const char *> compname_ptrs;
            compname_ptrs.reserve(compnames.size());
            for (const std::string &name : compnames)
              compname_ptrs.push_back(name.c_str());

            const std::vector<std::string> regionnames{"@level%d@n"};
            std::vector<const char *> regionname_ptrs;
            regionname_ptrs.reserve(regionnames.size());
            for (const std::string &name : regionnames)
              regionname_ptrs.push_back(name.c_str());

            std::array<std::vector<int>, ndims> data;
            for (int d = 0; d < ndims; ++d) {
              data[d].reserve(1);
              data[d].push_back(2);
            }
            std::array<const void *, ndims> data_ptrs;
            for (int d = 0; d < ndims; ++d)
              data_ptrs[d] = data[d].data();

            ierr = DBPutMrgvar(metafile.get(), levelrationame.c_str(),
                               "mrgTree", ndims, compname_ptrs.data(), nlevels,
                               regionname_ptrs.data(), DB_INT, data_ptrs.data(),
                               nullptr);
            assert(!ierr);
          }
        } // if !slice (MRG tree)

        // Extents feed the multimesh DBOPT_EXTENTS option, so this runs for
        // both slice and non-slice output; for a slice it is filtered to the
        // slab.
        typedef std::array<std::array<int, ndims>, 2> iextent_t;
        typedef std::array<std::array<CCTK_REAL, ndims>, 2> extent_t;
        std::vector<iextent_t> iextents;
        std::vector<extent_t> extents;

        for (const auto &patchdata : ghext->patchdata) {
          for (const auto &leveldata : patchdata.leveldata) {
            const auto &groupdata = *leveldata.groupdata.at(gi);
            const int tl = 0;
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const int ncomponents = mfab.size();
            if (patchdata.is_cartesian) {
              const amrex::Geometry &geom =
                  patchdata.amrcore->Geom(leveldata.level);
              const amrex::Real *const x0 = geom.ProbLo();
              const amrex::Real *const dx = geom.CellSize();
              for (int c = 0; c < ncomponents; ++c) {
                const amrex::Box &fabbox = mfab.fabbox(c); // exterior
                Arith::vect<int, 3> sub_lo{
                    fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)};
                Arith::vect<int, 3> sub_hi{fabbox.bigEnd(0) + 1,
                                           fabbox.bigEnd(1) + 1,
                                           fabbox.bigEnd(2) + 1};
                if (slice) {
                  const auto r = slice_slab(mfab, c, x0, dx, indextype, *slice);
                  if (!r)
                    continue; // skip components that miss the plane
                  sub_lo = r->first;
                  sub_hi = r->second;
                }
                iextent_t iextent;
                extent_t extent;
                for (int d = 0; d < ndims; ++d) {
                  iextent[0][d] = sub_lo[d];
                  iextent[1][d] = sub_hi[d] - 1;
                  extent[0][d] = x0[d] + sub_lo[d] * dx[d];
                  extent[1][d] = x0[d] + (sub_hi[d] - 1) * dx[d];
                }
                iextents.push_back(iextent);
                extents.push_back(extent);
              }
            } else { // if !patchdata.is_cartesian)
              // The function-entry guard rejects non-Cartesian slices, so this
              // branch is unreachable for slices and needs no slice_slab skip.
              assert(!slice);
              const auto &minima =
                  coordinate_minima.at(patchdata.patch).at(leveldata.level);
              const auto &maxima =
                  coordinate_maxima.at(patchdata.patch).at(leveldata.level);
              for (int c = 0; c < ncomponents; ++c) {
                const amrex::Box &fabbox = mfab.fabbox(c); // exterior
                iextent_t iextent;
                extent_t extent;
                for (int d = 0; d < ndims; ++d) {
                  iextent[0][d] = fabbox.smallEnd(d);
                  iextent[1][d] = fabbox.bigEnd(d);
                  extent[0][d] = minima.at(c)[d];
                  extent[1][d] = maxima.at(c)[d];
                }
                iextents.push_back(iextent);
                extents.push_back(extent);
              }
            } // if !patchdata.is_cartesian)
          }
        }

        // Write extents (MRG block; omitted for slices). Non-slice:
        // iextents.size() == ncomps_total.
        if (!slice) {
          const int ncomps = int(iextents.size());
          const std::string iextentsname = multimeshname + "_wmrgtree_ijkExts";
          const std::string extentsname = multimeshname + "_wmrgtree_xyzExts";

          const std::vector<std::string> icompnames{"iMin", "iMax", "jMin",
                                                    "jMax", "kMin", "kMax"};
          // TODO: These are local coordinates -- are they useful? Should they
          // be omitted? Look at `patchdata.is_cartesian`.
          const std::vector<std::string> compnames{"xMin", "xMax", "yMin",
                                                   "yMax", "zMin", "zMax"};
          std::vector<const char *> icompname_ptrs;
          icompname_ptrs.reserve(icompnames.size());
          for (const std::string &name : icompnames)
            icompname_ptrs.push_back(name.c_str());
          std::vector<const char *> compname_ptrs;
          compname_ptrs.reserve(compnames.size());
          for (const std::string &name : compnames)
            compname_ptrs.push_back(name.c_str());

          const std::vector<std::string> regionnames{"@patch%d@n"};
          std::vector<const char *> regionname_ptrs;
          regionname_ptrs.reserve(regionnames.size());
          for (const std::string &name : regionnames)
            regionname_ptrs.push_back(name.c_str());

          std::array<std::array<std::vector<int>, 2>, ndims> idata;
          std::array<std::array<std::vector<CCTK_REAL>, 2>, ndims> data;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata[d][f].reserve(ncomps);
              data[d][f].reserve(ncomps);
            }
          }
          for (int c = 0; c < ncomps; ++c) {
            for (int d = 0; d < ndims; ++d) {
              for (int f = 0; f < 2; ++f) {
                idata[d][f].push_back(iextents[c][f][d]);
                data[d][f].push_back(extents[c][f][d]);
              }
            }
          }
          std::array<std::array<const void *, 2>, ndims> idata_ptrs;
          std::array<std::array<const void *, 2>, ndims> data_ptrs;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata_ptrs[d][f] = idata[d][f].data();
              data_ptrs[d][f] = data[d][f].data();
            }
          }

          ierr = DBPutMrgvar(metafile.get(), iextentsname.c_str(), "mrgTree",
                             2 * ndims, icompname_ptrs.data(), ncomps,
                             regionname_ptrs.data(), DB_INT, idata_ptrs.data(),
                             nullptr);
          assert(!ierr);

          ierr = DBPutMrgvar(metafile.get(), extentsname.c_str(), "mrgTree",
                             2 * ndims, compname_ptrs.data(), ncomps,
                             regionname_ptrs.data(), db_datatype_v<CCTK_REAL>,
                             data_ptrs.data(), nullptr);
          assert(!ierr);

          // Write rank
          const std::vector<int> ranks(ncomps, ndims);
          const std::vector<const void *> rank_ptrs{ranks.data()};
          ierr = DBPutMrgvar(metafile.get(), "rank", "mrgTree", 1, nullptr,
                             ncomps, regionname_ptrs.data(), DB_INT,
                             rank_ptrs.data(), nullptr);
          assert(!ierr);
        }

        // Write multimesh

        std::vector<std::string> meshnames;
        for (const auto &patchdata : ghext->patchdata) {
          for (const auto &leveldata : patchdata.leveldata) {
            const auto &groupdata = *leveldata.groupdata.at(gi);
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const amrex::DistributionMapping &dm = mfab.DistributionMap();
            const int ncomponents = dm.size();
            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const amrex::Real *const x0 = geom.ProbLo();
            const amrex::Real *const dx = geom.CellSize();
            for (int c = 0; c < ncomponents; ++c) {
              if (slice && !slice_slab(mfab, c, x0, dx, indextype, *slice))
                continue; // skip components that miss the plane
              const int proc = dm[c];
              const std::string proc_filename =
                  make_subdirname(output_file, cctk_iteration) + "/" +
                  make_filename(output_file, cctk_iteration,
                                proc / ioproc_every);
              const std::string meshname =
                  proc_filename + ":" +
                  make_meshname(nghosts, patchdata.patch, leveldata.level, c);
              meshnames.push_back(meshname);
            }
          }
        }
        std::vector<const char *> meshname_ptrs;
        meshname_ptrs.reserve(meshnames.size());
        for (const auto &meshname : meshnames)
          meshname_ptrs.push_back(meshname.c_str());

        const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
        assert(optlist);

        int cycle = cctk_iteration;
        ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
        assert(!ierr);

        // float time = cctk_time;
        // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
        // assert(!ierr);
        double dtime = cctk_time;
        ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
        assert(!ierr);

        int quadmesh = DB_QUADMESH;
        ierr = DBAddOption(optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadmesh);
        assert(!ierr);

        int extents_size = 2 * ndims;
        // This needs to have type `double`, even if everything else is
        // `float`
        typedef std::array<std::array<double, ndims>, 2> dextent_t;
#ifdef CCTK_REAL_PRECISION_8
        std::vector<dextent_t> &dextents = extents;
#else
        std::vector<dextent_t> dextents;
        dextents.resize(extents.size());
        for (size_t n = 0; n < dextents.size(); ++n) {
          const auto &ext = extents[n];
          dextent_t dext;
          for (int f = 0; f < 2; ++f)
            for (int d = 0; d < ndims; ++d)
              dext[f][d] = ext[f][d];
          dextents[n] = dext;
        }
#endif
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS_SIZE, &extents_size);
        assert(!ierr);
        assert(extents.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_EXTENTS, dextents.data());
        assert(!ierr);

        std::vector<int> zonecounts;
        zonecounts.reserve(meshnames.size());
        for (const auto &patchdata : ghext->patchdata) {
          for (const auto &leveldata : patchdata.leveldata) {
            const auto &groupdata = *leveldata.groupdata.at(gi);
            const int tl = 0;
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const int ncomponents = mfab.size();
            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const amrex::Real *const x0 = geom.ProbLo();
            const amrex::Real *const dx = geom.CellSize();
            for (int c = 0; c < ncomponents; ++c) {
              const amrex::Box &fabbox = mfab.fabbox(c); // exterior
              Arith::vect<int, 3> sub_lo{fabbox.smallEnd(0), fabbox.smallEnd(1),
                                         fabbox.smallEnd(2)};
              Arith::vect<int, 3> sub_hi{fabbox.bigEnd(0) + 1,
                                         fabbox.bigEnd(1) + 1,
                                         fabbox.bigEnd(2) + 1};
              if (slice) {
                const auto r = slice_slab(mfab, c, x0, dx, indextype, *slice);
                if (!r)
                  continue; // skip components that miss the plane
                sub_lo = r->first;
                sub_hi = r->second;
              }
              std::array<int, ndims> dims_vc;
              for (int d = 0; d < ndims; ++d)
                dims_vc[d] =
                    (sub_hi[d] - sub_lo[d]) + int(indextype.cellCentered(d));
              int zonecount = 1;
              for (int d = 0; d < ndims; ++d)
                zonecount *= dims_vc[d];
              zonecounts.push_back(zonecount);
            }
          }
        }
        assert(zonecounts.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_ZONECOUNTS, zonecounts.data());
        assert(!ierr);

        // The multimesh references the MRG tree only when one was written.
        const std::string mrgtreename = "mrgtree";
        if (!slice) {
          ierr = DBAddOption(optlist.get(), DBOPT_MRGTREE_NAME,
                             const_cast<char *>(mrgtreename.c_str()));
          assert(!ierr);
        }

        ierr = DBPutMultimesh(metafile.get(), multimeshname.c_str(),
                              meshname_ptrs.size(), meshname_ptrs.data(),
                              nullptr, optlist.get());
        assert(!ierr);

        mesh_nblocks[mesh_props] = meshname_ptrs.size();
      } // if write multimesh

      // Write multivar
      {
        const std::string multimeshname = make_meshname(nghosts);

        const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
        assert(optlist);

        int cycle = cctk_iteration;
        ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
        assert(!ierr);

        // float time = cctk_time;
        // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
        // assert(!ierr);
        double dtime = cctk_time;
        ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
        assert(!ierr);

        ierr = DBAddOption(optlist.get(), DBOPT_MMESH_NAME,
                           const_cast<char *>(multimeshname.c_str()));
        assert(!ierr);

        int vartype_scalar = DB_VARTYPE_SCALAR;
        ierr = DBAddOption(optlist.get(), DBOPT_TENSOR_RANK, &vartype_scalar);
        assert(!ierr);

        int quadvar = DB_QUADVAR;
        ierr = DBAddOption(optlist.get(), DBOPT_MB_BLOCK_TYPE, &quadvar);
        assert(!ierr);

        for (int vi = 0; vi < numvars; ++vi) {
          const std::string multivarname = make_varname(gi, vi);

          std::vector<std::string> varnames;
          for (const auto &patchdata : ghext->patchdata) {
            for (const auto &leveldata : patchdata.leveldata) {
              const auto &groupdata = *leveldata.groupdata.at(gi);
              const int tl = 0;
              const amrex::MultiFab &mfab = *groupdata.mfab[tl];
              const amrex::DistributionMapping &dm = mfab.DistributionMap();
              const int ncomponents = dm.size();
              const amrex::Geometry &geom =
                  patchdata.amrcore->Geom(leveldata.level);
              const amrex::Real *const x0 = geom.ProbLo();
              const amrex::Real *const dx = geom.CellSize();
              for (int c = 0; c < ncomponents; ++c) {
                // Same guard as the multimesh loops, keeping multivar aligned.
                if (slice && !slice_slab(mfab, c, x0, dx, indextype, *slice))
                  continue; // skip components that miss the plane
                const int proc = dm[c];
                const std::string proc_filename =
                    make_subdirname(output_file, cctk_iteration) + "/" +
                    make_filename(output_file, cctk_iteration,
                                  proc / ioproc_every);
                const std::string varname =
                    proc_filename + ":" +
                    make_varname(gi, vi, patchdata.patch, leveldata.level, c);
                varnames.push_back(varname);
              }
            }
          }
          std::vector<const char *> varname_ptrs;
          varname_ptrs.reserve(varnames.size());
          for (const auto &varname : varnames)
            varname_ptrs.push_back(varname.c_str());

          // A slice multimesh and every multivar referencing it must enumerate
          // the same block set; a divergence means block membership has become
          // centering-dependent. Abort rather than emit a file VisIt cannot
          // read. Collective-safe (identical on every rank).
          if (slice) {
            const std::size_t nmesh = mesh_nblocks.at(mesh_props);
            if (varname_ptrs.size() != nmesh)
              CCTK_VERROR(
                  "2D Silo slice multivar \"%s\" enumerates %zu blocks but its "
                  "multimesh \"%s\" enumerates %zu. The slice plane likely "
                  "lies "
                  "on an inter-box boundary and block membership has become "
                  "centering-dependent (see restrict_box_interior in "
                  "io_slice.hxx / slice_slab in io_silo.cxx).",
                  multivarname.c_str(), varname_ptrs.size(),
                  multimeshname.c_str(), nmesh);
          }

          ierr = DBPutMultivar(metafile.get(), multivarname.c_str(),
                               varname_ptrs.size(), varname_ptrs.data(),
                               nullptr, optlist.get());
          assert(!ierr);
        } // for vi
      } // write multivar

    } // for gi

    // Write internal driver state
    {
      const std::string dirname = DB::legalize_name(driver_name);
      ierr = DBMkDir(metafile.get(), dirname.c_str());
      assert(!ierr);

      // Write number of levels
      {
        const auto num_patches = ghext->num_patches();

        std::vector<int> values(num_patches);

        for (const auto &patchdata : ghext->patchdata) {
          values.at(patchdata.patch) = ghext->num_levels(patchdata.patch);
        }

        const std::string varname =
            dirname + "/" + DB::legalize_name("nlevels");

        ierr = DBWrite(metafile.get(), varname.c_str(), values.data(),
                       &num_patches, 1, DB_INT);
        assert(!ierr);
      }

      // Write FabArrayBase (component positions and shapes)
      for (const auto &patchdata : ghext->patchdata) {
        for (const auto &leveldata : patchdata.leveldata) {
          const amrex::FabArrayBase &fab = *leveldata.fab;
          const int ncomponents = fab.size();
          std::vector<int> boxes(2 * ndims * ncomponents);
          for (int component = 0; component < ncomponents; ++component) {
            const amrex::Box &fabbox = fab.box(component); // valid region
            for (int d = 0; d < ndims; ++d)
              boxes[d + 2 * ndims * component] = fabbox.smallEnd(d);
            for (int d = 0; d < ndims; ++d)
              boxes[d + ndims + 2 * ndims * component] = fabbox.bigEnd(d);
          }
          const int dims[2] = {ncomponents, 2 * ndims};
          const std::string varname =
              dirname + "/" +
              make_fabarraybasename(patchdata.patch, leveldata.level);
          ierr = DBWrite(metafile.get(), varname.c_str(), boxes.data(), dims, 2,
                         DB_INT);
          assert(!ierr);
        }
      }

      // Write per-level iteration as rational (numerator/denominator)
      for (const auto &patchdata : ghext->patchdata) {
        for (const auto &leveldata : patchdata.leveldata) {
          {
            const std::string varname =
                dirname + "/" +
                DB::legalize_name("iteration_num.m" +
                                  std::to_string(patchdata.patch) + ".rl" +
                                  std::to_string(leveldata.level));
            long long val = leveldata.iteration.num;
            int dims = 1;
            ierr = DBWrite(metafile.get(), varname.c_str(), &val, &dims, 1,
                           DB_LONG_LONG);
            assert(!ierr);
          }
          {
            const std::string varname =
                dirname + "/" +
                DB::legalize_name("iteration_den.m" +
                                  std::to_string(patchdata.patch) + ".rl" +
                                  std::to_string(leveldata.level));
            long long val = leveldata.iteration.den;
            int dims = 1;
            ierr = DBWrite(metafile.get(), varname.c_str(), &val, &dims, 1,
                           DB_LONG_LONG);
            assert(!ierr);
          }
        }
      }
    }

    {
      const std::string visitname = [&]() {
        std::ostringstream buf;
        buf << output_dir << "/" << output_file << ".silo.visit";
        return buf.str();
      }();
      ofstream visit(visitname, ios::app);
      assert(visit.good());
      visit << make_filename(output_file, cctk_iteration) << "\n";
    }

    {
      output_file_description_t ofd;
      ofd.filename = metafilename;
      ofd.description = slice ? "2D CarpetX HDF5 Silo slice output"
                              : "3D CarpetX HDF5 Silo output";
      ofd.writer_thorn = CCTK_THORNSTRING;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!output_group.at(gi))
          continue;
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;
        const int numvars = CCTK_NumVarsInGroupI(gi);
        assert(numvars >= 0);
        const int firstvar = CCTK_FirstVarIndexI(gi);
        assert(numvars == 0 || firstvar >= 0);
        for (int vi = 0; vi < numvars; ++vi) {
          ofd.variables.push_back(CCTK_FullVarName(firstvar + vi));
        }
      }
      ofd.iterations = {cctk_iteration};
      if (slice)
        for (const int d : slice->inplane_dirs())
          ofd.output_directions.push_back(d);
      else
        for (int d = 0; d < dim; ++d)
          ofd.output_directions.push_back(d);
      ofd.format_name = "CarpetX/Silo/HDF5";
      ofd.format_version = {1, 0, 0};

      OutputMeta_RegisterOutputFile(std::move(ofd));
    }

  } // if write metadata

  interval_meta = nullptr;

  if (io_verbose)
    timer.print();
}

} // namespace CarpetX

#endif
