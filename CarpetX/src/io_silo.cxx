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
using namespace std::experimental;
#else
#include <filesystem>
#endif
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace std;

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
template <> struct db_datatype<char> {
  static constexpr int value = DB_CHAR;
};
template <> struct db_datatype<short> {
  static constexpr int value = DB_SHORT;
};
template <> struct db_datatype<int> {
  static constexpr int value = DB_INT;
};
template <> struct db_datatype<long> {
  static constexpr int value = DB_LONG;
};
template <> struct db_datatype<long long> {
  static constexpr int value = DB_LONG_LONG;
};
template <> struct db_datatype<float> {
  static constexpr int value = DB_FLOAT;
};
template <> struct db_datatype<double> {
  static constexpr int value = DB_DOUBLE;
};

struct mesh_props_t {
  amrex::IntVect ngrow;

  auto to_tuple() const { return make_tuple(ngrow); }
  friend bool operator==(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() == q.to_tuple();
  }
  friend bool operator<(const mesh_props_t &p, const mesh_props_t &q) {
    return p.to_tuple() < q.to_tuple();
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
  const regex filename_regex("[.]it(\\d+)[.]silo$",
                             regex_constants::ECMAScript);
  smatch sm;
  regex_search(file_name, sm, filename_regex);
  if (sm.empty())
    return -1;
  return stoi(sm.str(1));
}

std::string make_meshname(const int reflevel = -1, const int component = -1) {
  assert((reflevel == -1) == (component == -1));
  std::ostringstream buf;
  if (reflevel < 0)
    buf << "gh";
  else
    buf << "box"                                        //
        << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  return DB::legalize_name(buf.str());
}

std::string make_varname(const int gi, const int vi, const int reflevel = -1,
                         const int component = -1) {
  assert((reflevel == -1) == (component == -1));
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
    buf << ".rl" << setw(2) << setfill('0') << reflevel //
        << ".c" << setw(8) << setfill('0') << component;
  return DB::legalize_name(buf.str());
}

const std::string driver_name = "CarpetX";

std::string make_fabarraybasename(const int reflevel) {
  std::ostringstream buf;
  buf << "FabArrayBase.rl" << setw(2) << setfill('0') << reflevel;
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
  auto interval_setup = make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(nullptr);

  const int metafile_ioproc = 0;
  const bool read_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);

  // TODO: directories instead of carefully chosen names

  static Timer timer_parameters("InputSiloParameters.parameters");
  auto interval_parameters = make_unique<Interval>(timer_parameters);

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
  auto interval_setup = make_unique<Interval>(timer_setup);

  const MPI_Comm mpi_comm = MPI_COMM_WORLD;
  const int myproc = CCTK_MyProc(cctkGH);

  const int metafile_ioproc = 0;
  const bool read_metafile = myproc == metafile_ioproc;

  // Configure Silo library
  DBShowErrors(DB_ALL_AND_DRVR, nullptr);

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_meta("InputSiloGridStructure.meta");
  auto interval_meta = make_unique<Interval>(timer_meta);

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

  // TODOPATCH: Handle multiple patches
  assert(ghext->num_patches() == 1);
  const int patch = 0;
  auto &patchdata = ghext->patchdata.at(patch);

  int nlevels;
  if (read_metafile) {
    // Read number of levels
    const std::string varname = dirname + "/" + DB::legalize_name("nlevels");
    const int vartype = DBGetVarType(metafile.get(), varname.c_str());
    assert(vartype == DB_INT);
    const int varlength = DBGetVarLength(metafile.get(), varname.c_str());
    assert(varlength == 1);
    ierr = DBReadVar(metafile.get(), varname.c_str(), &nlevels);
    assert(!ierr);
  }
  MPI_Bcast(&nlevels, 1, MPI_INT, metafile_ioproc, mpi_comm);
  CCTK_VINFO("Found %d levels", nlevels);
  patchdata.amrcore->SetFinestLevel(nlevels - 1);

  // Read FabArrayBase (component positions and shapes)
  for (int level = 0; level < nlevels; ++level) {
    CCTK_VINFO("Reading level %d...", level);

    const std::string varname = dirname + "/" + make_fabarraybasename(level);

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
  }

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
  auto interval_setup = make_unique<Interval>(timer_setup);

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
  auto interval_data = make_unique<Interval>(timer_data);

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

    // TODOPATCH: Handle multiple patches
    assert(ghext->num_patches() == 1);
    const int patch = 0;
    auto &patchdata = ghext->patchdata.at(patch);

    // Loop over levels
    for (const auto &leveldata : patchdata.leveldata) {
      if (io_verbose)
        CCTK_VINFO("Reading patch %d level %d", patchdata.patch,
                   leveldata.level);

      // Loop over groups
      // set<mesh_props_t> have_meshes;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!input_group.at(gi))
          continue;
#warning "TODO: read grid arrays"
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;
        if (io_verbose)
          CCTK_VINFO("  Reading group %s", CCTK_FullGroupName(gi));

        auto &groupdata = *leveldata.groupdata.at(gi);
        const int numvars = groupdata.numvars;
        const int tl = 0;
        amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        // const amrex::IntVect &ngrow = mfab.nGrowVect();
        const amrex::DistributionMapping &dm = mfab.DistributionMap();

        // const mesh_props_t mesh_props{ngrow};
        // const bool have_mesh = have_meshes.count(mesh_props);

        // Loop over components (AMReX boxes)
        const int nfabs = dm.size();
        for (int component = 0; component < nfabs; ++component) {
          if (io_verbose)
            CCTK_VINFO("    Reading component %d", component);

          const int proc = dm[component];
          const int ioproc = proc / ioproc_every * ioproc_every;
          const bool recv_this_fab = proc == myproc;
          const bool read_this_fab = ioproc == myproc;
          if (!(recv_this_fab || read_this_fab))
            continue;

          const amrex::Box &fabbox = mfab.fabbox(component); // exterior

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);
          ptrdiff_t zonecount = 1;
          for (int d = 0; d < ndims; ++d)
            zonecount *= dims[d];
          assert(zonecount >= 0 && zonecount <= INT_MAX);

          // Communicate variable, part 1
          static Timer timer_mpi("InputSilo.mpi");
          auto interval_mpi = make_unique<Interval>(timer_mpi);
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
                      mpi_datatype<CCTK_REAL>::value, ioproc, mpi_tag, mpi_comm,
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

            const std::string meshname =
                make_meshname(leveldata.level, component);

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
              const std::string varname =
                  make_varname(gi, vi, leveldata.level, component);
              if (io_verbose)
                CCTK_VINFO("      Reading variable %s", varname.c_str());

              const DB::ptr<DBquadvar> quadvar =
                  DB::make(DBGetQuadvar(file.get(), varname.c_str()));
              assert(quadvar);

              assert(quadvar->ndims == ndims);
              assert(ndims <= 3);
              for (int d = 0; d < ndims; ++d)
                assert(quadvar->dims[d] == dims[d]);
              assert(quadvar->datatype == db_datatype<CCTK_REAL>::value);
              assert(quadvar->centering == centering);
              assert(quadvar->nvals == 1);

              // TODO: check DBOPT_COORDSYS: int cartesian = DB_CARTESIAN;
              const int column_major = 0;
              assert(quadvar->major_order == column_major);

              const void *const read_ptr = quadvar->vals[0];

              void *const data_ptr = data + vi * zonecount;
              memcpy(data_ptr, read_ptr, zonecount * sizeof(CCTK_REAL));
            } // for vi
          }   // if read_file

          // Communicate variable, part 2
          static Timer timer_wait("InputSilo.wait");
          auto interval_wait = make_unique<Interval>(timer_wait);
          if (recv_this_fab && read_this_fab) {
            // do nothing
          } else if (recv_this_fab) {
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
          } else {
            buffer.resize(numvars * zonecount);
            assert(numvars * zonecount <= INT_MAX);
            MPI_Send(buffer.data(), numvars * zonecount,
                     mpi_datatype<CCTK_REAL>::value, proc, mpi_tag, mpi_comm);
          }
          if (recv_this_fab)
            for (int vi = 0; vi < numvars; ++vi)
              groupdata.valid.at(tl).at(vi).set_all(
                  make_valid_all(), []() { return "read from Silo file"; });
          interval_wait = nullptr;

        } // for component

      } // for gi
    }   // for leveldata
  }     // write data

  interval_data = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void OutputSilo(const cGH *restrict const cctkGH,
                const std::vector<bool> &output_group,
                const std::string &output_dir, const std::string &output_file) {
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
  auto interval_setup = make_unique<Interval>(timer_setup);

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
  DBSetCompression("METHOD=GZIP");
  DBSetEnableChecksums(1);

  // TODO: directories instead of carefully chosen names

  constexpr int ndims = dim;

  interval_setup = nullptr;

  static Timer timer_data("OutputSilo.data");
  auto interval_data = make_unique<Interval>(timer_data);

  // Write data
  {
    static Timer timer_mkdir("OutputSilo.mkdir");
    auto interval_mkdir = make_unique<Interval>(timer_mkdir);
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

    // TODOPATCH: Handle multiple patches
    assert(ghext->num_patches() == 1);
    const int patch = 0;
    auto &patchdata = ghext->patchdata.at(patch);

    // Loop over levels
    for (const auto &leveldata : patchdata.leveldata) {

      // Loop over groups
      set<mesh_props_t> have_meshes;
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
        if (!output_group.at(gi))
          continue;
#warning "TODO: Output grid arrays"
        if (CCTK_GroupTypeI(gi) != CCTK_GF)
          continue;

        const auto &groupdata = *leveldata.groupdata.at(gi);
        const int numvars = groupdata.numvars;
        const int tl = 0;
        const amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        const amrex::IntVect &ngrow = mfab.nGrowVect();
        const amrex::DistributionMapping &dm = mfab.DistributionMap();

        const mesh_props_t mesh_props{ngrow};
        const bool have_mesh = have_meshes.count(mesh_props);

        // Loop over components (AMReX boxes)
        const int nfabs = dm.size();
        for (int component = 0; component < nfabs; ++component) {
          const int proc = dm[component];
          const int ioproc = proc / ioproc_every * ioproc_every;
          const bool send_this_fab = proc == myproc;
          const bool write_this_fab = ioproc == myproc;
          if (!(send_this_fab || write_this_fab))
            continue;

          // TODO: Check whether data are valid
          const amrex::Box &fabbox = mfab.fabbox(component); // exterior

          array<int, ndims> dims;
          for (int d = 0; d < ndims; ++d)
            dims[d] = fabbox.length(d);
          ptrdiff_t zonecount = 1;
          for (int d = 0; d < ndims; ++d)
            zonecount *= dims[d];
          assert(zonecount >= 0 && zonecount <= INT_MAX);

          if (write_file && !have_mesh) {
            static Timer timer_mesh("OutputSilo.mesh");
            Interval interval_mesh(timer_mesh);

            const std::string meshname =
                make_meshname(leveldata.level, component);

            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = dims[d] + int(indextype.cellCentered(d));

            const amrex::Geometry &geom =
                patchdata.amrcore->Geom(leveldata.level);
            const amrex::Real *const x0 = geom.ProbLo();
            const amrex::Real *const dx = geom.CellSize();
            array<vector<CCTK_REAL>, ndims> coords;
            for (int d = 0; d < ndims; ++d) {
              coords[d].resize(dims_vc[d]);
              for (int i = 0; i < dims_vc[d]; ++i)
                coords[d][i] = x0[d] + (fabbox.smallEnd(d) + i) * dx[d];
            }
            array<void *, ndims> coord_ptrs;
            for (int d = 0; d < ndims; ++d)
              coord_ptrs[d] = coords[d].data();

            const DB::ptr<DBoptlist> optlist = DB::make(DBMakeOptlist(10));
            assert(optlist);

            int cartesian = DB_CARTESIAN;
            ierr = DBAddOption(optlist.get(), DBOPT_COORDSYS, &cartesian);
            assert(!ierr);

            int cycle = cctk_iteration;
            ierr = DBAddOption(optlist.get(), DBOPT_CYCLE, &cycle);
            assert(!ierr);

            array<int, ndims> min_index, max_index;
            for (int d = 0; d < ndims; ++d) {
              min_index[d] = ngrow[d];
              max_index[d] = ngrow[d];
            }
            ierr =
                DBAddOption(optlist.get(), DBOPT_LO_OFFSET, min_index.data());
            assert(!ierr);
            ierr =
                DBAddOption(optlist.get(), DBOPT_HI_OFFSET, max_index.data());
            assert(!ierr);

            int column_major = 0;
            ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
            assert(!ierr);

            // float time = cctk_time;
            // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
            // assert(!ierr);
            double dtime = cctk_time;
            ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
            assert(!ierr);

            int hide_from_gui = 1;
            ierr =
                DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI, &hide_from_gui);
            assert(!ierr);

            ierr = DBPutQuadmesh(file.get(), meshname.c_str(), nullptr,
                                 coord_ptrs.data(), dims_vc.data(), ndims,
                                 db_datatype<CCTK_REAL>::value, DB_COLLINEAR,
                                 optlist.get());
            assert(!ierr);
          } // if write mesh

          // Communicate variable
          static Timer timer_mpi("OutputSilo.mpi");
          auto interval_mpi = make_unique<Interval>(timer_mpi);
          const int mpi_tag = 22900; // randomly chosen
          std::vector<CCTK_REAL> buffer;
          const CCTK_REAL *data = nullptr;
          if (send_this_fab && write_this_fab) {
            const amrex::FArrayBox &fab = mfab[component];
            data = fab.dataPtr();
          } else if (send_this_fab) {
            const amrex::FArrayBox &fab = mfab[component];
            assert(numvars * zonecount <= INT_MAX);
            MPI_Send(fab.dataPtr(), numvars * zonecount,
                     mpi_datatype<CCTK_REAL>::value, ioproc, mpi_tag, mpi_comm);
          } else {
            buffer.resize(numvars * zonecount);
            assert(numvars * zonecount <= INT_MAX);
            MPI_Recv(buffer.data(), numvars * zonecount,
                     mpi_datatype<CCTK_REAL>::value, proc, mpi_tag, mpi_comm,
                     MPI_STATUS_IGNORE);
            data = buffer.data();
          }
          interval_mpi = nullptr;

          // Write variable
          if (write_file) {
            static Timer timer_var("OutputSilo.var");
            Interval interval_var(timer_var);

            const std::string meshname =
                make_meshname(leveldata.level, component);

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
            ierr = DBAddOption(optlist.get(), DBOPT_MAJORORDER, &column_major);
            assert(!ierr);

            // float time = cctk_time;
            // ierr = DBAddOption(optlist.get(), DBOPT_TIME, &time);
            // assert(!ierr);
            double dtime = cctk_time;
            ierr = DBAddOption(optlist.get(), DBOPT_DTIME, &dtime);
            assert(!ierr);

            int hide_from_gui = 1;
            ierr =
                DBAddOption(optlist.get(), DBOPT_HIDE_FROM_GUI, &hide_from_gui);
            assert(!ierr);

            if (centering == DB_EDGECENT || centering == DB_FACECENT) {
              // Need to find the other 2 edge- or face-centered
              // variables, and output them as well. Maybe output all
              // 3 when the x- or xy-centered value is output? Maybe
              // add a "sibling" tag to such grid functions to find
              // these other components?
              assert(0);
            }

            for (int vi = 0; vi < numvars; ++vi) {
              const std::string varname =
                  make_varname(gi, vi, leveldata.level, component);

              const void *const data_ptr = data + vi * zonecount;

              ierr = DBPutQuadvar1(
                  file.get(), varname.c_str(), meshname.c_str(), data_ptr,
                  dims.data(), ndims, nullptr, 0, db_datatype<CCTK_REAL>::value,
                  centering, optlist.get());
            } // for vi
          }   // if write_file

        } // for component

      } // for gi
    }   // for leveldata
  }     // write data

  interval_data = nullptr;

  static Timer timer_meta("OutputSilo.meta");
  auto interval_meta = make_unique<Interval>(timer_meta);

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
    set<mesh_props_t> have_meshes;
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi) {
      if (!output_group.at(gi))
        continue;
#warning "TODO: Output grid arrays"
      if (CCTK_GroupTypeI(gi) != CCTK_GF)
        continue;

      const auto &patchdata0 = ghext->patchdata.at(0);
      const auto &leveldata0 = patchdata0.leveldata.at(0);
      const auto &groupdata0 = *leveldata0.groupdata.at(gi);
      const int numvars = groupdata0.numvars;
      const int tl = 0;
      const amrex::MultiFab &mfab0 = *groupdata0.mfab[tl];
      const amrex::IndexType &indextype = mfab0.ixType();
      const amrex::IntVect &ngrow = mfab0.nGrowVect();

      const mesh_props_t mesh_props{ngrow};
      const bool have_mesh = have_meshes.count(mesh_props);

      if (!have_mesh) {

        const std::string multimeshname = make_meshname();

        // TODOPATCH: Handle multiple patches
        assert(ghext->num_patches() == 1);
        const int patch = 0;
        auto &patchdata = ghext->patchdata.at(patch);

        // Count components per level
        const int nlevels = patchdata.leveldata.size();
        std::vector<int> ncomps_level;
        for (const auto &leveldata : patchdata.leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::DistributionMapping &dm = mfab.DistributionMap();
          const int nfabs = dm.size();
          ncomps_level.push_back(nfabs);
        }
        std::vector<int> firstcomp_level;
        int ncomps_total = 0;
        for (const int ncomps : ncomps_level) {
          firstcomp_level.push_back(ncomps_total);
          ncomps_total += ncomps;
        }

        // Describe which components belong to which level
        // Question: Can this name be changed?
        const std::string levelmaps_name = multimeshname + "_wmrgtree_lvlMaps";
        {
          std::vector<int> segment_types;
          std::vector<vector<int> > segment_data;
          segment_types.reserve(nlevels);
          segment_data.reserve(nlevels);
          for (int l = 0; l < nlevels; ++l) {
            const int comp0 = firstcomp_level.at(l);
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

          ierr = DBPutGroupelmap(metafile.get(), levelmaps_name.c_str(),
                                 nlevels, segment_types.data(),
                                 segment_lengths.data(), nullptr,
                                 segment_data_ptrs.data(), nullptr, 0, nullptr);
          assert(!ierr);
        }

        // Describe which components are children of which other components
        // Question: Can this name be changed?
        const std::string childmaps_name = multimeshname + "_wmrgtree_chldMaps";
        std::vector<int> num_children;
        {
          std::vector<int> segment_types;
          std::vector<vector<int> > segment_data;
          segment_types.reserve(ncomps_total);
          segment_data.reserve(ncomps_total);
          for (const auto &leveldata : patchdata.leveldata) {
            const int level = leveldata.level;
            const int fine_level = level + 1;
            if (fine_level < nlevels) {
              const auto &groupdata = *leveldata.groupdata.at(gi);
              const amrex::MultiFab &mfab = *groupdata.mfab[tl];
              const auto &fine_leveldata = patchdata.leveldata.at(fine_level);
              const auto &fine_groupdata = *fine_leveldata.groupdata.at(gi);
              const amrex::MultiFab &fine_mfab = *fine_groupdata.mfab[tl];

              const int ncomps = ncomps_level.at(level);
              const int fine_comp0 = firstcomp_level.at(fine_level);
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
          }

          std::vector<int> &segment_lengths = num_children;
          std::vector<const int *> segment_data_ptrs;
          segment_lengths.reserve(segment_data.size());
          segment_data_ptrs.reserve(segment_data.size());
          for (const auto &data : segment_data) {
            segment_lengths.push_back(data.size());
            segment_data_ptrs.push_back(data.data());
          }

          ierr = DBPutGroupelmap(metafile.get(), childmaps_name.c_str(),
                                 ncomps_total, segment_types.data(),
                                 segment_lengths.data(), nullptr,
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
            ierr = DBAddRegion(mrgtree.get(), "levels", 0, nlevels, nullptr, 0,
                               nullptr, nullptr, nullptr, nullptr);
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
            ierr = DBAddRegion(mrgtree.get(), "patches", 0, ncomps_total,
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

          array<vector<int>, ndims> data;
          for (int d = 0; d < ndims; ++d) {
            data[d].reserve(1);
            data[d].push_back(2);
          }
          array<const void *, ndims> data_ptrs;
          for (int d = 0; d < ndims; ++d)
            data_ptrs[d] = data[d].data();

          ierr = DBPutMrgvar(metafile.get(), levelrationame.c_str(), "mrgTree",
                             ndims, compname_ptrs.data(), nlevels,
                             regionname_ptrs.data(), DB_INT, data_ptrs.data(),
                             nullptr);
          assert(!ierr);
        }

        typedef array<array<int, ndims>, 2> iextent_t;
        typedef array<array<CCTK_REAL, ndims>, 2> extent_t;
        std::vector<iextent_t> iextents;
        std::vector<extent_t> extents;
        iextents.reserve(ncomps_total);
        extents.reserve(ncomps_total);
        for (const auto &leveldata : patchdata.leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::Geometry &geom =
              patchdata.amrcore->Geom(leveldata.level);
          const amrex::Real *const x0 = geom.ProbLo();
          const amrex::Real *const dx = geom.CellSize();
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const amrex::Box &fabbox = mfab.fabbox(c); // exterior
            iextent_t iextent;
            extent_t extent;
            for (int d = 0; d < ndims; ++d) {
              iextent[0][d] = fabbox.smallEnd(d);
              iextent[1][d] = fabbox.bigEnd(d);
              extent[0][d] = x0[d] + fabbox.smallEnd(d) * dx[d];
              extent[1][d] = x0[d] + fabbox.bigEnd(d) * dx[d];
            }
            iextents.push_back(iextent);
            extents.push_back(extent);
          }
        }

        // Write extents
        {
          const std::string iextentsname = multimeshname + "_wmrgtree_ijkExts";
          const std::string extentsname = multimeshname + "_wmrgtree_xyzExts";

          const std::vector<std::string> icompnames{"iMin", "iMax", "jMin",
                                                    "jMax", "kMin", "kMax"};
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

          array<array<vector<int>, 2>, ndims> idata;
          array<array<vector<CCTK_REAL>, 2>, ndims> data;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata[d][f].reserve(ncomps_total);
              data[d][f].reserve(ncomps_total);
            }
          }
          for (int c = 0; c < ncomps_total; ++c) {
            for (int d = 0; d < ndims; ++d) {
              for (int f = 0; f < 2; ++f) {
                idata[d][f].push_back(iextents[c][f][d]);
                data[d][f].push_back(extents[c][f][d]);
              }
            }
          }
          array<array<const void *, 2>, ndims> idata_ptrs;
          array<array<const void *, 2>, ndims> data_ptrs;
          for (int d = 0; d < ndims; ++d) {
            for (int f = 0; f < 2; ++f) {
              idata_ptrs[d][f] = idata[d][f].data();
              data_ptrs[d][f] = data[d][f].data();
            }
          }

          ierr = DBPutMrgvar(metafile.get(), iextentsname.c_str(), "mrgTree",
                             2 * ndims, icompname_ptrs.data(), ncomps_total,
                             regionname_ptrs.data(), DB_INT, idata_ptrs.data(),
                             nullptr);
          assert(!ierr);

          ierr = DBPutMrgvar(
              metafile.get(), extentsname.c_str(), "mrgTree", 2 * ndims,
              compname_ptrs.data(), ncomps_total, regionname_ptrs.data(),
              db_datatype<CCTK_REAL>::value, data_ptrs.data(), nullptr);
          assert(!ierr);

          // Write rank
          const std::vector<int> ranks(ncomps_total, ndims);
          const std::vector<const void *> rank_ptrs{ranks.data()};
          ierr = DBPutMrgvar(metafile.get(), "rank", "mrgTree", 1, nullptr,
                             ncomps_total, regionname_ptrs.data(), DB_INT,
                             rank_ptrs.data(), nullptr);
          assert(!ierr);
        }

        // Write multimesh

        std::vector<std::string> meshnames;
        for (const auto &leveldata : patchdata.leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const amrex::DistributionMapping &dm = mfab.DistributionMap();
          const int nfabs = dm.size();
          for (int c = 0; c < nfabs; ++c) {
            const int proc = dm[c];
            const std::string proc_filename =
                make_subdirname(output_file, cctk_iteration) + "/" +
                make_filename(output_file, cctk_iteration, proc / ioproc_every);
            const std::string meshname =
                proc_filename + ":" + make_meshname(leveldata.level, c);
            meshnames.push_back(meshname);
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
        // This needs to have type `double`, even if everything else is `float`
        typedef array<array<double, ndims>, 2> dextent_t;
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
        for (const auto &leveldata : patchdata.leveldata) {
          const auto &groupdata = *leveldata.groupdata.at(gi);
          const int tl = 0;
          const amrex::MultiFab &mfab = *groupdata.mfab[tl];
          const int nfabs = mfab.size();
          for (int c = 0; c < nfabs; ++c) {
            const amrex::Box &fabbox = mfab.fabbox(c); // exterior
            array<int, ndims> dims_vc;
            for (int d = 0; d < ndims; ++d)
              dims_vc[d] = fabbox.length(d) + int(indextype.cellCentered(d));
            int zonecount = 1;
            for (int d = 0; d < ndims; ++d)
              zonecount *= dims_vc[d];
            zonecounts.push_back(zonecount);
          }
        }
        assert(zonecounts.size() == meshname_ptrs.size());
        ierr = DBAddOption(optlist.get(), DBOPT_ZONECOUNTS, zonecounts.data());
        assert(!ierr);

        const std::string mrgtreename = "mrgtree";
        ierr = DBAddOption(optlist.get(), DBOPT_MRGTREE_NAME,
                           const_cast<char *>(mrgtreename.c_str()));
        assert(!ierr);

        ierr = DBPutMultimesh(metafile.get(), multimeshname.c_str(),
                              meshname_ptrs.size(), meshname_ptrs.data(),
                              nullptr, optlist.get());
        assert(!ierr);

        have_meshes.insert(mesh_props);
      } // if write multimesh

      // Write multivar
      {
        const std::string multimeshname = make_meshname();

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
          const int patch = 0;
          const auto &patchdata = ghext->patchdata.at(patch);
          for (const auto &leveldata : patchdata.leveldata) {
            const auto &groupdata = *leveldata.groupdata.at(gi);
            const int tl = 0;
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const amrex::DistributionMapping &dm = mfab.DistributionMap();
            const int nfabs = dm.size();
            for (int c = 0; c < nfabs; ++c) {
              const int proc = dm[c];
              const std::string proc_filename =
                  make_subdirname(output_file, cctk_iteration) + "/" +
                  make_filename(output_file, cctk_iteration,
                                proc / ioproc_every);
              const std::string varname =
                  proc_filename + ":" +
                  make_varname(gi, vi, leveldata.level, c);
              varnames.push_back(varname);
            }
          }
          std::vector<const char *> varname_ptrs;
          varname_ptrs.reserve(varnames.size());
          for (const auto &varname : varnames)
            varname_ptrs.push_back(varname.c_str());

          ierr = DBPutMultivar(metafile.get(), multivarname.c_str(),
                               varname_ptrs.size(), varname_ptrs.data(),
                               nullptr, optlist.get());
          assert(!ierr);
        } // for vi
      }   // write multivar

    } // for gi

    // Write internal driver state
    {
      const std::string dirname = DB::legalize_name(driver_name);
      ierr = DBMkDir(metafile.get(), dirname.c_str());
      assert(!ierr);

      const int patch = 0;
      const auto &patchdata = ghext->patchdata.at(patch);

      // Write number of levels
      {
        const int dims = 1;
        const int value = patchdata.leveldata.size();
        const std::string varname =
            dirname + "/" + DB::legalize_name("nlevels");
        ierr =
            DBWrite(metafile.get(), varname.c_str(), &value, &dims, 1, DB_INT);
        assert(!ierr);
      }

      // Write FabArrayBase (component positions and shapes)
      for (const auto &leveldata : patchdata.leveldata) {
        const amrex::FabArrayBase &fab = *leveldata.fab;
        const int nfabs = fab.size();
        std::vector<int> boxes(2 * ndims * nfabs);
        for (int component = 0; component < nfabs; ++component) {
          const amrex::Box &fabbox = fab.box(component); // valid region
          for (int d = 0; d < ndims; ++d)
            boxes[d + 2 * ndims * component] = fabbox.smallEnd(d);
          for (int d = 0; d < ndims; ++d)
            boxes[d + ndims + 2 * ndims * component] = fabbox.bigEnd(d);
        }
        const int dims[2] = {nfabs, 2 * ndims};
        const std::string varname =
            dirname + "/" + make_fabarraybasename(leveldata.level);
        ierr = DBWrite(metafile.get(), varname.c_str(), boxes.data(), dims, 2,
                       DB_INT);
        assert(!ierr);
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
      ofd.description = "3D CarpetX HDF5 Silo output";
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
