#include "io_adios2.hxx"

#include "driver.hxx"
#include "timer.hxx"

#include <div.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef HAVE_CAPABILITY_ADIOS2

#include <adios2.h>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
static inline int omp_in_parallel() { return 0; }
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <ios>
#include <iostream>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace CarpetX {

struct carpetx_adios2_t {
  carpetx_adios2_t() = default;

  carpetx_adios2_t(const carpetx_adios2_t &) = delete;
  carpetx_adios2_t(carpetx_adios2_t &&) = default;
  carpetx_adios2_t &operator=(const carpetx_adios2_t &) = delete;
  carpetx_adios2_t &operator=(carpetx_adios2_t &&) = default;

  static std::optional<carpetx_adios2_t> self;

  ////////////////////////////////////////////////////////////////////////////////

  static constexpr bool io_verbose = true;

  static constexpr bool combine_components = true;
  static constexpr bool combine_via_span = true;

  ////////////////////////////////////////////////////////////////////////////////

  template <typename T, std::size_t D> struct box_t {
    std::array<T, D> lo, hi;
    constexpr std::array<T, D> shape() const {
      std::array<T, D> sh;
      for (std::size_t d = 0; d < D; ++d)
        sh[d] = hi[d] < lo[d] ? T{0} : hi[d] - lo[d];
      return sh;
    }
    constexpr T size() const {
      const auto sh = shape();
      T sz{1};
      for (std::size_t d = 0; d < D; ++d)
        sz *= sh[d];
      return sz;
    }
  };

  template <typename T, typename I, std::size_t D> struct level_t {
    box_t<T, D> rdomain;
    std::array<bool, D> is_cell_centred;
    box_t<I, D> idomain;
    constexpr std::array<T, D> rcoord(const std::array<I, D> &icoord) const {
      std::array<T, D> r;
      for (std::size_t d = 0; d < D; ++d) {
        const T rlo = rdomain.lo[d];
        const T rhi = rdomain.hi[d];
        const I ilo2 = 2 * idomain.lo[d] - is_cell_centred;
        const I ihi2 = 2 * idomain.hi[d] + is_cell_centred;
        const I i2 = 2 * icoord[d];
        // The expression below has been carefully constructed to
        // avoid round-off errors at i=ilo and i=ihi. Do not rearrange
        // the terms without preserving this property.
        r[d] = (T(ihi2 - i2) / T(ihi2 - ilo2)) * rlo[d] +
               (T(i2 - ilo2) / T(ihi2 - ilo2)) * rhi[d];
      }
      return r;
    }
    std::vector<box_t<I, D> > grids;
    std::array<std::vector<I>, 2> offsets_sizes() const {
      std::vector<I> offsets(grids.size() + 1), sizes(grids.size());
      I offset{0};
      for (std::size_t n = 0; n < grids.size(); ++n) {
        const T size = grids[n].size();
        offsets[n] = offset;
        sizes[n] = size;
        offset += size;
      }
      offsets[grids.size()] = offset;
      return {offsets, sizes};
    }
  };

  template <typename T, typename I, std::size_t D> struct patch_t {
    box_t<T, D> rdomain;
    std::vector<level_t<T, I, D> > levels;
  };

  template <typename T, typename I, std::size_t D> struct grid_structure_t {
    std::vector<patch_t<T, I, D> > patches;
  };

  ////////////////////////////////////////////////////////////////////////////////

  static std::string make_varname(const int gi, const int vi,
                                  const int patch = -1, const int reflevel = -1,
                                  const int component = -1) {
    std::string varname;
    if (vi < 0) {
      assert(0);
      varname = CCTK_FullGroupName(gi);
    } else {
      const int v0 = CCTK_FirstVarIndexI(gi);
      varname = CCTK_FullVarName(v0 + vi);
      // varname = regex_replace(varname, regex("::"), "-");
      // for (auto &ch : varname)
      //   ch = tolower(ch);
    }
    std::ostringstream buf;
    buf << varname;
    if (patch >= 0)
      buf << ".m" << setw(2) << setfill('0') << patch;
    if (reflevel >= 0)
      buf << ".rl" << setw(2) << setfill('0') << reflevel;
    if (component >= 0)
      buf << ".c" << setw(8) << setfill('0') << component;
    return buf.str();
  }

  ////////////////////////////////////////////////////////////////////////////////

  adios2::ADIOS adios;
  adios2::IO io;
  adios2::Engine engine;

  void OutputADIOS2(const cGH *const cctkGH,
                    const std::vector<bool> &output_group,
                    const std::string &output_dir,
                    const std::string &output_file);
  ~carpetx_adios2_t();
};

////////////////////////////////////////////////////////////////////////////////

void OutputADIOS2(const cGH *const cctkGH,
                  const std::vector<bool> &output_group,
                  const std::string &output_dir,
                  const std::string &output_file) {
  if (!carpetx_adios2_t::self)
    carpetx_adios2_t::self = std::make_optional<carpetx_adios2_t>();
  carpetx_adios2_t::self->OutputADIOS2(cctkGH, output_group, output_dir,
                                       output_file);
}

void ShutdownADIOS2() { carpetx_adios2_t::self.reset(); }

////////////////////////////////////////////////////////////////////////////////

std::optional<carpetx_adios2_t> carpetx_adios2_t::self;

void carpetx_adios2_t::OutputADIOS2(const cGH *const cctkGH,
                                    const std::vector<bool> &output_group,
                                    const std::string &output_dir,
                                    const std::string &output_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr;

  // Set up timers
  static Timer timer("OutputADIOS2");
  Interval interval(timer);

  if (io_verbose)
    CCTK_VINFO("OutputADIOS2...");

  if (std::count(output_group.begin(), output_group.end(), true) == 0)
    return;

  try {

    // const MPI_Comm mpi_comm = MPI_COMM_WORLD;
    // const int myproc = CCTK_MyProc(cctkGH);
    const int nprocs = CCTK_nProcs(cctkGH);

    if (!io) {

      if (io_verbose)
        CCTK_VINFO("  Creating ADIOS object...");
      // "Fortran" enables column-major mode
      adios = adios2::ADIOS(MPI_COMM_WORLD, "Fortran");

      if (io_verbose)
        CCTK_VINFO("  Creating IO object...");
      // The name "IO" must be unique in the ADIOS object
      io = adios.DeclareIO("IO");
      const int num_aggregators = [&]() {
        if (CCTK_EQUALS(out_mode, "proc"))
          return nprocs;
        if (CCTK_EQUALS(out_mode, "np"))
          return div_ceil(nprocs, out_proc_every);
        if (CCTK_EQUALS(out_mode, "onefile"))
          return 1;
        assert(0);
      }();
      assert(num_aggregators > 0);
      io.SetParameters({
          {"NumAggregators", std::to_string(num_aggregators)},
          // {"Threads", std::to_string(omp_get_max_threads())},
      });

      if (io_verbose)
        CCTK_VINFO("  Creating engine...");
      // Create output engine
      std::ostringstream buf;
      const int mode = 0755;
      static once_flag create_directory;
      call_once(create_directory, [&]() {
        const int ierr = CCTK_CreateDirectory(mode, output_dir.c_str());
        assert(ierr >= 0);
      });
      buf << output_dir << "/" << output_file << ".bp";
      const std::string filename = buf.str();
      // This just confirms the default
      // io.SetEngine("BP4");
      engine = io.Open(filename, adios2::Mode::Write);

      if (io_verbose)
        CCTK_VINFO("  Defining variables...");

      // Loop over patches
      for (const auto &patchdata : ghext->patchdata) {

        // Loop over levels
        for (const auto &leveldata : patchdata.leveldata) {

          const int numgroups = CCTK_NumGroups();
          for (int gi = 0; gi < numgroups; ++gi) {
            if (output_group.at(gi)) {

              cGroup cgroup;
              ierr = CCTK_GroupData(gi, &cgroup);
              assert(!ierr);
              assert(cgroup.grouptype == CCTK_GF);
              assert(cgroup.vartype == CCTK_VARIABLE_REAL);
              assert(cgroup.dim == 3);
              // cGroupDynamicData cgroupdynamicdata;
              // ierr = CCTK_GroupDynamicData(cctkGH, gi, &cgroupdynamicdata);
              // assert(!ierr);
              // TODO: Check whether group has storage
              // TODO: Check whether data are valid

              const auto &groupdata = *leveldata.groupdata.at(gi);
              // const int firstvarindex = groupdata.firstvarindex;
              const int numvars = groupdata.numvars;

              const int tl = 0;
              const amrex::MultiFab &mfab = *groupdata.mfab[tl];
              const int num_local_components = mfab.local_size();

              // Loop over variables
              for (int vi = 0; vi < numvars; ++vi) {

                if (!combine_components) {

                  // Loop over components (AMReX boxes)
                  for (int local_component = 0;
                       local_component < num_local_components;
                       ++local_component) {
                    const int component = mfab.IndexArray().at(local_component);
                    const std::string varname = make_varname(
                        gi, vi, patchdata.patch, leveldata.level, component);
                    if (io_verbose)
                      CCTK_VINFO("      Defining variable %s...",
                                 varname.c_str());
                    const adios2::Variable<CCTK_REAL> var =
                        io.DefineVariable<CCTK_REAL>(varname, {}, {},
                                                     {1, 1, 1});
                  } // for local_component

                } else { // if combine_components

                  const std::string varname =
                      make_varname(gi, vi, patchdata.patch, leveldata.level);
                  if (io_verbose)
                    CCTK_VINFO("      Defining variable %s...",
                               varname.c_str());
                  const adios2::Variable<CCTK_REAL> var =
                      io.DefineVariable<CCTK_REAL>(varname, {}, {}, {1});

                } // if combine_components

              } // for vi
            }
          } // for gi

        } // for leveldata
      }   // for patchdata

    } // if !adios2_state

    if (io_verbose)
      CCTK_VINFO("  Beginning step...");
    engine.BeginStep();

    // Loop over patches
    std::vector<patch_t<CCTK_REAL, int, 3> > patches(ghext->patchdata.size());
    for (const auto &patchdata : ghext->patchdata) {

      // Loop over levels
      std::vector<level_t<CCTK_REAL, int, 3> > levels(
          patchdata.leveldata.size());
      for (const auto &leveldata : patchdata.leveldata) {
        if (io_verbose)
          CCTK_VINFO("  Writing patch %d level %d...", patchdata.patch,
                     leveldata.level);

        const int numgroups = CCTK_NumGroups();
        for (int gi = 0; gi < numgroups; ++gi) {
          if (output_group.at(gi)) {
            if (io_verbose)
              CCTK_VINFO("    Writing group %s...", CCTK_FullGroupName(gi));

            cGroup cgroup;
            ierr = CCTK_GroupData(gi, &cgroup);
            assert(!ierr);
            assert(cgroup.grouptype == CCTK_GF);
            assert(cgroup.vartype == CCTK_VARIABLE_REAL);
            assert(cgroup.dim == 3);
            // cGroupDynamicData cgroupdynamicdata;
            // ierr = CCTK_GroupDynamicData(cctkGH, gi, &cgroupdynamicdata);
            // assert(!ierr);
            // TODO: Check whether group has storage
            // TODO: Check whether data are valid

            const auto &groupdata = *leveldata.groupdata.at(gi);
            // const int firstvarindex = groupdata.firstvarindex;
            const int numvars = groupdata.numvars;
            const int tl = 0;
            const amrex::MultiFab &mfab = *groupdata.mfab[tl];
            const amrex::IndexType &indextype = mfab.ixType();
            // const amrex::IntVect &ngrow = mfab.nGrowVect();
            // const amrex::DistributionMapping &dm = mfab.DistributionMap();

            const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
            const double *const xlo = geom.ProbLo();
            const double *const xhi = geom.ProbHi();
            const box_t<CCTK_REAL, 3> rdomain{
              lo : {xlo[0], xlo[1], xlo[2]},
              hi : {xhi[0], xhi[1], xhi[2]}
            };
            const amrex::Box &dom = geom.Domain();
            const box_t<int, 3> idomain{
              lo : {dom.smallEnd(0), dom.smallEnd(1), dom.smallEnd(2)},
              hi : {dom.bigEnd(0) + 1, dom.bigEnd(1) + 1, dom.bigEnd(2) + 1}
            };
            const std::array<bool, 3> is_cell_centred{
                indextype.cellCentered(0), indextype.cellCentered(1),
                indextype.cellCentered(2)};

            const int num_local_components = mfab.local_size();
            std::vector<box_t<int, 3> > grids(num_local_components);
            for (int local_component = 0;
                 local_component < num_local_components; ++local_component) {
              const int component = mfab.IndexArray().at(local_component);
              const amrex::Box &fabbox = mfab.fabbox(component); // exterior
              grids.at(local_component) = box_t<int, 3>{
                lo : {fabbox.smallEnd(0), fabbox.smallEnd(1),
                      fabbox.smallEnd(2)},
                hi : {fabbox.bigEnd(0) + 1, fabbox.bigEnd(1) + 1,
                      fabbox.bigEnd(2) + 1}
              };
            } // for local_component
            levels.at(leveldata.level) = level_t<CCTK_REAL, int, 3>{
              rdomain : rdomain,
              is_cell_centred : is_cell_centred,
              idomain : idomain,
              grids : std::move(grids)
            };
            const level_t<CCTK_REAL, int, 3> &level =
                levels.at(leveldata.level);
            const auto offsets_sizes = level.offsets_sizes();
            const std::vector<int> &offsets = offsets_sizes[0];
            const std::vector<int> &sizes = offsets_sizes[1];

            // Loop over variables
            for (int vi = 0; vi < numvars; ++vi) {
              if (!combine_components) {

                // Loop over components (AMReX boxes)
                for (int local_component = 0;
                     local_component < num_local_components;
                     ++local_component) {
                  const int component = mfab.IndexArray().at(local_component);

                  adios2::Dims lsh(3);
                  for (int d = 0; d < 3; ++d)
                    lsh.at(d) = level.grids.at(local_component).shape()[d];
                  const int np = level.grids.at(local_component).size();
                  const amrex::FArrayBox &fab = mfab[component];

                  const std::string varname = make_varname(
                      gi, vi, patchdata.patch, leveldata.level, component);
                  if (io_verbose)
                    CCTK_VINFO("      Writing variable %s...", varname.c_str());
                  adios2::Variable<CCTK_REAL> var =
                      io.InquireVariable<CCTK_REAL>(varname);
                  assert(var);
                  var.SetSelection({{}, lsh});
                  const CCTK_REAL *const ptr = fab.dataPtr() + vi * np;
                  assert(ptr);
                  engine.Put(var, ptr);
                } // for local_component

              } else { // if combine_components

                const std::string varname =
                    make_varname(gi, vi, patchdata.patch, leveldata.level);
                if (io_verbose)
                  CCTK_VINFO("      Writing variable %s...", varname.c_str());

                const size_t total_np = offsets.back();
                adios2::Variable<CCTK_REAL> var =
                    io.InquireVariable<CCTK_REAL>(varname);
                var.SetSelection({{}, {total_np}});

                if (!combine_via_span) {

                  std::vector<CCTK_REAL> alldata(offsets.back());

                  // Loop over components (AMReX boxes)
                  for (int local_component = 0;
                       local_component < num_local_components;
                       ++local_component) {
                    const int component = mfab.IndexArray().at(local_component);

                    const size_t np = sizes.at(local_component);
                    const amrex::FArrayBox &fab = mfab[component];
                    const CCTK_REAL *const ptr = fab.dataPtr() + vi * np;
                    assert(ptr);

                    std::memcpy(alldata.data() + offsets.at(local_component),
                                ptr, sizeof(CCTK_REAL) * np);

                  } // for local_component

                  engine.Put(var, alldata.data());

                } else { // if combine_via_span

                  const adios2::Variable<CCTK_REAL>::Span span =
                      engine.Put(var);

                  for (int local_component = 0;
                       local_component < num_local_components;
                       ++local_component) {
                    const int component = mfab.IndexArray().at(local_component);
                    const size_t np = sizes.at(local_component);
                    const amrex::FArrayBox &fab = mfab[component];
                    const CCTK_REAL *const ptr = fab.dataPtr() + vi * np;
                    assert(ptr);

                    std::memcpy(span.data() + offsets.at(local_component), ptr,
                                sizeof(CCTK_REAL) * np);
                  } // for local_component

                } // if combine_via_span

              } // if combine_components

            } // for vi
          }
        } // for gi

      } // for leveldata
      assert(!levels.empty());
      const auto rdomain = levels.front().rdomain;
      patches.at(patchdata.patch) = patch_t<CCTK_REAL, int, 3>{
        rdomain : rdomain,
        levels : std::move(levels)
      };
    } // for patchdata
    assert(!patches.empty());
    const grid_structure_t<CCTK_REAL, int, 3>
    grid_structure{patches : std::move(patches)};

    if (io_verbose)
      CCTK_VINFO("  Performing puts...");
    engine.PerformPuts();

    if (io_verbose)
      CCTK_VINFO("  Ending step...");
    engine.EndStep();

  } catch (std::invalid_argument &e) {
    std::cerr << "Invalid argument exception: " << e.what() << "\n";
    CCTK_Abort(nullptr, 1);
  } catch (std::ios_base::failure &e) {
    std::cerr << "IO System base failure exception: " << e.what() << "\n";
    CCTK_Abort(nullptr, 1);
  } catch (std::exception &e) {
    std::cerr << "Exception: " << e.what() << "\n";
    CCTK_Abort(nullptr, 1);
  }

  if (io_verbose)
    CCTK_VINFO("OutputADIOS2 done.");

  if (io_verbose)
    timer.print();
}

carpetx_adios2_t::~carpetx_adios2_t() {
  if (engine) {
    if (io_verbose)
      CCTK_VINFO("ADIOS2: Closing engine...");
    engine.Close();
  }
}

} // namespace CarpetX

#else

namespace CarpetX {
void ShutdownADIOS2() {}
} // namespace CarpetX

#endif // #ifdef HAVE_CAPABILITY_ADIOS2
