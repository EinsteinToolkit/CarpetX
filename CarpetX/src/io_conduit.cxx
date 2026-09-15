#include "io_conduit.hxx"

#include "driver.hxx"
#include "io_meta.hxx"
#include "mpi_types.hxx"
#include "timer.hxx"

#include <CactusBase/IOUtil/src/ioutil_CheckpointRecovery.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef HAVE_CAPABILITY_Conduit

#include <conduit/conduit.hpp>
// #include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_blueprint_mpi.hpp>
#include <conduit/conduit_relay.hpp>

#include <conduit/conduit_relay_io_csv.hpp>

#ifdef CONDUIT_RELAY_IO_HDF5_ENABLED
#include <conduit/conduit_relay_io_hdf5.hpp>
// #include <conduit/conduit_relay_mpi_io_hdf5.hpp>
#endif

#ifdef CONDUIT_RELAY_IO_SILO_ENABLED
#include <conduit/conduit_relay_io_silo.hpp>
#endif

// #include <mpi.h>

#include <algorithm>
#include <array>
#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
using namespace std::experimental;
#else
#include <filesystem>
#endif
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <type_traits>

namespace CarpetX {

namespace {
constexpr bool io_verbose = true;

////////////////////////////////////////////////////////////////////////////////
//
// Convert from/to Conduit data types

template <std::size_t> struct conduit_uint_datatype_id;
template <>
struct conduit_uint_datatype_id<8>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::UINT8_ID> {};
template <>
struct conduit_uint_datatype_id<16>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::UINT16_ID> {};
template <>
struct conduit_uint_datatype_id<32>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::UINT32_ID> {};
template <>
struct conduit_uint_datatype_id<64>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::UINT64_ID> {};
template <std::size_t S>
constexpr conduit::DataType::TypeID conduit_uint_datatype_id_v =
    conduit_uint_datatype_id<S>::value;

template <std::size_t> struct conduit_int_datatype_id;
template <>
struct conduit_int_datatype_id<8>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::INT8_ID> {};
template <>
struct conduit_int_datatype_id<16>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::INT16_ID> {};
template <>
struct conduit_int_datatype_id<32>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::INT32_ID> {};
template <>
struct conduit_int_datatype_id<64>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::INT64_ID> {};
template <std::size_t S>
constexpr conduit::DataType::TypeID conduit_int_datatype_id_v =
    conduit_int_datatype_id<S>::value;

template <std::size_t> struct conduit_float_datatype_id;
template <>
struct conduit_float_datatype_id<32>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::FLOAT32_ID> {};
template <>
struct conduit_float_datatype_id<64>
    : std::integral_constant<conduit::DataType::TypeID,
                             conduit::DataType::FLOAT64_ID> {};
template <std::size_t S>
constexpr conduit::DataType::TypeID conduit_float_datatype_id_v =
    conduit_float_datatype_id<S>::value;

template <typename T>
constexpr std::enable_if_t<(std::is_integral_v<T> && std::is_unsigned_v<T>),
                           conduit::DataType::TypeID>
conduit_datatype_id() {
  return conduit_uint_datatype_id_v<8 * sizeof(T)>;
}

template <typename T>
constexpr std::enable_if_t<(std::is_integral_v<T> && std::is_signed_v<T>),
                           conduit::DataType::TypeID>
conduit_datatype_id() {
  return conduit_int_datatype_id_v<8 * sizeof(T)>;
}

template <typename T>
constexpr
    std::enable_if_t<std::is_floating_point_v<T>, conduit::DataType::TypeID>
    conduit_datatype_id() {
  return conduit_float_datatype_id_v<8 * sizeof(T)>;
}

template <typename T> conduit::DataType conduit_datatype() {
  return conduit::DataType::default_dtype(conduit_datatype_id<T>());
}

////////////////////////////////////////////////////////////////////////////////
//
// Helper functions

std::string clean_varname(std::string varname) {
  varname = regex_replace(varname, std::regex("::"), "-");
  varname = regex_replace(varname, std::regex("\\["), "-");
  varname = regex_replace(varname, std::regex("\\]"), "");
  for (auto &ch : varname)
    ch = tolower(ch);
  return varname;
}

} // namespace

////////////////////////////////////////////////////////////////////////////////
//
// Output

void OutputConduit(const cGH *cctkGH, const std::vector<bool> &output_group,
                   const std::string &output_dir,
                   const std::string &output_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Set up timers
  static Timer timer("OutputConduit");
  Interval interval(timer);

  if (std::count(output_group.begin(), output_group.end(), true) == 0)
    return;

  if (io_verbose)
    CCTK_VINFO("OutputConduit...");

  // Check parameters
  if (CCTK_EQUALS(out_conduit_protocol, "HDF5")) {
#ifndef CONDUIT_RELAY_IO_HDF5_ENABLED
    CCTK_ERROR("Unsupported Conduit output protocol HDF5 is not supported in "
               "this Conduit installation");
#endif
  }
  if (CCTK_EQUALS(out_conduit_protocol, "Silo")) {
#ifndef CONDUIT_RELAY_IO_HDF5_ENABLED
    CCTK_ERROR("Unsupported Conduit output protocol Silo is not supported in "
               "this Conduit installation");
#endif
  }

  const MPI_Comm mpi_comm = amrex::ParallelDescriptor::Communicator();
  const int mpi_tag = 0;
  const int myproc = CCTK_MyProc(nullptr);
  const int nprocs = CCTK_nProcs(nullptr);

  // I/O is performed in groups. In each group, only one process
  // performs I/O, called the "leader". The other processes
  // communicate with the leader via MPI.
  int out_proc_every1;
  if (CCTK_EQUALS(out_mode, "proc"))
    out_proc_every1 = 1;
  else if (CCTK_EQUALS(out_mode, "np"))
    out_proc_every1 = out_proc_every;
  else if (CCTK_EQUALS(out_mode, "onefile"))
    out_proc_every1 = nprocs;
  else
    CCTK_ERROR("internal error");
  // Calculate the I/O group number and its rank within that I/O group
  const auto io_group = [&](const int proc) { return proc / out_proc_every1; };
  const auto io_rank = [&](const int proc) { return proc % out_proc_every1; };
  const int my_io_group = io_group(myproc);
  const int my_io_rank = io_rank(myproc);

  // Create communicators for each I/O process and its group of
  // processes
  MPI_Comm mpi_io_group_comm;
  MPI_Comm_split(mpi_comm, my_io_group, my_io_rank, &mpi_io_group_comm);
  int io_group_myproc, io_group_nprocs;
  MPI_Comm_rank(mpi_io_group_comm, &io_group_myproc);
  MPI_Comm_size(mpi_io_group_comm, &io_group_nprocs);
  assert(io_group_myproc == my_io_rank);
  // The first processe in an I/O group is the "leader", which
  // actually performs the I/O
  const bool is_io_leader = io_group_myproc == 0;

  // Create a communicator for just the I/O leaders
  MPI_Comm mpi_io_leader_comm;
  MPI_Comm_split(mpi_comm, (io_group_myproc == 0 ? 0 : MPI_UNDEFINED), myproc,
                 &mpi_io_leader_comm);
  int io_leader_myproc = -1, io_leader_nprocs = -1;
  if (mpi_io_leader_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(mpi_io_leader_comm, &io_leader_myproc);
    MPI_Comm_size(mpi_io_leader_comm, &io_leader_nprocs);
    assert(my_io_group == io_leader_myproc);
  }
  // The first I/O leader is the meta I/O process
  const bool is_io_meta = io_leader_myproc == 0;

  conduit::Node node_gh;

  if (false) {
    node_gh["tensortypes/scalar3d/dims"] = std::vector<int>{};
    node_gh["tensortypes/scalar3d/stored/value"] = std::vector<int>{};
    node_gh["tensortypes/scalar3d/components/stored"] = "value";
    node_gh["tensortypes/scalar3d/components/parity"] = +1;

    node_gh["tensortypes/vector3d/dims"] = std::vector<int>{3};
    node_gh["tensortypes/vector3d/stored/x"] = std::vector<int>{0};
    node_gh["tensortypes/vector3d/stored/y"] = std::vector<int>{1};
    node_gh["tensortypes/vector3d/stored/z"] = std::vector<int>{2};
    node_gh["tensortypes/vector3d/components/0/stored"] = "x";
    node_gh["tensortypes/vector3d/components/0/parity"] = +1;
    node_gh["tensortypes/vector3d/components/1/stored"] = "y";
    node_gh["tensortypes/vector3d/components/1/parity"] = +1;
    node_gh["tensortypes/vector3d/components/2/stored"] = "z";
    node_gh["tensortypes/vector3d/components/2/parity"] = +1;

    node_gh["tensortypes/symtensor3d/dims"] = std::vector<int>{3, 3};
    node_gh["tensortypes/symtensor3d/stored/xx"] = std::vector<int>{0, 0};
    node_gh["tensortypes/symtensor3d/stored/xy"] = std::vector<int>{0, 1};
    node_gh["tensortypes/symtensor3d/stored/xz"] = std::vector<int>{0, 2};
    node_gh["tensortypes/symtensor3d/stored/yy"] = std::vector<int>{1, 1};
    node_gh["tensortypes/symtensor3d/stored/yz"] = std::vector<int>{1, 2};
    node_gh["tensortypes/symtensor3d/stored/zz"] = std::vector<int>{2, 2};
    node_gh["tensortypes/symtensor3d/components/0/0/stored"] = "xx";
    node_gh["tensortypes/symtensor3d/components/0/0/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/0/1/stored"] = "xy";
    node_gh["tensortypes/symtensor3d/components/0/1/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/0/2/stored"] = "xz";
    node_gh["tensortypes/symtensor3d/components/0/2/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/1/0/stored"] = "xy";
    node_gh["tensortypes/symtensor3d/components/1/0/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/1/1/stored"] = "yy";
    node_gh["tensortypes/symtensor3d/components/1/1/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/1/2/stored"] = "yz";
    node_gh["tensortypes/symtensor3d/components/1/2/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/2/0/stored"] = "xz";
    node_gh["tensortypes/symtensor3d/components/2/0/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/2/1/stored"] = "yz";
    node_gh["tensortypes/symtensor3d/components/2/1/parity"] = +1;
    node_gh["tensortypes/symtensor3d/components/2/2/stored"] = "zz";
    node_gh["tensortypes/symtensor3d/components/2/2/parity"] = +1;

    node_gh["tensortypes/antisymtensor3d/dims"] = std::vector<int>{3, 3};
    node_gh["tensortypes/antisymtensor3d/stored/xy"] = std::vector<int>{0, 1};
    node_gh["tensortypes/antisymtensor3d/stored/xz"] = std::vector<int>{0, 2};
    node_gh["tensortypes/antisymtensor3d/stored/yz"] = std::vector<int>{1, 2};
    node_gh["tensortypes/antisymtensor3d/components/0/0/parity"] = 0;
    node_gh["tensortypes/antisymtensor3d/components/0/1/stored"] = "xy";
    node_gh["tensortypes/antisymtensor3d/components/0/1/parity"] = +1;
    node_gh["tensortypes/antisymtensor3d/components/0/2/stored"] = "xz";
    node_gh["tensortypes/antisymtensor3d/components/0/2/parity"] = +1;
    node_gh["tensortypes/antisymtensor3d/components/1/0/stored"] = "xy";
    node_gh["tensortypes/antisymtensor3d/components/1/0/parity"] = -1;
    node_gh["tensortypes/antisymtensor3d/components/1/1/parity"] = 0;
    node_gh["tensortypes/antisymtensor3d/components/1/2/stored"] = "yz";
    node_gh["tensortypes/antisymtensor3d/components/1/2/parity"] = +1;
    node_gh["tensortypes/antisymtensor3d/components/2/0/stored"] = "xz";
    node_gh["tensortypes/antisymtensor3d/components/2/0/parity"] = -1;
    node_gh["tensortypes/antisymtensor3d/components/2/1/stored"] = "yz";
    node_gh["tensortypes/antisymtensor3d/components/2/1/parity"] = -1;
    node_gh["tensortypes/antisymtensor3d/components/2/2/parity"] = 0;

    node_gh["tensorfields/rho/tensortype"] = "scalar3d";
    node_gh["tensorfields/rho/stored/value"] = "rho";
    node_gh["tensorfields/vel/tensortype"] = "vector3d";
    node_gh["tensorfields/vel/stored/x"] = "velx";
    node_gh["tensorfields/vel/stored/y"] = "vely";
    node_gh["tensorfields/vel/stored/z"] = "velz";
  } // if false

  // Count domains
  using domain_id = std::array<int, 3>; // (patch, level, component)
  std::map<domain_id, int> domain_indices;
  int numdomains = 0;
  // The partition maps are valid only on the meta I/O process
  std::vector<int> partition_map_files;
  std::vector<int> partition_map_domains;
  for (const auto &patchdata : ghext->patchdata) {
    for (const auto &leveldata : patchdata.leveldata) {
      const amrex::DistributionMapping &dm = leveldata.fab->DistributionMap();
      const int numcomponents = dm.size();
      for (int component = 0; component < numcomponents; ++component) {
        const int proc = dm.ProcessorMap().at(component);
        const domain_id domid{patchdata.patch, leveldata.level, component};
        domain_indices[domid] = numdomains;
        if (is_io_meta) {
          partition_map_files.push_back(io_group(proc));
          partition_map_domains.push_back(numdomains);
        }
        ++numdomains;
      }
    }
  }

  const int numgroups = CCTK_NumGroups();

  // The display names valid only on the meta I/O process
  std::map<std::string, std::string> display_names;

  // Loop over patches and levels
  for (const auto &patchdata : ghext->patchdata) {
    for (const auto &leveldata : patchdata.leveldata) {

      const amrex::Geometry &geom = patchdata.amrcore->Geom(leveldata.level);
      const amrex::Real *restrict const geom_x0 = geom.ProbLo();
      const amrex::Real *restrict const geom_dx = geom.CellSize();

      const amrex::DistributionMapping &dm = leveldata.fab->DistributionMap();
      const int numcomponents = dm.size();

      // Loop over local components (AMReX boxes)
      for (int component = 0; component < numcomponents; ++component) {
        const int proc = dm.ProcessorMap().at(component);
        const int component_io_group = io_group(proc);
        const int component_io_rank = io_rank(proc);
        if (component_io_group != my_io_group)
          continue;

        const bool we_own_this_component = io_group_myproc == component_io_rank;
        const bool we_io_this_component = is_io_leader;
        if (!we_own_this_component && !we_io_this_component)
          continue;

        const domain_id domid{patchdata.patch, leveldata.level, component};
        const int domainindex = domain_indices.at(domid);

        // We create a multi-domain mesh, which is just a set of
        // single-domain meshes. Each domain needs to have a different
        // name (`domainname`).
        const std::string domainname = [&]() {
          std::ostringstream buf;
          buf << "domain_" << std::setw(8) << std::setfill('0') << domainindex;
          return buf.str();
        }();
        conduit::Node &node_domain = node_gh[domainname];

        const amrex::Box box = leveldata.fab->box(component);

        const int coordgroup = CCTK_GroupIndex("CoordinatesX::vertex_coords");
        if (coordgroup < 0)
          CCTK_ERROR("Group CoordinatesX::vertex_coords not found. Is the "
                     "thorn CoordinatesX active?");
        assert(coordgroup >= 0);
        const auto &coordgroupdata = *leveldata.groupdata.at(coordgroup);

        // Vertex-centred shape without ghosts
        std::array<int, dim> lsh, lbnd;
        for (int d = 0; d < dim; ++d) {
          lsh[d] = box.length(d) + 1;
          lbnd[d] = box.smallEnd(d);
        }
        int npoints = 1;
        for (int d = 0; d < dim; ++d)
          npoints *= lsh[d];
        std::vector<CCTK_REAL> values(npoints);

        // Define coordinates
        conduit::Node &node_coordsets = node_domain["coordsets"];
        const std::string coordname = "coords";
        conduit::Node &node_coords = node_coordsets[coordname];
        if (patchdata.is_cartesian) {
          // Uniform coordinates
          if (we_io_this_component) {
            node_coords["type"] = "uniform";

            std::array<CCTK_REAL, dim> x0, dx;
            for (int d = 0; d < dim; ++d) {
              dx[d] = geom_dx[d];
              x0[d] = geom_x0[d] + dx[d] * lbnd[d];
            }

            node_coords["dims/i"] = lsh[0];
            node_coords["dims/j"] = lsh[1];
            node_coords["dims/k"] = lsh[2];

            node_coords["origin/x"] = x0[0];
            node_coords["origin/y"] = x0[1];
            node_coords["origin/z"] = x0[2];

            node_coords["spacing/dx"] = dx[0];
            node_coords["spacing/dy"] = dx[1];
            node_coords["spacing/dz"] = dx[2];
          }
        } else {
          // 3D coordinate arrays

          const std::array<std::string, dim> coordnames{"x", "y", "z"};
          for (int dir = 0; dir < dim; ++dir) {

            if (we_own_this_component) {
              constexpr int timelevel = 0;
              const auto &coordmfab = *coordgroupdata.mfab.at(timelevel);
              const auto &coordfab = coordmfab[component];

              // Copy coordinate values without ghosts
              for (int k = 0; k < lsh[2]; ++k) {
                for (int j = 0; j < lsh[1]; ++j) {
                  for (int i = 0; i < lsh[0]; ++i) {
                    const int idx = i + lsh[0] * (j + lsh[1] * k);
                    const amrex::IntVect ivect(lbnd[0] + i, lbnd[1] + j,
                                               lbnd[2] + k);
                    assert(coordfab.box().contains(ivect));
                    values.at(idx) = coordfab(ivect, dir);
                  }
                }
              }
            }
            if (we_own_this_component && !we_io_this_component)
              MPI_Send(values.data(), values.size(),
                       mpi_datatype_v<decltype(values)::value_type>, 0, mpi_tag,
                       mpi_io_group_comm);
            if (!we_own_this_component && we_io_this_component)
              MPI_Recv(values.data(), values.size(),
                       mpi_datatype_v<decltype(values)::value_type>,
                       component_io_rank, mpi_tag, mpi_io_group_comm,
                       MPI_STATUS_IGNORE);

            if (we_io_this_component) {
              node_coords["type"] = "explicit";
              const std::string &coordname = coordnames[dir];
              // Note: We should avoid the double copy (use `set_external`)
              node_coords["values"][coordname].set(values);
            }
          } // for dir
        }

        // Define topologies
        conduit::Node &node_topologies = node_domain["topologies"];
        const std::string toponame = "topo";
        conduit::Node &node_topo = node_topologies[toponame];

        if (we_io_this_component) {
          if (patchdata.is_cartesian) {
            // A uniform topology
            node_topo["type"] = "uniform";
            node_topo["coordset"] = coordname;

            node_topo["elements/origin/i"] = lbnd[0];
            node_topo["elements/origin/j"] = lbnd[1];
            node_topo["elements/origin/k"] = lbnd[2];
          } else {
            // A structured topology
            node_topo["type"] = "structured";
            node_topo["coordset"] = coordname;

            // There one fewer elements (cells) than vertices
            node_topo["elements/dims/i"] = lsh[0] - 1;
            node_topo["elements/dims/j"] = lsh[1] - 1;
            node_topo["elements/dims/k"] = lsh[2] - 1;

            node_topo["elements/dims/i0"] = lbnd[0];
            node_topo["elements/dims/j0"] = lbnd[1];
            node_topo["elements/dims/k0"] = lbnd[2];
          }

          if (patchdata.leveldata.size() > 1) {
            // Define nestsets (AMR hierarchy)
            conduit::Node &node_nestsets = node_domain["nestsets"];
            const std::string nestname = "nest";
            conduit::Node &node_nest = node_nestsets[nestname];
            node_nest["association"] = "vertex";
            node_nest["topology"] = toponame;

            const auto emit_window =
                [&](const GHExt::PatchData::LevelData &cleveldata,
                    const int ccomponent,
                    const GHExt::PatchData::LevelData &fleveldata,
                    const int fcomponent, const bool emit_parent,
                    const bool emit_child) {
                  assert(cleveldata.level + 1 == fleveldata.level);
                  const amrex::Box cbox = cleveldata.fab->box(ccomponent);
                  const amrex::Box fbox = fleveldata.fab->box(fcomponent);

                  // Vertex-centred shapes without ghosts
                  std::array<int, dim> clsh, clbnd;
                  for (int d = 0; d < dim; ++d) {
                    clsh[d] = cbox.length(d) + 1;
                    clbnd[d] = cbox.smallEnd(d);
                  }
                  std::array<int, dim> flsh, flbnd;
                  for (int d = 0; d < dim; ++d) {
                    flsh[d] = fbox.length(d) + 1;
                    flbnd[d] = fbox.smallEnd(d);
                  }

                  // Calculate intersection (on fine grid)
                  // Note: We should use FPInfo and BoxArray, this would
                  // probably be more efficient
                  using std::max, std::min;
                  std::array<int, dim> wmin, wmax;
                  for (int d = 0; d < dim; ++d) {
                    wmin[d] = max(2 * clbnd[d], flbnd[d]);
                    wmax[d] = min(2 * (clbnd[d] + clsh[d]), flbnd[d] + flsh[d]);
                  }

                  // Do the boxes intersect?
                  bool do_intersect = true;
                  for (int d = 0; d < dim; ++d)
                    do_intersect &= wmax[d] >= wmin[d];
                  if (!do_intersect)
                    return;

                  const domain_id cdomid{patchdata.patch, cleveldata.level,
                                         ccomponent};
                  const domain_id fdomid{patchdata.patch, fleveldata.level,
                                         fcomponent};
                  const int cdomainindex = domain_indices.at(cdomid);
                  const int fdomainindex = domain_indices.at(fdomid);
                  const std::string windowname = [&]() {
                    std::ostringstream buf;
                    buf << "window_" << std::setw(8) << std::setfill('0')
                        << cdomainindex << "_" << std::setw(8)
                        << std::setfill('0') << fdomainindex;
                    return buf.str();
                  }();

                  if (emit_parent) {
                    // Parent
                    const std::string cdomainname = [&]() {
                      std::ostringstream buf;
                      buf << "domain_" << std::setw(8) << std::setfill('0')
                          << cdomainindex;
                      return buf.str();
                    }();
                    conduit::Node &node_cdomain = node_gh[cdomainname];
                    conduit::Node &node_cnestsets = node_cdomain["nestsets"];
                    conduit::Node &node_cnest = node_cnestsets[nestname];
                    conduit::Node &node_cwindows = node_cnest["windows"];
                    conduit::Node &node_cwindow = node_cwindows[windowname];

                    node_cwindow["domain_id"] = fdomainindex;
                    node_cwindow["domain_type"] = "parent";
                    node_cwindow["ratio/i"] = 2;
                    node_cwindow["ratio/j"] = 2;
                    node_cwindow["ratio/k"] = 2;
                    node_cwindow["origin/i"] = wmin[0] / 2 - 0 * clbnd[0];
                    node_cwindow["origin/j"] = wmin[1] / 2 - 0 * clbnd[1];
                    node_cwindow["origin/k"] = wmin[2] / 2 - 0 * clbnd[2];
                    node_cwindow["dims/i"] = (wmax[0] - wmin[0] + 1) / 2;
                    node_cwindow["dims/j"] = (wmax[1] - wmin[1] + 1) / 2;
                    node_cwindow["dims/k"] = (wmax[2] - wmin[2] + 1) / 2;
                  }

                  if (emit_child) {
                    // Child
                    const std::string fdomainname = [&]() {
                      std::ostringstream buf;
                      buf << "domain_" << std::setw(8) << std::setfill('0')
                          << fdomainindex;
                      return buf.str();
                    }();
                    conduit::Node &node_fdomain = node_gh[fdomainname];
                    conduit::Node &node_fnestsets = node_fdomain["nestsets"];
                    conduit::Node &node_fnest = node_fnestsets[nestname];
                    conduit::Node &node_fwindows = node_fnest["windows"];
                    conduit::Node &node_fwindow = node_fwindows[windowname];

                    node_fwindow["domain_id"] = cdomainindex;
                    node_fwindow["domain_type"] = "child";
                    node_fwindow["ratio/i"] = 2;
                    node_fwindow["ratio/j"] = 2;
                    node_fwindow["ratio/k"] = 2;
                    node_fwindow["origin/i"] = wmin[0] - 0 * flbnd[0];
                    node_fwindow["origin/j"] = wmin[1] - 0 * flbnd[1];
                    node_fwindow["origin/k"] = wmin[2] - 0 * flbnd[2];
                    node_fwindow["dims/i"] = wmax[0] - wmin[0];
                    node_fwindow["dims/j"] = wmax[1] - wmin[1];
                    node_fwindow["dims/k"] = wmax[2] - wmin[2];
                  }
                };

            // Look for all parents (where we are a child)
            if (leveldata.level > 0) {
              const auto &cleveldata =
                  patchdata.leveldata.at(leveldata.level - 1);
              const amrex::DistributionMapping &cdm =
                  cleveldata.fab->DistributionMap();
              const int cnumcomponents = cdm.size();
              for (int ccomponent = 0; ccomponent < cnumcomponents;
                   ++ccomponent)
                emit_window(cleveldata, ccomponent, leveldata, component, false,
                            true);
            }

            // Look for all children (where we are a parent)
            if (leveldata.level < int(patchdata.leveldata.size()) - 1) {
              const auto &fleveldata =
                  patchdata.leveldata.at(leveldata.level + 1);
              const amrex::DistributionMapping &fdm =
                  fleveldata.fab->DistributionMap();
              const int fnumcomponents = fdm.size();
              for (int fcomponent = 0; fcomponent < fnumcomponents;
                   ++fcomponent)
                emit_window(leveldata, component, fleveldata, fcomponent, true,
                            false);
            }
          } // if numlevels > 1

          // TODO: Define adjacency sets (multi-domain connectivity) for
          // multi-patch systems
        } // if we_io_this_component

        // Define fields (variables)
        conduit::Node &node_fields = node_domain["fields"];

        // Loop over groups
        for (int group = 0; group < numgroups; ++group) {
          if (!output_group.at(group))
            continue;
          cGroup cgroup;
          {
            const int ierr = CCTK_GroupData(group, &cgroup);
            assert(!ierr);
          }
          if (cgroup.grouptype != CCTK_GF)
            continue;
          const auto &groupdata = *leveldata.groupdata.at(group);

          // See above -- we set up only vertex-centred coordinates and
          // topologies
          for (int d = 0; d < dim; ++d)
            assert(groupdata.indextype[d] == 0);

          const int timelevel = 0;

          // Loop over variables
          for (int var = 0; var < groupdata.numvars; ++var) {

            if (we_own_this_component) {
              assert(cgroup.grouptype == CCTK_GF);
              const amrex::MultiFab &mfab = *groupdata.mfab[timelevel];
              const amrex::FArrayBox &fab = mfab[component];
              assert(cgroup.vartype == CCTK_VARIABLE_REAL);

              // Copy values without ghosts
              for (int k = 0; k < lsh[2]; ++k) {
                for (int j = 0; j < lsh[1]; ++j) {
                  for (int i = 0; i < lsh[0]; ++i) {
                    const int idx = i + lsh[0] * (j + lsh[1] * k);
                    const amrex::IntVect ivect(lbnd[0] + i, lbnd[1] + j,
                                               lbnd[2] + k);
                    assert(fab.box().contains(ivect));
                    values.at(idx) = fab(ivect, var);
                  }
                }
              }
            }
            if (we_own_this_component && !we_io_this_component)
              MPI_Send(values.data(), values.size(),
                       mpi_datatype_v<decltype(values)::value_type>, 0, mpi_tag,
                       mpi_io_group_comm);
            if (!we_own_this_component && we_io_this_component)
              MPI_Recv(values.data(), values.size(),
                       mpi_datatype_v<decltype(values)::value_type>,
                       component_io_rank, mpi_tag, mpi_io_group_comm,
                       MPI_STATUS_IGNORE);

            if (we_io_this_component) {
              const std::string varname = clean_varname(
                  CCTK_FullVarName(groupdata.firstvarindex + var));
              conduit::Node &node_var = node_fields[varname];

              // TODO: Check centering
              node_var["topology"] = toponame;
              node_var["association"] = "vertex";

              // Note: We should avoid the double copy (use `set_external`)
              node_var["values"].set(values);

              // Choose a nice name for this field
              if (is_io_meta)
                // TODO: Insert each name only once
                display_names[varname] = varname;
            }
          } // for var
        } // for group

        if (we_io_this_component) {
          conduit::Node &node_state = node_domain["state"];
          node_state["domain_id"] = domainindex;
          node_state["level_id"] = leveldata.level;
          node_state["cycle"] = cctk_iteration;
          node_state["time"] = cctk_time;
        }
      } // for component
    } // for level
  } // for patch

  MPI_Comm_free(&mpi_io_group_comm);

  // There is nothing else to do for regular processes. Only I/O processes
  // continue.
  if (!is_io_leader)
    return;

  // TODO: clean mesh

  for (int p = 0; p < io_leader_nprocs; ++p) {
    MPI_Barrier(mpi_io_leader_comm);
    if (io_leader_myproc == p) {
      CCTK_VINFO("Verifying meshes on process %d (I/O process %d)", myproc,
                 io_group_myproc);
      conduit::Node info;
      const bool isvalid = conduit::blueprint::mesh::verify(node_gh, info);
      if (!isvalid) {
        info.print();
        CCTK_ERROR("node_gh: Invalid blueprint mesh format");
      }
    }
  }

  const std::string filesuffix = [&]() {
    if (CCTK_EQUALS(out_conduit_protocol, "CSV"))
      return "csv";
    else if (CCTK_EQUALS(out_conduit_protocol, "HDF5"))
      return "h5";
    else if (CCTK_EQUALS(out_conduit_protocol, "Silo"))
      return "silo";
    else
      CCTK_VERROR("Unsupported Conduit output protocol \"%s\"",
                  out_conduit_protocol);
  }();

  const std::string metafilesuffix = [&]() {
    if (CCTK_EQUALS(out_conduit_protocol, "HDF5"))
      return "blueprint_root_hdf5";
    else
      return "blueprint_root";
  }();

  const std::string arrayfilename = [&]() {
    std::ostringstream buf;
    buf << output_dir << "/"                                            //
        << output_file                                                  //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
        << ".grid_arrays." << filesuffix;
    return buf.str();
  }();

  const std::string metafilename = [&]() {
    std::ostringstream buf;
    // VisIt recognizes this suffix
    // Should be one of
    // - `.blueprint_root`
    // - `.blueprint_root_hdf5`
    // - `.root`
    buf << output_dir << "/"                                            //
        << output_file                                                  //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
        << "." << metafilesuffix;
    return buf.str();
  }();

  const std::string procsubdirname = [&]() {
    std::ostringstream buf;
    buf << output_file                                                  //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
        << ".blueprint.dir";
    return buf.str();
  }();

  const std::string procfilename = [&]() {
    std::ostringstream buf;
    buf << output_dir << "/" << procsubdirname << "/"                    //
        << output_file                                                   //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration  //
        << ".p" << std::setw(6) << std::setfill('0') << io_leader_myproc //
        << "." << filesuffix;
    return buf.str();
  }();

  const std::string procfilepattern = [&]() {
    std::ostringstream buf;
    // relative to root file
    buf << procsubdirname << "/"                                        //
        << output_file                                                  //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
        << ".p%06d"                                                     //
        << "." << filesuffix;
    return buf.str();
  }();

  const std::string partitionpattern = [&]() {
    std::ostringstream buf;
    buf << procsubdirname << "/"                                        //
        << output_file                                                  //
        << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
        << ".p{:06}"                                                    //
        << "." << filesuffix << ":/"                                    //
        << "domain_{:08}";
    return buf.str();
  }();

  conduit::Node node_meta;
  conduit::Node &node_mesh = node_meta["blueprint_index"];

  conduit::Node &node_index = node_mesh["cctkGH"];
  assert(node_gh.number_of_children() > 0);
  // conduit::blueprint::mesh::generate_index(node_gh[0], "", numdomains,
  //                                          node_index);
  conduit::blueprint::mpi::mesh::generate_index(node_gh, "", node_index,
                                                mpi_io_leader_comm);

  if (is_io_meta) {
    // Add nice display names for the fields
    for (const auto &[path, name] : display_names)
      node_index["fields"][path]["display_name"] = name;

    // Provide a partition map (mapping domains to files)
    node_index["state/partition_pattern"] = partitionpattern;
    std::vector<int> partition_map_file(numdomains);
    std::vector<int> partition_map_domain(numdomains);
    node_index["state/partition_map/file"].set_external(partition_map_files);
    node_index["state/partition_map/domain"].set_external(
        partition_map_domains);
  }

  {
    conduit::Node info;
    const bool isvalid =
        conduit::blueprint::mesh::index::verify(node_index, info);
    if (!isvalid) {
      info.print();
      CCTK_ERROR("node_index: Invalid blueprint mesh index format");
    }
  }

  if (CCTK_EQUALS(out_conduit_protocol, "CSV")) {
    node_meta["protocol/name"] = "csv";
    node_meta["protocol/version"] = CONDUIT_VERSION;
  } else if (CCTK_EQUALS(out_conduit_protocol, "HDF5")) {
    node_meta["protocol/name"] = "hdf5";
    node_meta["protocol/version"] = CONDUIT_VERSION;
  } else if (CCTK_EQUALS(out_conduit_protocol, "Silo")) {
    node_meta["protocol/name"] = "silo";
    node_meta["protocol/version"] = CONDUIT_VERSION;
  } else {
    CCTK_VERROR("Unsupported Conduit output protocol \"%s\"",
                out_conduit_protocol);
  }

  node_meta["number_of_files"] = nprocs;
  node_meta["number_of_trees"] = numdomains;

  node_meta["file_pattern"] = procfilepattern;
  // Note: VisIt only supports up to `%08d`
  node_meta["tree_pattern"] = "domain_%08d";

  // const std::string summaryfilename = [&]() {
  //   std::ostringstream buf;
  //   buf << output_dir << "/" << output_file                   //
  //       << ".it" << std::setw(8) << std::setfill('0') << cctk_iteration //
  //       << ".txt";
  //   return buf.str();
  // }();

  // Scalars and arrays
  conduit::Node node_arrays;
  if (is_io_meta) {
    // Loop over groups
    for (int group = 0; group < numgroups; ++group) {
      if (!output_group.at(group))
        continue;
      cGroup cgroup;
      {
        const int ierr = CCTK_GroupData(group, &cgroup);
        assert(!ierr);
      }
      if (cgroup.grouptype == CCTK_GF)
        continue;
      const auto &groupdata = *ghext->globaldata.arraygroupdata.at(group);

      const std::string domainname = [&]() {
        std::ostringstream buf;
        buf << "group." << clean_varname(CCTK_FullGroupName(group));
        return buf.str();
      }();
      conduit::Node &node_domain = node_arrays[domainname];

      // Define coordinates
      conduit::Node &node_coordsets = node_domain["coordsets"];

      // 1D coordinates for each direction
      node_coordsets["coords/type"] = "uniform";
      // Scalars are not supported -- create a one-element one-dimensional
      // array instead
      if (groupdata.dimension == 0) {
        node_coordsets["coords/dims/i"] = 1;
      } else {
        for (int dir = 0; dir < groupdata.dimension; ++dir) {
          const std::string coordname = [&]() {
            std::ostringstream buf;
            buf << "coords/dims/" << char('i' + dir);
            return buf.str();
          }();
          node_coordsets[coordname] = groupdata.lsh[dir];
        } // for dir
      }

      // Define topologies
      conduit::Node &node_topologies = node_domain["topologies"];
      node_topologies["topo/type"] = "uniform";
      node_topologies["topo/coordset"] = "coords";

      // Define fields (variables)
      conduit::Node &node_fields = node_domain["fields"];

      const int timelevel = 0;

      // Loop over variables
      for (int var = 0; var < groupdata.numvars; ++var) {

        const std::string varname =
            clean_varname(CCTK_FullVarName(groupdata.firstvarindex + var));
        conduit::Node &node_var = node_fields[varname];

        assert(cgroup.grouptype != CCTK_GF);
        assert(cgroup.disttype == CCTK_DISTRIB_CONSTANT);
        assert(cgroup.vartype == CCTK_VARIABLE_REAL);
        const CCTK_REAL *const dataptr = static_cast<const CCTK_REAL *>(
            CCTK_VarDataPtrI(cctkGH, timelevel, groupdata.firstvarindex + var));
        std::size_t numpoints = 1;
        for (int dir = 0; dir < groupdata.dimension; ++dir)
          numpoints *= groupdata.lsh[dir];

        // TODO: Check centering
        node_var["topology"] = "topo";
        node_var["association"] = "vertex";
        // TODO: Avoid ghosts
        node_var["values"].set_external(const_cast<CCTK_REAL *>(dataptr),
                                        numpoints);
      } // for var

      conduit::Node &node_state = node_domain["state"];
      node_state["cycle"] = cctk_iteration;
      node_state["time"] = cctk_time;

      {
        conduit::Node info;
        const bool isvalid =
            conduit::blueprint::verify("mesh", node_arrays[domainname], info);
        if (!isvalid) {
          info.print();
          CCTK_VERROR("node_arrays[%s]: Invalid blueprint mesh format",
                      domainname.c_str());
        }
      }

    } // for group
  }

  {
    if (is_io_meta) {
      const int mode = 0755; // u=rwx g=r-x o=r-x
      static bool did_create_output_directory = false;
      if (!did_create_output_directory) {
        const int ierr = CCTK_CreateDirectory(mode, output_dir.c_str());
        assert(ierr >= 0);
        did_create_output_directory = true;
      }
      const std::string procdirname = output_dir + "/" + procsubdirname;
      const int ierr = CCTK_CreateDirectory(mode, procdirname.c_str());
      assert(ierr >= 0);
    }
    MPI_Barrier(mpi_io_leader_comm);
  }

  if (CCTK_EQUALS(out_conduit_protocol, "HDF5")) {
    conduit::Node hdf5_opts;
    hdf5_opts["chunking/compression/method"] = "gzip";
    hdf5_opts["chunking/compression/level"] = 9;

    conduit::relay::io::hdf5_save(node_gh, procfilename, hdf5_opts);

    if (is_io_meta) {
      conduit::relay::io::hdf5_save(node_arrays, arrayfilename, hdf5_opts);

      conduit::relay::io::hdf5_save(node_meta, metafilename, hdf5_opts);
      // node_index.to_summary_string_stream(summaryfilename,
      // conduit::Node()); node_gh.to_yaml_stream(output_dir + "/output.yaml");
    }
  } else {
    // Don't know yet how to save the other file formats (CSV and Silo)
    CCTK_VERROR("Unsupported Conduit output protocol \"%s\"",
                out_conduit_protocol);
  }

  MPI_Comm_free(&mpi_io_leader_comm);

#if 0
  if (is_io_meta) {
    const int nx = 50, ny = 50;
    const double x_min = -2.0, x_max = 2.0;
    const double y_min = -2.0, y_max = 2.0;
    const double c_re = 0.285, c_im = 0.001;
    const int levels = 3;
    conduit::Node julia;
    conduit::blueprint::mesh::examples::julia_nestsets_complex(
        nx, ny, x_min, x_max, y_min, y_max, c_re, c_im, levels, julia);
    {
      conduit::Node info;
      const bool isvalid = conduit::blueprint::mesh::verify(julia, info);
      if (!isvalid) {
        info.print();
        CCTK_ERROR("julia: Invalid blueprint mesh format");
      }
    }
    conduit::relay::io::hdf5_save(julia, "julia.h5", hdf5_opts);

    conduit::Node root;
    conduit::Node &index = root["blueprint_index/julia_nestset_complex"];
    const int ndomains = conduit::blueprint::mesh::number_of_domains(julia);
    conduit::blueprint::mesh::generate_index(julia["domain_000000"], "",
                                             ndomains, index);
    {
      conduit::Node info;
      const bool isvalid = conduit::blueprint::mesh::index::verify(index, info);
      if (!isvalid) {
        info.print();
        CCTK_ERROR("index: Invalid blueprint mesh index format");
      }
    }
    root["protocol/name"] = "hdf5";
    root["protocol/version"] = CONDUIT_VERSION;
    root["number_of_files"] = 1;
    root["number_of_trees"] = ndomains;
    root["file_pattern"] = "julia.h5";
    root["tree_pattern"] = "domain_%06d";
    conduit::relay::io::hdf5_save(root, "root.h5", hdf5_opts);
  }
#endif
}

} // namespace CarpetX

#endif
