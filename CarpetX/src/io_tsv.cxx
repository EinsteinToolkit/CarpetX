#include "io_tsv.hxx"

#include "driver.hxx"
#include "mpi_types.hxx"
#include "timer.hxx"

#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Misc.h>
#include <util_Network.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {

// Which information to output in the header
enum out_fileinfo_t {
  fileinfo_none = 0x0,
  fileinfo_create_date = 0x1,
  fileinfo_parameter_filename = 0x2,
  fileinfo_axis_labels = 0x4
};
static out_fileinfo_t get_out_fileinfo() {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_EQUALS(out_fileinfo, "none"))
    return fileinfo_none;
  if (CCTK_EQUALS(out_fileinfo, "creation date"))
    return fileinfo_create_date;
  if (CCTK_EQUALS(out_fileinfo, "parameter filename"))
    return fileinfo_parameter_filename;
  if (CCTK_EQUALS(out_fileinfo, "axis labels"))
    return fileinfo_axis_labels;
  if (CCTK_EQUALS(out_fileinfo, "all"))
    return out_fileinfo_t(fileinfo_create_date | fileinfo_parameter_filename |
                          fileinfo_axis_labels);
  CCTK_ERROR("Internal error");
}

// Collect the information we may want to put into headers
struct fileinfo_t {
  std::string current_date, current_time;
  std::string hostname;
  std::string parameter_filename;
};
static fileinfo_t get_fileinfo() {
  fileinfo_t fileinfo;

  char current_date[1000];
  Util_CurrentDate(sizeof current_date - 1, current_date);
  fileinfo.current_date = current_date;

  char current_time[1000];
  Util_CurrentTime(sizeof current_time - 1, current_time);
  fileinfo.current_time = current_time;

  char hostname[1000];
  Util_GetHostName(hostname, sizeof hostname - 1);
  fileinfo.hostname = hostname;

  char parameter_filename[10000];
  CCTK_ParameterFilename(sizeof parameter_filename - 1, parameter_filename);
  fileinfo.parameter_filename = parameter_filename;

  return fileinfo;
}

// How to output the headers
enum out_header_t { header_none, header_comment, header_plain };
static out_header_t get_out_header() {
  DECLARE_CCTK_PARAMETERS;
  if (CCTK_EQUALS(out_tsv_header, "none"))
    return header_none;
  if (CCTK_EQUALS(out_tsv_header, "comment"))
    return header_comment;
  if (CCTK_EQUALS(out_tsv_header, "plain"))
    return header_plain;
  CCTK_ERROR("Internal error");
}

void WriteTSVold(const cGH *restrict cctkGH, const std::string &filename,
                 int gi, const std::vector<std::string> &varnames,
                 const out_fileinfo_t out_fileinfo,
                 const out_header_t out_header) {
  std::ostringstream buf;
  buf << filename << ".tsv";
  std::ofstream file(buf.str());
  const std::string sep = "\t";

  // get more precision for floats, could also use
  // https://stackoverflow.com/a/30968371
  file << std::setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1)
       << std::scientific;

  // Output header
  switch (out_header) {
  case header_none:
    // no header
    break;
  case header_comment: {
    fileinfo_t fileinfo = get_fileinfo();
    if (out_fileinfo & fileinfo_create_date) {
      file << "# Date: " << fileinfo.current_date
           << "   time: " << fileinfo.current_time
           << "   hostname: " << fileinfo.hostname << "\n";
    }
    if (out_fileinfo & fileinfo_parameter_filename) {
      file << "# parameter file: \"" << fileinfo.parameter_filename << "\"\n";
    }
    if (out_fileinfo & fileinfo_axis_labels) {
      file << "# 1:iteration" << sep << "2:time" << sep << "3:patch" << sep
           << "4:level" << sep << "5:component" << sep << "6:i" << sep << "7:j"
           << sep << "8:k" << sep << "9:x" << sep << "10:y" << sep << "11:z";
      int col = 12;
      for (const auto &varname : varnames)
        file << sep << col++ << ":" << varname;
      file << "\n";
    }
    break;
  }
  case header_plain: {
    // We can only output axis labels when using the "plain" format; everything
    // else would break TSV "standard" (and confuse TSV readers)
    if (out_fileinfo & fileinfo_axis_labels) {
      file << "iteration" << sep << "time" << sep << "patch" << sep << "level"
           << sep << "component" << sep << "i" << sep << "j" << sep << "k"
           << sep << "x" << sep << "y" << sep << "z";
      for (const auto &varname : varnames)
        file << sep << varname;
      file << "\n";
    }
    break;
  }
  }

  for (const auto &patchdata : ghext->patchdata) {
    for (const auto &leveldata : patchdata.leveldata) {
      const auto &groupdata = *leveldata.groupdata.at(gi);
      const int tl = 0;
      const auto &geom = patchdata.amrcore->Geom(leveldata.level);
      const auto &mfab = *groupdata.mfab.at(tl);
      for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const amrex::Array4<const CCTK_REAL> &vars = mfab.array(mfi);
        const auto &imin = vars.begin;
        const auto &imax = vars.end;
        for (int k = imin.z; k < imax.z; ++k) {
          for (int j = imin.y; j < imax.y; ++j) {
            for (int i = imin.x; i < imax.x; ++i) {
              const std::array<int, dim> I{i, j, k};
              std::array<CCTK_REAL, dim> x;
              for (int d = 0; d < dim; ++d)
                x[d] = geom.ProbLo(d) +
                       (I[d] + CCTK_REAL(0.5) * groupdata.indextype[d]) *
                           geom.CellSize(d);
              file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time << sep
                   << patchdata.patch << sep << leveldata.level << sep
                   << mfi.index() << sep << I[0] << sep << I[1] << sep << I[2]
                   << sep << x[0] << sep << x[1] << sep << x[2];
              for (int n = 0; n < groupdata.numvars; ++n)
                file << sep << vars(i, j, k, n);
              file << "\n";
            }
          }
        }
      }
    }
  }

  file.close();
}

void OutputTSVold(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!out_tsv)
    return;

  static Timer timer("OutputTSV");
  Interval interval(timer);

  const out_fileinfo_t out_fileinfo_val = get_out_fileinfo();
  const out_header_t out_header = get_out_header();

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata0 =
        *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(gi);
    if (groupdata0.mfab.size() > 0) {
      const int tl = 0;

      std::string groupname = CCTK_FullGroupName(gi);
      groupname = regex_replace(groupname, std::regex("::"), "-");
      for (auto &c : groupname)
        c = tolower(c);
      std::ostringstream buf;
      buf << out_dir << "/" << groupname;
      buf << ".it" << std::setw(6) << std::setfill('0') << cctk_iteration;
      buf << ".p" << std::setw(4) << std::setfill('0') << CCTK_MyProc(nullptr);
      const std::string filename = buf.str();

      amrex::Vector<std::string> varnames(groupdata0.numvars);
      for (int vi = 0; vi < groupdata0.numvars; ++vi) {
        std::ostringstream buf;
        buf << CCTK_VarName(groupdata0.firstvarindex + vi);
        for (int i = 0; i < tl; ++i)
          buf << "_p";
        varnames.at(vi) = buf.str();
      }

      WriteTSVold(cctkGH, filename, gi, varnames, out_fileinfo_val, out_header);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

void WriteTSVScalars(const cGH *restrict cctkGH, const std::string &filename,
                     const int gi, const bool truncate_files,
                     const out_fileinfo_t out_fileinfo,
                     const out_header_t out_header) {
  // Output only on root process
  if (CCTK_MyProc(nullptr) > 0)
    return;

  const auto &arraygroupdata = *ghext->globaldata.arraygroupdata.at(gi);
  cGroup cgroup;
  const int ierr = CCTK_GroupData(gi, &cgroup);
  assert(!ierr);

  std::vector<std::string> varnames;
  for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
    varnames.push_back(CCTK_VarName(arraygroupdata.firstvarindex + vi));

  const std::string sep = "\t";
  static std::set<std::string> did_write_file;
  // Truncate only if (1) we should, and (2) this is the first time we're
  // writing to that file
  const bool truncate_this_file =
      truncate_files && did_write_file.count(filename) == 0;
  did_write_file.insert(filename);
  const std::ios::openmode mode =
      truncate_this_file ? std::ios::out : std::ios::app;
  std::ofstream file(filename, mode);
  // get more precision for floats, could also use
  // https://stackoverflow.com/a/30968371
  file << std::setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1)
       << std::scientific;

  if (file.tellp() == 0) {
    // Output header if the file is empty, i.e. the first time we write to this
    // file
    switch (out_header) {
    case header_none:
      // no header
      break;
    case header_comment: {
      fileinfo_t fileinfo = get_fileinfo();
      if (out_fileinfo & fileinfo_create_date) {
        file << "# Date: " << fileinfo.current_date
             << "   time: " << fileinfo.current_time
             << "   hostname: " << fileinfo.hostname << "\n";
      }
      if (out_fileinfo & fileinfo_parameter_filename) {
        file << "# parameter file: \"" << fileinfo.parameter_filename << "\"\n";
      }
      if (out_fileinfo & fileinfo_axis_labels) {
        file << "# 1:iteration" << sep << "2:time";
        int col = 3;
        for (const auto &varname : varnames)
          if (cgroup.vartype == CCTK_VARIABLE_REAL ||
              cgroup.vartype == CCTK_VARIABLE_INT) {
            file << sep << col++ << ":" << varname;
          } else if (cgroup.vartype == CCTK_VARIABLE_COMPLEX) {
            file << sep << col++ << ":" << varname << ".real";
            file << sep << col++ << ":" << varname << ".imag";
          } else {
            assert(0 && "Unexpected variable type");
          }
        file << "\n";
      }
      break;
    }
    case header_plain: {
      // We can only output axis labels when using the "plain" format;
      // everything else would break TSV "standard" (and confuse TSV readers)
      if (out_fileinfo & fileinfo_axis_labels) {
        file << "iteration" << sep << "time";
        for (const auto &varname : varnames)
          if (cgroup.vartype == CCTK_VARIABLE_REAL ||
              cgroup.vartype == CCTK_VARIABLE_INT) {
            file << sep << varname;
          } else if (cgroup.vartype == CCTK_VARIABLE_COMPLEX) {
            file << sep << varname << "_real";
            file << sep << varname << "_imag";
          } else {
            assert(0 && "Unexpected variable type");
          }
        file << "\n";
      }
      break;
    }
    }
  }

  // Output data
  file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
  const int tl = 0;
  const auto &data = arraygroupdata.data.at(tl);
  for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
    file << sep << data[vi];
  file << "\n";
}

void WriteTSVArrays(const cGH *restrict cctkGH, const std::string &filename,
                    const int gi, const int out_dir, const bool truncate_files,
                    const out_fileinfo_t out_fileinfo,
                    const out_header_t out_header) {
  // Output only on root process
  if (CCTK_MyProc(nullptr) > 0)
    return;

  const auto &arraygroupdata = *ghext->globaldata.arraygroupdata.at(gi);
  cGroup cgroup;
  const int ierr = CCTK_GroupData(gi, &cgroup);
  assert(!ierr);

  if (out_dir >= arraygroupdata.dimension)
    return;

  std::vector<std::string> varnames;
  for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
    varnames.push_back(CCTK_VarName(arraygroupdata.firstvarindex + vi));

  const std::string sep = "\t";
  static std::set<std::string> did_write_file;
  // Truncate only if (1) we should, and (2) this is the first time we're
  // writing to that file
  const bool truncate_this_file =
      truncate_files && did_write_file.count(filename) == 0;
  did_write_file.insert(filename);
  const std::ios::openmode mode =
      truncate_this_file ? std::ios::out : std::ios::app;
  std::ofstream file(filename, mode);
  // get more precision for floats, could also use
  // https://stackoverflow.com/a/30968371
  file << std::setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1)
       << std::scientific;

  if (file.tellp() == 0) {
    // Output header if the file is empty, i.e. the first time we write to this
    // file
    switch (out_header) {
    case header_none:
      // no header
      break;
    case header_comment: {
      fileinfo_t fileinfo = get_fileinfo();
      if (out_fileinfo & fileinfo_create_date) {
        file << "# Date: " << fileinfo.current_date
             << "   time: " << fileinfo.current_time
             << "   hostname: " << fileinfo.hostname << "\n";
      }
      if (out_fileinfo & fileinfo_parameter_filename) {
        file << "# parameter file: \"" << fileinfo.parameter_filename << "\"\n";
      }
      if (out_fileinfo & fileinfo_axis_labels) {
        int col = 1;
        file << "# " << col++ << ":iteration";
        file << sep << col++ << ":time";
        for (int dir = 0; dir < arraygroupdata.dimension; ++dir)
          file << sep << col++ << ":"
               << "ijk"[dir];
        for (const auto &varname : varnames)
          if (cgroup.vartype == CCTK_VARIABLE_REAL ||
              cgroup.vartype == CCTK_VARIABLE_INT) {
            file << sep << col++ << ":" << varname;
          } else if (cgroup.vartype == CCTK_VARIABLE_COMPLEX) {
            file << sep << col++ << ":" << varname << ".real";
            file << sep << col++ << ":" << varname << ".imag";
          } else {
            assert(0 && "Unexpected variable type");
          }
        file << "\n";
      }
      break;
    }
    case header_plain: {
      // We can only output axis labels when using the "plain" format;
      // everything else would break TSV "standard" (and confuse TSV readers)
      if (out_fileinfo & fileinfo_axis_labels) {
        file << "iteration";
        file << sep << "time";
        for (int dir = 0; dir < arraygroupdata.dimension; ++dir)
          file << sep << "ijk"[dir];
        for (const auto &varname : varnames)
          if (cgroup.vartype == CCTK_VARIABLE_REAL ||
              cgroup.vartype == CCTK_VARIABLE_INT) {
            file << sep << varname;
          } else if (cgroup.vartype == CCTK_VARIABLE_COMPLEX) {
            file << sep << varname << "_real";
            file << sep << varname << "_imag";
          } else {
            assert(0 && "Unexpected variable type");
          }
        file << "\n";
      }
      break;
    }
    }
  }

  constexpr int di = 1;
  const int dj = di * arraygroupdata.lsh[0];
  const int dk = dj * arraygroupdata.lsh[1];
  const int np = dk * arraygroupdata.lsh[2];
  const std::array<int, dim> DI{di, dj, dk};

  // Output data
  for (int i = 0; i < arraygroupdata.lsh[out_dir]; ++i) {
    file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
    for (int dir = 0; dir < arraygroupdata.dimension; ++dir)
      file << sep << (dir == out_dir ? i : 0);
    const int tl = 0;
    const auto &data = arraygroupdata.data.at(tl);
    for (int vi = 0; vi < arraygroupdata.numvars; ++vi)
      file << sep << data[np * vi + DI[out_dir] * i];
    file << "\n";
  }
}

void WriteTSVGFs(const cGH *restrict cctkGH, const std::string &filename,
                 const int gi, const vect<bool, dim> &outdirs,
                 const vect<CCTK_REAL, dim> &outcoords,
                 const bool truncate_files, const out_fileinfo_t out_fileinfo,
                 const out_header_t out_header,
                 const bool output_boundary_points) {
  const auto &groupdata0 =
      *ghext->patchdata.at(0).leveldata.at(0).groupdata.at(gi);

  // Number of values transmitted per grid point
  const int nintvalues = 1                  // patch
                         + 1                // level
                         + dim              // grid point index
                         + 1;               // isghost
  const int nvalues = nintvalues            // integer values
                      + dim                 // coordinates
                      + groupdata0.numvars; // grid function values

  // Data transmitted from this process
  std::vector<CCTK_REAL> data;
  data.reserve(10000);
  for (const auto &patchdata : ghext->patchdata) {
    for (const auto &leveldata : patchdata.leveldata) {
      const auto &groupdata = *leveldata.groupdata.at(gi);

      const int tl = 0;

      // Convert a (direction, face) pair to an AMReX Orientation
      const auto orient = [&](int d, int f) {
        return amrex::Orientation(d, amrex::Orientation::Side(f));
      };

      const auto &geom = patchdata.amrcore->Geom(leveldata.level);
      const amrex::Box &domain =
          patchdata.amrcore->Geom(leveldata.level).Domain();
      vect<CCTK_REAL, dim> x0, dx;
      for (int d = 0; d < dim; ++d) {
        dx[d] = geom.CellSize(d);
        x0[d] = geom.ProbLo(d) + groupdata.indextype[d] * dx[d] / 2;
      }
      vect<int, dim> icoord;
      for (int d = 0; d < dim; ++d)
        if (outdirs[d])
          icoord[d] = INT_MIN;
        else
          icoord[d] = lrint((outcoords[d] - x0[d]) / dx[d]);

      const auto &symmetries = ghext->patchdata.at(leveldata.patch).symmetries;

      const auto &mfab = *groupdata.mfab.at(tl);

      const vect<int, dim> nghosts = {mfab.nGrow(0), mfab.nGrow(1),
                                      mfab.nGrow(2)};

      for (amrex::MFIter mfi(mfab); mfi.isValid(); ++mfi) {
        const amrex::Array4<const CCTK_REAL> &vars = mfab.array(mfi);

        const amrex::Box &vbx =
            mfi.validbox(); // interior region (without ghosts)
        vect<vect<bool, dim>, 2> bbox;
        for (int d = 0; d < dim; ++d)
          for (int f = 0; f < 2; ++f)
            bbox[f][d] = vbx[orient(d, f)] ==
                             domain[orient(d, f)] +
                                 (groupdata.indextype[d] == 0 && f == 1) &&
                         symmetries[f][d] != symmetry_t::none;

        const vect<int, dim> varmin = {vars.begin.x, vars.begin.y,
                                       vars.begin.z};
        const vect<int, dim> varmax = {vars.end.x, vars.end.y, vars.end.z};

        // Skip ghost points but keep boundary points
        const vect<int, dim> intmin = varmin + !bbox[0] * nghosts;
        const vect<int, dim> intmax = varmax - !bbox[1] * nghosts;

        bool output_something = true;
        vect<int, dim> imin, imax;
        for (int d = 0; d < dim; ++d) {
          if (outdirs[d]) {
            // output everything
            imin[d] = varmin[d] + !output_boundary_points * nghosts[d];
            imax[d] = varmax[d] - !output_boundary_points * nghosts[d];
          } else if (icoord[d] >= varmin[d] && icoord[d] < varmax[d]) {
            // output one point
            imin[d] = icoord[d];
            imax[d] = icoord[d] + 1;
          } else {
            // output nothing
            output_something = false;
          }
        }

        if (output_something) {
          for (int k = imin[2]; k < imax[2]; ++k) {
            for (int j = imin[1]; j < imax[1]; ++j) {
              for (int i = imin[0]; i < imax[0]; ++i) {
                const vect<int, dim> I{i, j, k};
                const auto old_size = data.size();
                data.push_back(patchdata.patch);
                data.push_back(leveldata.level);
                for (int d = 0; d < dim; ++d)
                  data.push_back(I[d]);
                const bool isghost = any(I < intmin || I >= intmax);
                data.push_back(isghost);
                for (int d = 0; d < dim; ++d)
                  data.push_back(x0[d] + I[d] * dx[d]);
                for (int vi = 0; vi < groupdata.numvars; ++vi)
                  data.push_back(vars(i, j, k, vi));
                assert(data.size() == old_size + nvalues);
              }
            }
          }
        } // if output_something
      } // for mfi
    } // for leveldata
  } // for patchdata
  assert(data.size() % nvalues == 0);

  const MPI_Comm comm = amrex::ParallelDescriptor::Communicator();
  const int myproc = amrex::ParallelDescriptor::MyProc();
  const int nprocs = amrex::ParallelDescriptor::NProcs();
  const int ioproc = 0;

  assert(data.size() <= INT_MAX);
  const int npoints = data.size();

  std::vector<int> all_npoints;
  if (myproc == ioproc)
    all_npoints.resize(nprocs);
  MPI_Gather(&npoints, 1, MPI_INT, all_npoints.data(), 1, MPI_INT, ioproc,
             comm);

  int total_npoints = 0;
  std::vector<int> all_offsets;
  if (myproc == ioproc) {
    all_offsets.resize(nprocs);
    for (int p = 0; p < nprocs; ++p) {
      all_offsets.at(p) = total_npoints;
      assert(total_npoints <= INT_MAX - all_npoints.at(p));
      total_npoints += all_npoints.at(p);
    }
  }
  std::vector<CCTK_REAL> all_data;
  if (myproc == ioproc)
    all_data.resize(total_npoints);
  MPI_Gatherv(data.data(), npoints, mpi_datatype<CCTK_REAL>::value,
              all_data.data(), all_npoints.data(), all_offsets.data(),
              mpi_datatype<CCTK_REAL>::value, ioproc, comm);

  if (myproc == ioproc) {
    assert(total_npoints % nvalues == 0);
    std::vector<int> iptr(total_npoints / nvalues);
    iota(iptr.begin(), iptr.end(), 0);
    const auto compare_eq = [&](const int i, const int j) {
      std::array<int, nintvalues> pi, pj;
      for (int d = 0; d < nintvalues; ++d)
        pi[d] = int(all_data.at(i * nvalues + d));
      // Ignore `isghost` field
      pi[dim + 2] = 0;
      for (int d = 0; d < nintvalues; ++d)
        pj[d] = int(all_data.at(j * nvalues + d));
      // Ignore `isghost` field
      pj[dim + 2] = 0;
      return pi == pj;
    };
    const auto compare_lt = [&](const int i, const int j) {
      std::array<int, nintvalues> pi, pj;
      for (int d = 0; d < nintvalues; ++d)
        pi[d] = int(all_data.at(i * nvalues + d));
      for (int d = 0; d < nintvalues; ++d)
        pj[d] = int(all_data.at(j * nvalues + d));
      return pi < pj;
    };
    sort(iptr.begin(), iptr.end(), compare_lt);
    const auto last = unique(iptr.begin(), iptr.end(), compare_eq);
    iptr.erase(last, iptr.end());

    std::vector<std::string> varnames;
    for (int vi = 0; vi < groupdata0.numvars; ++vi)
      varnames.push_back(CCTK_VarName(groupdata0.firstvarindex + vi));

    const std::string sep = "\t";
    static std::set<std::string> did_write_file;
    // Truncate only if (1) we should, and (2) this is the first time we're
    // writing to that file
    const bool truncate_this_file =
        truncate_files && did_write_file.count(filename) == 0;
    did_write_file.insert(filename);
    const std::ios::openmode mode =
        truncate_this_file ? std::ios::out : std::ios::app;
    std::ofstream file(filename, mode);
    // get more precision for floats, could also use
    // https://stackoverflow.com/a/30968371
    file << std::setprecision(std::numeric_limits<CCTK_REAL>::digits10 + 1)
         << std::scientific;

    if (file.tellp() == 0) {
      // Output header if the file is empty, i.e. the first time we write to
      // this file
      switch (out_header) {
      case header_none:
        // no header
        break;
      case header_comment: {
        int col = 0;
        file << "# " << ++col << ":iteration";
        file << sep << ++col << ":time";
        file << sep << ++col << ":patch";
        file << sep << ++col << ":level";
        for (int d = 0; d < dim; ++d)
          file << sep << ++col << ":"
               << "ijk"[d];
        for (int d = 0; d < dim; ++d)
          file << sep << ++col << ":"
               << "xyz"[d];
        for (const auto &varname : varnames)
          file << sep << ++col << ":" << varname;
        file << "\n";
        break;
      }
      case header_plain: {
        file << "iteration";
        file << sep << "time";
        file << sep << "patch";
        file << sep << "level";
        for (int d = 0; d < dim; ++d)
          file << sep << "ijk"[d];
        for (int d = 0; d < dim; ++d)
          file << sep << "xyz"[d];
        for (const auto &varname : varnames)
          file << sep << varname;
        file << "\n";
        break;
      }
      }
    }

    // Output data
    for (const auto i : iptr) {
      int pos = nvalues * i;
      file << cctkGH->cctk_iteration << sep << cctkGH->cctk_time;
      for (int v = 0; v < nintvalues; ++v)
        if (v != dim + 2) // skip `isghost` marker
          file << sep << int(all_data.at(pos++));
        else
          pos++;
      for (int v = nintvalues; v < nvalues; ++v)
        file << sep << all_data.at(pos++);
      file << "\n";
      assert(pos % nvalues == 0);
    }
  }
}

void OutputTSV(const cGH *restrict cctkGH) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_tsv_vars[0] == '\0')
    return;

  static Timer timer("OutputTSVUni");
  Interval interval(timer);

  const out_fileinfo_t out_fileinfo_val = get_out_fileinfo();
  const out_header_t out_header = get_out_header();
  const bool combine_iterations = out_tsv_combine_iterations;
  const bool output_boundary_points = out_tsv_output_boundary_points;

  // `IO_TruncateOutputFiles` says whether files should be truncated
  // at the beginning of the simulation. Typically, we should truncate
  // if we do not recover, and vice versa. This is the decision that
  // `IO_TruncateOutputFiles` makes. We still need to figure out
  // ourselves whether this is the beginning of a simulation.
  const bool truncate_files = IO_TruncateOutputFiles(cctkGH);

  // Find output groups
  const std::vector<bool> group_enabled = [&] {
    std::vector<bool> enabled(CCTK_NumGroups(), false);
    const auto callback{
        [](const int index, const char *const optstring, void *const arg) {
          std::vector<bool> &enabled = *static_cast<std::vector<bool> *>(arg);
          enabled.at(CCTK_GroupIndexFromVarI(index)) = true;
        }};
    CCTK_TraverseString(out_tsv_vars, callback, &enabled, CCTK_GROUP_OR_VAR);
    if (verbose) {
      CCTK_VINFO("TSV output for groups:");
      for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
        if (group_enabled.at(gi))
          CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
    }
    return enabled;
  }();
  const auto num_out_groups =
      count(group_enabled.begin(), group_enabled.end(), true);
  if (num_out_groups == 0)
    return;

  const int numgroups = CCTK_NumGroups();
  for (int gi = 0; gi < numgroups; ++gi) {
    if (group_enabled.at(gi)) {
      std::string groupname = CCTK_FullGroupName(gi);
      groupname = regex_replace(groupname, std::regex("::"), "-");
      for (auto &ch : groupname)
        ch = tolower(ch);
      std::ostringstream buf;
      buf << out_dir << "/" << groupname;
      if (!combine_iterations)
        buf << ".it" << std::setw(6) << std::setfill('0') << cctk_iteration;
      const std::string basename = buf.str();
      switch (CCTK_GroupTypeI(gi)) {
      case CCTK_SCALAR:
        WriteTSVScalars(cctkGH, basename + ".tsv", gi, truncate_files,
                        out_fileinfo_val, out_header);
        break;
      case CCTK_ARRAY:
        WriteTSVArrays(cctkGH, basename + ".x.tsv", gi, 0, truncate_files,
                       out_fileinfo_val, out_header);
        WriteTSVArrays(cctkGH, basename + ".y.tsv", gi, 1, truncate_files,
                       out_fileinfo_val, out_header);
        WriteTSVArrays(cctkGH, basename + ".z.tsv", gi, 2, truncate_files,
                       out_fileinfo_val, out_header);
        break;
      case CCTK_GF:
        WriteTSVGFs(cctkGH, basename + ".x.tsv", gi, {true, false, false},
                    {0, out_xline_y, out_xline_z}, truncate_files,
                    out_fileinfo_val, out_header, output_boundary_points);
        WriteTSVGFs(cctkGH, basename + ".y.tsv", gi, {false, true, false},
                    {out_yline_x, 0, out_yline_z}, truncate_files,
                    out_fileinfo_val, out_header, output_boundary_points);
        WriteTSVGFs(cctkGH, basename + ".z.tsv", gi, {false, false, true},
                    {out_zline_x, out_zline_y, 0}, truncate_files,
                    out_fileinfo_val, out_header, output_boundary_points);
        break;
      default:
        assert(0);
      }
    }
  }
}

} // namespace CarpetX
