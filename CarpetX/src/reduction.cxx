#include "driver.hxx"
#include "mpi_types.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <cctk_Parameters.h>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Orientation.H>

#include <bitset>
#include <cstdlib>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace std;
template <typename T, int D> MPI_Datatype reduction_mpi_datatype() {
  static MPI_Datatype datatype = MPI_DATATYPE_NULL;
  if (datatype == MPI_DATATYPE_NULL) {
    MPI_Type_contiguous(sizeof(reduction<T, D>) / sizeof(T),
                        mpi_datatype<T>::value, &datatype);
    char name[MPI_MAX_OBJECT_NAME];
    int namelen;
    MPI_Type_get_name(mpi_datatype<T>::value, name, &namelen);
    ostringstream buf;
    buf << "reduction<" << name << "," << D << ">";
    string newname = buf.str();
    MPI_Type_set_name(datatype, newname.c_str());
    MPI_Type_commit(&datatype);
  }
  return datatype;
}

namespace {
template <typename T>
void mpi_reduce_typed(const void *restrict x0, void *restrict y0,
                      const int *restrict length) {
  const T *restrict x = static_cast<const T *>(x0);
  T *restrict y = static_cast<T *>(y0);
#pragma omp simd
  for (int i = 0; i < *length; ++i)
    y[i] += x[i];
}

void mpi_reduce(void *restrict x, void *restrict y, int *restrict length,
                MPI_Datatype *restrict datatype) {
  // Analyze MPI datatype
  int num_integers, num_addresses, num_datatypes, combiner;
  MPI_Type_get_envelope(*datatype, &num_integers, &num_addresses,
                        &num_datatypes, &combiner);
  assert(combiner == MPI_COMBINER_CONTIGUOUS);
  assert(num_integers == 1);
  assert(num_datatypes == 1);
  vector<int> integers(num_integers);
  vector<MPI_Aint> addresses(num_addresses);
  vector<MPI_Datatype> datatypes(num_datatypes);
  MPI_Type_get_contents(*datatype, num_integers, num_addresses, num_datatypes,
                        integers.data(), addresses.data(), datatypes.data());
  MPI_Datatype inner_datatype = datatypes.at(0);
  if (inner_datatype == MPI_FLOAT)
    return mpi_reduce_typed<reduction<float, dim> >(x, y, length);
  if (inner_datatype == MPI_DOUBLE)
    return mpi_reduce_typed<reduction<double, dim> >(x, y, length);
  if (inner_datatype == MPI_LONG_DOUBLE)
    return mpi_reduce_typed<reduction<long double, dim> >(x, y, length);
  CCTK_ERROR("Unsupported MPI datatype");
}
} // namespace

MPI_Op reduction_mpi_op() {
  static MPI_Op op = MPI_OP_NULL;
  if (op == MPI_OP_NULL)
    MPI_Op_create(mpi_reduce, 1 /*commutes*/, &op);
  return op;
}

////////////////////////////////////////////////////////////////////////////////

namespace {
template <typename T>
reduction<T, dim>
reduce_array(const amrex::Array4<const T> &restrict vars, const int n,
             const vect<int, dim> &tmin, const vect<int, dim> &tmax,
             const vect<int, dim> &indextype, const vect<int, dim> &imin,
             const vect<int, dim> &imax,
             const amrex::Array4<const int> *restrict const finemask,
             const vect<T, dim> &x0, const vect<T, dim> &dx) {
  constexpr vect<vect<int, dim>, dim> di = {vect<int, dim>::unit(0),
                                            vect<int, dim>::unit(1),
                                            vect<int, dim>::unit(2)};
  constexpr vect<vect<vect<int, dim>, dim>, 2> dirs = {-di, +di};

  constexpr vect<vect<std::bitset<8>, dim>, 2> faces = {
      {0b01010101, 0b00110011, 0b00001111},
      {0b10101010, 0b11001100, 0b11110000}};

  const vect<vect<int, dim>, 2> ibnd = {imin, imax - 1};

  const auto masked = [&](const vect<int, dim> &ipos) {
    return finemask && (*finemask)(ipos[0], ipos[1], ipos[2]);
  };

  const CCTK_REAL dV = prod(dx);

  // Use per-loop reduction objects to reduce round-off error
  reduction<T, dim> redk;
  // TODO: use loop.hxx code to loop over grid
  for (int k = tmin[2]; k < tmax[2]; ++k) {
    reduction<T, dim> redj;
    for (int j = tmin[1]; j < tmax[1]; ++j) {
      reduction<T, dim> redi;
      for (int i = tmin[0]; i < tmax[0]; ++i) {
        const vect<int, dim> ipos = {i, j, k};

        // For vertex-centred grids, ensure that points at the outer boundary
        // are counted with a weight of 1/2
        std::bitset<8> outer_active = 0b11111111;
        for (int f = 0; f < 2; ++f)
          for (int d = 0; d < dim; ++d)
            if (indextype[d] == 0 && ipos[d] == ibnd[f][d])
              outer_active &= ~faces[f][d];

        // For vertex-centred grids, ensure that points at refinement boundaries
        // are counted with a weight of 1/2
        std::bitset<8> inner_active;
        if (!masked(ipos)) {
          inner_active = 0b11111111;
        } else {
          inner_active = 0b00000000;
          for (int f = 0; f < 2; ++f)
            for (int d = 0; d < dim; ++d)
              if (indextype[d] == 0 &&
                  !(ipos[d] != ibnd[f][d] && masked(ipos + dirs[f][d])))
                inner_active |= faces[f][d];
        }

        assert((~outer_active & ~inner_active).none());

        const std::bitset<8> active = outer_active & inner_active;
        if (active.any()) {
          const T W = active.count() / T(active.size());

          const vect<T, dim> x = x0 + ipos * dx;
          redi += reduction<T, dim>(x, W * dV, vars(i, j, k, n));
        }
      }
      redj += redi;
    }
    redk += redj;
  }
  return redk;
}
} // namespace

reduction<CCTK_REAL, dim> reduce(int gi, int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;

  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  assert(group.grouptype == CCTK_GF);

  reduction<CCTK_REAL, dim> red;
  // TODO: Parallelize over patches and levels
  for (auto &restrict patchdata : ghext->patchdata) {
    for (auto &restrict leveldata : patchdata.leveldata) {
      const auto &restrict groupdata = *leveldata.groupdata.at(gi);
      const amrex::MultiFab &mfab = *groupdata.mfab.at(tl);
      unique_ptr<amrex::iMultiFab> finemask_imfab;

      warn_if_invalid(groupdata, vi, tl, make_valid_int(),
                      []() { return "Before reduction"; });

      const vect<int, dim> indextype = groupdata.indextype;

      const auto &restrict geom = patchdata.amrcore->Geom(leveldata.level);
      const CCTK_REAL *restrict const x01 = geom.ProbLo();
      const CCTK_REAL *restrict const dx1 = geom.CellSize();
      const vect<CCTK_REAL, dim> dx = {dx1[0], dx1[1], dx1[2]};
      const vect<CCTK_REAL, dim> x0v = {x01[0], x01[1], x01[2]};
      const auto x0 = x0v + indextype * dx / 2;

      const int fine_level = leveldata.level + 1;
      if (fine_level < int(patchdata.leveldata.size())) {
        const auto &restrict fine_leveldata =
            patchdata.leveldata.at(fine_level);
        const auto &restrict fine_groupdata = *fine_leveldata.groupdata.at(gi);
        const amrex::MultiFab &fine_mfab = *fine_groupdata.mfab.at(tl);

        const amrex::IntVect reffact{2, 2, 2};

        finemask_imfab = make_unique<amrex::iMultiFab>(makeFineMask(
            mfab, fine_mfab.boxArray(), reffact, geom.periodicity(),
            /*coarse value*/ 0, /* fine value */ 1));
      }

      auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
      // TODO: check that multi-threading actually helps (and we are
      // not dominated by memory latency)
      // TODO: document required version of OpenMP to use custom reductions
#pragma omp parallel reduction(reduction : red)
      for (amrex::MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
        const amrex::Box &bx = mfi.tilebox(); // current tile (without ghosts)
        const vect<int, dim> tmin{bx.smallEnd(0), bx.smallEnd(1),
                                  bx.smallEnd(2)};
        const vect<int, dim> tmax{bx.bigEnd(0) + 1, bx.bigEnd(1) + 1,
                                  bx.bigEnd(2) + 1};
        const amrex::Box &vbx =
            mfi.validbox(); // interior region (without ghosts)
        const vect<int, dim> imin{vbx.smallEnd(0), vbx.smallEnd(1),
                                  vbx.smallEnd(2)};
        const vect<int, dim> imax{vbx.bigEnd(0) + 1, vbx.bigEnd(1) + 1,
                                  vbx.bigEnd(2) + 1};

        const amrex::Array4<const CCTK_REAL> &vars = mfab.array(mfi);

        unique_ptr<amrex::Array4<const int> > finemask;
        if (finemask_imfab) {
          finemask = make_unique<amrex::Array4<const int> >(
              finemask_imfab->array(mfi));
          // Ensure the mask has the correct size
          assert(finemask->begin.x == vars.begin.x);
          assert(finemask->begin.y == vars.begin.y);
          assert(finemask->begin.z == vars.begin.z);
          assert(finemask->end.x == vars.end.x);
          assert(finemask->end.y == vars.end.y);
          assert(finemask->end.z == vars.end.z);
        }

        red += reduce_array(vars, vi, tmin, tmax, indextype, imin, imax,
                            finemask.get(), x0, dx);
      }
    }
  }

  MPI_Datatype datatype = reduction_mpi_datatype<CCTK_REAL, dim>();
  MPI_Op op = reduction_mpi_op();
  MPI_Allreduce(MPI_IN_PLACE, &red, 1, datatype, op, MPI_COMM_WORLD);

  return red;
}

} // namespace CarpetX
