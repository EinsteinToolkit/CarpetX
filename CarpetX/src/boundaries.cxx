#include "boundaries.hxx"
#include "boundaries_impl.hxx"
// #include "timer.hxx"

#include <vect.hxx>

#include <AMReX.H>

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <type_traits>
#include <vector>

namespace CarpetX {

// TODO: The algorithm below applies symmetry and boundary conditions
// in all directions "at once", and thus has to be able to apply
// multiple boundary conditions simultaneously. This adds complexity
// and makes the code slow to build.
//
// An alternative would be to apply symmetry and boundary conditions
// sequentially, first on faces, then edges, then corners. On edges
// and corners we would choose one (preferred) symmetry or boundary
// condition to be applied there. This would simplify the code, reduce
// build times, and also improve run time. It would not be necessary
// any more to run on edges or corners; instead, one could extend some
// of the faces to cover edges and corners.

////////////////////////////////////////////////////////////////////////////////

BoundaryCondition::BoundaryCondition(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::FArrayBox &dest)
    : groupdata(groupdata), patchdata(ghext->patchdata.at(groupdata.patch)),
      geom(patchdata.amrcore->Geom(groupdata.level)), dest(dest),
      imin{geom.Domain().smallEnd(0), geom.Domain().smallEnd(1),
           geom.Domain().smallEnd(2)},
      imax{geom.Domain().bigEnd(0) + 1 + !groupdata.indextype[0],
           geom.Domain().bigEnd(1) + 1 + !groupdata.indextype[1],
           geom.Domain().bigEnd(2) + 1 + !groupdata.indextype[2]},
      // Vertices are shifted outwards by 1/2 grid spacing
      xmin{geom.ProbLo(0) + CCTK_REAL(0.5) * groupdata.indextype[0] - 1,
           geom.ProbLo(1) + CCTK_REAL(0.5) * groupdata.indextype[1] - 1,
           geom.ProbLo(2) + CCTK_REAL(0.5) * groupdata.indextype[2] - 1},
      xmax{geom.ProbHi(0) - CCTK_REAL(0.5) * groupdata.indextype[0] + 1,
           geom.ProbHi(1) - CCTK_REAL(0.5) * groupdata.indextype[1] + 1,
           geom.ProbHi(2) - CCTK_REAL(0.5) * groupdata.indextype[2] + 1},
      dx{geom.CellSize(0), geom.CellSize(1), geom.CellSize(2)},
      dmin{dest.box().smallEnd(0), dest.box().smallEnd(1),
           dest.box().smallEnd(2)},
      dmax{dest.box().bigEnd(0) + 1, dest.box().bigEnd(1) + 1,
           dest.box().bigEnd(2) + 1},
      layout(dmin, dmax), destptr(dest.dataPtr()) {
  // static std::vector<Timer> timers = [&]() {
  //   std::vector<Timer> timers;
  //   timers.reserve(omp_get_max_threads());
  //   for (int t = 0; t < omp_get_max_threads(); ++t) {
  //     std::ostringstream buf;
  //     buf << "BoundaryCondition::BoundaryCondition#" << t;
  //     timers.emplace_back(buf.str());
  //   }
  //   return timers;
  // }();
  // Interval interval(timers.at(omp_get_thread_num()));

  // Check centering
  for (int d = 0; d < dim; ++d)
    assert((dest.box().type(d) == amrex::IndexType::CELL) ==
           groupdata.indextype[d]);

  // Ensure we have enough ghost zones
  // TODO: Implement this
  // for (int d = 0; d < dim; ++d)
  //   assert(... <= groupdata.nghostzones[d]);

  // Ensure the destination array is large enough
  for (int d = 0; d < dim; ++d) {
    assert(dest.box().smallEnd(d) <= dmin[d]);
    assert(dest.box().bigEnd(d) + 1 >= dmax[d]);
  }
}

void BoundaryCondition::apply() const {
  apply_on_face<NEG, NEG, NEG>();
  apply_on_face<INT, NEG, NEG>();
  apply_on_face<POS, NEG, NEG>();
  apply_on_face<NEG, INT, NEG>();
  apply_on_face<INT, INT, NEG>();
  apply_on_face<POS, INT, NEG>();
  apply_on_face<NEG, POS, NEG>();
  apply_on_face<INT, POS, NEG>();
  apply_on_face<POS, POS, NEG>();

  apply_on_face<NEG, NEG, INT>();
  apply_on_face<INT, NEG, INT>();
  apply_on_face<POS, NEG, INT>();
  apply_on_face<NEG, INT, INT>();
  // apply_on_face<INT, INT, INT>();
  apply_on_face<POS, INT, INT>();
  apply_on_face<NEG, POS, INT>();
  apply_on_face<INT, POS, INT>();
  apply_on_face<POS, POS, INT>();

  apply_on_face<NEG, NEG, POS>();
  apply_on_face<INT, NEG, POS>();
  apply_on_face<POS, NEG, POS>();
  apply_on_face<NEG, INT, POS>();
  apply_on_face<INT, INT, POS>();
  apply_on_face<POS, INT, POS>();
  apply_on_face<NEG, POS, POS>();
  apply_on_face<INT, POS, POS>();
  apply_on_face<POS, POS, POS>();
}

} // namespace CarpetX
