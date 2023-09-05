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

#if 1

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

#else

CCTK_DEVICE void BoundaryKernel::operator()(const Arith::vect<int, dim> &dst,
                                            const int cmin,
                                            const int cmax) const {
  // `src` is the point at which we are looking to determine the boundary value
  Arith::vect<int, dim> src = dst;
  // `delta` (if nonzero) describes a second point at which we are looking, so
  // that we can calculate a gradient for the boundary value
  Arith::vect<int, dim> delta{0, 0, 0};
  for (int d = 0; d < dim; ++d) {
    if (boundaries[d] == boundary_t::dirichlet) {
      // do nothing
    } else if (boundaries[d] == boundary_t::linear_extrapolation) {
      // Same slope:
      //   f'(0)       = f'(h)
      //   f(h) - f(0) = f(2h) - f(h)
      //          f(0) = 2 f(h) - f(2h)
      // f(0) is the boundary point
      src[d] = linear_extrapolation_source[d];
      delta[d] = -inormal[d];
    } else if (boundaries[d] == boundary_t::neumann) {
      // Same value:
      //   f(0) = f(h)
      // f(0) is the boundary point
      src[d] = neumann_source[d];
    } else if (boundaries[d] == boundary_t::robin) {
      // Robin condition, specialized to 1/r fall-off:
      //   f(r) = finf + C/r
      // Determine fall-off constant `C`:
      //   C = r * (f(r) - finf)
      // Solve for value at boundary:
      //   f(r+h) = finf + C / (r + h)
      //          = finf + r / (r + h) * (f(r) - finf)
      // Rewrite using Cartesian coordinates:
      //   C = |x| * (f(x) - finf)
      //   f(x') = finf + C / |x'|
      //         = finf + |x| / |x'| * (f(x) - finf)
      // f(x') is the boundary point
      src[d] = robin_source[d];
    } else if (symmetries[d] == symmetry_t::reflection) {
      src[d] = reflection_offset[d] - dst[d];
    } else if (symmetries[d] == symmetry_t::none &&
               boundaries[d] == boundary_t::none) {
      // this direction is not a boundary; do nothing
    } else {
      // std::cerr << " dst=" << dst << " d=" << d
      //           << " boundaries=" << boundaries
      //           << " symmetries=" << symmetries <<
      //           "\n";
      assert(0);
    }
  }
#ifdef CCTK_DEBUG
  for (int d = 0; d < dim; ++d)
    assert(src[d] >= dmin[d] && src[d] < dmax[d]);
  assert(!all(src == dst));
#endif

  for (int comp = cmin; comp < cmax; ++comp) {
    const CCTK_REAL dirichlet_value = dirichlet_values[comp];
    const CCTK_REAL robin_value = robin_values[comp];
    const CCTK_REAL reflection_parity = reflection_parities[comp];
    const Loop::GF3D2<CCTK_REAL> var(layout, destptr + comp * layout.np);

    CCTK_REAL val;
    if (any(boundaries == boundary_t::dirichlet)) {
      val = dirichlet_value;
    } else {
      val = var(src);
      if (any(boundaries == boundary_t::robin)) {
        for (int d = 0; d < dim; ++d) {
          if (boundaries[d] == boundary_t::robin) {
            using std::sqrt;
            // boundary point
            const auto xb = xmin + dst * dx;
            const auto rb = sqrt(sum(pow2(xb)));
            // interior point
            const auto xi = xmin + src * dx;
            const auto ri = sqrt(sum(pow2(xi)));
            const auto q = ri / rb;
            val = robin_value + q * (val - robin_value);
          }
        }
      }
      if (any(boundaries == boundary_t::linear_extrapolation)) {
        // Calculate gradient
        const CCTK_REAL grad = val - var(src + delta);
        using std::sqrt;
        val += sqrt(sum(pow2(dst - src)) / sum(pow2(delta))) * grad;
      }
#ifdef CCTK_DEBUG
      for (int d = 0; d < dim; ++d)
        assert(dst[d] >= dmin[d] && dst[d] < dmax[d]);
#endif
      if (any(symmetries == symmetry_t::reflection))
        val *= reflection_parity;
    }
    var.store(dst, val);
  }
}

void BoundaryCondition::apply() const {
  // static std::vector<Timer> timers = [&]() {
  //   std::vector<Timer> timers;
  //   timers.reserve(omp_get_max_threads());
  //   for (int t = 0; t < omp_get_max_threads(); ++t) {
  //     std::ostringstream buf;
  //     buf << "BoundaryCondition::apply#" << t;
  //     timers.emplace_back(buf.str());
  //   }
  //   return timers;
  // }();
  // static std::vector<Timer> kernel_timers = [&]() {
  //   std::vector<Timer> timers;
  //   timers.reserve(omp_get_max_threads());
  //   for (int t = 0; t < omp_get_max_threads(); ++t) {
  //     std::ostringstream buf;
  //     buf << "BoundaryKernel::operator()#" << t;
  //     timers.emplace_back(buf.str());
  //   }
  //   return timers;
  // }();
  // Interval interval(timers.at(omp_get_thread_num()));

  assert(tasks.empty());

  for (int nk = -1; nk <= +1; ++nk) {
    for (int nj = -1; nj <= +1; ++nj) {
      for (int ni = -1; ni <= +1; ++ni) {
        const Arith::vect<int, dim> inormal{ni, nj, nk};
        if (all(inormal == 0))
          continue;

        using std::max, std::min;
        Arith::vect<int, dim> bmin, bmax;
        for (int d = 0; d < dim; ++d) {
          if (inormal[d] < 0) {
            // We have to fill the lower boundary
            bmin[d] = dmin[d];
            bmax[d] = min(dmax[d], imin[d]);
          } else if (inormal[d] > 0) {
            // We have to fill the upper boundary
            bmin[d] = max(dmin[d], imax[d]);
            bmax[d] = dmax[d];
          } else {
            // This direction is not a boundary
            bmin[d] = max(dmin[d], imin[d]);
            bmax[d] = min(dmax[d], imax[d]);
          }
        }

        // Do we actually own part of this boundary?
        if (any(bmax <= bmin))
          continue;

        // Find which symmetry or boundary conditions apply to us. On edges
        // or in corners, multiple conditions will apply.
        Arith::vect<symmetry_t, dim> symmetries;
        Arith::vect<boundary_t, dim> boundaries;
        for (int dir = 0; dir < dim; ++dir) {
          if (inormal[dir] == 0) {
            // interior
            symmetries[dir] = symmetry_t::none;
            boundaries[dir] = boundary_t::none;
          } else {
            // face
            const int face = inormal[dir] > 0;
            const symmetry_t symmetry = patchdata.symmetries[face][dir];
            const boundary_t boundary = groupdata.boundaries[face][dir];
            if ((symmetry == symmetry_t::none &&
                 boundary == boundary_t::none) ||
                symmetry == symmetry_t::interpatch ||
                symmetry == symmetry_t::periodic) {
              symmetries[dir] = symmetry_t::none;
              boundaries[dir] = boundary_t::none;
            } else if (symmetry == symmetry_t::reflection) {
              symmetries[dir] = symmetry;
              boundaries[dir] = boundary_t::none;
            } else if (boundary == boundary_t::dirichlet ||
                       boundary == boundary_t::linear_extrapolation ||
                       boundary == boundary_t::neumann ||
                       boundary == boundary_t::robin) {
              symmetries[dir] = symmetry_t::none;
              boundaries[dir] = boundary;
            } else {
#pragma omp critical
              CCTK_ERROR("internal error");
            }
          }
        }

        if (!all(symmetries == symmetry_t::none &&
                 boundaries == boundary_t::none)) {

          const int ncomps = dest.nComp();
          assert(ncomps <= BoundaryKernel::maxncomps);
          const int cmin = 0;
          const int cmax = ncomps;

          std::array<CCTK_REAL, BoundaryKernel::maxncomps> dirichlet_values;
          for (int comp = 0; comp < ncomps; ++comp)
            dirichlet_values[comp] = groupdata.dirichlet_values.at(comp);

          Arith::vect<int, dim> neumann_source;
          for (int d = 0; d < dim; ++d) {
            if (boundaries[d] == boundary_t::neumann) {
              if (inormal[d] != 0) {
                neumann_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
                if (inormal[d] < 0)
                  assert(neumann_source[d] < dmax[d]);
                else
                  assert(neumann_source[d] >= dmin[d]);
              }
            } else {
              neumann_source[d] = 666666666; // go for a segfault
            }
          }

          Arith::vect<int, dim> linear_extrapolation_source;
          for (int d = 0; d < dim; ++d) {
            if (boundaries[d] == boundary_t::linear_extrapolation) {
              assert(inormal[d] != 0);
              linear_extrapolation_source[d] =
                  inormal[d] < 0 ? imin[d] : imax[d] - 1;
              if (inormal[d] < 0) {
                assert(linear_extrapolation_source[d] < dmax[d]);
                assert(linear_extrapolation_source[d] - inormal[d] < dmax[d]);
              } else {
                assert(linear_extrapolation_source[d] >= dmin[d]);
                assert(linear_extrapolation_source[d] - inormal[d] >= dmin[d]);
              }
            } else {
              linear_extrapolation_source[d] = 666666666; // go for a segfault
            }
          }

          Arith::vect<int, dim> robin_source;
          for (int d = 0; d < dim; ++d) {
            if (boundaries[d] == boundary_t::robin) {
              assert(inormal[d] != 0);
              robin_source[d] = inormal[d] < 0 ? imin[d] : imax[d] - 1;
              if (inormal[d] < 0)
                assert(robin_source[d] < dmax[d]);
              else
                assert(robin_source[d] >= dmin[d]);
            } else {
              robin_source[d] = 666666666; // go for a segfault
            }
          }

          std::array<CCTK_REAL, BoundaryKernel::maxncomps> robin_values;
          for (int comp = 0; comp < ncomps; ++comp)
            robin_values[comp] = groupdata.robin_values.at(comp);

          Arith::vect<int, dim> reflection_offset;
          for (int d = 0; d < dim; ++d) {
            if (symmetries[d] == symmetry_t::reflection) {
              assert(inormal[d] != 0);
              reflection_offset[d] =
                  inormal[d] < 0
                      ? 2 * imin[d] - groupdata.indextype.at(d)
                      : 2 * (imax[d] - 1) + groupdata.indextype.at(d);
              if (inormal[d] < 0)
                assert(reflection_offset[d] - bmin[d] < dmax[d]);
              else
                assert(reflection_offset[d] - (bmax[d] - 1) >= dmin[d]);
            } else {
              reflection_offset[d] = 666666666; // go for a segfault
            }
          }

          std::array<CCTK_REAL, BoundaryKernel::maxncomps> reflection_parities;
          for (int comp = 0; comp < ncomps; ++comp) {
            CCTK_REAL reflection_parity = +1;
            for (int d = 0; d < dim; ++d)
              if (symmetries[d] == symmetry_t::reflection)
                reflection_parity *= groupdata.parities.at(comp).at(d);
            using std::fabs;
            assert(fabs(reflection_parity) == 1);
            reflection_parities[comp] = reflection_parity;
          }

          const BoundaryKernel kernel = boundary_kernel(
              inormal, symmetries, boundaries,
              //
              dirichlet_values, neumann_source, linear_extrapolation_source,
              robin_source, robin_values, reflection_offset,
              reflection_parities);

          {
            // Interval kernel_interval(kernel_timers.at(omp_get_thread_num()));
            const amrex::Box box(
                amrex::IntVect(bmin[0], bmin[1], bmin[2]),
                amrex::IntVect(bmax[0] - 1, bmax[1] - 1, bmax[2] - 1));
            for (int comp = cmin; comp < cmax; ++comp)
              amrex::ParallelFor(box,
                                 [kernel, comp] CCTK_DEVICE(
                                     const int i, const int j, const int k)
                                     CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                       const Arith::vect<int, dim> p{i, j, k};
                                       kernel(p, comp, comp + 1);
                                     });
            // amrex::ParallelFor(box,
            //                    [kernel, cmin, cmax] CCTK_DEVICE(
            //                        const int i, const int j, const int k)
            //                        CCTK_ATTRIBUTE_ALWAYS_INLINE {
            //                          const Arith::vect<int, dim> p{i, j, k};
            //                          kernel(p, cmin, cmax);
            //                        });
          }
          // const task_t<BoundaryKernel> task =
          //     make_task(std::move(kernel), bmin, bmax, cmin, cmax);
          // run_task(task);
          // tasks.push_back(std::move(task));

        } // if !all(symmetries == none && boundaries == none)
      }   // for ni, nj, nk
    }
  }

  // run_tasks(tasks);
  // tasks.clear();
}

#endif

} // namespace CarpetX
