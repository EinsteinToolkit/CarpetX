#include "pdesolvers.hxx"

// TODO: Don't include files from other thorns; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/fillpatch.hxx"
#include "../../CarpetX/src/reduction.hxx"
#include "../../CarpetX/src/schedule.hxx"
#include "../../CarpetX/src/task_manager.hxx"

#include <arr.hxx>
#include <div.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_MultiFabUtil.H>

#include <mpi.h>

#include <petscdm.h>
#include <petscsnes.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace PDESolvers {

const int tl = 0;

// TODO: Generalize this
const Arith::arr<int, 3> indextype{0, 0, 0}; // vertex centred

////////////////////////////////////////////////////////////////////////////////

// Erik Schnetter the blocks for vertex centred grids overlap by one
// point at the block boundaries (see
// https://amrex-codes.github.io/amrex/docs_html/Basics.html#id4).
// does FillPatchSingleLevel ensure these points are consistent? or is
// the domain double-valued there, and only the ghost points (not
// shown in this figure) are set?

// weiqun If the original source data are consistent, the results of
// FillPatchSingleLevel are also consistent.  Otherwise, it's undefined which
// original values will used for which boxes.  It depends on how boxes are
// distributed.

// If you want to the data consistent, you can use MultiFab's OverrideSync
// functions.
//     void OverrideSync (const Periodicity& period =
//         Periodicity::NonPeriodic());
//     void OverrideSync (const iMultiFab& msk,
//         const Periodicity& period = Periodicity::NonPeriodic());

// The firs version will use data from the grid with the lowest grid number to
// override others.

// In the second version, you can pass your own mask.

// Erik Schnetter thanks!

// weiqun Even if you want to use the first version's mask, you might still want
// to call
//     //! Owner is the grid with the lowest grid number containing the data.
//     std::unique_ptr<iMultiFab> OwnerMask (const Periodicity& period =
//         Periodicity::NonPeriodic()) const;
// and then use the second function.

// If you need to call OverrideSync repeatedly, you could save the mask.

// Erik Schnetter thanks!

////////////////////////////////////////////////////////////////////////////////

// 1: interior (defined by user)
// 2: outer boundary (also defined by users)
// 3: defined via synchronization (copied from another box on same level)
// 4: defined via prolongation (interpolated from next coarser level)
// 5: defined via restriction (injected from next finer level)
enum class point_type_t {
  undf = 0,
  intr = 1,
  bdry = 2,
  sync = 3,
  prol = 4,
  rest = 5
};

void define_point_type() {
  // Decode Cactus variables
  const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);
  const int vn_ind = CCTK_VarIndex("PDESolvers::communication_indicator");
  assert(vn_ind >= 0);
  const int gi_ind = CCTK_GroupIndexFromVarI(vn_ind);
  assert(gi_ind >= 0);
  const int v0_ind = CCTK_FirstVarIndexI(gi_ind);
  assert(v0_ind >= 0);
  const int vi_ind = vn_ind - v0_ind;
  assert(vi_ind >= 0);

  // Initialize point type everywhere, assuming there is no synchronization,
  // prolongation, or restriction
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             gf_pt(p.I) = int(point_type_t::undf);
                           });
    grid.loop_int<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             assert(gf_pt(p.I) == int(point_type_t::undf));
                             gf_pt(p.I) = int(point_type_t::intr);
                           });
    grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             assert(gf_pt(p.I) == int(point_type_t::undf));
                             gf_pt(p.I) = int(point_type_t::bdry);
                           });
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    leveldata.groupdata.at(gi_pt)->valid.at(tl).at(vi_pt).set_all(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::define_point_type"; });
  });

  // Set indicator to level
  CarpetX::active_levels->loop_parallel(
      [&](const int patch, const int level, const int index,
          const int component, const cGH *restrict const cctkGH) {
        const Loop::GF3D2layout layout1(cctkGH, indextype);
        const Loop::GF3D2<CCTK_REAL> gf_ind(
            layout1,
            static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
        const CarpetX::GridDescBase grid(cctkGH);
        grid.loop_all<0, 0, 0>(grid.nghostzones,
                               [&](const Loop::PointDesc &p)
                                   ARITH_INLINE { gf_ind(p.I) = level; });
        const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
        const auto &leveldata = patchdata.leveldata.at(level);
        leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set_all(
            CarpetX::make_valid_all(),
            []() { return "PDESolver::define_point_type before restricting"; });
      });

  // Restrict
  for (const auto &patchdata : CarpetX::ghext->patchdata) {
    for (int level = int(patchdata.leveldata.size()) - 2; level >= 0; --level) {
      const auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      const auto &groupdata = *leveldata.groupdata.at(gi_ind);
      const auto &finegroupdata = *fineleveldata.groupdata.at(gi_ind);
      amrex::MultiFab &mfab_ind = *groupdata.mfab.at(tl);
      const amrex::MultiFab &finemfab_ind = *finegroupdata.mfab.at(tl);
      const amrex::IntVect reffact{2, 2, 2};
      const int rank = sum(Arith::vect<int, 3>(indextype));
      switch (rank) {
      case 0:
        average_down_nodal(finemfab_ind, mfab_ind, reffact);
        break;
      case 1:
        average_down_edges(finemfab_ind, mfab_ind, reffact);
        break;
      case 2:
        average_down_faces(finemfab_ind, mfab_ind, reffact);
        break;
      case 3:
        average_down(finemfab_ind, mfab_ind, 0, 1 /*nvars*/, reffact);
        break;
      default:
        assert(0);
      }
    }
  }

  // Check where the indicator changed; these are the restricted points
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             if (gf_ind(p.I) != level) {
                               assert(gf_pt(p.I) != int(point_type_t::bdry));
                               gf_pt(p.I) = int(point_type_t::rest);
                             }
                           });
  });

  // Set indicator to index
  CarpetX::active_levels->loop_parallel(
      [&](const int patch, const int level, const int index,
          const int component, const cGH *restrict const cctkGH) {
        const Loop::GF3D2layout layout1(cctkGH, indextype);
        const Loop::GF3D2<CCTK_REAL> gf_ind(
            layout1,
            static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
        const CarpetX::GridDescBase grid(cctkGH);
        grid.loop_all<0, 0, 0>(grid.nghostzones,
                               [&](const Loop::PointDesc &p)
                                   ARITH_INLINE { gf_ind(p.I) = index; });
        const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
        const auto &leveldata = patchdata.leveldata.at(level);
        leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set_all(
            CarpetX::make_valid_all(), []() {
              return "PDESolver::define_point_type before synchronizing";
            });
      });

  // Synchronize (as if there was no prolongation)
  {
    CarpetX::task_manager tasks1;
    CarpetX::task_manager tasks2;
    for (const auto &patchdata : CarpetX::ghext->patchdata) {
      for (const auto &leveldata : patchdata.leveldata) {
        const int level = leveldata.level;
        const auto &groupdata = *leveldata.groupdata.at(gi_ind);
        amrex::MultiFab &mfab_ind = *groupdata.mfab.at(tl);
        tasks1.submit_serially(
            [&tasks2, &groupdata, &mfab_ind, &patchdata, level]() {
              FillPatch_Sync(tasks2, groupdata, mfab_ind,
                             patchdata.amrcore->Geom(level));
            });
      }
    }
    tasks1.run_tasks_serially();
    tasks2.run_tasks_serially();
  }

  // Check where the indicator changed; these are the synchronized points.
  // If points are both restricted and synchronized, then we count them as
  // synchronized.
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             if (gf_ind(p.I) != index) {
                               assert(gf_pt(p.I) == int(point_type_t::undf) ||
                                      gf_pt(p.I) == int(point_type_t::bdry));
                               gf_pt(p.I) = int(point_type_t::sync);
                             }
                           });
  });

  // Set indicator to level
  CarpetX::active_levels->loop_parallel(
      [&](const int patch, const int level, const int index,
          const int component, const cGH *restrict const cctkGH) {
        const Loop::GF3D2layout layout1(cctkGH, indextype);
        const Loop::GF3D2<CCTK_REAL> gf_ind(
            layout1,
            static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
        const CarpetX::GridDescBase grid(cctkGH);
        grid.loop_all<0, 0, 0>(grid.nghostzones,
                               [&](const Loop::PointDesc &p)
                                   ARITH_INLINE { gf_ind(p.I) = level; });
        const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
        const auto &leveldata = patchdata.leveldata.at(level);
        leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set_all(
            CarpetX::make_valid_all(), []() {
              return "PDESolver::define_point_type before prolongating";
            });
      });

  // Prolongate and synchronize (we cannot just prolongate)
  {
    CarpetX::task_manager tasks1;
    CarpetX::task_manager tasks2;
    CarpetX::task_manager tasks3;
    for (const auto &patchdata : CarpetX::ghext->patchdata) {
      for (const auto &leveldata : patchdata.leveldata) {
        const int level = leveldata.level;
        if (level == 0)
          continue;
        const auto &coarseleveldata = patchdata.leveldata.at(level - 1);
        const auto &groupdata = *leveldata.groupdata.at(gi_ind);
        const auto &coarsegroupdata = *coarseleveldata.groupdata.at(gi_ind);
        amrex::MultiFab &mfab_ind = *groupdata.mfab.at(tl);
        amrex::MultiFab &coarsemfab_ind = *coarsegroupdata.mfab.at(tl);
        amrex::Interpolater *const interpolator = groupdata.interpolator;
        tasks1.submit_serially([&tasks2, &tasks3, &groupdata, &coarsegroupdata,
                                &mfab_ind, &coarsemfab_ind, &patchdata,
                                &interpolator, level]() {
          FillPatch_ProlongateGhosts(tasks2, tasks3, groupdata, coarsegroupdata,
                                     mfab_ind, coarsemfab_ind,
                                     patchdata.amrcore->Geom(level),
                                     patchdata.amrcore->Geom(level - 1),
                                     interpolator, groupdata.bcrecs);
        });
      }
    }
    tasks1.run_tasks_serially();
    tasks2.run_tasks_serially();
    tasks3.run_tasks_serially();
  }

  // Check where the indicator changed; these are the prolongated points. If
  // points are both restricted and prolongated, then we count them as
  // prolongated.
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          if (gf_ind(p.I) != level && gf_pt(p.I) != int(point_type_t::sync)) {
            assert(gf_pt(p.I) == int(point_type_t::undf) ||
                   gf_pt(p.I) == int(point_type_t::bdry));
            // Points cannot be both restricted and prolongated
            assert(gf_pt(p.I) != int(point_type_t::rest));
            gf_pt(p.I) = int(point_type_t::prol);
          }
        });
  });

  // Invalidate indicator
  CarpetX::active_levels->loop_serially([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set_all(
        CarpetX::valid_t(), []() { return "PDESolver::define_point_type"; });
  });

  // Collect some statistics
  Arith::vect<int, 6> npoints{0, 0, 0, 0, 0, 0};
  CarpetX::active_levels->loop_parallel(
      [&](const int patch, const int level, const int index,
          const int component, const cGH *restrict const cctkGH) {
        const Loop::GF3D2layout layout1(cctkGH, indextype);
        const Loop::GF3D2<const CCTK_REAL> gf_pt(
            layout1, static_cast<const CCTK_REAL *>(
                         CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
        const CarpetX::GridDescBase grid(cctkGH);
        Arith::vect<int, 6> npoints1{0, 0, 0, 0, 0, 0};
        grid.loop_all<0, 0, 0>(grid.nghostzones,
                               [&](const Loop::PointDesc &p) ARITH_INLINE {
                                 const int pt = int(gf_pt(p.I));
                                 assert(pt >= 0 && pt < 6);
                                 ++npoints1[pt];
                               });
        for (int n = 0; n < 6; ++n)
#pragma omp atomic
          npoints[n] += npoints1[n];
      });
  CCTK_VINFO("Point type counts:");
  CCTK_VINFO("  undefined:    %d", npoints[int(point_type_t::undf)]);
  CCTK_VINFO("  interior:     %d", npoints[int(point_type_t::intr)]);
  CCTK_VINFO("  boundary:     %d", npoints[int(point_type_t::bdry)]);
  CCTK_VINFO("  synchronized: %d", npoints[int(point_type_t::sync)]);
  CCTK_VINFO("  prolongated:  %d", npoints[int(point_type_t::prol)]);
  CCTK_VINFO("  restricted:   %d", npoints[int(point_type_t::rest)]);
  CCTK_VINFO("  total:        %d", sum(npoints));
  assert(npoints[int(point_type_t::undf)] == 0);
}

////////////////////////////////////////////////////////////////////////////////

void enumerate_points(
    int &restrict npoints_local, int &restrict npoints_global,
    std::vector<std::vector<int> > &restrict component_offsets,
    std::vector<std::vector<int> > &restrict component_sizes,
    int &restrict npoints_prolongated_local,
    int &restrict npoints_prolongated_global,
    std::vector<std::vector<int> > &restrict component_prolongated_offsets,
    std::vector<std::vector<int> > &restrict component_prolongated_sizes,
    csr_t &Jp) {
  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  static const int vn_idx = CCTK_VarIndex("PDESolvers::idx");
  assert(vn_idx >= 0);
  static const int gi_idx = CCTK_GroupIndexFromVarI(vn_idx);
  assert(gi_idx >= 0);
  static const int v0_idx = CCTK_FirstVarIndexI(gi_idx);
  assert(v0_idx >= 0);
  static const int vi_idx = vn_idx - v0_idx;
  assert(vi_idx >= 0);
  const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);
  // const int vn_vcx = CCTK_VarIndex("CoordinatesX::vcoordx");
  // assert(vn_vcx >= 0);
  // const int gi_vcx = CCTK_GroupIndexFromVarI(vn_vcx);
  // assert(gi_vcx >= 0);
  // const int v0_vcx = CCTK_FirstVarIndexI(gi_vcx);
  // assert(v0_vcx >= 0);
  // const int vi_vcx = vn_vcx - v0_vcx;
  // assert(vi_vcx >= 0);

  // Determine number of components per level
  assert(CarpetX::ghext->num_patches() == 1);
  const auto &patchdata = CarpetX::ghext->patchdata.at(0);
  std::vector<int> level_sizes(patchdata.leveldata.size(), 0);
  std::vector<int> level_maxcomponents(patchdata.leveldata.size(), -1);
  CarpetX::active_levels->loop_serially(
      [&](const int patch, const int level, const int index,
          const int component, const cGH *restrict const cctkGH) {
        // #pragma omp atomic
        ++level_sizes.at(level);
        int &maxcomponent = level_maxcomponents.at(level);
        using std::max;
        // #pragma omp critical
        maxcomponent = max(maxcomponent, component);
      });

  // Allocate data structure
  component_offsets.resize(patchdata.leveldata.size()); // process local
  component_sizes.resize(patchdata.leveldata.size());
  component_prolongated_offsets.resize(
      patchdata.leveldata.size()); // process local
  component_prolongated_sizes.resize(patchdata.leveldata.size());
  for (std::size_t level = 0; level < patchdata.leveldata.size(); ++level) {
    assert(level_maxcomponents.at(level) + 1 == level_sizes.at(level));
    component_offsets.at(level).resize(level_sizes.at(level));
    component_sizes.at(level).resize(level_sizes.at(level));
    component_prolongated_offsets.at(level).resize(level_sizes.at(level));
    component_prolongated_sizes.at(level).resize(level_sizes.at(level));
  }

  // Enumerate and count points
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    int npoints = 0;
    int npoints_prolongated = 0;
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             switch (int(gf_pt(p.I))) {
                             case int(point_type_t::intr):
                               ++npoints;
                               break;
                             case int(point_type_t::bdry):
                             case int(point_type_t::rest):
                             case int(point_type_t::sync):
                               // ignore this point
                               break;
                             case int(point_type_t::prol):
                               ++npoints_prolongated;
                               break;
                             default:
                               assert(0);
                             }
                           });
    component_sizes.at(level).at(component) = npoints;
    component_prolongated_sizes.at(level).at(component) = npoints_prolongated;
  });

  // Local exclusive prefix sum
  int npoints = 0;
  int npoints_prolongated = 0;
  for (std::size_t level = 0; level < component_offsets.size(); ++level) {
    for (std::size_t component = 0;
         component < component_offsets.at(level).size(); ++component) {
      component_offsets.at(level).at(component) = npoints;
      npoints += component_sizes.at(level).at(component);
      component_prolongated_offsets.at(level).at(component) =
          npoints_prolongated;
      npoints_prolongated +=
          component_prolongated_sizes.at(level).at(component);
    }
  }
  npoints_local = npoints;
  npoints_prolongated_local = npoints_prolongated;

  // Global exclusive prefix sum
  int npoints_offset = 0;
  MPI_Exscan(&npoints_local, &npoints_offset, 1, MPI_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Allreduce(&npoints_local, &npoints_global, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);
  int npoints_prolongated_offset = 0;
  MPI_Exscan(&npoints_prolongated_local, &npoints_prolongated_offset, 1,
             MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&npoints_prolongated_local, &npoints_prolongated_global, 1,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Set indices
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_idx(
        layout1,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    int idx = npoints_offset + component_offsets.at(level).at(component);
    int idx_prolongated = npoints_prolongated_offset +
                          component_prolongated_offsets.at(level).at(component);
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
            // solve for this point
            gf_idx(p.I) = idx++;
            break;
          case int(point_type_t::bdry):
            gf_idx(p.I) = -1;
            break;
          case int(point_type_t::rest):
          case int(point_type_t::sync):
            // undefined (will be set later)
            gf_idx(p.I) = -2;
            break;
          case int(point_type_t::prol):
            gf_idx(p.I) = prolongation_index_offset + idx_prolongated++;
            break;
          default:
            assert(0);
          }
        });
    assert(idx == npoints_offset + component_offsets.at(level).at(component) +
                      component_sizes.at(level).at(component));
    if (!(idx_prolongated ==
          npoints_prolongated_offset +
              component_prolongated_offsets.at(level).at(component) +
              component_prolongated_sizes.at(level).at(component))) {
      std::cout << "level=" << level << "\n";
      std::cout << "component=" << component << "\n";
      std::cout << "idx_prolongated=" << idx_prolongated << "\n";
      std::cout << "npoints_prolongated_offset=" << npoints_prolongated_offset
                << "\n";
      std::cout << "component_prolongated_offsets="
                << component_prolongated_offsets.at(level).at(component)
                << "\n";
      std::cout << "component_prolongated_sizes="
                << component_prolongated_sizes.at(level).at(component) << "\n";
    }
    assert(idx_prolongated ==
           npoints_prolongated_offset +
               component_prolongated_offsets.at(level).at(component) +
               component_prolongated_sizes.at(level).at(component));
    const auto &leveldata = patchdata.leveldata.at(level);
    leveldata.groupdata.at(gi_idx)->valid.at(tl).at(vi_idx).set_all(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::enumerate_points"; });
  });

  // Restrict index
  for (const auto &patchdata : CarpetX::ghext->patchdata) {
    for (int level = int(patchdata.leveldata.size()) - 2; level >= 0; --level) {
      const auto &leveldata = patchdata.leveldata.at(level);
      const auto &fineleveldata = patchdata.leveldata.at(level + 1);
      amrex::MultiFab &mfab_idx = *leveldata.groupdata.at(gi_idx)->mfab.at(tl);
      const amrex::MultiFab &finemfab_idx =
          *fineleveldata.groupdata.at(gi_idx)->mfab.at(tl);
      const amrex::IntVect reffact{2, 2, 2};
      const int rank = sum(Arith::vect<int, 3>(indextype));
      switch (rank) {
      case 0:
        average_down_nodal(finemfab_idx, mfab_idx, reffact);
        break;
      case 1:
        average_down_edges(finemfab_idx, mfab_idx, reffact);
        break;
      case 2:
        average_down_faces(finemfab_idx, mfab_idx, reffact);
        break;
      case 3:
        average_down(finemfab_idx, mfab_idx, 0, 1 /*nvars*/, reffact);
        break;
      default:
        assert(0);
      }
    }
  }

  // Synchronize index
  {
    CarpetX::task_manager tasks1;
    CarpetX::task_manager tasks2;
    for (const auto &patchdata : CarpetX::ghext->patchdata) {
      for (const auto &leveldata : patchdata.leveldata) {
        const int level = leveldata.level;
        const auto &groupdata = *leveldata.groupdata.at(gi_idx);
        amrex::MultiFab &mfab_idx = *groupdata.mfab.at(tl);
        tasks1.submit_serially(
            [&tasks2, &patchdata, &groupdata, &mfab_idx, level]() {
              FillPatch_Sync(tasks2, groupdata, mfab_idx,
                             patchdata.amrcore->Geom(level));
            });
      }
    }
    tasks1.run_tasks_serially();
    tasks2.run_tasks_serially();
  }

  // Check that restriction and synchronization worked
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_idx(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
          case int(point_type_t::sync):
          case int(point_type_t::rest):
            assert(gf_idx(p.I) >= 0 && gf_idx(p.I) < prolongation_index_offset);
            break;
          case int(point_type_t::bdry):
            assert(gf_idx(p.I) == -1);
            break;
          case int(point_type_t::prol):
            assert(gf_idx(p.I) >= prolongation_index_offset);
            break;
          default:
            assert(0);
          }
        });
    const auto &leveldata = patchdata.leveldata.at(level);
    leveldata.groupdata.at(gi_idx)->valid.at(tl).at(vi_idx).set_all(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::enumerate_points"; });
  });

  // Calculate prolongation Jacobian
  {
    // Find all coarse indices which are sources of prolongation
    std::vector<
        std::vector<std::vector<std::tuple<int, Arith::vect<int, 3> > > > >
        locations(patchdata.leveldata.size());
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      locations.at(level).resize(level_sizes.at(level));
    CarpetX::active_levels->loop_parallel(
        [&](const int patch, const int level, const int index,
            const int component, const cGH *restrict const cctkGH) {
          // Level 0 has no prolongated points
          if (level == 0)
            assert(component_prolongated_sizes.at(level).at(component) == 0);
          if (component_prolongated_sizes.at(level).at(component) == 0)
            return;
          const auto &mfab =
              *patchdata.leveldata.at(level).groupdata.at(gi_idx)->mfab.at(tl);
          const auto &fabbox = mfab.fabbox(index);
          const Arith::vect<int, 3> amrex_origin{
              fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)};
          const Loop::GF3D2layout layout1(cctkGH, indextype);
          const Loop::GF3D2<const CCTK_REAL> gf_pt(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
          const Loop::GF3D2<const CCTK_REAL> gf_idx(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
          const CarpetX::GridDescBase grid(cctkGH);
          // Arith::vect<int, 3> imin, imax;
          // grid.box_all<0, 0, 0>(grid.nghostzones, imin, imax);
          const auto offset = amrex_origin;

          std::vector<std::tuple<int, Arith::vect<int, 3> > > locs;
          grid.loop_all<0, 0, 0>(
              grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
                if (int(gf_pt(p.I)) == int(point_type_t::prol)) {
                  const int fL = level;
                  const Arith::vect<int, 3> fI = p.I + offset;
                  const int cL = fL - 1;
                  const Arith::vect<int, 3> cI0 = div_floor(fI, 2);
                  const Arith::vect<int, 3> cIm = mod_floor(fI, 2);
                  // linear interpolation
                  for (int k = 0; k < cIm[2] + 1; ++k) {
                    for (int j = 0; j < cIm[1] + 1; ++j) {
                      for (int i = 0; i < cIm[0] + 1; ++i) {
                        const Arith::vect<int, 3> cI =
                            cI0 + Arith::vect<int, 3>{i, j, k};
                        locs.emplace_back(cL, cI);
                      }
                    }
                  }
                }
              });
          locations.at(level).at(component) = std::move(locs);
        });

    // Find the components where the prolongation sources live, and determine
    // their Jacobian indices
    std::vector<std::vector<std::vector<int> > > indices(
        patchdata.leveldata.size());
    // std::vector<std::vector<std::vector<CCTK_REAL> > > vcoordxs(
    //     patchdata.leveldata.size());
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level) {
      indices.at(level).resize(level_sizes.at(level));
      // vcoordxs.at(level).resize(level_sizes.at(level));
      for (int component = 0; component < level_sizes.at(level); ++component) {
        indices.at(level).at(component).resize(
            locations.at(level).at(component).size(), 0x80000000U);
        // vcoordxs.at(level).at(component).resize(
        //     locations.at(level).at(component).size(), -1.0 / 0.0);
      }
    }
    CarpetX::active_levels->loop_parallel(
        [&](const int patch, const int level, const int index,
            const int component, const cGH *restrict const cctkGH) {
          const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
          // The finest level has no prolongation sources
          if (level == int(patchdata.leveldata.size()) - 1)
            return;
          const auto &mfab =
              *patchdata.leveldata.at(level).groupdata.at(gi_idx)->mfab.at(tl);
          const auto &fabbox = mfab.fabbox(index);
          const Arith::vect<int, 3> amrex_origin{
              fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)};
          const Loop::GF3D2layout layout1(cctkGH, indextype);
          const Loop::GF3D2<const CCTK_REAL> gf_pt(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
          const Loop::GF3D2<const CCTK_REAL> gf_idx(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
          // const Loop::GF3D2<const CCTK_REAL> gf_vcoordx(
          //     layout1,
          //     static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl,
          //     vn_vcx)));
          const CarpetX::GridDescBase grid(cctkGH);
          Arith::vect<int, 3> imin, imax;
          grid.box_all<0, 0, 0>(grid.nghostzones, imin, imax);
          const auto offset = amrex_origin;

          for (int level1 = 0; level1 < int(patchdata.leveldata.size());
               ++level1) {
            for (int component1 = 0; component1 < level_sizes.at(level1);
                 ++component1) {
              const auto &locs = locations.at(level1).at(component1);
              auto &idxs = indices.at(level1).at(component1);
              // auto &vcxs = vcoordxs.at(level1).at(component1);
              for (int n = 0; n < int(locs.size()); ++n) {
                const auto &loc = locs.at(n);
                const auto &cL = std::get<0>(loc);
                const auto &cI = std::get<1>(loc);
                if (cL == level) {
                  const Arith::vect<int, 3> I = cI - offset;
                  if (all(imin <= I && I < imax)) {
                    const int idx = int(gf_idx(I));
                    // We require that a prolongation source is not prolongated
                    // itself. (We could extend the algorithm to handle this.)
                    assert(idx < prolongation_index_offset);
                    int &idxn = idxs.at(n);
#pragma omp atomic write
                    idxn = idx;
                    // CCTK_REAL &vcxn = vcxs.at(n);
                    // #pragma omp atomic write
                    // vcxn = gf_vcoordx(I);
                  }
                }
              }
            }
          }
        });
    // Note: Prolongation boundaries do not yet work with multiple processes
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level) {
      for (int component = 0; component < level_sizes.at(level); ++component) {
        for (const auto &idx : indices.at(level).at(component))
          assert(idx != int(0x80000000U));
        // for (const auto &vcx : vcoordxs.at(level).at(component))
        //   assert(vcx != -1.0 / 0.0);
      }
    }

    // Calculate the Jacobian stencils
    std::vector<std::vector<std::vector<std::tuple<int, int, CCTK_REAL> > > >
        Jpvalss(patchdata.leveldata.size());
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      Jpvalss.at(level).resize(level_sizes.at(level));
    CarpetX::active_levels->loop_parallel(
        [&](const int patch, const int level, const int index,
            const int component, const cGH *restrict const cctkGH) {
          // Level 0 has no prolongated points
          if (level == 0)
            assert(component_prolongated_sizes.at(level).at(component) == 0);
          if (component_prolongated_sizes.at(level).at(component) == 0)
            return;
          const auto &mfab = *CarpetX::ghext->patchdata.at(patch)
                                  .leveldata.at(level)
                                  .groupdata.at(gi_idx)
                                  ->mfab.at(tl);
          const auto &fabbox = mfab.fabbox(index);
          const Arith::vect<int, 3> amrex_origin{
              fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)};
          const Loop::GF3D2layout layout1(cctkGH, indextype);
          const Loop::GF3D2<const CCTK_REAL> gf_pt(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
          const Loop::GF3D2<const CCTK_REAL> gf_idx(
              layout1, static_cast<const CCTK_REAL *>(
                           CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
          // const Loop::GF3D2<const CCTK_REAL> gf_vcoordx(
          //     layout1,
          //     static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl,
          //     vn_vcx)));
          const CarpetX::GridDescBase grid(cctkGH);
          const auto offset = amrex_origin;

          const auto &locs = locations.at(level).at(component);
          const auto &idxs = indices.at(level).at(component);
          // const auto &vcxs = vcoordxs.at(level).at(component);
          assert(idxs.size() == locs.size());
          // assert(vcxs.size() == locs.size());
          int n = 0;
          std::vector<std::tuple<int, int, CCTK_REAL> > Jpvals;
          Jpvals.reserve(locs.size());
          grid.loop_all<0, 0, 0>(
              grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
                if (int(gf_pt(p.I)) == int(point_type_t::prol)) {
                  assert(n < int(idxs.size()));
                  const int idx = int(gf_idx(p.I)) - prolongation_index_offset;
                  const int cL = level - 1;
                  const Arith::vect<int, 3> cIm = mod_floor(p.I + offset, 2);
                  const Arith::vect<int, 3> cI0 = div_floor(p.I + offset, 2);
                  for (int k = 0; k < cIm[2] + 1; ++k) {
                    for (int j = 0; j < cIm[1] + 1; ++j) {
                      for (int i = 0; i < cIm[0] + 1; ++i) {
                        const Arith::vect<int, 3> cI =
                            cI0 + Arith::vect<int, 3>{i, j, k};
                        assert(std::get<0>(locs.at(n)) == cL);
                        assert(all(std::get<1>(locs.at(n)) == cI));
                        const int row = idx;
                        const int col = idxs.at(n);
                        // Ignore boundary points
                        if (col >= 0) {
                          // We assume that a prolongation source is not
                          // prolongated itself. (We could extend the algorithm
                          // to handle this.)
                          assert(col < prolongation_index_offset);
                          const CCTK_REAL val =
                              1.0 /
                              ((cIm[0] + 1) * (cIm[1] + 1) * (cIm[2] + 1));
                          Jpvals.emplace_back(row, col, val);
                        }
                        ++n;
                      }
                    }
                  }

                  // if (all(cIm == Arith::vect<int, 3>::pure(0))) {
                  //   const CCTK_REAL vcx = gf_vcoordx(p.I);
                  //   assert(vcxs.at(n) == vcx);
                  // }
                }
              });
          assert(n == int(locs.size()));
          Jpvalss.at(level).at(component) = std::move(Jpvals);
        });

    // Convert Jpvals into sparse matrix
    Jp = csr_t(npoints_prolongated_global, npoints_global, Jpvalss);
  }
}

void copy_Cactus_to_PETSc(
    Vec vec, const std::vector<int> &varinds,
    const std::vector<std::vector<int> > &component_offsets,
    const std::vector<std::vector<int> > &component_sizes) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }
  static const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  static const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  static const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  static const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);

  // Get vector from PETSc
  CCTK_REAL *vec_ptr;
  ierr = VecGetArray(vec, &vec_ptr);
  assert(!ierr);
  CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy Cactus vector to PETSc
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const CarpetX::GridDescBase grid(cctkGH);
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    for (int n = 0; n < nvars; ++n)
      CarpetX::error_if_invalid(*leveldata.groupdata.at(gis.at(n)), vis.at(n),
                                tl, CarpetX::make_valid_int(), []() {
                                  return "PDESolver::copy_Cactus_to_PETSc";
                                });
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    std::vector<Loop::GF3D2<const CCTK_REAL> > gfs;
    gfs.reserve(nvars);
    for (int n = 0; n < nvars; ++n)
      gfs.emplace_back(layout1, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                                    cctkGH, tl, varinds.at(n))));
    const int component_offset = component_offsets.at(level).at(component);
    int nelems = 0;
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
            for (int n = 0; n < nvars; ++n)
              petsc_ptr[nvars * component_offset + nelems++] = gfs.at(n)(p.I);
            break;
          case int(point_type_t::bdry):
          case int(point_type_t::rest):
          case int(point_type_t::sync):
          case int(point_type_t::prol):
            // ignore this point
            break;
          default:
            assert(0);
          }
        });
    assert(nelems == nvars * component_sizes.at(level).at(component));
  });

  ierr = VecRestoreArray(vec, &vec_ptr);
  assert(!ierr);
}

void copy_PETSc_to_Cactus(
    Vec vec, const std::vector<int> &varinds,
    const std::vector<std::vector<int> > &component_offsets,
    const std::vector<std::vector<int> > &component_sizes) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }
  static const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  static const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  static const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  static const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);

  // Get vector from PETSc
  const CCTK_REAL *vec_ptr;
  ierr = VecGetArrayRead(vec, &vec_ptr);
  assert(!ierr);
  const CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy PETSc vector to Cactus
  CarpetX::active_levels->loop_parallel([&](const int patch, const int level,
                                            const int index,
                                            const int component,
                                            const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const CarpetX::GridDescBase grid(cctkGH);
    const auto &patchdata = CarpetX::ghext->patchdata.at(patch);
    const auto &leveldata = patchdata.leveldata.at(level);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    std::vector<Loop::GF3D2<CCTK_REAL> > gfs;
    gfs.reserve(nvars);
    for (int n = 0; n < nvars; ++n)
      gfs.emplace_back(layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                                    cctkGH, tl, varinds.at(n))));
    const int component_offset = component_offsets.at(level).at(component);
    int nelems = 0;
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
            for (int n = 0; n < nvars; ++n)
              gfs.at(n)(p.I) = petsc_ptr[nvars * component_offset + nelems++];
            break;
          case int(point_type_t::bdry):
          case int(point_type_t::rest):
          case int(point_type_t::sync):
          case int(point_type_t::prol):
            // ignore this point
            break;
          default:
            assert(0);
          }
        });
    assert(nelems == nvars * component_sizes.at(level).at(component));
    for (int n = 0; n < nvars; ++n)
      leveldata.groupdata.at(gis.at(n))->valid.at(tl).at(vis.at(n)).set_all(
          CarpetX::make_valid_int(),
          []() { return "PDESolver::copy_PETSc_to_Cactus"; });
  });

  ierr = VecRestoreArrayRead(vec, &vec_ptr);
  assert(!ierr);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void PDESolvers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Setup;
  DECLARE_CCTK_PARAMETERS;

  // Parse Cactus parameter into PETSc command line options
  /*static*/ std::vector<std::string> args;
  args.push_back("cactus");
  std::istringstream buf(petsc_options);
  std::copy(std::istream_iterator<std::string>(buf),
            std::istream_iterator<std::string>(), std::back_inserter(args));

  /*static*/ std::vector<char *> argptrs;
  for (auto &arg : args)
    argptrs.push_back(arg.data());
  argptrs.push_back(nullptr);

  /*static*/ int argc = args.size();
  /*static*/ char **argv = argptrs.data();
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, NULL);
  assert(!ierr);
}

namespace {
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *fun) {
  return (
      *static_cast<std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> *>(
          fun))(snes, x, f);
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat B, void *fun) {
  return (*static_cast<
          std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> *>(
      fun))(snes, x, J, B);
}
} // namespace

extern "C" void PDESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("PDESolvers_Solve: Solving on %d levels",
             CarpetX::ghext->num_levels());

  PetscErrorCode ierr;

  // Grid layout

  define_point_type();

  int npoints_local, npoints_global;
  int npoints_prolongated_local, npoints_prolongated_global;
  // component_offsets are process local
  std::vector<std::vector<int> > component_offsets, component_sizes;
  std::vector<std::vector<int> > component_prolongated_offsets,
      component_prolongated_sizes;
  csr_t Jp;
  enumerate_points(npoints_local, npoints_global, component_offsets,
                   component_sizes, npoints_prolongated_local,
                   npoints_prolongated_global, component_prolongated_offsets,
                   component_prolongated_sizes, Jp);

  // TODO: fix this
  const std::vector<int> solinds{CCTK_VarIndex("Poisson::sol")};
  const std::vector<int> resinds{CCTK_VarIndex("Poisson::res")};
  const int nvars = solinds.size();
  assert(int(resinds.size()) == nvars);

  // Nonlinear solver

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
  assert(!ierr);
  ierr = SNESSetFromOptions(snes);
  assert(!ierr);
  ierr = SNESSetAlwaysComputesFinalResidual(snes, PETSC_TRUE);
  assert(!ierr);

  // State vector and evaluation function

  Vec r;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, nvars * npoints_local,
                      nvars * npoints_global, &r);
  assert(!ierr);
  ierr = VecSetFromOptions(r);
  assert(!ierr);

  std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> evalf =
      [&](SNES snes, Vec x, Vec f) {
        copy_PETSc_to_Cactus(x, solinds, component_offsets, component_sizes);
        CallScheduleGroup(cctkGH, "PDESolvers_Residual");
        copy_Cactus_to_PETSc(f, resinds, component_offsets, component_sizes);
        return 0;
      };
  ierr = SNESSetFunction(snes, r, FormFunction, &evalf);
  assert(!ierr);

  // Matrix and Jacobian evaluation function

  Mat J;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, nvars * npoints_local,
                      nvars * npoints_local, nvars * npoints_global,
                      nvars * npoints_global, dnz, NULL, onz, NULL, &J);
  assert(!ierr);
  // ierr = MatCreate(PETSC_COMM_WORLD, &J);
  // assert(!ierr);
  // ierr = MatSetSizes(J, nvars * npoints_local, nvars * npoints_local,
  //                    nvars * npoints_global, nvars * npoints_global);
  // assert(!ierr);
  ierr = MatSetFromOptions(J);
  assert(!ierr);
  jacobians = std::make_optional<jacobians_t>();

  std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> evalJ =
      [&](SNES snes, Vec x, Mat J, Mat B) {
        copy_PETSc_to_Cactus(x, solinds, component_offsets, component_sizes);
        CallScheduleGroup(cctkGH, "PDESolvers_Jacobian");
        jacobians->define_matrix(Jp, J);
        jacobians->clear();
        if (false) {
          MatInfo info;
          PetscErrorCode ierr;
          ierr = MatGetInfo(J, MAT_GLOBAL_SUM, &info);
          assert(!ierr);
          CCTK_VINFO("Jacobian info: nz_allocated=%g nz_used=%g nz_unneeded=%g "
                     "memory=%g",
                     double(info.nz_allocated), double(info.nz_used),
                     double(info.nz_unneeded), double(info.memory));
        }
        return 0;
      };
  ierr = SNESSetJacobian(snes, J, J, FormJacobian, &evalJ);
  assert(!ierr);

  {
    double atol, rtol, stol;
    int maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);
    assert(!ierr);
    CCTK_VINFO("Settings: atol=%g, rtol=%g, stol=%g, maxit=%d, maxf=%d", atol,
               rtol, stol, maxit, maxf);
  }

  // Initialize solution

  Vec x;
  ierr = VecDuplicate(r, &x);
  assert(!ierr);
  copy_Cactus_to_PETSc(x, solinds, component_offsets, component_sizes);

  // Solve

  CCTK_INFO("Solve");
  ierr = SNESSolve(snes, NULL, x);
  assert(!ierr);

  {
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(snes, &reason);
    assert(!ierr);
    switch (reason) {
    case SNES_CONVERGED_FNORM_ABS:
    case SNES_CONVERGED_FNORM_RELATIVE:
    case SNES_CONVERGED_SNORM_RELATIVE:
    case SNES_CONVERGED_ITS:
      // SNES converged; do nothing
      CCTK_INFO("SNES iterations converged");
      break;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_LINEAR_SOLVE:
    case SNES_DIVERGED_FNORM_NAN:
    case SNES_DIVERGED_MAX_IT:
    case SNES_DIVERGED_LINE_SEARCH:
    case SNES_DIVERGED_INNER:
    case SNES_DIVERGED_LOCAL_MIN:
    case SNES_DIVERGED_DTOL:
    case SNES_DIVERGED_JACOBIAN_DOMAIN:
    case SNES_DIVERGED_TR_DELTA:
    case SNES_CONVERGED_ITERATING:
      CCTK_WARN(CCTK_WARN_ALERT, "****************************************");
      CCTK_WARN(CCTK_WARN_ALERT, "SNES iterations diverged");
      CCTK_WARN(CCTK_WARN_ALERT, "****************************************");
      break;
    default:
      CCTK_ERROR("Unknown SNES converged reason");
    }
  }

  {
    int niters;
    ierr = SNESGetIterationNumber(snes, &niters);
    assert(!ierr);
    int liters;
    ierr = SNESGetLinearSolveIterations(snes, &liters);
    assert(!ierr);
    CCTK_VINFO("Iterations: %d Newton, %d linear", niters, liters);
  }

  // Extract solution

  copy_PETSc_to_Cactus(x, solinds, component_offsets, component_sizes);
  CallScheduleGroup(cctkGH, "PDESolvers_Residual");
  {
    const int nvars = resinds.size();
    std::vector<int> gis, vis;
    gis.reserve(nvars);
    vis.reserve(nvars);
    for (int vn : resinds) {
      assert(vn >= 0);
      const int gi = CCTK_GroupIndexFromVarI(vn);
      assert(gi >= 0);
      gis.push_back(gi);
      const int v0 = CCTK_FirstVarIndexI(gi);
      assert(v0 >= 0);
      const int vi = vn - v0;
      assert(vi >= 0);
      vis.push_back(vi);
    }
    const CarpetX::reduction<CCTK_REAL, 3> res_norm =
        CarpetX::reduce(gis.at(0), vis.at(0), tl);
    std::ostringstream buf;
    buf << res_norm;
    CCTK_VINFO("Residual norm: %s", buf.str().c_str());
  }

  // Clean up

  jacobians.reset();
  VecDestroy(&x);
  VecDestroy(&r);
  MatDestroy(&J);
  SNESDestroy(&snes);
  CCTK_INFO("Done.");
}

extern "C" void PDESolvers_Shutdown(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Shutdown;
  PetscFinalize();
}

} // namespace PDESolvers
