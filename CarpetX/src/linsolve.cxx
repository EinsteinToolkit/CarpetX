#include "driver.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <cmath>
#include <vector>

namespace CarpetX {

extern "C" void CarpetX_SolvePoisson(const CCTK_INT gi_sol,
                                     const CCTK_INT gi_rhs,
                                     const CCTK_INT gi_res,
                                     const CCTK_REAL reltol,
                                     const CCTK_REAL abstol,
                                     CCTK_REAL *restrict const res_initial,
                                     CCTK_REAL *restrict const res_final) {
  assert(gi_rhs >= 0);
  assert(gi_sol >= 0);
  const bool have_res = gi_res >= 0;
  if (have_res)
    assert(gi_res >= 0);
  const int tl = 0;
  const int vi = 0;

  // Create operator

  assert(ghext->patchdata.size() == 1);
  const int patch = 0;
  const auto &patchdata = ghext->patchdata.at(patch);
  amrex::Vector<amrex::Geometry> geoms(patchdata.leveldata.size());
  amrex::Vector<amrex::BoxArray> grids(patchdata.leveldata.size());
  amrex::Vector<amrex::DistributionMapping> dmaps(patchdata.leveldata.size());
  for (int level = 0; level < int(patchdata.leveldata.size()); ++level) {
    geoms.at(level) = patchdata.amrcore->Geom(level);
    grids.at(level) = patchdata.amrcore->boxArray(level);
    dmaps.at(level) = patchdata.amrcore->DistributionMap(level);
  }

  amrex::MLNodeLaplacian mlnodelaplacian(geoms, grids, dmaps);

  // TODO: Look at Geometry for this
  mlnodelaplacian.setDomainBC(
      {amrex::LinOpBCType::Dirichlet, amrex::LinOpBCType::Dirichlet,
       amrex::LinOpBCType::Dirichlet},
      {amrex::LinOpBCType::Dirichlet, amrex::LinOpBCType::Dirichlet,
       amrex::LinOpBCType::Dirichlet});

  vector<amrex::MultiFab> sigmas(patchdata.leveldata.size());
  for (int level = 0; level < int(patchdata.leveldata.size()); ++level) {
    auto &sigma = sigmas.at(level);
    sigma.define(patchdata.amrcore->boxArray(level),
                 patchdata.amrcore->DistributionMap(level), 1, 0);
    sigma.setVal(1.0);
    mlnodelaplacian.setSigma(level, sigma);
  }

  mlnodelaplacian.setVerbose(10);

  // Create solver

  amrex::MLMG mlmg(mlnodelaplacian);

  amrex::Vector<amrex::MultiFab *> ress(patchdata.leveldata.size());
  amrex::Vector<amrex::MultiFab *> sols(patchdata.leveldata.size());
  amrex::Vector<const amrex::MultiFab *> rhss(patchdata.leveldata.size());
  for (int level = 0; level < int(patchdata.leveldata.size()); ++level) {
    const auto &restrict leveldata = patchdata.leveldata.at(level);
    const auto &restrict groupdata_rhs = *leveldata.groupdata.at(gi_rhs);
    rhss.at(level) = groupdata_rhs.mfab.at(tl).get();
    const auto &restrict groupdata_sol = *leveldata.groupdata.at(gi_sol);
    sols.at(level) = groupdata_sol.mfab.at(tl).get();
    if (have_res) {
      const auto &restrict groupdata_res = *leveldata.groupdata.at(gi_res);
      ress.at(level) = groupdata_res.mfab.at(tl).get();
    }
  }

  mlmg.setVerbose(10);
  mlmg.setBottomVerbose(10);

  // Solve

  if (have_res) {
    mlmg.compResidual(ress, sols, rhss);
    *res_initial = 0;
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      *res_initial =
          fmax(*res_initial, ress.at(level)->norminf(vi, 0, false, true));
  } else {
    *res_initial = NAN;
  }

#pragma omp critical
  {
    CCTK_VINFO("Before solving:");
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      CCTK_VINFO("norm_inf rhs[%d]: %g", level,
                 double(rhss.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      CCTK_VINFO("norm_inf sol[%d]: %g", level,
                 double(sols.at(level)->norminf(vi, 0, false, true)));
    if (have_res)
      for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
        CCTK_VINFO("norm_inf res[%d]: %g", level,
                   double(ress.at(level)->norminf(vi, 0, false, true)));
  }

  // const CCTK_REAL rtol = 0.0;
  // const CCTK_REAL atol = 1.0e-12;
  const CCTK_REAL maxerr = mlmg.solve(sols, rhss, reltol, abstol);
#pragma omp critical
  CCTK_VINFO("Solution error (norm_inf): %g", double(maxerr));

  if (have_res) {
    mlmg.compResidual(ress, sols, rhss);
    *res_final = 0;
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      *res_final =
          fmax(*res_final, ress.at(level)->norminf(vi, 0, false, true));
  } else {
    *res_final = NAN;
  }

#pragma omp critical
  {
    CCTK_VINFO("After solving:");
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      CCTK_VINFO("norm_inf rhs[%d]: %g", level,
                 double(rhss.at(level)->norminf(vi, 0, false, true)));
    for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
      CCTK_VINFO("norm_inf sol[%d]: %g", level,
                 double(sols.at(level)->norminf(vi, 0, false, true)));
    if (have_res)
      for (int level = 0; level < int(patchdata.leveldata.size()); ++level)
        CCTK_VINFO("norm_inf res[%d]: %g", level,
                   double(ress.at(level)->norminf(vi, 0, false, true)));
  }
}

} // namespace CarpetX
