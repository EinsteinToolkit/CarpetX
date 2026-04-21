#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

const int dim = 3; // FIXME: used by reduction.hxx

#include "loop_device.hxx" //  FIXME: gives access to amrex:: used by reduction...
#include "reduction.hxx"

#include <cassert>
#include <cmath>

void CenterOfMass_ComputeIntegrands(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CenterOfMass_ComputeIntegrands;

  const CCTK_REAL x0 = 0.5;
  const CCTK_REAL r0 = 0.25;

  grid.loop_all_device<0, 0, 0>(grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // just a contant density slab of radius r0 at location x0 + cctk_time
    rho(p.I) = fabs(vcoordx(p.I) - (x0 + cctk_time)) < r0 ? 1. : 0.;
    x_rho(p.I) = rho(p.I) * vcoordx(p.I);
  });
}

void CenterOfMass_ComputeIntegrals(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CenterOfMass_ComputeIntegrals;

  const int vi_rho = CCTK_VarIndex("CenterOfMass::rho");
  assert(vi_rho >= 0);
  const int gi_rho = CCTK_GroupIndexFromVarI(vi_rho);
  assert(gi_rho >= 0);
  const int first_vi_rho = CCTK_FirstVarIndexI(gi_rho);
  assert(first_vi_rho >= 0);

  const int vi_x_rho = CCTK_VarIndex("CenterOfMass::x_rho");
  assert(vi_x_rho >= 0);
  const int gi_x_rho = CCTK_GroupIndexFromVarI(vi_x_rho);
  assert(gi_x_rho >= 0);
  const int first_vi_x_rho = CCTK_FirstVarIndexI(gi_x_rho);
  assert(first_vi_x_rho >= 0);

  const int tl = 0;
  const CarpetX::reduction<CCTK_REAL, dim> red_rho = CarpetX::reduce(gi_rho, vi_rho - first_vi_rho, tl);
  const CarpetX::reduction<CCTK_REAL, dim> red_x_rho = CarpetX::reduce(gi_x_rho, vi_x_rho - first_vi_x_rho, tl);

  CCTK_VINFO("tims: %g, total mass: %g, dipole moment: %g, center of mass: %g", cctk_time, red_rho.sum, red_x_rho.sum, red_x_rho.sum / red_rho.sum);
}
