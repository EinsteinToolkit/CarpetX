#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <cmath>
#include <cassert>

#include "loop_device.hxx"

// TODO: this should not call itself "grid" since it only provides spacings and
// coords and would alos be used with non-Cartesian grids but Cartesian point
// labels

extern "C"
void GridX_SetRanges(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_GridX_SetRanges;

  static int level_run_on = -1;
  static int patch_run_on = -1;

  #pragma omp critical
  if(level_run_on == -1 || level_run_on > cctk_level ||
      patch_run_on == -1 || patch_run_on > cctk_patch) {
    assert((level_run_on == -1 && patch_run_on == -1) ||
           (level_run_on != -1 || patch_run_on != -1));
    *coarse_dx = cctk_delta_space[0];
    *coarse_dy = cctk_delta_space[1];
    *coarse_dz = cctk_delta_space[2];
    level_run_on = cctk_level;
    patch_run_on = cctk_patch;
  }
}

extern "C"
void GridX_SetCoordinates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTSX_GridX_SetCoordinates;

  grid.loop_all_device<0,0,0>(
    grid.nghostzones,
    [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      x(p.I) = p.x;
      y(p.I) = p.x;
      z(p.I) = p.x;
      r(p.I) = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    });
}
