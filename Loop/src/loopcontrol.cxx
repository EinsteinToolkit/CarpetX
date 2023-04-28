#include <cctk.h>

#include "loop.hxx"
#include "loopcontrol.h"

extern "C" GridDescBase_t LC_CreateGridDesc(const cGH *cctkGH) {
  static_assert(LC_DIM == Loop::dim);
  const Loop::GridDescBase lgrid(cctkGH);
  GridDescBase_t grid;

  for (int d = 0; d < LC_DIM; ++d) {
    grid.gsh[d] = lgrid.gsh[d];
    grid.lbnd[d] = lgrid.lbnd[d];
    grid.ubnd[d] = lgrid.ubnd[d];
    grid.lsh[d] = lgrid.lsh[d];
    grid.ash[d] = lgrid.ash[d];
    grid.bbox[2 * d + 0] = lgrid.bbox[0][d];
    grid.bbox[2 * d + 1] = lgrid.bbox[1][d];
    grid.nghostzones[d] = lgrid.nghostzones[d];
    grid.tmin[d] = lgrid.tmin[d];
    grid.tmax[d] = lgrid.tmax[d];
    grid.x0[d] = lgrid.x0[d];
    grid.dx[d] = lgrid.dx[d];
  }

  return grid;
}
