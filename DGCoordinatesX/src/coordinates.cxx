#include "dg.hxx"

#include <div.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>

namespace DGCoordinatesX {
using namespace Arith;
using namespace DG;

using Loop::dim;

extern "C" void DGCoordinatesX_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_DGCoordinatesX_Setup;
  DECLARE_CCTK_PARAMETERS;

  // The `grid` variable stores vertex-centred metadata which are one grid point
  // larger than the cell-centred metadata we need.
  CCTK_VINFO("DG setup:");
  CCTK_VINFO("  dg_npoints=%d", dg_npoints);
  CCTK_VINFO("  nghostzones=[%d,%d,%d]", grid.nghostzones[0],
             grid.nghostzones[1], grid.nghostzones[2]);
  CCTK_VINFO("  gsh=[%d,%d,%d]", grid.gsh[0] - 1, grid.gsh[1] - 1,
             grid.gsh[2] - 1);
  CCTK_VINFO("  lbnd=[%d,%d,%d]", grid.lbnd[0], grid.lbnd[1], grid.lbnd[2]);
  CCTK_VINFO("  lsh=[%d,%d,%d]", grid.lsh[0] - 1, grid.lsh[1] - 1,
             grid.lsh[2] - 1);
  for (int d = 0; d < dim; ++d)
    assert(grid.nghostzones[d] == 1);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.gsh[d] - 1 - 2 * grid.nghostzones[d], dg_npoints) ==
           0);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.lbnd[d], dg_npoints) == 0);
  for (int d = 0; d < dim; ++d)
    assert(mod_floor(grid.lsh[d] - 1 - 2 * grid.nghostzones[d], dg_npoints) ==
           0);

  switch (dg_npoints) {

  case 2: {
    constexpr int N = 2;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx(p.I) = coords.coord[0];
      coordy(p.I) = coords.coord[1];
      coordz(p.I) = coords.coord[2];
      cvol(p.I) = coords.vol;
    });
    break;
  }

  case 3: {
    constexpr int N = 3;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx(p.I) = coords.coord[0];
      coordy(p.I) = coords.coord[1];
      coordz(p.I) = coords.coord[2];
      cvol(p.I) = coords.vol;
    });
    break;
  }

  case 4: {
    constexpr int N = 4;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx(p.I) = coords.coord[0];
      coordy(p.I) = coords.coord[1];
      coordz(p.I) = coords.coord[2];
      cvol(p.I) = coords.vol;
    });
    break;
  }

  case 8: {
    constexpr int N = 8;
    grid.loop_all<1, 1, 1>(grid.nghostzones, [&](const Loop::PointDesc &p) {
      const get_coords<CCTK_REAL, N> coords(grid, p);
      coordx(p.I) = coords.coord[0];
      coordy(p.I) = coords.coord[1];
      coordz(p.I) = coords.coord[2];
      cvol(p.I) = coords.vol;
    });
    break;
  }

  default:
    CCTK_ERROR("Unsupported DG element size");
  }
}

} // namespace DGCoordinatesX
