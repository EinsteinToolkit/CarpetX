#include "loop.hxx"

namespace Loop {

std::ostream &operator<<(std::ostream &os, const where_t where) {
  switch (where) {
  case where_t::everywhere:
    return os << "everywhere";
  case where_t::interior:
    return os << "interior";
  case where_t::boundary:
    return os << "boundary";
#if 0
  case where_t::ghosts_inclusive:
    return os << "ghosts_inclusive";
#endif
  case where_t::ghosts:
    return os << "ghosts";
  default:
    assert(0);
  }
}

std::ostream &operator<<(std::ostream &os, const PointDesc &p) {
  return os << "PointDesc{"
            << "level:" << p.level << ", "
            << "patch:" << p.patch << ", "
            << "block:" << p.block << ", "
            << "I:" << p.I << ", "
            << "iter:" << p.iter << ", "
            << "NI:" << p.NI << ", "
            << "I0:" << p.I0 << ", "
            << "BI:" << p.BI << ", "
            << "X:" << p.X << ", "
            << "DX:" << p.DX << ", "
            << "imin,imax:{" << p.imin << "," << p.imax << "}}";
}

GridDescBase::GridDescBase() {}

GridDescBase::GridDescBase(const cGH *restrict cctkGH) {
  level = cctkGH->cctk_level;
  patch = cctkGH->cctk_patch;
  block = cctkGH->cctk_block;

  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_gsh[d] != undefined);
    gsh[d] = cctkGH->cctk_gsh[d];
  }
  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_lbnd[d] != undefined);
    assert(cctkGH->cctk_ubnd[d] != undefined);
    lbnd[d] = cctkGH->cctk_lbnd[d];
    ubnd[d] = cctkGH->cctk_ubnd[d];
  }
  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_tile_min[d] != undefined);
    assert(cctkGH->cctk_tile_max[d] != undefined);
    tmin[d] = cctkGH->cctk_tile_min[d];
    tmax[d] = cctkGH->cctk_tile_max[d];
  }
  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_lsh[d] != undefined);
    lsh[d] = cctkGH->cctk_lsh[d];
  }
  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_ash[d] != undefined);
    ash[d] = cctkGH->cctk_ash[d];
  }
  for (int f = 0; f < 2; ++f) {
    for (int d = 0; d < dim; ++d) {
      assert(cctkGH->cctk_bbox[2 * d + f] != undefined);
      bbox[f][d] = cctkGH->cctk_bbox[2 * d + f];
    }
  }
  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_nghostzones[d] != undefined);
    nghostzones[d] = cctkGH->cctk_nghostzones[d];
  }

  for (int d = 0; d < dim; ++d) {
    assert(cctkGH->cctk_levfac[d] != undefined);
    assert(cctkGH->cctk_levoff[d] != undefined);
    assert(cctkGH->cctk_levoffdenom[d] != 0);
    dx[d] = cctkGH->cctk_delta_space[d] / cctkGH->cctk_levfac[d];
    x0[d] = cctkGH->cctk_origin_space[d] +
            dx[d] * cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d];
  }
}

std::ostream &operator<<(std::ostream &os, const GridDescBase &grid) {
  // Convert to vertex-centred boundaries
  const auto x0 = grid.x0 - grid.dx / 2;
  const auto x1 = x0 + (grid.gsh - 1) * grid.dx;
  return os << "GridDescBase{"
            << "level=" << grid.level << ",patch=" << grid.patch
            << ",block=" << grid.block << ",gsh=" << grid.gsh
            << ",lbnd=" << grid.lbnd << ",ubnd=" << grid.ubnd
            << ",lsh=" << grid.lsh << ",ash=" << grid.ash
            << ",bbox=" << grid.bbox << ",nghostzones=" << grid.nghostzones
            << ",tmin=" << grid.tmin << ",tmax=" << grid.tmax
            << ",x0[vc]=" << x0 << ",x1[vc]=" << x1 << ",dx=" << grid.dx << "}";
}

} // namespace Loop
