#include "multipatch.hxx"

#include <cctk_Parameters.h>

namespace MultiPatch {

namespace Cartesian {

Patch makePatch(const PatchTransformations &pt) {
  DECLARE_CCTK_PARAMETERS;
  const PatchFace outer_boundary{true, -1};
  Patch patch0;
  patch0.name = "cartesian";
  patch0.ncells = {cartesian_ncells_i, cartesian_ncells_j, cartesian_ncells_k};
  patch0.xmin = {cartesian_xmin, cartesian_ymin, cartesian_zmin};
  patch0.xmax = {cartesian_xmax, cartesian_ymax, cartesian_zmax};
  patch0.is_cartesian = true;
  patch0.faces = {{outer_boundary, outer_boundary, outer_boundary},
                  {outer_boundary, outer_boundary, outer_boundary}};
  return patch0;
}

// Implementations
CCTK_DEVICE CCTK_HOST std_tuple<int, vec<CCTK_REAL, dim> >
global2local_impl(const PatchTransformations &pt,
                  const vec<CCTK_REAL, dim> &x) {
  return std_make_tuple(0, x);
}

CCTK_DEVICE
    CCTK_HOST std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
                        vec<smat<CCTK_REAL, dim>, dim> >
    d2local_dglobal2_impl(const PatchTransformations &pt, int patch,
                          const vec<CCTK_REAL, dim> &a) {
  switch (patch) {
  case 0:
    return std_make_tuple(a, zero<vec<vec<CCTK_REAL, dim>, dim> >()(),
                          zero<vec<smat<CCTK_REAL, dim>, dim> >()());
  default:
    assert(0);
  }
}

CCTK_DEVICE
    CCTK_HOST std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> >
    dlocal_dglobal_impl(const PatchTransformations &pt, int patch,
                        const vec<CCTK_REAL, dim> &a) {
  const auto x_dx_ddx = d2local_dglobal2_impl(pt, patch, a);
  return std_make_tuple(std::get<0>(x_dx_ddx), std::get<1>(x_dx_ddx));
}

CCTK_DEVICE CCTK_HOST vec<CCTK_REAL, dim>
local2global_impl(const PatchTransformations &pt, int patch,
                  const vec<CCTK_REAL, dim> &a) {
  const auto x_dx = dlocal_dglobal_impl(pt, patch, a);
  return std::get<0>(x_dx);
}

// Host functions
std_tuple<int, vec<CCTK_REAL, dim> >
global2local(const PatchTransformations &pt, const vec<CCTK_REAL, dim> &x) {
  return global2local_impl(pt, x);
}
vec<CCTK_REAL, dim> local2global(const PatchTransformations &pt, int patch,
                                 const vec<CCTK_REAL, dim> &a) {
  return local2global_impl(pt, patch, a);
}
std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> >
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const vec<CCTK_REAL, dim> &a) {
  return dlocal_dglobal_impl(pt, patch, a);
}
std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
          vec<smat<CCTK_REAL, dim>, dim> >
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const vec<CCTK_REAL, dim> &a) {
  return d2local_dglobal2_impl(pt, patch, a);
}

// Device functions
CCTK_DEVICE std_tuple<int, vec<CCTK_REAL, dim> >
global2local_device(const PatchTransformations &pt,
                    const vec<CCTK_REAL, dim> &x) {
  return global2local_impl(pt, x);
}
CCTK_DEVICE vec<CCTK_REAL, dim>
local2global_device(const PatchTransformations &pt, int patch,
                    const vec<CCTK_REAL, dim> &a) {
  return local2global_impl(pt, patch, a);
}
CCTK_DEVICE
std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> >
dlocal_dglobal_device(const PatchTransformations &pt, int patch,
                      const vec<CCTK_REAL, dim> &a) {
  return dlocal_dglobal_impl(pt, patch, a);
}
CCTK_DEVICE
std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
          vec<smat<CCTK_REAL, dim>, dim> >
d2local_dglobal2_device(const PatchTransformations &pt, int patch,
                        const vec<CCTK_REAL, dim> &a) {
  return d2local_dglobal2_impl(pt, patch, a);
}

} // namespace Cartesian

PatchSystem SetupCartesian() {
  PatchTransformations pt;
  pt.global2local = &Cartesian::global2local;
  pt.local2global = &Cartesian::local2global;
  pt.dlocal_dglobal = &Cartesian::dlocal_dglobal;
  pt.d2local_dglobal2 = &Cartesian::d2local_dglobal2;
  pt.global2local_device = &Cartesian::global2local_device;
  pt.local2global_device = &Cartesian::local2global_device;
  pt.dlocal_dglobal_device = &Cartesian::dlocal_dglobal_device;
  pt.d2local_dglobal2_device = &Cartesian::d2local_dglobal2_device;

  return PatchSystem("Cartesian", std::vector<Patch>{Cartesian::makePatch(pt)},
                     std::move(pt));
}

} // namespace MultiPatch
