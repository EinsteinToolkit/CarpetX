#include "multipatch.hxx"

#include <cctk_Parameters.h>

#include <cmath>

namespace MultiPatch {

namespace CubedSphere {

constexpr int ncells_i = 10;
constexpr int ncells_j = 10;
constexpr int ncells_k = 10;

constexpr int ncells_rad = 10;

Patch makePatch0(const PatchTransformations &pt) {
  const CCTK_REAL rmin = pt.cubed_sphere_rmin;
  const CCTK_REAL rmax = pt.cubed_sphere_rmax;
  Patch patch0;
  patch0.name = "cubed_sphere 0";
  patch0.ncells = {ncells_i, ncells_j, ncells_k};
  patch0.xmin = {-rmin, -rmin, -rmin};
  patch0.xmax = {+rmax, +rmax, +rmax};
  patch0.is_cartesian = true;
  patch0.faces = {{{false, 1}, {false, 2}, {false, 3}},
                  {{false, 4}, {false, 5}, {false, 6}}};
  return patch0;
}

Patch makePatch1(const PatchTransformations &pt) {
  const CCTK_REAL rmin = pt.cubed_sphere_rmin;
  const CCTK_REAL rmax = pt.cubed_sphere_rmax;
  Patch patch1;
  patch1.name = "cubed_sphere 1";
  patch1.ncells = {ncells_rad, ncells_j, ncells_k};
  patch1.xmin = {rmin, -1.0, -1.0};
  patch1.xmax = {rmax, +1.0, +1.0};
  patch1.is_cartesian = false;
  patch1.faces = {{{false, 0}, {false, 2}, {false, 3}},
                  {{true, -1}, {false, 5}, {false, 6}}};
  return patch1;
}

// Implementations
CCTK_DEVICE CCTK_HOST std_tuple<int, vec<CCTK_REAL, dim> >
global2local_impl(const PatchTransformations &pt,
                  const vec<CCTK_REAL, dim> &x) {
  const CCTK_REAL rmin = pt.cubed_sphere_rmin;
  if (abs(x(0)) <= rmin && abs(x(1)) <= rmin && abs(x(2)) <= rmin)
    return std_make_tuple(0, x);
  if (x(0) >= abs(x(1)) && x(0) >= abs(x(2)))
    return std_make_tuple(1, Arith::nan<vec<CCTK_REAL, dim> >()());
  if (x(0) <= abs(x(1)) && x(0) <= abs(x(2)))
    return std_make_tuple(4, Arith::nan<vec<CCTK_REAL, dim> >()());
  if (x(1) >= abs(x(0)) && x(1) >= abs(x(2)))
    return std_make_tuple(2, Arith::nan<vec<CCTK_REAL, dim> >()());
  if (x(1) <= abs(x(0)) && x(1) <= abs(x(2)))
    return std_make_tuple(5, Arith::nan<vec<CCTK_REAL, dim> >()());
  if (x(2) >= abs(x(0)) && x(2) >= abs(x(1)))
    return std_make_tuple(3, Arith::nan<vec<CCTK_REAL, dim> >()());
  if (x(2) <= abs(x(0)) && x(2) <= abs(x(1)))
    return std_make_tuple(6, Arith::nan<vec<CCTK_REAL, dim> >()());
  assert(0);
}

CCTK_DEVICE
    CCTK_HOST std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
                        vec<smat<CCTK_REAL, dim>, dim> >
    d2local_dglobal2_impl(const PatchTransformations &pt, int patch,
                          const vec<CCTK_REAL, dim> &a) {
  const CCTK_REAL rmin = pt.cubed_sphere_rmin;
  const CCTK_REAL rmax = pt.cubed_sphere_rmax;
  switch (patch) {
  case 0:
    return std_make_tuple(a, zero<vec<vec<CCTK_REAL, dim>, dim> >()(),
                          zero<vec<smat<CCTK_REAL, dim>, dim> >()());
  case 1: {
    const auto rho = a(0);
    const auto mu_y = a(1);
    const auto mu_z = a(2);

    // Conditions:
    // - \mu \in [-1; +1]
    // - for rho = rmin, \mu lie on a cube at rmin
    // - for rho = rmax, \mu lie on a sphere at rmax

    // for rho = rmin:
    const auto x0 = rmin;
    const auto y0 = rmin * mu_y;
    const auto z0 = rmin * mu_z;
    const auto r0 = rmin * sqrt(pow2(x0) + pow2(y0) + pow2(z0));
    // y = r \cos \theta
    // z = r \cos \theta
    const auto theta_y = acos(y0 / r0);
    const auto theta_z = acos(z0 / r0);
    // for rho = rmax:
    const auto r1 = rmax;
    const auto r = lincom(rmin, r0, rmax, r1, rho);

    const auto y = r * cos(theta_y);
    const auto z = r * cos(theta_z);
    const auto x = sqrt(pow2(r) - pow2(y) - pow2(z));

    return std_make_tuple(vec<CCTK_REAL, dim>{x, y, z},
                          Arith::nan<vec<vec<CCTK_REAL, dim>, dim> >()(),
                          Arith::nan<vec<smat<CCTK_REAL, dim>, dim> >()());
  }
  // case 2:
  //   ...;
  // case 3:
  //   ...;
  // case 4:
  //   ...;
  // case 5:
  //   ...;
  // case 6:
  //   ...;
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

} // namespace CubedSphere

PatchSystem SetupCubedSphere() {
  PatchTransformations pt;
  pt.global2local = &CubedSphere::global2local;
  pt.local2global = &CubedSphere::local2global;
  pt.dlocal_dglobal = &CubedSphere::dlocal_dglobal;
  pt.d2local_dglobal2 = &CubedSphere::d2local_dglobal2;
  pt.global2local_device = &CubedSphere::global2local_device;
  pt.local2global_device = &CubedSphere::local2global_device;
  pt.dlocal_dglobal_device = &CubedSphere::dlocal_dglobal_device;
  pt.d2local_dglobal2_device = &CubedSphere::d2local_dglobal2_device;

  return PatchSystem("Cubed Sphere",
                     std::vector<Patch>{CubedSphere::makePatch0(pt),
                                        CubedSphere::makePatch1(pt)},
                     std::move(pt));
}

} // namespace MultiPatch
