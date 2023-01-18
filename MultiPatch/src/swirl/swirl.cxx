#include "multipatch.hxx"

#include <sum.hxx>

#include <cctk_Parameters.h>

#include <cmath>

namespace MultiPatch {

namespace Swirl {

CCTK_DEVICE CCTK_HOST vec<CCTK_REAL, dim> zero3() {
  return zero<vec<CCTK_REAL, dim> >()();
}
CCTK_DEVICE CCTK_HOST smat<CCTK_REAL, dim> zero33() {
  return zero<smat<CCTK_REAL, dim> >()();
}
CCTK_DEVICE CCTK_HOST smat<CCTK_REAL, dim> delta33() {
  return one<smat<CCTK_REAL, dim> >()();
}

template <typename T> static CCTK_DEVICE CCTK_HOST T pow4(const T &x) {
  return pow2(pow2(x));
}

Patch makePatch(const PatchTransformations &pt) {
  const int ncells_i = pt.swirl_ncells_i;
  const int ncells_j = pt.swirl_ncells_j;
  const int ncells_k = pt.swirl_ncells_k;
  const PatchFace outer_boundary{true, -1};
  Patch patch0;
  patch0.name = "swirl";
  patch0.ncells = {ncells_i, ncells_j, ncells_k};
  patch0.xmin = {-1, -1, -1};
  patch0.xmax = {+1, +1, +1};
  patch0.is_cartesian = false;
  patch0.faces = {{outer_boundary, outer_boundary, outer_boundary},
                  {outer_boundary, outer_boundary, outer_boundary}};
  return patch0;
}

// Implementations
CCTK_DEVICE CCTK_HOST std_tuple<int, vec<CCTK_REAL, dim> >
global2local_impl(const PatchTransformations &pt,
                  const vec<CCTK_REAL, dim> &x) {
  using std::cbrt, std::cos, std::sin, std::sqrt;
  // alpha=0 at origin, alpha=1 at boundary
  const CCTK_REAL alpha = 1 - sqrt((pow2(x(0)) + pow2(x(1)) + pow2(x(2))) / 3);
  const CCTK_REAL cos_alpha = cbrt(pow2(cos(M_PI / 2 * alpha)));
  const CCTK_REAL sin_alpha = cbrt(pow2(sin(M_PI / 2 * alpha)));

  const vec<vec<CCTK_REAL, dim>, dim> Rinvelts = {
      {cos_alpha * cos_alpha, -cos_alpha * sin_alpha, sin_alpha * sin_alpha},
      {sin_alpha * sin_alpha, cos_alpha * cos_alpha, -cos_alpha * sin_alpha},
      {-cos_alpha * sin_alpha, sin_alpha * sin_alpha, cos_alpha * cos_alpha},
  };
  const mat<CCTK_REAL, dim> Rinv([&](int i, int j) { return Rinvelts(i)(j); });
  // Note: det Rinv = 1, i.e. |x| = |a|

  const vec<CCTK_REAL, dim> a([&](int i) {
    return sum<dim>([&](int j) { return Rinv(i, j) * x(j); });
  });

  return std_make_tuple(0, a);
}

CCTK_DEVICE
    CCTK_HOST std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim>,
                        vec<smat<CCTK_REAL, dim>, dim> >
    d2local_dglobal2_impl(const PatchTransformations &pt, int patch,
                          const vec<CCTK_REAL, dim> &a) {
  switch (patch) {

  case 0: {
    using std::cbrt, std::cos, std::sin, std::sqrt;
    // alpha=1 at origin, alpha=0 at boundary
    const CCTK_REAL alpha =
        1 - sqrt((pow2(a(0)) + pow2(a(1)) + pow2(a(2))) / 3);
    const CCTK_REAL cos_alpha = cbrt(pow2(cos(M_PI / 2 * alpha)));
    const CCTK_REAL sin_alpha = cbrt(pow2(sin(M_PI / 2 * alpha)));
    // x=(rotated a) for alpha=1 (origin), x=a for alpha=0 (boundary)

    const vec<vec<CCTK_REAL, dim>, dim> Relts = {
        {cos_alpha, sin_alpha, 0},
        {0, cos_alpha, sin_alpha},
        {sin_alpha, 0, cos_alpha},
    };
    const mat<CCTK_REAL, dim> R([&](int i, int j) { return Relts(i)(j); });
    // Note: det R = 1, i.e. |x| = |a|

    const vec<CCTK_REAL, dim> x(
        [&](int i) { return sum<dim>([&](int j) { return R(i, j) * a(j); }); });

    const vec<CCTK_REAL, dim> dalpha_da = a / (3 * (alpha - 1));
    const vec<CCTK_REAL, dim> dcos_alpha_da =
        -(M_PI / 3) * sin(M_PI / 2 * alpha) / cbrt(cos(M_PI / 2 * alpha)) *
        dalpha_da;
    const vec<CCTK_REAL, dim> dsin_alpha_da =
        (M_PI / 3) * cos(M_PI / 2 * alpha) / cbrt(sin(M_PI / 2 * alpha)) *
        dalpha_da;

    const vec<vec<vec<CCTK_REAL, dim>, dim>, dim> dRelts_da = {
        {dcos_alpha_da, dsin_alpha_da, zero3()},
        {zero3(), dcos_alpha_da, dsin_alpha_da},
        {dsin_alpha_da, zero3(), dcos_alpha_da},
    };
    const mat<vec<CCTK_REAL, dim>, dim> dR_da([&](int i, int j) {
      return vec<CCTK_REAL, dim>([&](int k) { return dRelts_da(i)(j)(k); });
    });

    const vec<vec<CCTK_REAL, dim>, dim> dxda([&](int i) {
      return vec<CCTK_REAL, dim>([&](int j) {
        return sum<dim>([&](int k) { return dR_da(i, k)(j) * a(k); }) + R(i, j);
      });
    });

    const smat<CCTK_REAL, dim> ddalpha_dada([&](int i, int j) {
      return (delta33()(i, j) * (alpha - 1) - a(i) * dalpha_da(j)) /
             (3 * pow2(alpha - 1));
    });
    const smat<CCTK_REAL, dim> ddcos_alpha_dada([&](int i, int j) {
      return -(pow2(M_PI) / 18) *
                 (3 * pow2(cos(M_PI / 2 * alpha)) +
                  pow2(sin(M_PI / 2 * alpha))) /
                 cbrt(pow4(cos(M_PI / 2 * alpha))) * dalpha_da(i) *
                 dalpha_da(j) -
             (M_PI / 3) * sin(M_PI / 2 * alpha) / cbrt(cos(M_PI / 2 * alpha)) *
                 ddalpha_dada(i, j);
    });
    const smat<CCTK_REAL, dim> ddsin_alpha_dada([&](int i, int j) {
      return -(pow2(M_PI) / 18) *
                 (pow2(cos(M_PI / 2 * alpha)) +
                  3 * pow2(sin(M_PI / 2 * alpha))) /
                 cbrt(pow4(sin(M_PI / 2 * alpha))) * dalpha_da(i) *
                 dalpha_da(j) +
             (M_PI / 3) * cos(M_PI / 2 * alpha) / cbrt(sin(M_PI / 2 * alpha)) *
                 ddalpha_dada(i, j);
    });

    const vec<vec<smat<CCTK_REAL, dim>, dim>, dim> ddRelts_dada = {
        {ddcos_alpha_dada, ddsin_alpha_dada, zero33()},
        {zero33(), ddcos_alpha_dada, ddsin_alpha_dada},
        {ddsin_alpha_dada, zero33(), ddcos_alpha_dada},
    };
    const mat<smat<CCTK_REAL, dim>, dim> ddR_dada([&](int i, int j) {
      return smat<CCTK_REAL, dim>(
          [&](int k, int l) { return ddRelts_dada(i)(j)(k, l); });
    });

    const vec<smat<CCTK_REAL, dim>, dim> ddxdada([&](int i) {
      return smat<CCTK_REAL, dim>([&](int j, int k) {
        return sum<dim>([&](int l) { return ddR_dada(i, l)(j, k) * a(l); }) +
               dR_da(i, k)(j) + dR_da(i, j)(k);
      });
    });

    return std_make_tuple(x, dxda, ddxdada);
  }

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
} // namespace Swirl

PatchSystem SetupSwirl() {
  PatchTransformations pt;
  pt.global2local = &Swirl::global2local;
  pt.local2global = &Swirl::local2global;
  pt.dlocal_dglobal = &Swirl::dlocal_dglobal;
  pt.d2local_dglobal2 = &Swirl::d2local_dglobal2;
  pt.global2local_device = &Swirl::global2local_device;
  pt.local2global_device = &Swirl::local2global_device;
  pt.dlocal_dglobal_device = &Swirl::dlocal_dglobal_device;
  pt.d2local_dglobal2_device = &Swirl::d2local_dglobal2_device;

  return PatchSystem("Swirl", std::vector<Patch>{Swirl::makePatch(pt)},
                     std::move(pt));
}

} // namespace MultiPatch
