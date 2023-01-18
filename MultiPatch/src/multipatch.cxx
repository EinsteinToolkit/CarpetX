#include "multipatch.hxx"

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace MultiPatch {

PatchTransformations::PatchTransformations()
    : // Cartesian
      cartesian_xmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_xmax;
      }()),
      cartesian_xmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_xmin;
      }()),
      cartesian_ymax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ymax;
      }()),
      cartesian_ymin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ymin;
      }()),
      cartesian_zmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_zmax;
      }()),
      cartesian_zmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_zmin;
      }()),
      cartesian_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_i;
      }()),
      cartesian_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_j;
      }()),
      cartesian_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return cartesian_ncells_k;
      }()),

      // Cubed sphere
      cubed_sphere_rmin([] {
        DECLARE_CCTK_PARAMETERS;
        return cubed_sphere_rmin;
      }()),
      cubed_sphere_rmax([] {
        DECLARE_CCTK_PARAMETERS;
        return cubed_sphere_rmax;
      }()),

      // Swirl
      swirl_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_i;
      }()),
      swirl_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_j;
      }()),
      swirl_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return swirl_ncells_k;
      }()),

      // Cake
      cake_outer_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_outer_boundary_radius;
      }()),
      cake_inner_boundary_radius([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_inner_boundary_radius;
      }()),
      cake_cartesian_ncells_i([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_i;
      }()),
      cake_cartesian_ncells_j([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_j;
      }()),
      cake_cartesian_ncells_k([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_cartesian_ncells_k;
      }()),
      cake_angular_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_angular_cells;
      }()),
      cake_radial_cells([] {
        DECLARE_CCTK_PARAMETERS;
        return cake_radial_cells;
      }()) {}

std::unique_ptr<PatchSystem> the_patch_system;

namespace {
template <typename T, int D>
CCTK_DEVICE CCTK_HOST inline T det(const vec<vec<T, D>, D> &A) {
  return A(0)(0) * (A(1)(1) * A(2)(2) - A(1)(2) * A(2)(1)) +
         A(0)(1) * (A(1)(2) * A(2)(0) - A(1)(0) * A(2)(2)) +
         A(0)(2) * (A(1)(0) * A(2)(1) - A(1)(1) * A(2)(0));
}
} // namespace

// Aliased functions

extern "C" CCTK_INT
MultiPatch1_GetSystemSpecification(CCTK_INT *restrict const npatches) {
  *npatches = the_patch_system->num_patches();
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetPatchSpecification(
    const CCTK_INT ipatch, const CCTK_INT size, CCTK_INT *restrict const ncells,
    CCTK_REAL *restrict const xmin, CCTK_REAL *restrict const xmax) {
  assert(ipatch >= 0 && ipatch < the_patch_system->num_patches());
  assert(size == dim);
  const Patch &patch = the_patch_system->patches.at(ipatch);
  for (int d = 0; d < dim; ++d) {
    ncells[d] = patch.ncells[d];
    xmin[d] = patch.xmin[d];
    xmax[d] = patch.xmax[d];
  }
  return 0;
}

extern "C" CCTK_INT MultiPatch1_GetBoundarySpecification2(
    const CCTK_INT ipatch, const CCTK_INT size,
    CCTK_INT *restrict const is_interpatch_boundary) {
  assert(ipatch >= 0 && ipatch < the_patch_system->num_patches());
  assert(size == 2 * dim);
  const Patch &patch = the_patch_system->patches.at(ipatch);
  for (int f = 0; f < 2; ++f)
    for (int d = 0; d < dim; ++d)
      is_interpatch_boundary[2 * d + f] = !patch.faces[f][d].is_outer_boundary;
  return 0;
}

extern "C" void MultiPatch1_GlobalToLocal(
    const CCTK_INT npoints, const CCTK_REAL *restrict const globalsx,
    const CCTK_REAL *restrict const globalsy,
    const CCTK_REAL *restrict const globalsz, CCTK_INT *restrict const patches,
    CCTK_REAL *restrict const localsx, CCTK_REAL *restrict const localsy,
    CCTK_REAL *restrict const localsz) {
  const PatchTransformations &transformations =
      the_patch_system->transformations;
  const auto &global2local = *transformations.global2local;

  for (int n = 0; n < npoints; ++n) {
    const vec<CCTK_REAL, dim> x{globalsx[n], globalsy[n], globalsz[n]};
    const auto patch_a = global2local(transformations, x);
    const auto patch = std::get<0>(patch_a);
    const auto a = std::get<1>(patch_a);
    patches[n] = patch;
    localsx[n] = a(0);
    localsy[n] = a(1);
    localsz[n] = a(2);
  }
}

// Scheduled functions

extern "C" int MultiPatch_Setup() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(patch_system, "none"))
    the_patch_system = nullptr;
  else if (CCTK_EQUALS(patch_system, "Cartesian"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCartesian());
  else if (CCTK_EQUALS(patch_system, "Cubed sphere"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCubedSphere());
  else if (CCTK_EQUALS(patch_system, "Swirl"))
    the_patch_system = std::make_unique<PatchSystem>(SetupSwirl());
  else if (CCTK_EQUALS(patch_system, "Cake"))
    the_patch_system = std::make_unique<PatchSystem>(SetupCake());
  else
    CCTK_VERROR("Unknown patch system \"%s\"", patch_system);

  return 0;
}

extern "C" void MultiPatch_Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultiPatch_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  // const Loop::GridDescBase grid(cctkGH);
  const Loop::GridDescBaseDevice grid(cctkGH);

  const std::array<int, dim> indextype_vc = {0, 0, 0};
  const std::array<int, dim> indextype_cc = {1, 1, 1};
  const Loop::GF3D2layout layout_vc(cctkGH, indextype_vc);
  const Loop::GF3D2layout layout_cc(cctkGH, indextype_cc);

  const Loop::GF3D2<CCTK_REAL> gf_vcoordx(layout_vc, vcoordx);
  const Loop::GF3D2<CCTK_REAL> gf_vcoordy(layout_vc, vcoordy);
  const Loop::GF3D2<CCTK_REAL> gf_vcoordz(layout_vc, vcoordz);

  const Loop::GF3D2<CCTK_REAL> gf_ccoordx(layout_cc, ccoordx);
  const Loop::GF3D2<CCTK_REAL> gf_ccoordy(layout_cc, ccoordy);
  const Loop::GF3D2<CCTK_REAL> gf_ccoordz(layout_cc, ccoordz);

  const Loop::GF3D2<CCTK_REAL> gf_cvol(layout_cc, cvol);

  // Vertex centered Jacobian matrix
  const Loop::GF3D2<CCTK_REAL> gf_vJ_da_dx(layout_vc, vJ_da_dx);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_da_dy(layout_vc, vJ_da_dy);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_da_dz(layout_vc, vJ_da_dz);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_db_dx(layout_vc, vJ_db_dx);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_db_dy(layout_vc, vJ_db_dy);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_db_dz(layout_vc, vJ_db_dz);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_dc_dx(layout_vc, vJ_dc_dx);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_dc_dy(layout_vc, vJ_dc_dy);
  const Loop::GF3D2<CCTK_REAL> gf_vJ_dc_dz(layout_vc, vJ_dc_dz);

  // Vertex centered Jacobian matrix derivative
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dxdx(layout_vc, vdJ_d2a_dxdx);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dxdy(layout_vc, vdJ_d2a_dxdy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dxdz(layout_vc, vdJ_d2a_dxdz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dydy(layout_vc, vdJ_d2a_dydy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dydz(layout_vc, vdJ_d2a_dydz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2a_dzdz(layout_vc, vdJ_d2a_dzdz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dxdx(layout_vc, vdJ_d2b_dxdx);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dxdy(layout_vc, vdJ_d2b_dxdy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dxdz(layout_vc, vdJ_d2b_dxdz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dydy(layout_vc, vdJ_d2b_dydy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dydz(layout_vc, vdJ_d2b_dydz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2b_dzdz(layout_vc, vdJ_d2b_dzdz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dxdx(layout_vc, vdJ_d2c_dxdx);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dxdy(layout_vc, vdJ_d2c_dxdy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dxdz(layout_vc, vdJ_d2c_dxdz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dydy(layout_vc, vdJ_d2c_dydy);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dydz(layout_vc, vdJ_d2c_dydz);
  const Loop::GF3D2<CCTK_REAL> gf_vdJ_d2c_dzdz(layout_vc, vdJ_d2c_dzdz);

  const PatchTransformations pt = the_patch_system->transformations;
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones,
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index(layout_vc, p.I);
        const vec<CCTK_REAL, dim> a = {p.x, p.y, p.z};

        const auto d2J_tuple = pt.d2local_dglobal2(pt, cctk_patch, a);
        const auto x = std::get<0>(d2J_tuple);
        const auto J = std::get<1>(d2J_tuple);
        const auto dJ = std::get<2>(d2J_tuple);

        gf_vcoordx(index) = x(0);
        gf_vcoordy(index) = x(1);
        gf_vcoordz(index) = x(2);

        gf_vJ_da_dx(index) = J(0)(0);
        gf_vJ_da_dy(index) = J(0)(1);
        gf_vJ_da_dz(index) = J(0)(2);
        gf_vJ_db_dx(index) = J(1)(0);
        gf_vJ_db_dy(index) = J(1)(1);
        gf_vJ_db_dz(index) = J(1)(2);
        gf_vJ_dc_dx(index) = J(2)(0);
        gf_vJ_dc_dy(index) = J(2)(1);
        gf_vJ_dc_dz(index) = J(2)(2);

        gf_vdJ_d2a_dxdx(index) = dJ(0)(0, 0);
        gf_vdJ_d2a_dxdy(index) = dJ(0)(0, 1);
        gf_vdJ_d2a_dxdz(index) = dJ(0)(0, 2);
        gf_vdJ_d2a_dydy(index) = dJ(0)(1, 1);
        gf_vdJ_d2a_dydz(index) = dJ(0)(1, 2);
        gf_vdJ_d2a_dzdz(index) = dJ(0)(2, 2);
        gf_vdJ_d2b_dxdx(index) = dJ(1)(0, 0);
        gf_vdJ_d2b_dxdy(index) = dJ(1)(0, 1);
        gf_vdJ_d2b_dxdz(index) = dJ(1)(0, 2);
        gf_vdJ_d2b_dydy(index) = dJ(1)(1, 1);
        gf_vdJ_d2b_dydz(index) = dJ(1)(1, 2);
        gf_vdJ_d2b_dzdz(index) = dJ(1)(2, 2);
        gf_vdJ_d2c_dxdx(index) = dJ(2)(0, 0);
        gf_vdJ_d2c_dxdy(index) = dJ(2)(0, 1);
        gf_vdJ_d2c_dxdz(index) = dJ(2)(0, 2);
        gf_vdJ_d2c_dydy(index) = dJ(2)(1, 1);
        gf_vdJ_d2c_dydz(index) = dJ(2)(1, 2);
        gf_vdJ_d2c_dzdz(index) = dJ(2)(2, 2);
      });

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] ARITH_DEVICE(const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index(layout_cc, p.I);
        const vec<CCTK_REAL, dim> a = {p.x, p.y, p.z};
        const std_tuple<vec<CCTK_REAL, dim>, vec<vec<CCTK_REAL, dim>, dim> >
            x_dadx = pt.dlocal_dglobal(pt, cctk_patch, a);
        const vec<CCTK_REAL, dim> &x = std::get<0>(x_dadx);
        const vec<vec<CCTK_REAL, dim>, dim> &dadx = std::get<1>(x_dadx);
        const CCTK_REAL det_dadx = det(dadx);
        const CCTK_REAL vol =
            (p.dx * p.dy * p.dz) *
            det_dadx; // TODO: Is it det_dadx or sqrt(det_dadx)? As far as we
                      // know no sqrt is necessary. Note that det_dadx is
                      // negative in the cake patch, so that causes NaNs
        gf_ccoordx(index) = x(0);
        gf_ccoordy(index) = x(1);
        gf_ccoordz(index) = x(2);
        gf_cvol(index) = vol;
      });
}

/**
 * TODO: Fill with more parameter checks, if appropriate
 */
extern "C" void MultiPatch_Check_Parameters(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultiPatch_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  if (!(cake_outer_boundary_radius > 4 * cake_inner_boundary_radius))
    CCTK_PARAMWARN(
        "Make sure that outer boundary radius is larger than 4 times "
        "the inner boudary radius.");
}

} // namespace MultiPatch
