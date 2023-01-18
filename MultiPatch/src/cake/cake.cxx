#include "multipatch.hxx"
#include "tests.hxx"
#include "cake.hxx"
#include "cake_inverse_core_derivs.hxx"
#include "cake_jacobians.hxx"

#include <cassert>
#include <cmath>
#include <string>

namespace MultiPatch {
namespace Cake {

/**
 * @brief Equation (18) of https://arxiv.org/pdf/gr-qc/0512001v1.pdf
 *
 * @param pt The patch data
 * @param c The local c coordinate
 * @return The result of equation (18) of
 * https://arxiv.org/pdf/gr-qc/0512001v1.pdf.
 */
CCTK_DEVICE CCTK_HOST inline CCTK_REAL paper_r(const PatchTransformations &pt,
                                               CCTK_REAL c) {
  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};

  return (r0 * (1 - c) + r1 * (1 + c)) / 2;
}

/**
 * @brief Equation (19) of https://arxiv.org/pdf/gr-qc/0512001v1.pdf
 *
 * @param pt The patch data.
 * @param a The local a coordinate
 * @param b The local b coordinate
 * @param c The local c coordinate
 * @return The result of equation (19) of
 * https://arxiv.org/pdf/gr-qc/0512001v1.pdf
 */
CCTK_DEVICE CCTK_HOST inline CCTK_REAL
paper_F(const PatchTransformations &pt, CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};
  const auto E{1 + a * a + b * b};
  const auto r{paper_r(pt, c)};

  return sqrt(((r1 - r) + (r - r0) * E) / (r1 - r0));
}

/**
 * @brief Core function of the cake local -> global coordinate transformations.
 *
 * @param pt The cake patch data.
 * @param a The a local coordinate, ranging from (-1, 1)
 * @param b The b local coordiante, ranging from (-1, 1)
 * @param c the c local coordiante, ranging from (-1, 1)
 * @return The value of the core.
 */
CCTK_DEVICE CCTK_HOST inline CCTK_REAL
local_to_global_cake_core(const PatchTransformations &pt, CCTK_REAL a,
                          CCTK_REAL b, CCTK_REAL c) {
  const auto r{paper_r(pt, c)};
  const auto F{paper_F(pt, a, b, c)};
  return r / F;
}

/**
 * @brief Core function of the cake global -> local coordinate transformations
 *
 * @param pt The patch data.
 * @param a local a coordinate in terms of global variables
 * @param b local b coordinate in terms of global variables
 * @param var The signature variable of the patch (x for +-x, y for +-y, z for
 * +- z)
 * @return The value of the core.
 */
CCTK_DEVICE CCTK_HOST inline CCTK_REAL
global_to_local_cake_core(const PatchTransformations &pt, CCTK_REAL a,
                          CCTK_REAL b, CCTK_REAL var) {
  using std::sqrt;

  const auto r0{pt.cake_inner_boundary_radius};
  const auto r1{pt.cake_outer_boundary_radius};
  const auto E{1 + a * a + b * b};

  const auto part1{r0 * r0 - r1 * r1 + (E - 1) * var * var};
  const auto part2{4 * (r0 - r1) * (E * r0 - r1) * var * var +
                   (E - 1) * (E - 1) * var * var * var * var};
  const auto part3{(r0 - r1) * (r0 - r1)};

  return (part1 + sqrt(part2)) / part3;
}

/**
 * @brief The local to global coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable for passing to CarpetX.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A spatial vector containing the coordinate transformations.
 */
CCTK_DEVICE CCTK_HOST svec local2global(const PatchTransformations &pt,
                                        int patch, const svec &local_vars) {

  svec global_vars = {0.0, 0.0, 0.0};

  const auto a = local_vars(0);
  const auto b = local_vars(1);
  const auto c = local_vars(2);

  const auto core = local_to_global_cake_core(pt, a, b, c);

  switch (patch) {
  case static_cast<int>(patch_piece::cartesian):
    global_vars = local_vars;
    break;
  case static_cast<int>(patch_piece::plus_x):
    global_vars = {core, b * core, a * core};
    break;
  case static_cast<int>(patch_piece::minus_x):
    global_vars = {-core, -b * core, a * core};
    break;
  case static_cast<int>(patch_piece::plus_y):
    global_vars = {-b * core, core, a * core};
    break;
  case static_cast<int>(patch_piece::minus_y):
    global_vars = {b * core, -core, a * core};
    break;
  case static_cast<int>(patch_piece::plus_z):
    global_vars = {-a * core, b * core, core};
    break;
  case static_cast<int>(patch_piece::minus_z):
    global_vars = {a * core, b * core, -core};
    break;
  default:
#ifndef __CUDACC__
    CCTK_VERROR("No local -> global transformations available for patch %s",
                piece_name(static_cast<patch_piece>(patch)).c_str());
#else
    assert(0);
#endif
    break;
  }

  return global_vars;
}

/**
 * @brief The global to local coordinate transformation implementation.
 * This is where the actual coordinate transformations are performed but it is
 * also suitable to be passed to CarpetX without any wrappers.
 *
 * @param pt The patch data
 * @param global_vars The values of the local global (x, y, z)
 * @return A tuple containing the patch piece and the global coordinate triplet.
 */
CCTK_DEVICE CCTK_HOST std_tuple<int, svec>
global2local(const PatchTransformations &pt, const svec &global_vars) {
  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto piece = get_owner_patch(pt, global_vars);

  CCTK_REAL a{0}, b{0}, var{0};

  switch (static_cast<int>(piece)) {

  case static_cast<int>(patch_piece::cartesian):
    return std_make_tuple(static_cast<int>(piece), global_vars);

  case static_cast<int>(patch_piece::plus_x):
    a = z / x;
    b = y / x;
    var = x;
    break;

  case static_cast<int>(patch_piece::minus_x):
    a = -z / x;
    b = y / x;
    var = x;
    break;

  case static_cast<int>(patch_piece::plus_y):
    a = z / y;
    b = -x / y;
    var = y;
    break;

  case static_cast<int>(patch_piece::minus_y):
    a = -z / y;
    b = -x / y;
    var = y;
    break;

  case static_cast<int>(patch_piece::plus_z):
    a = -x / z;
    b = y / z;
    var = z;
    break;

  case static_cast<int>(patch_piece::minus_z):
    a = -x / z;
    b = -y / z;
    var = z;
    break;

  default:
#ifndef __CUDACC__
    throw std::runtime_error("error");
    CCTK_VERROR("At point (%f, %f, %f): No global -> local transformations "
                "available for patch %s.",
                x, y, z, piece_name(piece).c_str());
#else
    assert(0);
#endif
    break;
  }

  const CCTK_REAL c = global_to_local_cake_core(pt, a, b, var);

  return std_make_tuple(static_cast<int>(piece), svec{a, b, c});
}

/**
 * @brief This function computes the local to global coordinate
 * transformation, the jacobian and it's derivative. It can be passed directly
 * to CarpetX.
 *
 * @note The Jacobians are defined as
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 * dJ(i)(j,k) = $dJ^{i}_{j k} = \frac{d^2 a^i}{d x^j d x^k}
 * \right)$.
 *
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation,
 * the local to global jacobian matrix and it's derivative.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t, djac_t>
d2local_dglobal2(const PatchTransformations &pt, int patch,
                 const svec &local_vars) {

  auto local_to_global_result = pt.local2global(pt, patch, local_vars);
  auto jacobian_results = cake_jacs(pt, patch, local_to_global_result);

  return std_make_tuple(local_to_global_result, std::get<0>(jacobian_results),
                        std::get<1>(jacobian_results));
} // namespace Cake

/**
 * @brief This function computes the local to global coordinate transformation
 * and the jacobian. It can be passed directly to CarpetX.
 *
 * @note NThe Jacobians is defined as:
 * J(i)(j) = $J^{i}_{j} = \frac{d a^i}{d x^j}$.
 *
 * @param pt The patch data
 * @param patch The index of the patch to transform.
 * @param local_vars The values of the local variables (a,b,c)
 * @return A tuple containing the local to global coordinate transformation
 * and the local to global jacobian matrix.
 */
CCTK_DEVICE CCTK_HOST std_tuple<svec, jac_t>
dlocal_dglobal(const PatchTransformations &pt, int patch,
               const svec &local_vars) {

  const auto data = pt.d2local_dglobal2(pt, patch, local_vars);
  return std_make_tuple(std::get<0>(data), std::get<1>(data));
}

/**
 * Creates a cake patch piece
 *
 * @tparam p The piece of the patch to make.
 * @param pt The patch transformation object with patch data.
 * @return The constructed patch piece.
 */
template <patch_piece p> Patch make_patch(const PatchTransformations &pt) {
  Patch patch;

  patch.name = piece_name(p);

  // Basic configuration for a thornburg patch piece
  patch.ncells = {pt.cake_angular_cells, pt.cake_angular_cells,
                  pt.cake_radial_cells};

  patch.xmin = {-1.0, -1.0, -1.0};
  patch.xmax = {1.0, 1.0, 1.0};

  patch.is_cartesian = false;

  PatchFace co{false, static_cast<int>(patch_piece::cartesian)};
  PatchFace ex{true, static_cast<int>(patch_piece::exterior)};
  PatchFace px{false, static_cast<int>(patch_piece::plus_x)};
  PatchFace mx{false, static_cast<int>(patch_piece::minus_x)};
  PatchFace py{false, static_cast<int>(patch_piece::plus_y)};
  PatchFace my{false, static_cast<int>(patch_piece::minus_y)};
  PatchFace pz{false, static_cast<int>(patch_piece::plus_z)};
  PatchFace mz{false, static_cast<int>(patch_piece::minus_z)};

  if constexpr (p == patch_piece::cartesian) {
    patch.ncells = {pt.cake_cartesian_ncells_i, pt.cake_cartesian_ncells_j,
                    pt.cake_cartesian_ncells_k};

    patch.xmin = {-pt.cake_inner_boundary_radius,
                  -pt.cake_inner_boundary_radius,
                  -pt.cake_inner_boundary_radius};

    patch.xmax = {pt.cake_inner_boundary_radius, pt.cake_inner_boundary_radius,
                  pt.cake_inner_boundary_radius};

    patch.is_cartesian = true;

    patch.faces = {{mx, my, mz}, {px, py, pz}};

  } else if constexpr (p == patch_piece::plus_x) {
    patch.faces = {{mz, my, co}, {pz, py, ex}};

  } else if constexpr (p == patch_piece::minus_x) {
    patch.faces = {{mz, py, co}, {pz, my, ex}};

  } else if constexpr (p == patch_piece::plus_y) {
    patch.faces = {{mz, px, co}, {pz, mx, ex}};

  } else if constexpr (p == patch_piece::minus_y) {
    patch.faces = {{mz, mx, co}, {pz, px, ex}};

  } else if constexpr (p == patch_piece::plus_z) {
    patch.faces = {{px, my, co}, {mx, py, ex}};

  } else if constexpr (p == patch_piece::minus_z) {
    patch.faces = {{mx, my, co}, {px, py, ex}};
  }

  return patch;
}

} // namespace Cake

/**
 * TODO: Add correct host/device annotations for the functions. This work as is
 * if not compiling with the CUDA compiler Creates a Cake patch system
 *
 * @return A PatchSystem object with Cake data and functions
 */
PatchSystem SetupCake() {
  PatchTransformations pt;
  pt.global2local = &Cake::global2local;
  pt.local2global = &Cake::local2global;
  pt.dlocal_dglobal = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2 = &Cake::d2local_dglobal2;
  pt.global2local_device = &Cake::global2local;
  pt.local2global_device = &Cake::local2global;
  pt.dlocal_dglobal_device = &Cake::dlocal_dglobal;
  pt.d2local_dglobal2_device = &Cake::d2local_dglobal2;

  const auto patches =
      std::vector<Patch>{Cake::make_patch<Cake::patch_piece::cartesian>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_x>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_y>(pt),
                         Cake::make_patch<Cake::patch_piece::minus_z>(pt),
                         Cake::make_patch<Cake::patch_piece::plus_z>(pt)};

  return PatchSystem("Cake", std::move(patches), std::move(pt));
}

} // namespace MultiPatch
