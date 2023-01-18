#include "cake.hxx"

#include <cassert>

CCTK_DEVICE CCTK_HOST MultiPatch::Cake::patch_piece
MultiPatch::Cake::get_owner_patch(const PatchTransformations &pt,
                                  const svec &global_vars) {
  using MultiPatchTests::at_boundary;
  using MultiPatchTests::isapprox;
  using MultiPatchTests::within;
  using std::abs;
  using std::sqrt;

  const auto x = global_vars(0);
  const auto y = global_vars(1);
  const auto z = global_vars(2);

  const auto r0 = pt.cake_inner_boundary_radius;
  // const auto r1 = pt.cake_outer_boundary_radius;
  const auto r = sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2));

  // Are we on the exterior?
  // if (r > r1) {
  //  return patch_piece::exterior;
  //}

  // We could be on the interior or on one of the 6 patches
  if (r <= r0) {
    return patch_piece::cartesian;
  }

  /* We are on one of the 6 Thornburg patches. The patch direction is the
   * determined by the coordinate with the largest absolute value
   */
  const auto abs_x = abs(x);
  const auto abs_y = abs(y);
  const auto abs_z = abs(z);

  const auto max_pair_idx = abs_x > abs_y ? 0 : 1;
  const auto max_pair = abs(global_vars(max_pair_idx));

  const auto max_coord_idx = abs_z > max_pair ? 2 : max_pair_idx;

  /* After the direction is know, we must dtermine the sense, which is
   * represented by the sign of the coordinate with the largest absolute
   * value.
   */
  const int max_coord_sign = global_vars(max_coord_idx) > 0.0 ? 1 : -1;

  if (max_coord_sign < 0 && max_coord_idx == 0) {
    return patch_piece::minus_x;
  } else if (max_coord_sign > 0 && max_coord_idx == 0) {
    return patch_piece::plus_x;
  } else if (max_coord_sign < 0 && max_coord_idx == 1) {
    return patch_piece::minus_y;
  } else if (max_coord_sign > 0 && max_coord_idx == 1) {
    return patch_piece::plus_y;
  } else if (max_coord_sign < 0 && max_coord_idx == 2) {
    return patch_piece::minus_z;
  } else /*if (max_coord_sign > 0 && max_coord_idx == 2)*/ {
    return patch_piece::plus_z;
  }
}
