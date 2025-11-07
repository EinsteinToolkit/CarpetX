#include <defs.hxx>
#include <loop_device.hxx>
#include <sum.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace BoxInBox {
using namespace Loop;

constexpr int max_num_regions = 3;
constexpr int max_num_levels = 30;

enum class shape_t { sphere, cube };

template <typename T>
constexpr bool should_refine_sphere(const vect<T, dim> &X,
                                    const vect<T, dim> &Rmax) {
  return sum(pow2(X / Rmax)) <= 1;
}

template <typename T>
constexpr bool should_refine_cube(const vect<T, dim> &X,
                                  const vect<T, dim> &Rmax) {
  using std::abs;
  return maximum(abs(X / Rmax)) <= 1;
}

template <typename T>
constexpr bool should_refine(const shape_t shape, const vect<T, dim> &X,
                             const vect<T, dim> &Rmax) {
  switch (shape) {
  case shape_t::sphere:
    return should_refine_sphere(X, Rmax);
  case shape_t::cube:
    return should_refine_cube(X, Rmax);
  default:
    assert(0);
  };
}

extern "C" void BoxInBox_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_BoxInBox_Init;
  DECLARE_CCTK_PARAMETERS;

  // Copy parameters into grid scalars

  {
    const int region = 0;
    active[region] = active_1;
    num_levels[region] = num_levels_1;
    position_x[region] = position_x_1;
    position_y[region] = position_y_1;
    position_z[region] = position_z_1;
    for (int level = 0; level < max_num_levels; ++level) {
      radius[level + max_num_levels * region] = radius_1[level];
      radius_x[level + max_num_levels * region] = radius_x_1[level];
      radius_y[level + max_num_levels * region] = radius_y_1[level];
      radius_z[level + max_num_levels * region] = radius_z_1[level];
    }
  }

  {
    const int region = 1;
    active[region] = active_2;
    num_levels[region] = num_levels_2;
    position_x[region] = position_x_2;
    position_y[region] = position_y_2;
    position_z[region] = position_z_2;
    for (int level = 0; level < max_num_levels; ++level) {
      radius[level + max_num_levels * region] = radius_2[level];
      radius_x[level + max_num_levels * region] = radius_x_2[level];
      radius_y[level + max_num_levels * region] = radius_y_2[level];
      radius_z[level + max_num_levels * region] = radius_z_2[level];
    }
  }

  {
    const int region = 2;
    active[region] = active_3;
    num_levels[region] = num_levels_3;
    position_x[region] = position_x_3;
    position_y[region] = position_y_3;
    position_z[region] = position_z_3;
    for (int level = 0; level < max_num_levels; ++level) {
      radius[level + max_num_levels * region] = radius_3[level];
      radius_x[level + max_num_levels * region] = radius_x_3[level];
      radius_y[level + max_num_levels * region] = radius_y_3[level];
      radius_z[level + max_num_levels * region] = radius_z_3[level];
    }
  }
}

extern "C" void BoxInBox_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_BoxInBox_Setup;
  DECLARE_CCTK_PARAMETERS;

  // Refine only patch 0
  if (cctk_patch != 0)
    return;

  const auto get_shape = [](const auto &shape) {
    if (CCTK_EQUALS(shape, "sphere"))
      return shape_t::sphere;
    if (CCTK_EQUALS(shape, "cube"))
      return shape_t::cube;
    CCTK_ERROR("internal error");
  };
  const vect<shape_t, max_num_regions> shapes = {
      get_shape(shape_1), get_shape(shape_2), get_shape(shape_3)};

  // The level about which we're deciding on, i.e. the current level + 1
  using std::ilogb;
  const int level = ilogb(CCTK_REAL(cctk_levfac[0])) + 1;

  const auto active1 = vect<bool, max_num_regions>::make(
      [&](int region) { return bool(active[region]); });
  const auto num_levels1 = vect<int, max_num_regions>::make(
      [&](int region) { return int(num_levels[region]); });
  const auto positions1 =
      vect<vect<CCTK_REAL, dim>, max_num_regions>::make([&](int region) {
        return vect<CCTK_REAL, dim>{position_x[region], position_y[region],
                                    position_z[region]};
      });
  const auto radii1 =
      vect<vect<CCTK_REAL, dim>, max_num_regions>::make([&](int region) {
        const int level_region = level + max_num_levels * region;
        return vect<CCTK_REAL, dim>{
            radius_x[level_region] < 0 ? radius[level_region]
                                       : radius_x[level_region],
            radius_y[level_region] < 0 ? radius[level_region]
                                       : radius_y[level_region],
            radius_z[level_region] < 0 ? radius[level_region]
                                       : radius_z[level_region]};
      });

  for (int region = 0; region < max_num_regions; ++region)
    if (active1[region])
      if (level < num_levels1[region])
        if (any(radii1[region] <= 0))
          CCTK_VERROR("Non-positive radius found: Refined region %d, "
                      "refinement level %d, radii are [%.17g,%.17g,%.17g]",
                      region, level, double(radii1[region][0]),
                      double(radii1[region][1]), double(radii1[region][2]));

  const auto positions2 = positions1;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool do_refine = false;

        for (int region = 0; region < max_num_regions; ++region) {
          if (active1[region]) {
            if (level < num_levels1[region]) {
              const auto position = positions2[region];
              const auto radii = radii1[region];
              // Decrease the box radius by one grid spacing since AMReX seems
              // to expand the refined region by one cell
              const auto Rmax = radii - p.DX;
              const bool sr =
                  should_refine(shapes[region], p.X - position, Rmax);
              do_refine |= sr;
            }
          }
        }

        regrid_error(p.I) = do_refine ? 1 : 0;
      });
}

} // namespace BoxInBox
