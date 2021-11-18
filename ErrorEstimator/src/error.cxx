#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace ErrorEstimator {
using namespace Loop;
using namespace std;

template <typename T> CCTK_DEVICE constexpr T pow2(const T x) { return x * x; }

template <typename T>
CCTK_DEVICE constexpr T fmax3(const T x, const T y, const T z) {
  return fmax(fmax(x, y), z);
}

enum class shape_t { sphere, cube };

template <typename T>
constexpr T radius_sphere(const T x, const T y, const T z) {
  return sqrt(pow2(x) + pow2(y) + pow2(z));
}

template <typename T> constexpr T radius_cube(const T x, const T y, const T z) {
  DECLARE_CCTK_PARAMETERS;
  return fmax3(fabs(x)/error_factor_x, fabs(y)/error_factor_y, fabs(z)/error_factor_z);
}

template <typename T, int CI, int CJ, int CK>
T lap(const GF3D<const T, CI, CJ, CK> &var, const vect<int, dim> &I) {
  const auto DI = vect<int, dim>::unit(0);
  const auto DJ = vect<int, dim>::unit(1);
  const auto DK = vect<int, dim>::unit(2);
  return fabs(var(I - DI) - 2 * var(I) + var(I + DI)) +
         fabs(var(I - DJ) - 2 * var(I) + var(I + DJ)) +
         fabs(var(I - DK) - 2 * var(I) + var(I + DK));
}

extern "C" void ErrorEstimator_Estimate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ErrorEstimator_Estimate;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL scalefactor = scale_by_resolution ? cbrt(dx * dy * dz) : 1;

  const GF3D2<CCTK_REAL> regrid_error_(layout, regrid_error);

  const shape_t shape = [&]() {
    if (CCTK_EQUALS(region_shape, "sphere"))
      return shape_t::sphere;
    if (CCTK_EQUALS(region_shape, "cube"))
      return shape_t::cube;
    abort();
  }();

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL x = fabs(p.x) + fabs(dx);
        const CCTK_REAL y = fabs(p.y) + fabs(dy);
        const CCTK_REAL z = fabs(p.z) + fabs(dz);
        const CCTK_REAL r = [&]() {
          switch (shape) {
          case shape_t::sphere:
            return radius_sphere(x, y, z);
          case shape_t::cube:
            return radius_cube(x, y, z);
          default:
            assert(0);
          };
        }();
        const CCTK_REAL err = scalefactor / fmax(r, epsilon);
        regrid_error_(p.I) = err;
      });
}

} // namespace ErrorEstimator
