#ifndef MULTIPATCH_CAKE_HXX
#define MULTIPATCH_CAKE_HXX

#include "multipatch.hxx"
#include "tests.hxx"

#include <string>
#include <type_traits>

namespace MultiPatch {
namespace Cake {

/**
 * Stores a 3-vector
 */
using svec = vec<CCTK_REAL, dim>;

/**
 * Stores a 3-matrix
 */
using smat = smat<CCTK_REAL, dim>;

/**
 * Stores a Jacobian matrix of 3 dimensions
 */
using jac_t = vec<vec<CCTK_REAL, dim>, dim>;

/**
 * Stores a the derivatives of a Jacobian matrix of 3 dimensions
 */
using djac_t = vec<smat, dim>;

/**
 * Precondition assertion. If the precondition fails, the code is aborted.
 *
 * @param predicate The predicate to test.
 * @param msg The message to display while aborting the code.
 */
CCTK_DEVICE CCTK_HOST inline void expects(bool predicate, const char *msg) {
  if (!predicate) {
#ifndef __CUDACC__
    CCTK_ERROR(msg);
#else
    assert(0);
#endif
  }
}

/**
 * Compatibility wrapper for replacing Mathematica's Power function.
 *
 * @param n The exponent.
 * @param x The base.
 * @return The n-th power of x
 */
template <typename base_type, typename exponent_type>
CCTK_DEVICE CCTK_HOST static inline auto Power(base_type x, exponent_type n) {
  using std::pow;
  return pow(x, n);
}

/**
 * Compatibility wrapper for replacing Mathematica's Sqrt function.
 *
 * @param x The radicand.
 */
template <typename T> CCTK_DEVICE CCTK_HOST static inline auto Sqrt(T x) {
  using std::sqrt;
  return sqrt(x);
}

/**
 * Tags for each patch piece in the cake
 */
enum class patch_piece : int {
  cartesian = 0,

  plus_x = 1,
  minus_x = 2,

  plus_y = 3,
  minus_y = 4,

  plus_z = 5,
  minus_z = 6,

  inner_boundary = 7,
  outer_boundary = 8,

  exterior = -1
};

/**
 * Gets name from patch piece
 *
 * @tparam p A patch piece type.
 * @return A string representing the name of the piece.
 */
inline const std::string piece_name(const patch_piece &p) {
  switch (static_cast<int>(p)) {
  case static_cast<int>(patch_piece::cartesian):
    return "cartesian";
  case static_cast<int>(patch_piece::plus_x):
    return "plus x";
  case static_cast<int>(patch_piece::minus_x):
    return "minus x";
  case static_cast<int>(patch_piece::plus_y):
    return "plus y";
  case static_cast<int>(patch_piece::minus_y):
    return "minus y";
  case static_cast<int>(patch_piece::plus_z):
    return "plus z";
  case static_cast<int>(patch_piece::minus_z):
    return "minus z";
  case static_cast<int>(patch_piece::inner_boundary):
    return "inner boundary";
  case static_cast<int>(patch_piece::outer_boundary):
    return "outer boundary";
  default:
    return "exterior";
  }
}

/**
 * Get the patch piece that owns a global coordinate point.
 *
 * @param pt The PatchTransformations structure describing the patch system.
 * @param global_vars The global coordinate triplet to locate the owner for.
 * @return The patch piece owning the global coordinates.
 */
CCTK_DEVICE CCTK_HOST patch_piece
get_owner_patch(const PatchTransformations &pt, const svec &global_vars);

} // namespace Cake
} // namespace MultiPatch

#endif // MULTIPATCH_CAKE_HXX
