/// Generic, paper-agnostic SBP infrastructure. No coefficients live here. Every
/// concrete operator plugs its tables into the structs defined in this file.
#ifndef CARPETX_SBPOPERATORS_SBP_HXX
#define CARPETX_SBPOPERATORS_SBP_HXX

#include <cctk.h>
#include <loop_device.hxx>
#include <vec.hxx>
#include <vect.hxx>
#include <rational.hxx>

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

namespace SBPOperators {

using Rational = Arith::rational<CCTK_INT>;

/// Represents a derivative stencil: a fixed set of grid-point offsets,
/// relative to the point being evaluated, and their matching coefficients.
///
/// `C` is the coefficient type. Operator tables are written as `Stencil<N,
/// Rational>` so the published coefficients can be transcribed exactly, but the
/// stencils that reach a kernel must hold plain `CCTK_REAL`: `Arith::rational`
/// arithmetic is built on `Arith::checked_*`, which throws on overflow, and a
/// SYCL kernel cannot call into the exception machinery even along a branch
/// that is never taken. `to_real` does the conversion, at compile time.
template <std::size_t size, typename C = CCTK_REAL> struct Stencil {
  /// Operator grid point offsets, relative to the point being evaluated.
  std::array<std::ptrdiff_t, size> offsets{};

  /// Operator coefficients
  std::array<C, size> coefficients{};
};

/// Converts an exact rational stencil into the `CCTK_REAL` stencil a kernel can
/// use. Host-only, and meant to be evaluated at compile time.
template <std::size_t size>
constexpr CCTK_HOST Stencil<size> to_real(const Stencil<size, Rational> &s) {
  Stencil<size> r{};
  r.offsets = s.offsets;
  for (std::size_t i = 0; i < size; ++i) {
    r.coefficients[i] = static_cast<CCTK_REAL>(s.coefficients[i]);
  }
  return r;
}

/// Applies a single stencil (the interior stencil, or one boundary-closure row)
/// at point `p` along direction `dir`. Skips zero coefficients, so we don't
/// waste cycles accessing a grid function for no good reason.
template <std::size_t size, typename T, template <typename> typename GF3D_T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T
apply_stencil(const Stencil<size> &s, const Loop::PointDesc &p, const int dir,
              const GF3D_T<const T> &gf) {

  T coeff_sum{0};
  const auto inv_DX = T{1} / p.DX[dir];

  for (std::size_t i = 0; i < size; ++i) {
    // The tables spell structural zeros as exact rationals, so they convert to
    // exactly 0 and this comparison is safe.
    if (s.coefficients[i] == 0) {
      continue;
    }
    const auto I = p.I + s.offsets[i] * p.DI[dir];
    coeff_sum += gf(I) * static_cast<T>(s.coefficients[i]);
  }

  return coeff_sum * inv_DX;
}

/// `closures` is a `std::tuple<Stencil<N0>, Stencil<N1>, ...>`, one element
/// type per boundary row, since SBP boundary rows are not generally all the
/// same width. `row` (0 = the edge point itself) is only known at run time, so
/// this unrolls over every tuple index at compile time and applies the one that
/// matches.
template <typename Tuple, typename T, template <typename> typename GF3D_T,
          std::size_t... Is>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T apply_boundary_tuple(
    const Tuple &closures, const std::size_t row, const Loop::PointDesc &p,
    const int dir, const GF3D_T<const T> &gf, std::index_sequence<Is...>) {
  T result{0};
  ((row == Is
        ? (result = apply_stencil(std::get<Is>(closures), p, dir, gf), true)
        : false) ||
   ...);
  return result;
}

template <typename Tuple, typename T, template <typename> typename GF3D_T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE T apply_boundary_tuple(
    const Tuple &closures, const std::size_t row, const Loop::PointDesc &p,
    const int dir, const GF3D_T<const T> &gf) {
  return apply_boundary_tuple(
      closures, row, p, dir, gf,
      std::make_index_sequence<std::tuple_size_v<Tuple> >{});
}

/// Diagonal norm weights for the boundary block of an SBP operator.
/// Entry i is H_{ii}/h (dimensionless), where H is the SBP norm matrix and
/// h = dx.  The interior norm weight is always 1 and is not stored here.
/// These are needed to compute SAT penalty magnitudes: sigma = 1/(H_{ii} * h).
///
/// As with `Stencil`, `C` is the coefficient type: the tables are written as
/// exact rationals and converted to `CCTK_REAL` by `to_real` before they reach
/// an operator.
template <std::size_t P, typename C = CCTK_REAL> struct NormWeights {
  std::array<C, P> weights{};
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE C
  operator[](std::size_t i) const {
    return weights[i];
  }
};

/// Converts exact rational norm weights into `CCTK_REAL` ones. Host-only, and
/// meant to be evaluated at compile time.
template <std::size_t P>
constexpr CCTK_HOST NormWeights<P> to_real(const NormWeights<P, Rational> &nw) {
  NormWeights<P> r{};
  for (std::size_t i = 0; i < P; ++i) {
    r.weights[i] = static_cast<CCTK_REAL>(nw.weights[i]);
  }
  return r;
}

/// An SBP operator on a finite 1D direction: a single interior stencil, plus
/// one-sided boundary-closure rows for each of the lower and upper edges,
/// and the diagonal boundary-block norm weights needed for SAT penalties.
///
/// `LTuple`/`RTuple` are `std::tuple<Stencil<N0>, ..., Stencil<N_{P-1}>>`:
/// boundary rows are not required to share a single width so a homogeneous
/// `std::array` can't hold them, but a tuple can.
///
/// The two edges are stored independently here, rather than one being
/// derived from the other by reflection: this represents operators that
/// aren't mirror-symmetric the same way as ones that are, at the cost of
/// every operator family transcribing both edges by hand.
template <std::size_t int_size, typename LTuple, typename RTuple>
class SBPOperator {
  static_assert(
      std::tuple_size_v<LTuple> == std::tuple_size_v<RTuple>,
      "left and right boundary closures must have the same number of rows");

public:
  static constexpr std::size_t num_boundary_rows = std::tuple_size_v<LTuple>;

private:
  LTuple l_boundary_closures{};
  RTuple r_boundary_closures{};
  Stencil<int_size> interior_stencil{};
  NormWeights<num_boundary_rows> norm_weights_{};

public:
  /// Host-only: an operator is built once, at compile time, from the tables.
  /// Kernels receive a finished operator by value, so they never run this.
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST
  SBPOperator(LTuple l, RTuple r, Stencil<int_size> i,
              NormWeights<num_boundary_rows> nw)
      : l_boundary_closures(std::move(l)), r_boundary_closures(std::move(r)),
        interior_stencil(std::move(i)), norm_weights_(std::move(nw)) {}

  /// H_{ii}/h for boundary row i (0 = outermost).  Interior rows have weight 1.
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE CCTK_REAL
  boundary_h(std::size_t row) const {
    return norm_weights_[row];
  }

  /// Apply the operator in a given direction.
  template <std::size_t dir, template <typename> typename GF3D_T, typename T>
  inline auto CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE
  apply(const Loop::PointDesc &p, const GF3D_T<const T> &gf) const -> T {

    static_assert(dir == 0 || dir == 1 || dir == 2,
                  "Invalid derivative direction");

    // Distance from the lower boundary
    const auto dist_lo = p.I[dir] - p.bnd_min[dir];

    // Distance from the upper boundary
    const auto dist_hi = (p.bnd_max[dir] - 1) - p.I[dir];

    const bool is_l_bnd =
        0 <= dist_lo && dist_lo < static_cast<int>(num_boundary_rows);
    const bool is_r_bnd =
        0 <= dist_hi && dist_hi < static_cast<int>(num_boundary_rows);

    // If the domain is narrower than 2*num_boundary_rows, a point can
    // satisfy both conditions; In this case, we will choose to apply the lowwer
    // closure.
    if (is_l_bnd) {
      return apply_boundary_tuple(
          l_boundary_closures, static_cast<std::size_t>(dist_lo), p, dir, gf);
    } else if (is_r_bnd) {
      return apply_boundary_tuple(
          r_boundary_closures, static_cast<std::size_t>(dist_hi), p, dir, gf);
    } else {
      return apply_stencil(interior_stencil, p, dir, gf);
    }
  }
};

/// Builds an `SBPOperator` with its template arguments deduced from the
/// arguments. This is what class template argument deduction would give us for
/// free, but nvcc's frontend (`cudafe++`) asserts when CTAD is used inside a
/// `constexpr` function it also has to compile for the device, so we deduce
/// through a function template instead.
template <std::size_t int_size, typename LTuple, typename RTuple>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST
    SBPOperator<int_size, LTuple, RTuple>
    make_sbp_operator(LTuple l, RTuple r, Stencil<int_size> i,
                      NormWeights<std::tuple_size_v<LTuple> > nw) {
  return SBPOperator<int_size, LTuple, RTuple>(std::move(l), std::move(r),
                                               std::move(i), std::move(nw));
}

/// Diener, Dorband, Schnetter and Tiglio, 2007
/// https://arxiv.org/abs/gr-qc/0512001
namespace ddst2007 {

/// Boundary closure rows of the 4-2 operator: four rows, widths 4, 3, 5 and 6.
using op_42_closures =
    std::tuple<Stencil<4>, Stencil<3>, Stencil<5>, Stencil<6> >;

using op_42_t = SBPOperator<5, op_42_closures, op_42_closures>;

/// Host-only: this assembles the operator from its tables, which is a
/// compile-time job. Call it on the host and capture the result by value into a
/// kernel; `SBPOperator::apply` is the part that runs on the device.
constexpr inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST op_42_t get_op_42() {
  // Left closures
  constexpr Stencil<4, Rational> lb_0{{0, 1, 2, 3},
                                      {Rational(-24, 17), Rational(59, 34),
                                       Rational(-4, 17), Rational(-3, 34)}};

  constexpr Stencil<3, Rational> lb_1{
      {-1, 0, 1}, {Rational(-1, 2), Rational(0), Rational(1, 2)}};

  constexpr Stencil<5, Rational> lb_2{{-2, -1, 0, 1, 2},
                                      {Rational(4, 43), Rational(-59, 86),
                                       Rational(0), Rational(59, 86),
                                       Rational(-4, 43)}};

  constexpr Stencil<6, Rational> lb_3{
      {-3, -2, -1, 0, 1, 2},
      {Rational(3, 98), Rational(0), Rational(-59, 98), Rational(0),
       Rational(32, 49), Rational(-4, 49)}};

  // Right closures
  constexpr Stencil<4, Rational> rb_0{
      {0, -1, -2, -3},
      {Rational(24, 17), Rational(-59, 34), Rational(4, 17), Rational(3, 34)}};

  constexpr Stencil<3, Rational> rb_1{
      {1, 0, -1}, {Rational(1, 2), Rational(0), Rational(-1, 2)}};

  constexpr Stencil<5, Rational> rb_2{{2, 1, 0, -1, -2},
                                      {Rational(-4, 43), Rational(59, 86),
                                       Rational(0), Rational(-59, 86),
                                       Rational(4, 43)}};

  constexpr Stencil<6, Rational> rb_3{
      {3, 2, 1, 0, -1, -2},
      {Rational(-3, 98), Rational(0), Rational(59, 98), Rational(0),
       Rational(-32, 49), Rational(4, 49)}};

  // Interior
  constexpr Stencil<5, Rational> interior{{-2, -1, 0, 1, 2},
                                          {Rational(1, 12), Rational(-2, 3),
                                           Rational(0), Rational(2, 3),
                                           Rational(-1, 12)}};

  // Diagonal boundary-block norm weights H_{ii}/h from DDST07 Table 1
  constexpr NormWeights<4, Rational> nw{
      {Rational(17, 48), Rational(59, 48), Rational(43, 48), Rational(49, 48)}};

  // Operator. Everything below this line is in CCTK_REAL: the exact rationals
  // above exist only to make the tables checkable against the paper.
  constexpr auto lb =
      std::make_tuple(to_real(lb_0), to_real(lb_1), to_real(lb_2),
                      to_real(lb_3));
  constexpr auto rb =
      std::make_tuple(to_real(rb_0), to_real(rb_1), to_real(rb_2),
                      to_real(rb_3));

  return make_sbp_operator(lb, rb, to_real(interior), to_real(nw));
}

} // namespace ddst2007

} // namespace SBPOperators

#endif // #ifndef CARPETX_SBPOPERATORS_SBP_HXX
