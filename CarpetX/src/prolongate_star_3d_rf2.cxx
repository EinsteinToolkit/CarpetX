#include "prolongate_star_3d_rf2.hxx"

#include "driver.hxx"

#include <defs.hxx>
#include <vect.hxx>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace CarpetX {

using namespace Arith;

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////

template <typename T> static constexpr T minmod(const T x, const T y) {
  using std::abs, std::copysign, std::min, std::signbit;
  if (!(signbit(x) == signbit(y)))
    return 0;
  return copysign(min(abs(x), abs(y)), x);
}

template <typename T>
static constexpr T minmod(const T x, const T y, const T z) {
  using std::abs, std::copysign, std::min, std::signbit;
  if (!(signbit(x) == signbit(y) && signbit(x) == signbit(z)))
    return 0;
  return copysign(min({abs(x), abs(y), abs(z)}), x);
}

////////////////////////////////////////

// Vertex-centred interpolation
template <int ORDER> struct interp_vc_poly;

// Order 1 (linear)
template <> struct interp_vc_poly<1> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 1;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F,
            typename T = std::remove_reference_t<
                std::remove_reference_t<std::invoke_result_t<F, int> > > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return 1 / T(2) * (values(0) + values(1));
  }
};

// Order 3 (cubic)
template <> struct interp_vc_poly<3> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 3;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return -1 / T(16) * (values(-1) + values(2)) //
           + 9 / T(16) * (values(0) + values(1));
  }
};

// Order 5
template <> struct interp_vc_poly<5> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 5;
  static constexpr int required_ghosts = 3;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return 3 / T(256) * (values(-2) + values(3))    //
           - 25 / T(256) * (values(-1) + values(2)) //
           + 75 / T(128) * (values(0) + values(1));
  }
};

// Order 7
template <> struct interp_vc_poly<7> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 7;
  static constexpr int required_ghosts = 4;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return -5 / T(2048) * (values(-3) + values(4))    //
           + 49 / T(2048) * (values(-2) + values(3))  //
           - 245 / T(2048) * (values(-1) + values(2)) //
           + 1225 / T(2048) * (values(0) + values(1));
  }
};

// Order 9
template <> struct interp_vc_poly<9> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 9;
  static constexpr int required_ghosts = 5;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return 35 / T(65536) * (values(-4) + values(5))     //
           - 405 / T(65536) * (values(-3) + values(4))  //
           + 567 / T(16384) * (values(-2) + values(3))  //
           - 2205 / T(16384) * (values(-1) + values(2)) //
           + 19845 / T(32768) * (values(0) + values(1));
  }
};

// Order 11
template <> struct interp_vc_poly<11> {
  static inline const std::string name = "vc_poly";
  static constexpr int order = 11;
  static constexpr int required_ghosts = 6;
  static constexpr bool stencil_is_cc = false;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return -63 / T(524288) * (values(-5) + values(6))     //
           + 847 / T(524288) * (values(-4) + values(5))   //
           - 5445 / T(524288) * (values(-3) + values(4))  //
           + 22869 / T(524288) * (values(-2) + values(3)) //
           - 38115 / T(262144) * (values(-1) + values(2)) //
           + 160083 / T(262144) * (values(0) + values(1));
  }
};

////////////////////////////////////////

// Cell-centred polynomial (non-conservative) interpolation
template <int ORDER> struct interp_cc_poly;

// Order 0 (constant)
template <> struct interp_cc_poly<0> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 0;
  static constexpr int required_ghosts = 0;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return values(0);
  }
};

// Order 1 (linear)
template <> struct interp_cc_poly<1> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 1;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(1) although it is available
    return 1 / T(4) * values(-1) //
           + 3 / T(4) * values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 2 (quadratic)
template <> struct interp_cc_poly<2> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 2;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return 5 / T(32) * values(-1)   //
           + 15 / T(16) * values(0) //
           - 3 / T(32) * values(1);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 3 (cubic)
template <> struct interp_cc_poly<3> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 3;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(2) although it is available
    return -5 / T(128) * values(-2)   //
           + 35 / T(128) * values(-1) //
           + 105 / T(128) * values(0) //
           - 7 / T(128) * values(1);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 4
template <> struct interp_cc_poly<4> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 4;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return -45 / T(2048) * values(-2)  //
           + 105 / T(512) * values(-1) //
           + 945 / T(1024) * values(0) //
           - 63 / T(512) * values(1)   //
           + 35 / T(2048) * values(2);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 5
template <> struct interp_cc_poly<5> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 5;
  static constexpr int required_ghosts = 3;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(3) although it is available
    return 63 / T(8192) * values(-3)     //
           - 495 / T(8192) * values(-2)  //
           + 1155 / T(4096) * values(-1) //
           + 3465 / T(4096) * values(0)  //
           - 693 / T(8192) * values(1)   //
           + 77 / T(8192) * values(2);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 6
template <> struct interp_cc_poly<6> {
  static inline const std::string name = "cc_poly";
  static constexpr int order = 6;
  static constexpr int required_ghosts = 3;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return 273 / T(65536) * values(-3)     //
           - 1287 / T(32768) * values(-2)  //
           + 15015 / T(65536) * values(-1) //
           + 15015 / T(16384) * values(0)  //
           - 9009 / T(65536) * values(1)   //
           + 1001 / T(32768) * values(2)   //
           - 231 / T(65536) * values(3);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// ... more orders available upon request ...

////////////////////////////////////////

// Cell-centred conservative interpolation
template <int ORDER> struct interp_cc_cons;

// Note: The even orders here have a coefficient of 1 for values(0). This means
// they can be used as star-shaped stencil: a star-shaped stencil needs to be a
// combination of (1) a copy of the coarse-grid value and (2) a correction. This
// can also be used for ENO operators.

// Order 0 (constant)
// This is the same as the non-constant polynomial interpolation
// interp_cc_poly<0>.
template <> struct interp_cc_cons<0> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 0;
  static constexpr int required_ghosts = 0;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return values(0);
  }
};

// Order 1 (linear)
// This is the same as the non-constant polynomial interpolation
// interp_cc_poly<1>.
template <> struct interp_cc_cons<1> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 1;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(1) although it is available
    return 1 / T(4) * values(-1) + 3 / T(4) * values(0);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 2 (quadratic)
template <> struct interp_cc_cons<2> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 2;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return 1 / T(8) * values(-1)  //
           + 1 / T(1) * values(0) //
           - 1 / T(8) * values(1);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 3 (cubic)
template <> struct interp_cc_cons<3> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 3;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(2) although it is available
    return -3 / T(64) * values(-2)   //
           + 17 / T(64) * values(-1) //
           + 55 / T(64) * values(0)  //
           - 5 / T(64) * values(1);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 4
template <> struct interp_cc_cons<4> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 4;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return -3 / T(128) * values(-2)  //
           + 11 / T(64) * values(-1) //
           + 1 / T(1) * values(0)    //
           - 11 / T(64) * values(1)  //
           + 3 / T(128) * values(2);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 5
template <> struct interp_cc_cons<5> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 5;
  static constexpr int required_ghosts = 3;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // Skipping values(3) although it is available
    return +5 / T(512) * values(-3)   //
           - 37 / T(512) * values(-2) //
           + 69 / T(256) * values(-1) //
           + 231 / T(256) * values(0) //
           - 63 / T(512) * values(1)  //
           + 7 / T(512) * values(2);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 6
template <> struct interp_cc_cons<6> {
  static inline const std::string name = "cc_cons";
  static constexpr int order = 6;
  static constexpr int required_ghosts = 3;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = false;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    return +5 / T(1024) * values(-3)    //
           - 11 / T(256) * values(-2)   //
           + 201 / T(1024) * values(-1) //
           + 1 / T(1) * values(0)       //
           - 201 / T(1024) * values(1)  //
           + 11 / T(256) * values(2)    //
           - 5 / T(1024) * values(3);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// ... more orders available upon request ...

////////////////////////////////////////

// Cell-centred minmod interpolation. This is conservative.
template <int ORDER> struct interp_cc_cons_minmod;

// Order 1 (linear)
template <> struct interp_cc_cons_minmod<1> {
  static inline const std::string name = "cc_cons_minmod";
  static constexpr int order = 1;
  static constexpr int required_ghosts = 1;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = true;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    // First calculate the left and right derivatives on the coarse grid.
    // We skip dividing by and multiplying by the grid spacing.
    const T dx_lo = values(0) - values(-1);
    const T dx_up = values(1) - values(0);
    const T dx = minmod(dx_lo, dx_up);
    return values(0) - dx / 4;
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

////////////////////////////////////////

// Cell-centred (non-conservative) ENO interpolation
template <int ORDER> struct interp_cc_poly_eno;

// Order 1 (linear) is the same as interp_cc_poly<1> since the stencil is so
// small that it cannot be shifted without becoming an extrapolation.

// Order 2 (quadratic)
template <> struct interp_cc_poly_eno<2> {
  static inline const std::string name = "cc_poly_eno";
  static constexpr int order = 2;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = true;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    using std::abs;
    // Skipping values(2) although it is available
    // First calculate the undivided differences on the coarse grid (first
    // derivatives)
    const T dx_m1 = values(-1) - values(-2);
    const T dx_p0 = values(0) - values(-1);
    const T dx_p1 = values(1) - values(0);
    // Then calculate their undivided differences (second derivatives)
    const T ddx_m1 = dx_p0 - dx_m1;
    const T ddx_p0 = dx_p1 - dx_p0;
    // Choose the smoother one
    if (abs(ddx_m1) <= abs(ddx_p0)) {
      return -3 / T(32) * values(-2)  //
             + 7 / T(16) * values(-1) //
             + 21 / T(32) * values(0);
    } else {
      return 5 / T(32) * values(-1)   //
             + 15 / T(16) * values(0) //
             - 3 / T(32) * values(1);
    }
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 3 (cubic)
template <> struct interp_cc_poly_eno<3> {
  static inline const std::string name = "cc_poly_eno";
  static constexpr int order = 3;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = false;
  static constexpr bool stencil_is_star = true;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    using std::abs;
    // First calculate the undivided differences on the coarse grid (first
    // derivatives)
    const T dx_m1 = values(-1) - values(-2);
    const T dx_p0 = values(0) - values(-1);
    const T dx_p1 = values(1) - values(0);
    const T dx_p2 = values(2) - values(1);
    // Then calculate their undivided differences (second derivatives)
    const T ddx_m1 = dx_p0 - dx_m1;
    const T ddx_p0 = dx_p1 - dx_p0;
    const T ddx_p1 = dx_p2 - dx_p1;
    // Then calculate their undivided differences (third derivatives)
    const T dddx_m1 = ddx_p0 - ddx_m1;
    const T dddx_p0 = ddx_p1 - ddx_p0;
    // Choose the smoother one
    if (abs(dddx_m1) <= abs(dddx_p0)) {
      return -5 / T(128) * values(-2)   //
             + 35 / T(128) * values(-1) //
             + 105 / T(128) * values(0) //
             - 7 / T(128) * values(1);
    } else {
      return 15 / T(128) * values(-1)   //
             + 135 / T(128) * values(0) //
             - 27 / T(128) * values(1)  //
             + 5 / T(128) * values(2);
    }
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// ... more orders available upon request ...

////////////////////////////////////////

// Cell-centred (conservative) ENO interpolation
template <int ORDER> struct interp_cc_cons_eno;

// Order 2 (quadratic)
template <> struct interp_cc_cons_eno<2> {
  static inline const std::string name = "cc_cons_eno";
  static constexpr int order = 2;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = true;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    using std::abs, std::min;
    // First calculate the undivided differences on the coarse grid (first
    // derivatives)
    const T dx_m1 = values(-1) - values(-2);
    const T dx_p0 = values(0) - values(-1);
    const T dx_p1 = values(1) - values(0);
    const T dx_p2 = values(2) - values(1);
    // Then calculate their undivided differences (second derivatives)
    const T ddx_m1 = dx_p0 - dx_m1;
    const T ddx_p0 = dx_p1 - dx_p0;
    const T ddx_p1 = dx_p2 - dx_p1;
    // Choose the smoothest one
    if (abs(ddx_p0) <= min(abs(ddx_m1), abs(ddx_p1))) {
      // Centred stencil
      // xi = x0 - 1/4 [-1/2 0 1/2]
      return 1 / T(8) * values(-1)  //
             + 1 / T(1) * values(0) //
             - 1 / T(8) * values(1);
    } else if (abs(ddx_m1) <= min(abs(ddx_p0), abs(ddx_p1))) {
      // Left-shifted stencil
      // xi = x0 - 1/4 [1/2 -2 3/2]
      return -1 / T(8) * values(-2)  //
             + 1 / T(2) * values(-1) //
             + 5 / T(8) * values(0);
    } else {
      // Right-shifted stencil
      // xi = x0 - 1/4 [-3/2 2 -1/2]
      return 11 / T(8) * values(0)  //
             - 1 / T(2) * values(1) //
             + 1 / T(8) * values(2);
    }
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

// Order 3 (cubic)
template <> struct interp_cc_cons_eno<3> {
  static inline const std::string name = "cc_cons_eno";
  static constexpr int order = 3;
  static constexpr int required_ghosts = 2;
  static constexpr bool stencil_is_cc = true;
  static constexpr bool stencil_is_conservative = true;
  static constexpr bool stencil_is_star = true;
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_lo(const F &values) {
    using std::abs;
    // First calculate the undivided differences on the coarse grid (first
    // derivatives)
    const T dx_m1 = values(-1) - values(-2);
    const T dx_p0 = values(0) - values(-1);
    const T dx_p1 = values(1) - values(0);
    const T dx_p2 = values(2) - values(1);
    // Then calculate their undivided differences (second derivatives)
    const T ddx_m1 = dx_p0 - dx_m1;
    const T ddx_p0 = dx_p1 - dx_p0;
    const T ddx_p1 = dx_p2 - dx_p1;
    // Then calculate their undivided differences (third derivatives)
    const T dddx_m1 = ddx_p0 - ddx_m1;
    const T dddx_p0 = ddx_p1 - ddx_p0;
    // Choose the smoother one
    if (abs(dddx_m1) <= abs(dddx_p0))
      return -3 / T(64) * values(-2)   //
             + 17 / T(64) * values(-1) //
             + 55 / T(64) * values(0)  //
             - 5 / T(64) * values(1);
    else
      return 5 / T(64) * values(-1)   //
             + 73 / T(64) * values(0) //
             - 17 / T(64) * values(1) //
             + 3 / T(64) * values(2);
  }
  template <typename F, typename T = std::remove_reference_t<
                            std::invoke_result_t<F, int> > >
  static constexpr T interp_up(const F &values) {
    return interp_lo([&](int i) { return values(-i); });
  }
};

////////////////////////////////////////////////////////////////////////////////

// Test interpolation operators

template <typename Interp, typename T> struct test_interp {
  constexpr static int nghosts = Interp::required_ghosts;
  static_assert(nghosts >= 0);

  static void test_interp_poly(const int nforbid_lo, const int nforbid_up,
                               const int order) {
    using std::isfinite;

    assert(order >= 0);
    assert(order <= Interp::order);

    constexpr int radius = nghosts;
    constexpr int N = 2 * radius + 1;
    std::array<T, N> grid;
    const auto values = [&](int i) -> T & { return grid.at(i + radius); };

    bool found_error = false;

    // Test the stencil for all orders
    for (int p = 0; p <= order + 1; ++p) {
      // Our test polynomial is i^p
      const auto poly = [&](T x) { return pown(x, p); };

      // Set grid to polynomial
      for (int i = -nghosts; i <= +nghosts; ++i)
        values(i) = poly(2 + i);
      // Forbid some points
      for (int i = -nghosts; i <= -nghosts + nforbid_lo - 1; ++i)
        values(i) -= 1000;
      for (int i = +nghosts - nforbid_up + 1; i <= +nghosts; ++i)
        values(i) -= 1000;

      for (int offset = 0; offset < 2; ++offset) {
        // For cell-centred grids the interpolation points are different than
        // for vertex-centred grids
        constexpr T cc_offset = Interp::stencil_is_cc ? -1 / T(4) : 0;
        const T x = 2 + offset / T(2) + cc_offset;
        const T expected_result = poly(x);
        const T result =
            offset == 0 ? Interp::interp_lo(values) : Interp::interp_up(values);

        if (p == order + 1 && !(!Interp::stencil_is_cc && offset == 0)) {
          assert(result != expected_result);
          continue;
        }

        if (!isfinite(result) || result != expected_result) {
          found_error = true;
          std::cerr << "1d interpolation stencil failure"                   //
                    << ": name=" << Interp::name                            //
                    << ", original order=" << Interp::order                 //
                    << ", reduced order=" << order                          //
                    << ", centering=" << Interp::stencil_is_cc              //
                    << ", conservative=" << Interp::stencil_is_conservative //
                    << ", nforbid_lo=" << nforbid_lo                        //
                    << ", nforbid_up=" << nforbid_up                        //
                    << ", p=" << p                                          //
                    << ", offset=" << offset                                //
                    << ", expected=" << expected_result                     //
                    << ", result=" << result                                //
                    << "\n";
        }
      } // for offset
    } // for p

    if (found_error)
      std::abort();
  }

  static void test_interp_cons(const int nforbid_lo, const int nforbid_up,
                               const int order) {
    using std::isfinite;

    assert(order >= 0);
    assert(order <= Interp::order);

    constexpr int radius = nghosts;
    constexpr int N = 2 * radius + 1;
    std::array<T, N> grid;
    const auto values = [&](int i) -> T & { return grid.at(i + radius); };

    bool found_error = false;

    // Test the stencil for all orders
    for (int p = 0; p <= order + 1; ++p) {
      // Our test polynomial is (p+1) i^p
      const auto poly = [&](T x) { return (p + 1) * pown(x, p); };
      // Integrated polynomial
      const auto ipoly = [&](T x) { return pown(x, p + 1); };
      // Finite difference of the integrated polynomial for the coarse grid
      // (these will be the coarse grid values)
      const auto cpoly = [&](T x) {
        return ipoly(x + 1 / T(2)) - ipoly(x - 1 / T(2));
      };
      // Finite difference of the integrated polynomial for the coarse grid
      // (these will be the fine grid values)
      const auto fpoly = [&](T x) {
        return 2 * (ipoly(x + 1 / T(4)) - ipoly(x - 1 / T(4)));
      };

      // Set grid to polynomial
      for (int i = -nghosts; i <= +nghosts; ++i)
        values(i) = cpoly(2 + i);
      // Forbid some cells
      for (int i = -nghosts; i <= -nghosts + nforbid_lo - 1; ++i)
        values(i) -= 1000;
      for (int i = +nghosts - nforbid_up + 1; i <= +nghosts; ++i)
        values(i) -= 1000;

      for (int offset = 0; offset < 2; ++offset) {
        // For cell-centred grids the interpolation points are different than
        // for vertex-centred grids
        constexpr T cc_offset = Interp::stencil_is_cc ? -1 / T(4) : 0;
        const T x = 2 + offset / T(2) + cc_offset;
        const T expected_result = fpoly(x);
        const T result =
            offset == 0 ? Interp::interp_lo(values) : Interp::interp_up(values);

        if (p == order + 1) {
          // assert(result != expected_result);
          continue;
        }

        if (!isfinite(result) || result != expected_result) {
          found_error = true;
          std::cerr << "1d interpolation stencil failure"                   //
                    << ": name=" << Interp::name                            //
                    << ", original order=" << Interp::order                 //
                    << ", reduced order=" << order                          //
                    << ", centering=" << Interp::stencil_is_cc              //
                    << ", conservative=" << Interp::stencil_is_conservative //
                    << ", nforbid_lo=" << nforbid_lo                        //
                    << ", nforbid_up=" << nforbid_up                        //
                    << ", p=" << p                                          //
                    << ", offset=" << offset                                //
                    << ", expected=" << expected_result                     //
                    << ", result=" << result                                //
                    << "\n";
          std::cerr << "Values [";
          for (int i = -nghosts; i <= +nghosts; ++i)
            std::cerr << values(i) << ", ";
          std::cerr << "]\n";
        }
      } // for offset
    } // for p

    if (found_error)
      std::abort();
  }

  // Constructor
  test_interp(const int nforbid_lo = 0, const int nforbid_up = 0,
              const int order = Interp::order) {
    // std::cout << "Testing 1d interpolation stencil <" << Interp::name
    //           << ">...\n";

    if (Interp::stencil_is_conservative)
      test_interp_cons(nforbid_lo, nforbid_up, order);
    else
      test_interp_poly(nforbid_lo, nforbid_up, order);
  }
};

// Create static test object for all stencils. These tests will run at
// startup.

static test_interp<interp_vc_poly<1>, CCTK_REAL> test_interp_vc_poly_o1;
static test_interp<interp_vc_poly<3>, CCTK_REAL> test_interp_vc_poly_o3;
static test_interp<interp_vc_poly<5>, CCTK_REAL> test_interp_vc_poly_o5;
static test_interp<interp_vc_poly<7>, CCTK_REAL> test_interp_vc_poly_o7;
static test_interp<interp_vc_poly<9>, CCTK_REAL> test_interp_vc_poly_o9;
static test_interp<interp_vc_poly<11>, CCTK_REAL> test_interp_vc_poly_o11;

static test_interp<interp_cc_poly<0>, CCTK_REAL> test_interp_cc_poly_o0;
static test_interp<interp_cc_poly<1>, CCTK_REAL> test_interp_cc_poly_o1;
static test_interp<interp_cc_poly<2>, CCTK_REAL> test_interp_cc_poly_o2;
static test_interp<interp_cc_poly<3>, CCTK_REAL> test_interp_cc_poly_o3;
static test_interp<interp_cc_poly<4>, CCTK_REAL> test_interp_cc_poly_o4;
static test_interp<interp_cc_poly<5>, CCTK_REAL> test_interp_cc_poly_o5;
static test_interp<interp_cc_poly<6>, CCTK_REAL> test_interp_cc_poly_o6;

static test_interp<interp_cc_cons<0>, CCTK_REAL> test_interp_cc_cons_o0;
static test_interp<interp_cc_cons<1>, CCTK_REAL> test_interp_cc_cons_o1;
static test_interp<interp_cc_cons<2>, CCTK_REAL> test_interp_cc_cons_o2;
static test_interp<interp_cc_cons<3>, CCTK_REAL> test_interp_cc_cons_o3;
static test_interp<interp_cc_cons<4>, CCTK_REAL> test_interp_cc_cons_o4;
static test_interp<interp_cc_cons<5>, CCTK_REAL> test_interp_cc_cons_o5;
static test_interp<interp_cc_cons<6>, CCTK_REAL> test_interp_cc_cons_o6;

static test_interp<interp_cc_poly_eno<2>, CCTK_REAL> test_interp_cc_poly_eno_o2;
static test_interp<interp_cc_poly_eno<2>, CCTK_REAL>
    test_interp_cc_poly_eno_o2_fl(1, 0, 2);
static test_interp<interp_cc_poly_eno<2>, CCTK_REAL>
    test_interp_cc_poly_eno_o2_fu(0, 1, 2);
static test_interp<interp_cc_poly_eno<3>, CCTK_REAL> test_interp_cc_poly_eno_o3;
static test_interp<interp_cc_poly_eno<3>, CCTK_REAL>
    test_interp_cc_poly_eno_o3_fl(1, 0, 3);
static test_interp<interp_cc_poly_eno<3>, CCTK_REAL>
    test_interp_cc_poly_eno_o3_fu(0, 1, 3);

static test_interp<interp_cc_cons_eno<2>, CCTK_REAL> test_interp_cc_cons_eno_o2;
static test_interp<interp_cc_cons_eno<2>, CCTK_REAL>
    test_interp_cc_cons_eno_o2_fll(2, 0, 2);
static test_interp<interp_cc_cons_eno<2>, CCTK_REAL>
    test_interp_cc_cons_eno_o2_flu(1, 1, 2);
static test_interp<interp_cc_cons_eno<2>, CCTK_REAL>
    test_interp_cc_cons_eno_o2_fuu(0, 2, 2);
static test_interp<interp_cc_cons_eno<3>, CCTK_REAL> test_interp_cc_cons_eno_o3;
static test_interp<interp_cc_cons_eno<3>, CCTK_REAL>
    test_interp_cc_cons_eno_o3_fl(1, 0, 3);
static test_interp<interp_cc_cons_eno<3>, CCTK_REAL>
    test_interp_cc_cons_eno_o3_fu(0, 1, 3);

static test_interp<interp_cc_cons_minmod<1>, CCTK_REAL>
    test_interp_cc_cons_minmod_o1;
static test_interp<interp_cc_cons_minmod<1>, CCTK_REAL>
    test_interp_cc_cons_minmod_o1_fl(1, 0, 0);
static test_interp<interp_cc_cons_minmod<1>, CCTK_REAL>
    test_interp_cc_cons_minmod_o1_fu(0, 1, 0);

////////////////////////////////////////////////////////////////////////////////

// Decide how the 1d stencils should be combined into a 3d stencil. This
// depends on:
// - whether the 1d stencils want to combine into a star-shaped or a
//   box-shaped stencil
// - which offsets we have with respect to the coarse grid
// When combining box and star directions (which is a strange thing to do), then
// handle the box directions first. This is good when the star directions make a
// choice (e.g. for ENO or minmod), because these choices are then made last and
// their impact survives.
template <typename InterpI, typename InterpJ, typename InterpK, int offset_i,
          int offset_j, int offset_k, typename F,
          typename T =
              std::remove_reference_t<std::invoke_result_t<F, int, int, int> > >
static CCTK_DEVICE CCTK_HOST inline CCTK_ATTRIBUTE_ALWAYS_INLINE T
combine_1d_stencils(const F &get_crse) {
  constexpr bool star_i = InterpI::stencil_is_star;
  constexpr bool star_j = InterpJ::stencil_is_star;
  constexpr bool star_k = InterpK::stencil_is_star;
  // Combine stencil shape into a single variable
  // 0b000 ... 0b111, defined as 0bijk for the three directions i, j, k
  constexpr int star_ijk = (star_i << 2) | (star_j << 1) | (star_k << 0);

  // Shortcuts for interpolation stencils for all directions
  const auto interp_i = [&](const auto &values) {
    if constexpr (offset_i == 0)
      return InterpI::interp_lo(values);
    else
      return InterpI::interp_up(values);
  };
  const auto interp_j = [&](const auto &values) {
    if constexpr (offset_j == 0)
      return InterpJ::interp_lo(values);
    else
      return InterpJ::interp_up(values);
  };
  const auto interp_k = [&](const auto &values) {
    if constexpr (offset_k == 0)
      return InterpK::interp_lo(values);
    else
      return InterpK::interp_up(values);
  };

  T result;
  if constexpr (star_ijk == 0b000) {
    // Box stencil: tensor product
    result = interp_k([&](int dk) {
      return interp_j([&](int dj) {
        return interp_i([&](int di) { return get_crse(di, dj, dk); });
      });
    });
  } else if constexpr (star_ijk == 0b100) {
    // Star-shaped in the i direction
    const auto get_crse_i = [&](int di) {
      return interp_k([&](int dk) {
        return interp_j([&](int dj) { return get_crse(di, dj, dk); });
      });
    };
    result = interp_i([&](int di) { return get_crse_i(di); });
  } else if constexpr (star_ijk == 0b010) {
    // Star-shaped in the j direction
    const auto get_crse_j = [&](int dj) {
      return interp_k([&](int dk) {
        return interp_i([&](int di) { return get_crse(di, dj, dk); });
      });
    };
    result = interp_j([&](int dj) { return get_crse_j(dj); });
  } else if constexpr (star_ijk == 0b001) {
    // Star-shaped in the k direction
    const auto get_crse_k = [&](int dk) {
      return interp_j([&](int dj) {
        return interp_i([&](int di) { return get_crse(di, dj, dk); });
      });
    };
    result = interp_k([&](int dk) { return get_crse_k(dk); });
  } else if constexpr (star_ijk == 0b110) {
    // Star-shaped in the i and j direction
    const auto get_crse_ij = [&](int di, int dj) {
      return interp_k([&](int dk) { return get_crse(di, dj, dk); });
    };
    result = -get_crse_ij(0, 0)                                     //
             + interp_i([&](int di) { return get_crse_ij(di, 0); }) //
             + interp_j([&](int dj) { return get_crse_ij(0, dj); });
  } else if constexpr (star_ijk == 0b101) {
    // Star-shaped in the i and k direction
    const auto get_crse_ik = [&](int di, int dk) {
      return interp_j([&](int dj) { return get_crse(di, dj, dk); });
    };
    result = -get_crse_ik(0, 0)                                     //
             + interp_i([&](int di) { return get_crse_ik(di, 0); }) //
             + interp_k([&](int dk) { return get_crse_ik(0, dk); });
  } else if constexpr (star_ijk == 0b011) {
    // Star-shaped in the j and k direction
    const auto get_crse_jk = [&](int dj, int dk) {
      return interp_i([&](int di) { return get_crse(di, dj, dk); });
    };
    result = -get_crse_jk(0, 0)                                     //
             + interp_j([&](int dj) { return get_crse_jk(dj, 0); }) //
             + interp_k([&](int dk) { return get_crse_jk(0, dk); });
  } else if constexpr (star_ijk == 0b111) {
    // Star stencil: Coarse value plus 3 corrections
    result = -2 * get_crse(0, 0, 0)                                 //
             + interp_i([&](int di) { return get_crse(di, 0, 0); }) //
             + interp_j([&](int dj) { return get_crse(0, dj, 0); }) //
             + interp_k([&](int dk) { return get_crse(0, 0, dk); });
  }
#ifndef __CUDACC__
  // This does not work with CUDA although it should
  else {
    // Stencil shape not implemented
    static_assert(false);
  }
#endif

  return result;
}

////////////////////////////////////////////////////////////////////////////////

template <typename InterpI, typename InterpJ, typename InterpK, typename T>
struct test_interp_3d {
  constexpr static int order_i = InterpI::order;
  constexpr static int order_j = InterpJ::order;
  constexpr static int order_k = InterpK::order;
  static_assert(order_i >= 0);
  static_assert(order_j >= 0);
  static_assert(order_k >= 0);

  constexpr static int nghosts_i = InterpI::required_ghosts;
  constexpr static int nghosts_j = InterpJ::required_ghosts;
  constexpr static int nghosts_k = InterpK::required_ghosts;
  static_assert(nghosts_i >= 0);
  static_assert(nghosts_j >= 0);
  static_assert(nghosts_k >= 0);

  test_interp_3d() {
    using std::abs, std::isfinite, std::max;

    // std::cout << "Testing 3d interpolation stencil <" << InterpI::name << "["
    //           << InterpI::order << "]," << InterpJ::name << "["
    //           << InterpJ::order << "]," << InterpK::name << "["
    //           << InterpK::order << "]>...\n";

    constexpr int radius_i = nghosts_i;
    constexpr int radius_j = nghosts_j;
    constexpr int radius_k = nghosts_k;
    constexpr int NI = 2 * radius_i + 1;
    constexpr int NJ = 2 * radius_j + 1;
    constexpr int NK = 2 * radius_k + 1;
    std::vector<T> grid(NI * NJ * NK);
    const auto values = [&](int i, int j, int k) -> T & {
      assert(i >= -radius_i && i <= +radius_i);
      assert(j >= -radius_j && j <= +radius_j);
      assert(k >= -radius_k && k <= +radius_k);
      return grid.at((i + radius_i) +
                     NI * ((j + radius_j) + NJ * (k + radius_k)));
    };

    bool found_error = false;

    // Test the stencil for all orders
    for (int pk = 0; pk <= order_k; ++pk) {
      for (int pj = 0; pj <= order_j; ++pj) {
        for (int pi = 0; pi <= order_i; ++pi) {

          // Our test polynomial is (pi+1) i^pi
          const auto polyi = [&](T x) { return (pi + 1) * pown(x, pi); };
          // Integrated polynomial
          const auto ipolyi = [&](T x) { return pown(x, pi + 1); };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the coarse grid values)
          const auto cpolyi = [&](T x) {
            return ipolyi(x + 1 / T(2)) - ipolyi(x - 1 / T(2));
          };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the fine grid values)
          const auto fpolyi = [&](T x) {
            return 2 * (ipolyi(x + 1 / T(4)) - ipolyi(x - 1 / T(4)));
          };

          // Our test polynomial is (pj+1) j^pj
          const auto polyj = [&](T y) { return (pj + 1) * pown(y, pj); };
          // Integrated polynomial
          const auto ipolyj = [&](T y) { return pown(y, pj + 1); };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the coarse grid values)
          const auto cpolyj = [&](T y) {
            return ipolyj(y + 1 / T(2)) - ipolyj(y - 1 / T(2));
          };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the fine grid values)
          const auto fpolyj = [&](T y) {
            return 2 * (ipolyj(y + 1 / T(4)) - ipolyj(y - 1 / T(4)));
          };

          // Our test polynomial is (pk+1) k^pk
          const auto polyk = [&](T z) { return (pk + 1) * pown(z, pk); };
          // Integrated polynomial
          const auto ipolyk = [&](T z) { return pown(z, pk + 1); };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the coarse grid values)
          const auto cpolyk = [&](T z) {
            return ipolyk(z + 1 / T(2)) - ipolyk(z - 1 / T(2));
          };
          // Finite difference of the integrated polynomial for the coarse grid
          // (these will be the fine grid values)
          const auto fpolyk = [&](T z) {
            return 2 * (ipolyk(z + 1 / T(4)) - ipolyk(z - 1 / T(4)));
          };

          constexpr bool cons_i = InterpI::stencil_is_conservative;
          constexpr bool cons_j = InterpJ::stencil_is_conservative;
          constexpr bool cons_k = InterpK::stencil_is_conservative;
          const auto cpoly = [&](T x, T y, T z) {
            return (cons_i ? cpolyi(x) : polyi(x))   //
                   + (cons_j ? cpolyj(y) : polyj(y)) //
                   + (cons_k ? cpolyk(z) : polyk(z));
          };
          const auto fpoly = [&](T x, T y, T z) {
            return (cons_i ? fpolyi(x) : polyi(x))   //
                   + (cons_j ? fpolyj(y) : polyj(y)) //
                   + (cons_k ? fpolyk(z) : polyk(z));
          };

          // Set grid to polynomial
          for (int k = -nghosts_k; k <= +nghosts_k; ++k)
            for (int j = -nghosts_j; j <= +nghosts_j; ++j)
              for (int i = -nghosts_i; i <= +nghosts_i; ++i)
                values(i, j, k) = cpoly(2 + i, 2 + j, 2 + k);

          // For cell-centred grids the interpolation points are different than
          // for vertex-centred grids
          constexpr T cc_offset_i = InterpI::stencil_is_cc ? -1 / T(4) : 0;
          constexpr T cc_offset_j = InterpJ::stencil_is_cc ? -1 / T(4) : 0;
          constexpr T cc_offset_k = InterpK::stencil_is_cc ? -1 / T(4) : 0;
          for (int offset_k = 0; offset_k < 2; ++offset_k) {
            for (int offset_j = 0; offset_j < 2; ++offset_j) {
              for (int offset_i = 0; offset_i < 2; ++offset_i) {
                const T x = 2 + offset_i * 1 / T(2) + cc_offset_i;
                const T y = 2 + offset_j * 1 / T(2) + cc_offset_j;
                const T z = 2 + offset_k * 1 / T(2) + cc_offset_k;
                const T expected_result = fpoly(x, y, z);

                // Combine interpolation offsets
                // 0b000 ... 0b111, defined as 0ijk for the three directions i,
                // j, k
                const int offset_ijk =
                    (offset_i << 2) | (offset_j << 1) | (offset_k << 0);

                T result;
                switch (offset_ijk) {
                case 0b000:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 0, 0>(
                          values);
                  break;
                case 0b100:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 0, 0>(
                          values);
                  break;
                case 0b010:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 1, 0>(
                          values);
                  break;
                case 0b110:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 1, 0>(
                          values);
                  break;
                case 0b001:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 0, 1>(
                          values);
                  break;
                case 0b101:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 0, 1>(
                          values);
                  break;
                case 0b011:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 1, 1>(
                          values);
                  break;
                case 0b111:
                  result =
                      combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 1, 1>(
                          values);
                  break;
                default:
                  std::abort();
                }

                // How much error d we allow?
                constexpr T maxrelerr =
                    max({order_i, order_j, order_k}) <= 10 ? 0 : 1.0e-13;

                if (!isfinite(result) || abs(result - expected_result) >
                                             maxrelerr * abs(expected_result)) {
                  found_error = true;
                  std::ios old_state(nullptr);
                  old_state.copyfmt(std::cerr); // save format flags & precision
                  std::cerr
                      << std::setprecision(std::numeric_limits<T>::max_digits10)
                      << "3d interpolation stencil failure:\n"
                      << "  InterpI: name=" << InterpI::name << "\n"
                      << "           order=" << InterpI::order << "\n"
                      << "           centering=" << InterpI::stencil_is_cc
                      << "\n"
                      << "           conservative="
                      << InterpI::stencil_is_conservative << "\n"
                      << "           star=" << InterpI::stencil_is_star << "\n"
                      << "           p=" << pi << "\n"
                      << "           offset=" << offset_i << "\n"
                      << "  InterpJ: name=" << InterpJ::name << "\n"
                      << "           order=" << InterpJ::order << "\n"
                      << "           centering=" << InterpJ::stencil_is_cc
                      << "\n"
                      << "           conservative="
                      << InterpJ::stencil_is_conservative << "\n"
                      << "           star=" << InterpJ::stencil_is_star << "\n"
                      << "           p=" << pj << "\n"
                      << "           offset=" << offset_j << "\n"
                      << "  InterpK: name=" << InterpK::name << "\n"
                      << "           order=" << InterpK::order << "\n"
                      << "           centering=" << InterpK::stencil_is_cc
                      << "\n"
                      << "           conservative="
                      << InterpK::stencil_is_conservative << "\n"
                      << "           star=" << InterpK::stencil_is_star << "\n"
                      << "           p=" << pk << "\n"
                      << "           offset=" << offset_k << "\n"
                      << "  Result:  expected=" << expected_result << "\n"
                      << "           result=" << result << "\n";
                  std::cerr.copyfmt(old_state); // restore
                }

              } // for offset_i
            } // for offset_j
          } // for offset_k

        } // for pi
      } // for pj
    } // for pk

    if (found_error)
      std::abort();
  }
};

} // namespace

////////////////////////////////////////////////////////////////////////////////

template <typename InterpI, typename InterpJ, typename InterpK>
class prolongate_star_3d_rf2 final : public amrex::Interpolater {

  // // Create test objects for all 1d stencils. These tests will run at
  // startup.
  // test_interp<InterpI, CCTK_REAL> test_interpi_CCTK_REAL;
  // test_interp<InterpJ, CCTK_REAL> test_interpj_CCTK_REAL;
  // test_interp<InterpK, CCTK_REAL> test_interpk_CCTK_REAL;

  // Create a test object for these stencils. These tests will run at startup.
  test_interp_3d<InterpI, InterpJ, InterpK, CCTK_REAL> test_interp_3d_CCTK_REAL;

public:
  virtual ~prolongate_star_3d_rf2() override;

  virtual amrex::Box CoarseBox(const amrex::Box &fine, int ratio) override;
  virtual amrex::Box CoarseBox(const amrex::Box &fine,
                               const amrex::IntVect &ratio) override;

  virtual void interp(const amrex::FArrayBox &crse, int crse_comp,
                      amrex::FArrayBox &fine, int fine_comp, int ncomp,
                      const amrex::Box &fine_region,
                      const amrex::IntVect &ratio,
                      const amrex::Geometry &crse_geom,
                      const amrex::Geometry &fine_geom,
                      const amrex::Vector<amrex::BCRec> &bcr, int actual_comp,
                      int actual_state, amrex::RunOn gpu_or_cpu) override;

  virtual void interp_face(const amrex::FArrayBox &crse, int crse_comp,
                           amrex::FArrayBox &fine, int fine_comp, int ncomp,
                           const amrex::Box &fine_region,
                           const amrex::IntVect &ratio,
                           const amrex::IArrayBox &solve_mask,
                           const amrex::Geometry &crse_geom,
                           const amrex::Geometry &fine_geom,
                           const amrex::Vector<amrex::BCRec> &bcr,
                           int actual_bccomp, amrex::RunOn gpu_or_cpu) override;
};

template <typename InterpI, typename InterpJ, typename InterpK>
prolongate_star_3d_rf2<InterpI, InterpJ, InterpK>::~prolongate_star_3d_rf2() {}

template <typename InterpI, typename InterpJ, typename InterpK>
amrex::Box prolongate_star_3d_rf2<InterpI, InterpJ, InterpK>::CoarseBox(
    const amrex::Box &fine, const int ratio) {
  return CoarseBox(fine, amrex::IntVect(ratio));
}

template <typename InterpI, typename InterpJ, typename InterpK>
amrex::Box prolongate_star_3d_rf2<InterpI, InterpJ, InterpK>::CoarseBox(
    const amrex::Box &fine, const amrex::IntVect &ratio) {
  constexpr vect<int, dim> required_ghosts = {
      InterpI::required_ghosts,
      InterpJ::required_ghosts,
      InterpK::required_ghosts,
  };
  for (int d = 0; d < dim; ++d)
    assert(ratio.getVect()[d] == 2);
  amrex::Box crse = amrex::coarsen(fine, 2);
  for (int d = 0; d < dim; ++d)
    crse.grow(d, required_ghosts[d]);
  return crse;
}

template <typename InterpI, typename InterpJ, typename InterpK>
void prolongate_star_3d_rf2<InterpI, InterpJ, InterpK>::interp(
    const amrex::FArrayBox &crse, const int crse_comp, amrex::FArrayBox &fine,
    const int fine_comp, const int ncomp, const amrex::Box &fine_region,
    const amrex::IntVect &ratio, const amrex::Geometry &crse_geom,
    const amrex::Geometry &fine_geom, const amrex::Vector<amrex::BCRec> &bcr,
    const int actual_comp, const int actual_state,
    const amrex::RunOn gpu_or_cpu) {
  assert(actual_comp == 0);  // ???
  assert(actual_state == 0); // ???

  constexpr int NT = AMREX_GPU_MAX_THREADS;

  const amrex::Array4<CCTK_REAL> fine_array = fine.array();
  const amrex::Array4<const CCTK_REAL> crse_array = crse.const_array();

  amrex::ParallelFor<NT>(
      fine_region, ncomp,
      [=] AMREX_GPU_DEVICE(const int i, const int j, const int k,
                           const int comp) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ci = i / 2;
        const int cj = j / 2;
        const int ck = k / 2;
        const int oi = i % 2;
        const int oj = j % 2;
        const int ok = k % 2;

        // Access the coarse grid with various stencil offsets
        const auto get_crse = [&](const int di, const int dj, const int dk) {
          return crse_array(ci + di, cj + dj, ck + dk, crse_comp + comp);
        };

        // Combine interpolation offsets
        // 0b000 ... 0b111, defined as 0ijk for the three directions i, j, k
        const int crse_offset_ijk = (oi << 2) | (oj << 1) | (ok << 0);

        CCTK_REAL result;
        switch (crse_offset_ijk) {
        case 0b000:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 0, 0>(get_crse);
          break;
        case 0b100:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 0, 0>(get_crse);
          break;
        case 0b010:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 1, 0>(get_crse);
          break;
        case 0b110:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 1, 0>(get_crse);
          break;
        case 0b001:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 0, 1>(get_crse);
          break;
        case 0b101:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 0, 1>(get_crse);
          break;
        case 0b011:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 0, 1, 1>(get_crse);
          break;
        case 0b111:
          result =
              combine_1d_stencils<InterpI, InterpJ, InterpK, 1, 1, 1>(get_crse);
          break;
        }

        fine_array(i, j, k, fine_comp + comp) = result;
      });
}

template <typename InterpI, typename InterpJ, typename InterpK>
void prolongate_star_3d_rf2<InterpI, InterpJ, InterpK>::interp_face(
    const amrex::FArrayBox &crse, const int crse_comp, amrex::FArrayBox &fine,
    const int fine_comp, const int ncomp, const amrex::Box &fine_region,
    const amrex::IntVect &ratio, const amrex::IArrayBox &solve_mask,
    const amrex::Geometry &crse_geom, const amrex::Geometry &fine_geom,
    const amrex::Vector<amrex::BCRec> &bcr, const int actual_bccomp,
    const amrex::RunOn gpu_or_cpu) {
  assert(actual_bccomp == 0); // ???
  interp(crse, crse_comp, fine, fine_comp, ncomp, fine_region, ratio, crse_geom,
         fine_geom, bcr, 0, 0, gpu_or_cpu);
}

////////////////////////////////////////////////////////////////////////////////

// Geometry prolongation: higher orders, only vertex centred

static prolongate_star_3d_rf2<interp_vc_poly<1>, interp_vc_poly<1>,
                              interp_vc_poly<1> >
    prolongate_star_poly_3d_rf2_vvv_o1;
static prolongate_star_3d_rf2<interp_vc_poly<3>, interp_vc_poly<3>,
                              interp_vc_poly<3> >
    prolongate_star_poly_3d_rf2_vvv_o3;
static prolongate_star_3d_rf2<interp_vc_poly<5>, interp_vc_poly<5>,
                              interp_vc_poly<5> >
    prolongate_star_poly_3d_rf2_vvv_o5;
static prolongate_star_3d_rf2<interp_vc_poly<7>, interp_vc_poly<7>,
                              interp_vc_poly<7> >
    prolongate_star_poly_3d_rf2_vvv_o7;
static prolongate_star_3d_rf2<interp_vc_poly<9>, interp_vc_poly<9>,
                              interp_vc_poly<9> >
    prolongate_star_poly_3d_rf2_vvv_o9;
static prolongate_star_3d_rf2<interp_vc_poly<11>, interp_vc_poly<11>,
                              interp_vc_poly<11> >
    prolongate_star_poly_3d_rf2_vvv_o11;

const interpolator_family_t prolongate_star_poly_3d_rf2{
    {1,
     {
         &prolongate_star_poly_3d_rf2_vvv_o1,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
    {3,
     {
         &prolongate_star_poly_3d_rf2_vvv_o3,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
    {5,
     {
         &prolongate_star_poly_3d_rf2_vvv_o5,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
    {7,
     {
         &prolongate_star_poly_3d_rf2_vvv_o7,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
    {9,
     {
         &prolongate_star_poly_3d_rf2_vvv_o9,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
    {11,
     {
         &prolongate_star_poly_3d_rf2_vvv_o11,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
         nullptr,
     }},
};

// Hydrodynamics prolongation for regridding and for fluxes: Conservative,
// only cell or face centred

static prolongate_star_3d_rf2<interp_cc_cons<0>, interp_cc_cons<0>,
                              interp_cc_cons<0> >
    prolongate_star_cons_3d_rf2_ccc_o0;
static prolongate_star_3d_rf2<interp_vc_poly<1>, interp_cc_cons<0>,
                              interp_cc_cons<0> >
    prolongate_star_cons_3d_rf2_vcc_o0;
static prolongate_star_3d_rf2<interp_cc_cons<0>, interp_vc_poly<1>,
                              interp_cc_cons<0> >
    prolongate_star_cons_3d_rf2_cvc_o0;
static prolongate_star_3d_rf2<interp_cc_cons<0>, interp_cc_cons<0>,
                              interp_vc_poly<1> >
    prolongate_star_cons_3d_rf2_ccv_o0;

// We use minmod for order 1
static prolongate_star_3d_rf2<interp_cc_cons_minmod<1>,
                              interp_cc_cons_minmod<1>,
                              interp_cc_cons_minmod<1> >
    prolongate_star_cons_3d_rf2_ccc_o1;
static prolongate_star_3d_rf2<interp_vc_poly<1>, interp_cc_cons_minmod<1>,
                              interp_cc_cons_minmod<1> >
    prolongate_star_cons_3d_rf2_vcc_o1;
static prolongate_star_3d_rf2<interp_cc_cons_minmod<1>, interp_vc_poly<1>,
                              interp_cc_cons_minmod<1> >
    prolongate_star_cons_3d_rf2_cvc_o1;
static prolongate_star_3d_rf2<interp_cc_cons_minmod<1>,
                              interp_cc_cons_minmod<1>, interp_vc_poly<1> >
    prolongate_star_cons_3d_rf2_ccv_o1;

// We use ENO for orders 2 and 3. We probably need to introduce a fallback here.
static prolongate_star_3d_rf2<interp_cc_cons_eno<2>, interp_cc_cons_eno<2>,
                              interp_cc_cons_eno<2> >
    prolongate_star_cons_3d_rf2_ccc_o2;
static prolongate_star_3d_rf2<interp_vc_poly<3>, interp_cc_cons_eno<2>,
                              interp_cc_cons_eno<2> >
    prolongate_star_cons_3d_rf2_vcc_o2;
static prolongate_star_3d_rf2<interp_cc_cons_eno<2>, interp_vc_poly<3>,
                              interp_cc_cons_eno<2> >
    prolongate_star_cons_3d_rf2_cvc_o2;
static prolongate_star_3d_rf2<interp_cc_cons_eno<2>, interp_cc_cons_eno<2>,
                              interp_vc_poly<3> >
    prolongate_star_cons_3d_rf2_ccv_o2;

static prolongate_star_3d_rf2<interp_cc_cons_eno<3>, interp_cc_cons_eno<3>,
                              interp_cc_cons_eno<3> >
    prolongate_star_cons_3d_rf2_ccc_o3;
static prolongate_star_3d_rf2<interp_vc_poly<3>, interp_cc_cons_eno<3>,
                              interp_cc_cons_eno<3> >
    prolongate_star_cons_3d_rf2_vcc_o3;
static prolongate_star_3d_rf2<interp_cc_cons_eno<3>, interp_vc_poly<3>,
                              interp_cc_cons_eno<3> >
    prolongate_star_cons_3d_rf2_cvc_o3;
static prolongate_star_3d_rf2<interp_cc_cons_eno<3>, interp_cc_cons_eno<3>,
                              interp_vc_poly<3> >
    prolongate_star_cons_3d_rf2_ccv_o3;

const interpolator_family_t prolongate_star_cons_3d_rf2{
    {0,
     {
         nullptr,
         nullptr,
         nullptr,
         &prolongate_star_cons_3d_rf2_vcc_o0,
         nullptr,
         &prolongate_star_cons_3d_rf2_cvc_o0,
         &prolongate_star_cons_3d_rf2_ccv_o0,
         &prolongate_star_cons_3d_rf2_ccc_o0,
     }},
    {1,
     {
         nullptr,
         nullptr,
         nullptr,
         &prolongate_star_cons_3d_rf2_vcc_o1,
         nullptr,
         &prolongate_star_cons_3d_rf2_cvc_o1,
         &prolongate_star_cons_3d_rf2_ccv_o1,
         &prolongate_star_cons_3d_rf2_ccc_o1,
     }},
    {2,
     {
         nullptr,
         nullptr,
         nullptr,
         &prolongate_star_cons_3d_rf2_vcc_o2,
         nullptr,
         &prolongate_star_cons_3d_rf2_cvc_o2,
         &prolongate_star_cons_3d_rf2_ccv_o2,
         &prolongate_star_cons_3d_rf2_ccc_o2,
     }},
    {3,
     {
         nullptr,
         nullptr,
         nullptr,
         &prolongate_star_cons_3d_rf2_vcc_o3,
         nullptr,
         &prolongate_star_cons_3d_rf2_cvc_o3,
         &prolongate_star_cons_3d_rf2_ccv_o3,
         &prolongate_star_cons_3d_rf2_ccc_o3,
     }},
};

} // namespace CarpetX
