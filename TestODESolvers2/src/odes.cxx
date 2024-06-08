#include <loop.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <cmath>
#include <limits>

namespace TestODESolvers2 {
using namespace std;

#if 0
// pos(t) = sin(omega t)
// vel(t) = omega cos(omega t)

// d/dt pos = vel
// d/dt vel = -omega^2 pos

const CCTK_REAL omega = 1;
#endif

// u(t) = (1 + t)^p
// d/dt u = p (1 + t)^(p-1) = p u(t)^((p-1) / p)

// v(t) = exp alpha t
// d/dt v = alpha exp alpha t = alpha v

////////////////////////////////////////////////////////////////////////////////

extern "C" void TestODESolvers2_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_Initial;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_(cctkGH, time);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_(cctkGH, poly);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp1_(cctkGH, exp1);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp2_(cctkGH, exp2);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    time_(p.I) = cctk_time;
    poly_(p.I) = pow(1 + cctk_time, porder);
    exp1_(p.I) = exp(cctk_time);
    exp2_(p.I) = exp(cctk_time / 2);
  });
}

extern "C" void TestODESolvers2_Boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_Boundary;
  DECLARE_CCTK_PARAMETERS;

  // do nothing
}

extern "C" void TestODESolvers2_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_RHS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> time_(cctkGH, time);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp1_(cctkGH, exp1);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp2_(cctkGH, exp2);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_rhs_(cctkGH, time_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_rhs_(cctkGH, poly_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp1_rhs_(cctkGH, exp1_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp2_rhs_(cctkGH, exp2_rhs);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (porder > 0)
      if (abs(time_(p.I) - cctk_time) >
          10 * numeric_limits<CCTK_REAL>::epsilon())
        CCTK_VERROR("Time is incorrect: time=%.17g, cctk_time=%.17g",
                    double(time_(p.I)), double(cctk_time));
    time_rhs_(p.I) = 1;
    poly_rhs_(p.I) = porder == 0 ? 0 : porder * pow(1 + cctk_time, porder - 1);
    exp1_rhs_(p.I) = exp1_(p.I);
    exp2_rhs_(p.I) = exp2_(p.I) / 2;
  });
}

extern "C" void TestODESolvers2_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_Error;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> time_(cctkGH, time);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> poly_(cctkGH, poly);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp1_(cctkGH, exp1);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp2_(cctkGH, exp2);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_err_(cctkGH, time_err);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_err_(cctkGH, poly_err);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp_order_(cctkGH, exp_order);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    time_err_(p.I) = time_(p.I) - cctk_time;
    poly_err_(p.I) = poly_(p.I) - pow(1 + cctk_time, porder);
    if (cctk_iteration == 0) {
      exp_order_(p.I) = 0; // undefined
    } else {
      const CCTK_REAL exp1_err = exp1_(p.I) - exp(cctk_time);
      const CCTK_REAL exp2_err = exp2_(p.I) - exp(cctk_time / 2);
      const CCTK_REAL order = log(exp1_err / exp2_err) / log(CCTK_REAL(2)) - 1;
      // Re-map the order so that the error is of the same magnitude as the
      // floating-point epsilon instead of cctk_delta_time
      const CCTK_REAL order1 = pow(order - porder, porder + 2) + porder;
      exp_order_(p.I) = order1;
      CCTK_VINFO("exp1_err=%.17g exp2_err=%.17g order=%.17g", double(exp1_err),
                 double(exp2_err), double(exp_order_(p.I)));
    }
  });
}

} // namespace TestODESolvers2
