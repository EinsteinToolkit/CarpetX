#include <loop.hxx>
#include <subcycling_native_gate.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

namespace TestODESolvers2 {
using namespace std;

namespace {
optional<CarpetX::NativeGateSession> native_gate_session;
}

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

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_(cctkGH, time);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_(cctkGH, poly);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp1_(cctkGH, exp1);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp2_(cctkGH, exp2);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_dep_(cctkGH, time_dep);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_dep_(cctkGH, poly_dep);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp1_dep_(cctkGH, exp1_dep);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp2_dep_(cctkGH, exp2_dep);

  const auto set_boundary_state = [&](const Loop::PointDesc &p) {
    time_(p.I) = cctk_time;
    poly_(p.I) = pow(1 + cctk_time, porder);
    exp1_(p.I) = exp(cctk_time);
    exp2_(p.I) = exp(cctk_time / 2);
  };
  Loop::loop_bnd<1, 1, 1>(cctkGH, set_boundary_state);
  Loop::loop_ghosts<1, 1, 1>(cctkGH, set_boundary_state);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    time_dep_(p.I) = time_(p.I) + cctk_time;
    poly_dep_(p.I) = poly_(p.I) + cctk_time;
    exp1_dep_(p.I) = exp1_(p.I) + cctk_time;
    exp2_dep_(p.I) = exp2_(p.I) + cctk_time;
  });
}

extern "C" void TestODESolvers2_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_RHS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> time_(cctkGH, time);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> poly_(cctkGH, poly);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp1_(cctkGH, exp1);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp2_(cctkGH, exp2);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> time_dep_(cctkGH, time_dep);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> poly_dep_(cctkGH, poly_dep);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp1_dep_(cctkGH, exp1_dep);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> exp2_dep_(cctkGH, exp2_dep);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> time_rhs_(cctkGH, time_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> poly_rhs_(cctkGH, poly_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp1_rhs_(cctkGH, exp1_rhs);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> exp2_rhs_(cctkGH, exp2_rhs);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const auto tolerance =
        32 * numeric_limits<CCTK_REAL>::epsilon() *
        max({CCTK_REAL(1), abs(cctk_time), abs(time_(p.I)),
             abs(poly_(p.I)), abs(exp1_(p.I)), abs(exp2_(p.I))});
    if (abs(time_dep_(p.I) - (time_(p.I) + cctk_time)) > tolerance ||
        abs(poly_dep_(p.I) - (poly_(p.I) + cctk_time)) > tolerance ||
        abs(exp1_dep_(p.I) - (exp1_(p.I) + cctk_time)) > tolerance ||
        abs(exp2_dep_(p.I) - (exp2_(p.I) + cctk_time)) > tolerance)
      CCTK_VERROR("Native gate PostStep dependent is stale at t=%.17g",
                  double(cctk_time));
    const auto time_tolerance =
        32 * numeric_limits<CCTK_REAL>::epsilon() *
        max({CCTK_REAL(1), abs(cctk_time), abs(time_(p.I))});
    if (porder > 0)
      if (abs(time_(p.I) - cctk_time) > time_tolerance)
        CCTK_VERROR("Time is incorrect: time=%.17g, cctk_time=%.17g",
                    double(time_(p.I)), double(cctk_time));
    time_rhs_(p.I) = 1;
    poly_rhs_(p.I) = porder == 0 ? 0 : porder * pow(1 + cctk_time, porder - 1);
    exp1_rhs_(p.I) = exp1_(p.I);
    exp2_rhs_(p.I) = exp2_(p.I) / 2;
  });
}

extern "C" void TestODESolvers2_NativeGateInventory(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_NativeGateInventory;
  DECLARE_CCTK_PARAMETERS;
  if (!native_subcycling_gate || cctk_itlast != 0)
    return;
  try {
    CarpetX::write_native_gate_inventory(cctkGH, cctk_itlast);
    CCTK_INFO("CARPETX_NATIVE_GATE_INVENTORY_PASS");
  } catch (const exception &error) {
    CCTK_VERROR("Native gate inventory failed: %s", error.what());
  } catch (...) {
    CCTK_ERROR("Native gate inventory failed with an unknown exception");
  }
}

extern "C" void TestODESolvers2_NativeGateBegin(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_NativeGateBegin;
  DECLARE_CCTK_PARAMETERS;
  if (!native_subcycling_gate)
    return;
  try {
    if (native_gate_session.has_value())
      throw logic_error("native gate session was opened twice");
    const auto contract =
        CarpetX::native_gate_method_contract(cctk_iteration);
    const int status = CCTK_ParameterSet(
        "method", "ODESolvers", contract.ode_parameter_value);
    if (status != 0)
      throw runtime_error("CCTK_ParameterSet(ODESolvers::method) returned " +
                          to_string(status));
    native_gate_session.emplace(
        CarpetX::begin_native_gate(cctkGH, cctk_itlast));
  } catch (const exception &error) {
    native_gate_session.reset();
    CCTK_VERROR("Native gate begin failed: %s", error.what());
  } catch (...) {
    native_gate_session.reset();
    CCTK_ERROR("Native gate begin failed with an unknown exception");
  }
}

extern "C" void TestODESolvers2_NativeGateEnd(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers2_NativeGateEnd;
  DECLARE_CCTK_PARAMETERS;
  if (!native_subcycling_gate)
    return;
  try {
    if (!native_gate_session.has_value())
      throw logic_error("native gate session is absent");
    const auto receipt =
        CarpetX::end_native_gate(std::move(*native_gate_session));
    native_gate_session.reset();
    CCTK_VINFO("CARPETX_NATIVE_GATE_METHOD_PASS rhs=%zu controls=%zu",
               receipt.extra_rhs_evaluations, receipt.control_count);
    if (cctk_iteration == 3)
      CCTK_INFO("CARPETX_NATIVE_GATE_PASS");
  } catch (const exception &error) {
    native_gate_session.reset();
    CCTK_VERROR("Native gate end failed: %s", error.what());
  } catch (...) {
    native_gate_session.reset();
    CCTK_ERROR("Native gate end failed with an unknown exception");
  }
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
