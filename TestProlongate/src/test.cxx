#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>

namespace TestProlongate {
using namespace std;

template <typename T>
T fun1d(const T x, const T dx, const bool avgx, const int order) {
  if (!avgx)
    return pow(x, order);
  else
    return (pow(x + dx / 2, order + 1) - pow(x - dx / 2, order + 1)) /
           ((order + 1) * dx);
}

// The grid stores the average values of the underlying function (x*y*z)**n
// which amounts to storing differences of the anti-derivative 1/(n+1)**3 *
// (x*y*z)**(n+1)
template <typename T>
T fun(const T x, const T y, const T z, const T dx, const T dy, const T dz,
      const bool avgx, const bool avgy, const bool avgz, const int order) {
  return fun1d(x, dx, avgx, order) * fun1d(y, dy, avgy, order) *
         fun1d(z, dz, avgz, order);
}

extern "C" void TestProlongate_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // Get prolongation order from driver, the parmeter is private since really
  // there is normally no reason to depend on it
  int order_type;
  const void *const order_p =
      CCTK_ParameterGet("prolongation_order", "CarpetX", &order_type);
  assert(order_p);
  assert(order_type == PARAMETER_INT);
  const CCTK_INT operator_order = *static_cast<const CCTK_INT *>(order_p);

  int prolongation_type_type;
  const void *const prolongation_type_p = CCTK_ParameterGet(
      "prolongation_type", "CarpetX", &prolongation_type_type);
  assert(prolongation_type_p);
  assert(prolongation_type_type == PARAMETER_KEYWORD);
  const char *const prolongation_type =
      *static_cast<const char *const *>(prolongation_type_p);
  bool conservative_prolongation;
  if (CCTK_EQUALS(prolongation_type, "interpolate"))
    conservative_prolongation = false;
  else if (CCTK_EQUALS(prolongation_type, "conservative"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf-eno"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf-hermite"))
    conservative_prolongation = true;
  else
    assert(0);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(
        p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
        conservative_prolongation, conservative_prolongation, operator_order);
    data[p.idx] = good_data;
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false,
                                    false, false, operator_order);
    gf000[p.idx] = good_data;
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, false,
            conservative_prolongation, operator_order);
    gf001[p.idx] = good_data;
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, conservative_prolongation,
            false, operator_order);
    gf010[p.idx] = good_data;
  });
  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, conservative_prolongation,
            conservative_prolongation, operator_order);
    gf011[p.idx] = good_data;
  });
  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation, false,
            false, operator_order);
    gf100[p.idx] = good_data;
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation, false,
            conservative_prolongation, operator_order);
    gf101[p.idx] = good_data;
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
            conservative_prolongation, false, operator_order);
    gf110[p.idx] = good_data;
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(
        p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
        conservative_prolongation, conservative_prolongation, operator_order);
    gf111[p.idx] = good_data;
  });

#pragma omp critical
  *max_diff = 0;
}

extern "C" void TestProlongate_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  return; // do nothing
}

extern "C" void TestProlongate_Regrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Setting grid at iteration %d", cctk_iteration);
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (cctk_iteration < regrid_after) {
      regrid_error[p.idx] = 0;
    } else {
      if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
          fabs(p.z) <= refined_radius) {
        regrid_error[p.idx] = 1;
      } else {
        regrid_error[p.idx] = 0;
      }
    }
  });
}

extern "C" void TestProlongate_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // there is no reduction yet, so for now I require that only 1 MPI rank be in
  // use
  assert(CCTK_nProcs(cctkGH) == 1 && "Must be single rank");

  // get prolongation order from driver, the parmeter is private since really
  // there is normally no reason to depend on it
  int order_type;
  const void *const order_p =
      CCTK_ParameterGet("prolongation_order", "CarpetX", &order_type);
  assert(order_p);
  assert(order_type == PARAMETER_INT);
  const CCTK_INT operator_order = *static_cast<const CCTK_INT *>(order_p);

  int prolongation_type_type;
  const void *const prolongation_type_p = CCTK_ParameterGet(
      "prolongation_type", "CarpetX", &prolongation_type_type);
  assert(prolongation_type_p);
  assert(prolongation_type_type == PARAMETER_KEYWORD);
  const char *const prolongation_type =
      *static_cast<const char *const *>(prolongation_type_p);
  bool conservative_prolongation;
  if (CCTK_EQUALS(prolongation_type, "interpolate"))
    conservative_prolongation = false;
  else if (CCTK_EQUALS(prolongation_type, "conservative"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf-eno"))
    conservative_prolongation = true;
  else if (CCTK_EQUALS(prolongation_type, "ddf-hermite"))
    conservative_prolongation = true;
  else
    assert(0);

  CCTK_REAL my_max_diff = 0;

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(
        p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
        conservative_prolongation, conservative_prolongation, operator_order);
    const CCTK_REAL diff = fabs(data[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false,
                                    false, false, operator_order);
    const CCTK_REAL diff = fabs(gf000[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, false,
            conservative_prolongation, operator_order);
    const CCTK_REAL diff = fabs(gf001[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, conservative_prolongation,
            false, operator_order);
    const CCTK_REAL diff = fabs(gf010[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, false, conservative_prolongation,
            conservative_prolongation, operator_order);
    const CCTK_REAL diff = fabs(gf011[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation, false,
            false, operator_order);
    const CCTK_REAL diff = fabs(gf100[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation, false,
            conservative_prolongation, operator_order);
    const CCTK_REAL diff = fabs(gf101[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
            conservative_prolongation, false, operator_order);
    const CCTK_REAL diff = fabs(gf110[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data = fun(
        p.x, p.y, p.z, p.dx, p.dy, p.dz, conservative_prolongation,
        conservative_prolongation, conservative_prolongation, operator_order);
    const CCTK_REAL diff = fabs(gf111[p.idx] - good_data);
    my_max_diff = fmax(diff, my_max_diff);
  });

#pragma omp critical
  *max_diff = fmax(*max_diff, my_max_diff);
}

} // namespace TestProlongate
