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

template <int CI, int CJ, int CK>
bool is_near_outer_boundary(const cGH *const cctkGH, const Loop::PointDesc &p) {
  return (cctkGH->cctk_bbox[0] && p.i < 2 * cctkGH->cctk_nghostzones[0]) ||
         (cctkGH->cctk_bbox[2] && p.j < 2 * cctkGH->cctk_nghostzones[1]) ||
         (cctkGH->cctk_bbox[4] && p.k < 2 * cctkGH->cctk_nghostzones[2]) ||
         (cctkGH->cctk_bbox[1] &&
          p.i >= cctkGH->cctk_lsh[0] - 2 * cctkGH->cctk_nghostzones[0] - CI) ||
         (cctkGH->cctk_bbox[3] &&
          p.j >= cctkGH->cctk_lsh[1] - 2 * cctkGH->cctk_nghostzones[1] - CJ) ||
         (cctkGH->cctk_bbox[5] &&
          p.k >= cctkGH->cctk_lsh[2] - 2 * cctkGH->cctk_nghostzones[2] - CK);
}

bool is_near_jump(const Loop::PointDesc &p, const bool enox, const bool enoy,
                  const bool enoz) {
  // Stay away from the discontinuity.
  // `1 dx` on the coarse grid is `2 dx` on the fine grid. The coarse grid is
  // also offset by 1/2 grid spacing if it's cell centred.
  return (enox && fabs(p.x) < 2.1 * p.dx) || (enoy && fabs(p.y) < 2.1 * p.dy) ||
         (enoz && fabs(p.z) < 2.1 * p.dz);
}

template <typename T> T fun1d_base(const T x, const int order) {
  return order == 0 ? 1 : pow(x, order);
}

// The `DDF` operators have order `order` for VC, but only `order-1` for CC. The
// caller needs to handle this.
template <typename T>
T fun1d(const T x, const T dx, const bool consx, const int orderx) {
  T res = 0;
  for (int p = 0; p <= orderx; ++p)
    if (!consx)
      res += fun1d_base(x, p);
    else
      res += (fun1d_base(x + dx / 2, p + 1) - fun1d_base(x - dx / 2, p + 1)) /
             ((p + 1) * dx);
  return res;
}

// The grid stores the average values of the underlying function (x*y*z)**n
// which amounts to storing differences of the anti-derivative 1/(n+1)**3 *
// (x*y*z)**(n+1)
template <typename T>
T fun(const T x, const T y, const T z, const T dx, const T dy, const T dz,
      const bool consx, const bool consy, const bool consz, const bool enox,
      const bool enoy, const bool enoz, const int orderx, const int ordery,
      const int orderz) {
  return fun1d(x, dx, consx, orderx) * fun1d(y, dy, consy, ordery) *
             fun1d(z, dz, consz, orderz) +
         ((enox * signbit(x)) ^ (enoy * signbit(y)) ^ (enoz * signbit(z)));
}

extern "C" void TestProlongate_Set(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongate_Set;
  DECLARE_CCTK_PARAMETERS;

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
  bool vc_cons, cc_cons, cc_eno;
  constexpr bool vc_eno = false;
  int vc_order, cc_order;
  if (CCTK_EQUALS(prolongation_type, "interpolate")) {
    vc_cons = false;
    cc_cons = false;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order;
  } else if (CCTK_EQUALS(prolongation_type, "conservative")) {
    vc_cons = true;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order == 0 ? 0 : operator_order + 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf-eno")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = true;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf-hermite")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else {
    assert(0);
  }

  Loop::loop_int<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, vc_cons, vc_cons, vc_eno,
            vc_eno, vc_eno, vc_order, vc_order, vc_order);
    gf000(p.I) = good_data;
  });
  Loop::loop_int<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, vc_cons, cc_cons, vc_eno,
            vc_eno, cc_eno, vc_order, vc_order, cc_order);
    gf001(p.I) = good_data;
  });
  Loop::loop_int<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, cc_cons, vc_cons, vc_eno,
            cc_eno, vc_eno, vc_order, cc_order, vc_order);
    gf010(p.I) = good_data;
  });
  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, cc_cons, cc_cons, vc_eno,
            cc_eno, cc_eno, vc_order, cc_order, cc_order);
    gf011(p.I) = good_data;
  });
  Loop::loop_int<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, vc_cons, vc_cons, cc_eno,
            vc_eno, vc_eno, cc_order, vc_order, vc_order);
    gf100(p.I) = good_data;
  });
  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, vc_cons, cc_cons, cc_eno,
            vc_eno, cc_eno, cc_order, vc_order, cc_order);
    gf101(p.I) = good_data;
  });
  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, cc_cons, vc_cons, cc_eno,
            cc_eno, vc_eno, cc_order, cc_order, vc_order);
    gf110(p.I) = good_data;
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL good_data =
        fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, cc_cons, cc_cons, cc_eno,
            cc_eno, cc_eno, cc_order, cc_order, cc_order);
    gf111(p.I) = good_data;
  });

#pragma omp critical
  {
    *gf000_max_error = 0;
    *gf001_max_error = 0;
    *gf010_max_error = 0;
    *gf011_max_error = 0;
    *gf100_max_error = 0;
    *gf101_max_error = 0;
    *gf110_max_error = 0;
    *gf111_max_error = 0;

    *gf000_sum_count = 0;
    *gf001_sum_count = 0;
    *gf010_sum_count = 0;
    *gf011_sum_count = 0;
    *gf100_sum_count = 0;
    *gf101_sum_count = 0;
    *gf110_sum_count = 0;
    *gf111_sum_count = 0;
  }
}

extern "C" void TestProlongate_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongate_Sync;
  DECLARE_CCTK_PARAMETERS;

  // do nothing
}

extern "C" void TestProlongate_Regrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongate_Regrid;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Setting grid at iteration %d", cctk_iteration);
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (cctk_iteration < regrid_after) {
      regrid_error(p.I) = 0;
    } else {
      if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
          fabs(p.z) <= refined_radius) {
        regrid_error(p.I) = 1;
      } else {
        regrid_error(p.I) = 0;
      }
    }
  });
}

extern "C" void TestProlongate_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestProlongate_Check;
  DECLARE_CCTK_PARAMETERS;

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
  bool vc_cons, cc_cons, cc_eno;
  constexpr bool vc_eno = false;
  int vc_order, cc_order;
  if (CCTK_EQUALS(prolongation_type, "interpolate")) {
    vc_cons = false;
    cc_cons = false;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order;
  } else if (CCTK_EQUALS(prolongation_type, "conservative")) {
    vc_cons = true;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order == 0 ? 0 : operator_order + 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf-eno")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = true;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else if (CCTK_EQUALS(prolongation_type, "ddf-hermite")) {
    vc_cons = false;
    cc_cons = true;
    cc_eno = false;
    vc_order = operator_order;
    cc_order = operator_order - 1;
  } else {
    assert(0);
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<0, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<0, 0, 0>(cctkGH, p) ||
          is_near_jump(p, vc_eno, vc_eno, vc_eno)) {
        gf000_error(p.I) = 0;
        gf000_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, vc_cons, vc_cons,
                vc_eno, vc_eno, vc_eno, vc_order, vc_order, vc_order);
        gf000_error(p.I) = gf000(p.I) - good_data;
        gf000_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf000_error(p.I)));
        sum_count += gf000_count(p.I);
      }
    });
#pragma omp critical
    *gf000_max_error = fmax(*gf000_max_error, max_error);
#pragma omp atomic
    *gf000_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<0, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<0, 0, 1>(cctkGH, p) ||
          is_near_jump(p, vc_eno, vc_eno, cc_eno)) {
        gf001_error(p.I) = 0;
        gf001_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, vc_cons, cc_cons,
                vc_eno, vc_eno, cc_eno, vc_order, vc_order, cc_order);
        gf001_error(p.I) = gf001(p.I) - good_data;
        gf001_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf001_error(p.I)));
        sum_count += gf001_count(p.I);
      }
    });
#pragma omp critical
    *gf001_max_error = fmax(*gf001_max_error, max_error);
#pragma omp atomic
    *gf001_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<0, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<0, 1, 0>(cctkGH, p) ||
          is_near_jump(p, vc_eno, cc_eno, vc_eno)) {
        gf010_error(p.I) = 0;
        gf010_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, cc_cons, vc_cons,
                vc_eno, cc_eno, vc_eno, vc_order, cc_order, vc_order);
        gf010_error(p.I) = gf010(p.I) - good_data;
        gf010_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf010_error(p.I)));
        sum_count += gf010_count(p.I);
      }
    });
#pragma omp critical
    *gf010_max_error = fmax(*gf010_max_error, max_error);
#pragma omp atomic
    *gf010_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<0, 1, 1>(cctkGH, p) ||
          is_near_jump(p, vc_eno, cc_eno, cc_eno)) {
        gf011_error(p.I) = 0;
        gf011_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, vc_cons, cc_cons, cc_cons,
                vc_eno, cc_eno, cc_eno, vc_order, cc_order, cc_order);
        gf011_error(p.I) = gf011(p.I) - good_data;
        gf011_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf011_error(p.I)));
        sum_count += gf011_count(p.I);
      }
    });
#pragma omp critical
    *gf011_max_error = fmax(*gf011_max_error, max_error);
#pragma omp atomic
    *gf011_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<1, 0, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<1, 0, 0>(cctkGH, p) ||
          is_near_jump(p, cc_eno, vc_eno, vc_eno)) {
        gf100_error(p.I) = 0;
        gf100_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, vc_cons, vc_cons,
                cc_eno, vc_eno, vc_eno, cc_order, vc_order, vc_order);
        gf100_error(p.I) = gf100(p.I) - good_data;
        gf100_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf100_error(p.I)));
        sum_count += gf100_count(p.I);
      }
    });
#pragma omp critical
    *gf100_max_error = fmax(*gf100_max_error, max_error);
#pragma omp atomic
    *gf100_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<1, 0, 1>(cctkGH, p) ||
          is_near_jump(p, cc_eno, vc_eno, cc_eno)) {
        gf101_error(p.I) = 0;
        gf101_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, vc_cons, cc_cons,
                cc_eno, vc_eno, cc_eno, cc_order, vc_order, cc_order);
        gf101_error(p.I) = gf101(p.I) - good_data;
        gf101_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf101_error(p.I)));
        sum_count += gf101_count(p.I);
      }
    });
#pragma omp critical
    *gf101_max_error = fmax(*gf101_max_error, max_error);
#pragma omp atomic
    *gf101_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<1, 1, 0>(cctkGH, p) ||
          is_near_jump(p, cc_eno, cc_eno, vc_eno)) {
        gf110_error(p.I) = 0;
        gf110_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, cc_cons, vc_cons,
                cc_eno, cc_eno, vc_eno, cc_order, cc_order, vc_order);
        gf110_error(p.I) = gf110(p.I) - good_data;
        gf110_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf110_error(p.I)));
        sum_count += gf110_count(p.I);
      }
    });
#pragma omp critical
    *gf110_max_error = fmax(*gf110_max_error, max_error);
#pragma omp atomic
    *gf110_sum_count += sum_count;
  }

  {
    CCTK_REAL max_error = 0, sum_count = 0;
    Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (is_near_outer_boundary<1, 1, 1>(cctkGH, p) ||
          is_near_jump(p, cc_eno, cc_eno, cc_eno)) {
        gf111_error(p.I) = 0;
        gf111_count(p.I) = 0;
      } else {
        const CCTK_REAL good_data =
            fun(p.x, p.y, p.z, p.dx, p.dy, p.dz, cc_cons, cc_cons, cc_cons,
                cc_eno, cc_eno, cc_eno, cc_order, cc_order, cc_order);
        gf111_error(p.I) = gf111(p.I) - good_data;
        gf111_count(p.I) = 1;
        max_error = fmax(max_error, fabs(gf111_error(p.I)));
        sum_count += gf111_count(p.I);
      }
    });
#pragma omp critical
    *gf111_max_error = fmax(*gf111_max_error, max_error);
#pragma omp atomic
    *gf111_sum_count += sum_count;
  }
}

} // namespace TestProlongate
