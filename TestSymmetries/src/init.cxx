#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>

namespace TestSymmetries {

bool get_parameter(const std::string &name, const std::string &thorn) {
  int type;
  const void *const ptr = CCTK_ParameterGet(name.c_str(), thorn.c_str(), &type);
  assert(ptr);
  assert(type == PARAMETER_BOOLEAN);
  const bool val = bool(*static_cast<const int *>(ptr));
  return val;
}

template <int CI, int CJ, int CK, Loop::where_t where, typename F>
void map_centering_parity(const cGH *restrict const cctkGH,
                          const std::array<int, 3> &parity, const F &f) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const std::array<bool, 3> periodic{
      get_parameter("periodic_x", "CarpetX"),
      get_parameter("periodic_y", "CarpetX"),
      get_parameter("periodic_z", "CarpetX"),
  };
  const std::array<bool, 3> reflection_lower{
      get_parameter("reflection_x", "CarpetX"),
      get_parameter("reflection_y", "CarpetX"),
      get_parameter("reflection_z", "CarpetX"),
  };
  const std::array<bool, 3> reflection_upper{
      get_parameter("reflection_upper_x", "CarpetX"),
      get_parameter("reflection_upper_y", "CarpetX"),
      get_parameter("reflection_upper_z", "CarpetX"),
  };

  // We assume the domain is [-1; +1]

  const auto makevalue =
      [=] CCTK_DEVICE CCTK_HOST(const std::array<CCTK_REAL, Loop::dim> &x)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            using std::acos, std::cos, std::sin;
            const CCTK_REAL pi = acos(CCTK_REAL(-1));
            CCTK_REAL value = 1;
            for (int d = 0; d < Loop::dim; ++d) {
              CCTK_REAL val;
              if (periodic[d]) {
                val = cos(pi * (x[d] + CCTK_REAL(1) / 8));
              } else {
                if (parity[d] > 0) {
                  // even parity (scalar): f(x) = f(-x)
                  if (!reflection_lower[d] && !reflection_upper[d])
                    val = x[d] + 2;
                  else if (reflection_lower[d] && !reflection_upper[d])
                    val = cos(x[d] + 1);
                  else if (!reflection_lower[d] && reflection_upper[d])
                    val = cos(x[d] - 1);
                  else
                    val = cos(pi * (x[d] + 1));
                } else {
                  // odd parity (normal vector): f(x) = -f(-x)
                  if (!reflection_lower[d] && !reflection_upper[d])
                    val = x[d] + 2;
                  else if (reflection_lower[d] && !reflection_upper[d])
                    val = sin(x[d] + 1);
                  else if (!reflection_lower[d] && reflection_upper[d])
                    val = sin(x[d] - 1);
                  else
                    val = cos(pi / 2 * x[d]);
                }
              }
              value *= val;
            }
            return value;
          };

  // Test makevalue function
  {
    // Test symmetries
    for (int k = -3; k <= 3; ++k) {
      for (int j = -3; j <= 3; ++j) {
        for (int i = -3; i <= 3; ++i) {
          const std::array<CCTK_REAL, 3> x{
              CCTK_REAL(0.5) * i, CCTK_REAL(0.5) * j, CCTK_REAL(0.5) * k};
          std::array<CCTK_REAL, 3> x0 = x;
          CCTK_REAL f = makevalue(x);
          // map f back into the domain
          for (int d = 0; d < 3; ++d) {
            if (periodic[d]) {
              x0[d] = x0[d] < -1 ? x0[d] + 2 : x0[d] >= +1 ? x0[d] - 2 : x0[d];
              assert(x0[d] >= -1 && x0[d] < +1);
              f = f;
            } else {
              if (reflection_lower[d]) {
                if (x0[d] < -1) {
                  x0[d] = -(x0[d] + 1) - 1;
                  f *= parity[d];
                }
                assert(x0[d] >= -1);
              }
              if (reflection_upper[d]) {
                if (x0[d] > +1) {
                  x0[d] = -(x0[d] - 1) + 1;
                  f *= parity[d];
                }
                assert(x0[d] <= +1);
              }
            }
          } // for d
          {
            const CCTK_REAL f0 = makevalue(x0);
            using std::abs;
            if (!(abs(f - f0) <=
                  10 * std::numeric_limits<CCTK_REAL>::epsilon()))
              CCTK_VINFO("parity=[%d,%d,%d] idx=[%d,%d,%d] x0=[%g,%g,%g] "
                         "x=[%g,%g,%g] f0=%g f=%g",
                         parity[0], parity[1], parity[2], i, j, k, x0[0], x0[1],
                         x0[2], x[0], x[1], x[2], f0, f);
            assert(abs(f - f0) <=
                   10 * std::numeric_limits<CCTK_REAL>::epsilon());
          }
        } // for i
      }   // for j
    }     // for k
    // Test that the function is not just zero
    for (int k = -5; k <= 5; ++k) {
      for (int j = -5; j <= 5; ++j) {
        for (int i = -5; i <= 5; ++i) {
          if (abs(i) >= 3 && abs(j) >= 3 && abs(k) >= 3) {
            const std::array<CCTK_REAL, 3> x{
                CCTK_REAL(0.25) * i, CCTK_REAL(0.25) * j, CCTK_REAL(0.25) * k};
            const CCTK_REAL f = makevalue(x);
            const bool want_zero =
                (parity[0] < 0 && reflection_lower[0] && i == -4) ||
                (parity[1] < 0 && reflection_lower[1] && j == -4) ||
                (parity[2] < 0 && reflection_lower[2] && k == -4) ||
                (parity[0] < 0 && reflection_upper[0] && i == +4) ||
                (parity[1] < 0 && reflection_upper[1] && j == +4) ||
                (parity[2] < 0 && reflection_upper[2] && k == +4);
            using std::abs;
            const bool is_zero =
                abs(f) <= 10 * std::numeric_limits<CCTK_REAL>::epsilon();
            if (!(want_zero == is_zero))
              CCTK_VINFO("parity=[%d,%d,%d] idx=[%d,%d,%d] x=[%g,%g,%g] f=%g",
                         parity[0], parity[1], parity[2], i, j, k, x[0], x[1],
                         x[2], f);
            assert(want_zero == is_zero);
          }
        }
      }
    }
  }

  constexpr std::array<int, Loop::dim> centering{CI, CJ, CK};
  const Loop::GF3D2layout layout(cctkGH, centering);
  const std::array<std::array<CCTK_REAL *, 8>, 2> gfptrs{{
      {{
          var_vvv_mmm,
          var_vvc_mmm,
          var_vcv_mmm,
          var_vcc_mmm,
          var_cvv_mmm,
          var_cvc_mmm,
          var_ccv_mmm,
          var_ccc_mmm,
      }},
      {{
          var_vvv_ppp,
          var_vvc_ppp,
          var_vcv_ppp,
          var_vcc_ppp,
          var_cvv_ppp,
          var_cvc_ppp,
          var_ccv_ppp,
          var_ccc_ppp,
      }},
  }};
  assert(parity[1] == parity[0] && parity[2] == parity[0]);
  const Loop::GF3D2<CCTK_REAL> var(
      layout, gfptrs[parity[0] > 0][(CI << 2) + (CJ << 1) + (CK << 0)]);

  const Loop::GridDescBaseDevice grid(cctkGH);
  grid.loop_device<CI, CJ, CK, where>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const auto value = makevalue(p.X);
        f(p, centering, parity, var(p.I), value);
      });
}

template <int CI, int CJ, int CK, Loop::where_t where, typename F>
void map_centering_all_parities(const cGH *restrict const cctkGH, const F &f) {
  for (int pk = -1; pk <= +1; pk += 2) {
    for (int pj = -1; pj <= +1; pj += 2) {
      for (int pi = -1; pi <= +1; pi += 2) {
        const std::array<int, 3> parity{pi, pj, pk};

        // Other parities are not (yet?) implemented and tested
        if (!((pi < 0 && pj < 0 && pk < 0) || (pi > 0 && pj > 0 && pk > 0)))
          continue;

        map_centering_parity<CI, CJ, CK, where>(cctkGH, parity, f);
      }
    }
  }
}

template <Loop::where_t where, typename F>
void map_all_centerings_all_parities(const cGH *restrict const cctkGH,
                                     const F &f) {
  map_centering_all_parities<0, 0, 0, where>(cctkGH, f);
  map_centering_all_parities<0, 0, 1, where>(cctkGH, f);
  map_centering_all_parities<0, 1, 0, where>(cctkGH, f);
  map_centering_all_parities<0, 1, 1, where>(cctkGH, f);
  map_centering_all_parities<1, 0, 0, where>(cctkGH, f);
  map_centering_all_parities<1, 0, 1, where>(cctkGH, f);
  map_centering_all_parities<1, 1, 0, where>(cctkGH, f);
  map_centering_all_parities<1, 1, 1, where>(cctkGH, f);
}

extern "C" void TestSymmetries_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_Init;
  DECLARE_CCTK_PARAMETERS;

  map_all_centerings_all_parities<Loop::where_t::interior>(
      cctkGH, [] CCTK_DEVICE(const auto &p, const auto &centering,
                             const auto &parity, auto &var, const auto &value)
                  CCTK_ATTRIBUTE_ALWAYS_INLINE { var = value; });
}

extern "C" void TestSymmetries_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  map_all_centerings_all_parities<Loop::where_t::boundary>(
      cctkGH, [] CCTK_DEVICE(const auto &p, const auto &centering,
                             const auto &parity, auto &var, const auto &value)
                  CCTK_ATTRIBUTE_ALWAYS_INLINE { var = value; });
}

extern "C" void TestSymmetries_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_Sync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

int check;

extern "C" void TestSymmetries_CheckInit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_CheckInit;
  DECLARE_CCTK_PARAMETERS;

  check = 0;
}

extern "C" void TestSymmetries_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_Check;
  DECLARE_CCTK_PARAMETERS;

  int local_check = 0;

  map_all_centerings_all_parities<Loop::where_t::everywhere>(
      cctkGH, [=, local_check_ptr = &local_check] CCTK_DEVICE(
                  const auto &p, const auto &centering, const auto &parity,
                  auto &var, const auto &value) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        using std::abs;
        if (abs(var - value) > 10 * std::numeric_limits<CCTK_REAL>::epsilon()) {
#if !defined __CUDACC__ && !defined __HIPCC__ && !defined SYCL_LANGUAGE_VERSION
          CCTK_VERROR(
              "Grid function symmetry check failed: I=[%d,%d,%d] X=[%g,%g,%g] "
              "centering=[%d,%d,%d] parity=[%d,%d,%d] var=%.17g value=%.17g",
              p.I[0], p.I[1], p.I[2], p.X[0], p.X[1], p.X[2], centering[0],
              centering[1], centering[2], parity[0], parity[1], parity[2], var,
              value);
#endif
          *local_check_ptr = 1;
        }
      });

  if (local_check)
    CCTK_VERROR("Grid function symmetry check failed");

#pragma omp atomic update
  check += local_check;
}

extern "C" void TestSymmetries_CheckFinalize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestSymmetries_CheckFinalize;
  DECLARE_CCTK_PARAMETERS;

  if (check != 0)
    CCTK_VERROR("Grid function symmetry check failed");
}

} // namespace TestSymmetries
