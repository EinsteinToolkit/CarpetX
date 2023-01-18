#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

#include <cmath>

namespace TestODESolvers {
using namespace std;

extern "C" void TestODESolvers_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_Initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL u0 = pow(1 + cctk_time, order);

  const int di = 1;
  const int dj = di * (cctk_ash[0] - 1);
  const int dk = dj * (cctk_ash[1] - 1);

  for (int k = 0; k < cctk_lsh[2] - 1; ++k) {
    for (int j = 0; j < cctk_lsh[1] - 1; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0] - 1; ++i) {
        int ind = i * di + j * dj + k * dk;
        state[ind] = u0;
        state2[ind] = u0;
      }
    }
  }
}

extern "C" void TestODESolvers_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_RHS;
  DECLARE_CCTK_PARAMETERS;

  // I want to step state2 with twice the stepsize, so need to adjust the time
  // at the RK substeps to have t_(k) = t0 + k dt * 2
  // TODO: this will break when using more than 1 thread
  static int last_it = -1;
  static CCTK_REAL original_time;
#pragma omp critical
  if (last_it != cctk_iteration) {
    last_it = cctk_iteration;
    original_time = cctk_time;
  }
  const CCTK_REAL cctk_time2 = original_time + 2 * (cctk_time - original_time);

  // u(t) = (1+t)^p
  // d/dt u = p (1+t)^(p-1)

  const int di = 1;
  const int dj = di * (cctk_ash[0] - 1);
  const int dk = dj * (cctk_ash[1] - 1);

  for (int k = 0; k < cctk_lsh[2] - 1; ++k) {
    for (int j = 0; j < cctk_lsh[1] - 1; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0] - 1; ++i) {
        int ind = i * di + j * dj + k * dk;
        // solving u(t) for (1 + t)^(order - 1) = u_inverse gives:
        const CCTK_REAL u_inverse =
            pow(state[ind], (order - 1) / CCTK_REAL(order));
        // mix in some of the state vector and some of the analytic RHS
        rhs[ind] = 0.5 * order * (u_inverse + pow(1 + cctk_time, order - 1));

        if (cctk_iteration % 2) {
          const CCTK_REAL u_inverse2 =
              pow(state2[ind], (order - 1) / CCTK_REAL(order));
          rhs2[ind] =
              2 * 0.5 * order * (u_inverse2 + pow(1 + cctk_time2, order - 1));
        } else {
          rhs2[ind] = 0;
        }
      }
    }
  }
}

extern "C" void TestODESolvers_PostStep(CCTK_ARGUMENTS) {
  // Do nothing
}

extern "C" void TestODESolvers_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_Error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL u0 = pow(1 + cctk_time, order);

  const int di = 1;
  const int dj = di * (cctk_ash[0] - 1);
  const int dk = dj * (cctk_ash[1] - 1);

  for (int k = 0; k < cctk_lsh[2] - 1; ++k) {
    for (int j = 0; j < cctk_lsh[1] - 1; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0] - 1; ++i) {
        int ind = i * di + j * dj + k * dk;
        error[ind] = state[ind] - u0;
        error2[ind] = state2[ind] - u0;

        if (error[ind] == 0 || error2[ind] == 0) {
          // happens reliably at cctk_time == cctk_initial_time
          corder[ind] = 0;
        } else {
          corder[ind] = log(abs(error2[ind] / error[ind])) / log(2.0);
        }
      }
    }
  }
}

} // namespace TestODESolvers
