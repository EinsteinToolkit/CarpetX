#include <loop.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>

namespace Punctures {
using namespace Arith;
using namespace Loop;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

template <typename T> T pown(T x, int n) {
  if (n < 0) {
    x = 1 / x;
    n = -n;
  }
  T r = 1;
  while (n > 0) {
    if (n % 2)
      r *= x;
    x *= x;
    n /= 2;
  }
  return r;
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL delta(int a, int b) {
  return a == b;
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_REAL epsilon(int a, int b, int c) {
  if (a == 0 && b == 1 && c == 2)
    return 1;
  if (a == 0 && b == 2 && c == 1)
    return -1;
  if (a == 1 && b == 0 && c == 2)
    return -1;
  if (a == 1 && b == 2 && c == 0)
    return 1;
  if (a == 2 && b == 0 && c == 1)
    return 1;
  if (a == 2 && b == 1 && c == 0)
    return -1;
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE void
fcalc(const PointDesc &p, const T &u, T &alpha, mat<T, 3> &K, T &rhs, T &psi) {
  DECLARE_CCTK_PARAMETERS;

  T alpha1 = 0;
  const mat<CCTK_REAL, 3> g([](int a, int b) { return a == b; });
  K = mat<T, 3>();
  for (int i = 0; i < npunctures; ++i) {
    const vec<T, 3> x{posx[i], posy[i], posz[i]};
    const vec<T, 3> P{momx[i], momy[i], momz[i]};
    const vec<T, 3> S{amomx[i], amomy[i], amomz[i]};
    const T r = sqrt(sum<3>([&](int a) { return pow(p.X[a] - x(a), 2); }));
    const vec<T, 3> n([&](int a) { return r < rmin ? 0 : x(a) / r; });

    alpha1 += mass[i] / (2 * fmax(rmin, r));

    K += mat<T, 3>([&](int a, int b) {
      return 3 / (2 * pow(fmax(rmin, r), 2)) *
                 (P(a) * n(b) + P(b) * n(a) -
                  (g(a, b) - n(a) * n(b)) *
                      sum<3>([&](int c) { return P(c) * n(c); })) +
             3 / pow(fmax(rmin, r), 3) * sum<3>([&](int c, int d) {
               return epsilon(a, c, d) * S(c) * n(d) * n(b) +
                      epsilon(b, c, d) * S(c) * n(d) * n(a);
             });
    });
  }
  alpha = 1 / alpha1;

  if (isinf(alpha)) {
    // Infinitely far away, or there are no black holes
    rhs = 0;
    psi = 1;
    return;
  }

  if (alpha == 0) {
    assert(0); // handled by rmin
    // At a puncture
    K = mat<T, 3>();
    rhs = 0;
    psi = INFINITY;
    return;
  }

  const T beta = 1 / T(8) * pow(alpha, 7) *
                 sum<3>([&](int a, int b) { return K(a, b) * K(a, b); });

  rhs = -beta * pow(1 + alpha * u, -7);
  psi = 1 / alpha + u;
}

template <typename T> T frhs(const PointDesc &p, const T &u) {
  T alpha;
  mat<T, 3> K;
  T rhs;
  T psi;
  fcalc(p, u, alpha, K, rhs, psi);
  return rhs;
}

// u = 1 + C / r + ...
template <typename T> T fbnd(const PointDesc &p) { return 1; }

template <typename T> T fpsi(const PointDesc &p, const T &u) {
  T alpha;
  mat<T, 3> K;
  T rhs;
  T psi;
  fcalc(p, u, alpha, K, rhs, psi);
  return psi;
}

template <typename T> mat<T, 3> fK(const PointDesc &p, const T &u) {
  T alpha;
  mat<T, 3> K;
  T rhs;
  T psi;
  fcalc(p, u, alpha, K, rhs, psi);
  return K;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void Punctures_init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Punctures_init;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> usol_(cctkGH, usol);

  // Initialize all points (including ghosts)
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { usol_(p.I) = 1.0; });

  // Note: The boundary conditions are applied on the outermost layer
  // of the interior points, not on the boundary points.

  // Set boundary conditions
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    if (cctk_lbnd[0] + p.i <= 1 || cctk_lbnd[0] + p.i >= cctk_gsh[0] - 1 ||
        cctk_lbnd[1] + p.j <= 1 || cctk_lbnd[1] + p.j >= cctk_gsh[1] - 1 ||
        cctk_lbnd[2] + p.k <= 1 || cctk_lbnd[2] + p.k >= cctk_gsh[2] - 1)
      usol_(p.I) = fbnd<CCTK_REAL>(p);
  });
}

extern "C" void Punctures_solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Punctures_solve;
  DECLARE_CCTK_PARAMETERS;

  const int max_iters = 100;
  for (int iter = 1;; ++iter) {
    CCTK_VINFO("Nonlinear iteration #%d", iter);

    // Set up RHS
    CallScheduleGroup(cctkGH, "Punctures_solve1");

    const int gi_sol = CCTK_GroupIndex("Punctures::usol");
    assert(gi_sol >= 0);
    const int gi_rhs = CCTK_GroupIndex("Punctures::urhs");
    assert(gi_rhs >= 0);
    const int gi_res = CCTK_GroupIndex("Punctures::ures");
    assert(gi_res >= 0);

    // linear solver accuracy
    const CCTK_REAL reltol = 0.0;
    const CCTK_REAL abstol = 1.0e-10;
    CCTK_REAL res_initial, res_final;
    SolvePoisson(gi_sol, gi_rhs, gi_res, reltol, abstol, &res_initial,
                 &res_final);
    CCTK_VINFO("Linear residual before solve: %g", double(res_initial));
    CCTK_VINFO("Linear residual after solve:  %g", double(res_final));

    // Correct boundaries
    CallScheduleGroup(cctkGH, "Punctures_solve2");

    if (res_initial <= abstol) {
      CCTK_VINFO("Nonlinear iterations finished after %d iterations; accuracy "
                 "goal reached",
                 iter);
      break;
    }
    if (iter == max_iters) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Nonlinear iterations finished after %d iterations; accuracy "
                 "goal NOT reached",
                 iter);
      break;
    }
  }
}

extern "C" void Punctures_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Punctures_rhs;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> usol_(cctkGH, usol);
  const GF3D<CCTK_REAL, 0, 0, 0> urhs_(cctkGH, urhs);

  // Set RHS
  loop_all<0, 0, 0>(
      cctkGH, [&](const PointDesc &p) { urhs_(p.I) = frhs(p, usol_(p.I)); });
}

extern "C" void Punctures_boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Punctures_boundary;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> usol_(cctkGH, usol);

  // Set boundary conditions
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    if (cctk_lbnd[0] + p.i <= 1 || cctk_lbnd[0] + p.i >= cctk_gsh[0] - 1 ||
        cctk_lbnd[1] + p.j <= 1 || cctk_lbnd[1] + p.j >= cctk_gsh[1] - 1 ||
        cctk_lbnd[2] + p.k <= 1 || cctk_lbnd[2] + p.k >= cctk_gsh[2] - 1)
      usol_(p.I) = fbnd<CCTK_REAL>(p);
  });
}

extern "C" void Punctures_ADMBase(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Punctures_ADMBase;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> usol_(cctkGH, usol);

  const GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

  const GF3D<CCTK_REAL, 0, 0, 0> Kxx_(cctkGH, kxx);
  const GF3D<CCTK_REAL, 0, 0, 0> Kxy_(cctkGH, kxy);
  const GF3D<CCTK_REAL, 0, 0, 0> Kxz_(cctkGH, kxz);
  const GF3D<CCTK_REAL, 0, 0, 0> Kyy_(cctkGH, kyy);
  const GF3D<CCTK_REAL, 0, 0, 0> Kyz_(cctkGH, kyz);
  const GF3D<CCTK_REAL, 0, 0, 0> Kzz_(cctkGH, kzz);

  const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  const GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
  const GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
  const GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

  const GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

  const GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    const CCTK_REAL psi = fpsi(p, usol_(p.I));
    const mat<CCTK_REAL, 3> g([](int a, int b) { return a == b; });
    const mat<CCTK_REAL, 3> K = fK(p, usol_(p.I));
    const mat<CCTK_REAL, 3> gph = pow(psi, 4) * g;
    const mat<CCTK_REAL, 3> Kph = pow(psi, -2) * K;
    gxx_(p.I) = gph(0, 0);
    gxy_(p.I) = gph(0, 1);
    gxz_(p.I) = gph(0, 2);
    gyy_(p.I) = gph(1, 1);
    gyz_(p.I) = gph(1, 2);
    gzz_(p.I) = gph(2, 2);
    Kxx_(p.I) = Kph(0, 0);
    Kxy_(p.I) = Kph(0, 1);
    Kxz_(p.I) = Kph(0, 2);
    Kyy_(p.I) = Kph(1, 1);
    Kyz_(p.I) = Kph(1, 2);
    Kzz_(p.I) = Kph(2, 2);
    alp_(p.I) = 1; // TODO
    betax_(p.I) = 0;
    betay_(p.I) = 0;
    betaz_(p.I) = 0;
    dtalp_(p.I) = 0;
    dtbetax_(p.I) = 0;
    dtbetay_(p.I) = 0;
    dtbetaz_(p.I) = 0;
  });
}

} // namespace Punctures
