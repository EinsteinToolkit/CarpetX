#include "discretization.hxx"
#include "physics.hxx"

#include <dual.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <ssht/ssht.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace AHFinder {
using namespace Arith;

template <typename T> using vec3 = vec<T, 3>;
template <typename T> using mat3 = mat<T, 3>;
template <typename T> using smat3 = smat<T, 3>;

template <typename F> auto sum3(F &&f) { return sum<3>(std::forward<F>(f)); }

template <typename T> struct coords_t {
  geom_t geom;
  vec3<scalar_aij_t<T> > x;

  coords_t() = delete;
  coords_t(const geom_t &geom)
      : geom(geom), x{
                        scalar_aij_t<T>(geom),
                        scalar_aij_t<T>(geom),
                        scalar_aij_t<T>(geom),
                    } {}
};

template <typename T>
coords_t<T> coords_from_shape(const vec3<T> &pos, const scalar_aij_t<T> &h) {
  DECLARE_CCTK_PARAMETERS;
  using std::cos, std::sin;
  const geom_t &geom = h.geom;
  coords_t<T> coords(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      const T r = h()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      coords.x(0)()(i, j) = pos(0) + r * sin(theta) * cos(phi);
      coords.x(1)()(i, j) = pos(1) + r * sin(theta) * sin(phi);
      coords.x(2)()(i, j) = pos(2) + r * cos(theta);
    }
  }
  return coords;
}

template <typename T> struct metric_t {
  geom_t geom;
  smat3<scalar_aij_t<T> > g, K;
  smat3<vec3<scalar_aij_t<T> > > dg;

  metric_t() = delete;
  metric_t(const geom_t &geom)
      : geom(geom),
        g{
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
        },
        K{
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
        },
        dg{
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
        } {}
};

template <typename T>
metric_t<T> brill_lindquist_metric(const cGH *const cctkGH,
                                   const coords_t<T> &coords) {
  DECLARE_CCTK_PARAMETERS;

  const geom_t &geom = coords.geom;
  metric_t<T> metric(geom);

  const T M = Brill_Lindquist_mass;
  const vec3<T> x0{Brill_Lindquist_x, Brill_Lindquist_y, Brill_Lindquist_z};
  const mat3<T> fx{
      Brill_Lindquist_fx,
      0.0,
      0.0, //
      0.0,
      Brill_Lindquist_fy,
      0.0, //
      0.0,
      0.0,
      Brill_Lindquist_fz, //
  };

  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";
      using U = vect<T, 3>;
      using DT = dual<T, U>; // spatial derivatives (x, y, z)
      const vec3<DT> x{
          {coords.x(0)()(i, j), {1, 0, 0}},
          {coords.x(1)()(i, j), {0, 1, 0}},
          {coords.x(2)()(i, j), {0, 0, 1}},
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const vec3<DT> X = matmul(fx, x - x0);
      if (0)
        std::cout << "  X:" << X << "\n";
      using std::pow, std::sqrt;
      const DT R = sqrt(pow(X(0), 2) + pow(X(1), 2) + pow(X(2), 2));
      const DT Psi = pow(1 + M / (2 * R), 4);
      const smat3<T> eta([&](int a, int b) { return pow(X(a).eps[b], 2); });
      if (0)
        std::cout << "  eta:" << eta << "\n";
      const smat3<DT> gg([&](int a, int b) { return Psi * eta(a, b); });
      const smat3<T> g([&](int a, int b) { return gg(a, b).val; });
      const smat3<vec3<T> > dg([&](int a, int b) {
        return vec3<T>([&](int c) { return gg(a, b).eps[c]; });
      });
      if (0)
        std::cout << "  g:" << g << "\n";
      if (0)
        std::cout << "  dg:" << dg << "\n";
      const smat3<T> K([&](int a, int b) { return 0; });
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          metric.g(a, b)()(i, j) = g(a, b);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          for (int c = 0; c < 3; ++c)
            metric.dg(a, b)(c)()(i, j) = dg(a, b)(c);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          metric.K(a, b)()(i, j) = K(a, b);
    }
  }

  return metric;
}

template <typename T>
metric_t<T> interpolate_metric(const cGH *const cctkGH,
                               const coords_t<T> &coords) {
  const int gxx_ind = CCTK_VarIndex("ADMBase::gxx");
  const int gxy_ind = CCTK_VarIndex("ADMBase::gxy");
  const int gxz_ind = CCTK_VarIndex("ADMBase::gxz");
  const int gyy_ind = CCTK_VarIndex("ADMBase::gyy");
  const int gyz_ind = CCTK_VarIndex("ADMBase::gyz");
  const int gzz_ind = CCTK_VarIndex("ADMBase::gzz");
  const int kxx_ind = CCTK_VarIndex("ADMBase::kxx");
  const int kxy_ind = CCTK_VarIndex("ADMBase::kxy");
  const int kxz_ind = CCTK_VarIndex("ADMBase::kxz");
  const int kyy_ind = CCTK_VarIndex("ADMBase::kyy");
  const int kyz_ind = CCTK_VarIndex("ADMBase::kyz");
  const int kzz_ind = CCTK_VarIndex("ADMBase::kzz");

  constexpr int nvars = 6 * (1 + 3 + 1);
  const array<CCTK_INT, nvars> varinds{
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      kxx_ind, kxy_ind, kxz_ind, kyy_ind, kyz_ind, kzz_ind, //
  };
  const array<CCTK_INT, nvars> operations{
      0, 0, 0, 0, 0, 0, //
      1, 1, 1, 1, 1, 1, //
      2, 2, 2, 2, 2, 2, //
      3, 3, 3, 3, 3, 3, //
      0, 0, 0, 0, 0, 0, //
  };

  const geom_t &geom = coords.geom;
  metric_t<T> metric(geom);
  array<T *, nvars> ptrs{
      metric.g(0, 0)().data(),     metric.g(0, 1)().data(),
      metric.g(0, 2)().data(),     metric.g(1, 1)().data(),
      metric.g(1, 2)().data(),     metric.g(2, 2)().data(),
      metric.dg(0, 0)(0)().data(), metric.dg(0, 1)(0)().data(),
      metric.dg(0, 2)(0)().data(), metric.dg(1, 1)(0)().data(),
      metric.dg(1, 2)(0)().data(), metric.dg(2, 2)(0)().data(),
      metric.dg(0, 0)(1)().data(), metric.dg(0, 1)(1)().data(),
      metric.dg(0, 2)(1)().data(), metric.dg(1, 1)(1)().data(),
      metric.dg(1, 2)(1)().data(), metric.dg(2, 2)(1)().data(),
      metric.dg(0, 0)(2)().data(), metric.dg(0, 1)(2)().data(),
      metric.dg(0, 2)(2)().data(), metric.dg(1, 1)(2)().data(),
      metric.dg(1, 2)(2)().data(), metric.dg(2, 2)(2)().data(),
      metric.K(0, 0)().data(),     metric.K(0, 1)().data(),
      metric.K(0, 2)().data(),     metric.K(1, 1)().data(),
      metric.K(1, 2)().data(),     metric.K(2, 2)().data()};

  Interpolate(cctkGH, geom.npoints, coords.x(0)().data(), coords.x(1)().data(),
              coords.x(2)().data(), nvars, varinds.data(), operations.data(),
              ptrs.data());

  metric.g(1, 0)() = metric.g(0, 1)();
  metric.g(2, 0)() = metric.g(0, 2)();
  metric.g(2, 1)() = metric.g(1, 2)();
  metric.dg(1, 0)() = metric.dg(0, 1)();
  metric.dg(2, 0)() = metric.dg(0, 2)();
  metric.dg(2, 1)() = metric.dg(1, 2)();
  metric.K(1, 0)() = metric.K(0, 1)();
  metric.K(2, 0)() = metric.K(0, 2)();
  metric.K(2, 1)() = metric.K(1, 2)();

  return metric;
}

////////////////////////////////////////////////////////////////////////////////

enum class which_Theta_t { Theta_l, Theta_n };

template <typename T> struct Theta_t {
  which_Theta_t which_Theta;
  scalar_alm_t<std::complex<T> > Thetalm;
  scalar_alm_t<std::complex<T> > rhoThetalm;
};

template <typename T>
Theta_t<T> expansion(const metric_t<T> &metric, const vec3<T> &pos,
                     const scalar_alm_t<std::complex<T> > &hlm,
                     const which_Theta_t which_Theta = which_Theta_t::Theta_l) {
  DECLARE_CCTK_PARAMETERS;

  const geom_t &geom = hlm.geom;

  // Cartesian coordinates:
  //     x = r \cos\theta \cos\phi
  //     y = r \cos\theta \sin\phi
  //     z = r \sin\theta

  // Surface coordinates:
  //     F(r,\theta,\phi) = r - h(\theta,\phi)
  //     F   = r - h(\theta,\phi)
  //     \mu = \theta
  //     \nu = \phi

  // Derivatives Cartesian/spherical:
  //     dx/dr      =     \cos\theta \cos\phi
  //     dy/dr      =     \cos\theta \sin\phi
  //     dz/dr      =     \sin\theta
  //     dx/d\theta = - r \sin\theta \cos\phi
  //     dy/d\theta = - r \sin\theta \sin\phi
  //     dz/d\theta =   r \cos\theta
  //     dx/d\phi   = - r \cos\theta \sin\phi
  //     dy/d\phi   =   r \cos\theta \cos\phi
  //     dz/d\phi   =   0

  // Derivatives surface/spherical:
  //     dF  /dr      = 1
  //     d\mu/dr      = 0
  //     d\nu/dr      = 0
  //     df  /d\theta = - dh/d\theta(\theta,\phi)
  //     d\mu/d\theta = 1
  //     d\nu/d\theta = 0
  //     df  /d\phi   = - dh/d\phi(\theta,\phi)
  //     d\mu/d\phi   = 0
  //     d\nu/d\phi   = 1

  // Evaluate h and its derivatives
  const scalar_aij_t<T> hij = evaluate(hlm);
  const vector_alm_t<std::complex<T> > dhlm = gradient(hlm);
  const vector_aij_t<T> dhij = evaluate(dhlm);
  if (0)
    std::cout << "dhij:" << dhij << "\n";

  // s^a
  const vec3<scalar_aij_t<T> > sij{
      scalar_aij_t<T>(geom),
      scalar_aij_t<T>(geom),
      scalar_aij_t<T>(geom),
  };
  scalar_aij_t<T> rhoij(geom); // (28)
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";

      // Coordinates
      const T r = hij()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      if (0)
        std::cout << "  r:" << r << " theta:" << theta << " phi:" << phi
                  << "\n";
      using std::cos, std::sin;
      const vec3<T> x{
          pos(0) + r * sin(theta) * cos(phi),
          pos(1) + r * sin(theta) * sin(phi),
          pos(2) + r * cos(theta),
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const mat3<T> dxdr{
          sin(theta) * cos(phi),
          r * cos(theta) * cos(phi),
          -r * sin(phi), // dx/dphi / sin(theta)

          sin(theta) * sin(phi),
          r * cos(theta) * sin(phi),
          r * cos(phi), // dy/dphi / sin(theta)

          cos(theta),
          -r * sin(theta),
          0, // dz/dphi / sin(theta)
      };
      if (0)
        std::cout << "  dxdr:" << dxdr << "\n";
      const mat3<T> drdx = inv(dxdr);
      if (0)
        std::cout << "  drdx:" << drdx << "\n";
      const vec3<T> grad_r([&](int a) { return drdx(0, a); });

      // Level set function
      // F(r, theta, phi) = r - h(theta, phi)
      // const T F = 0
      const mat3<T> dFdr{
          1,              // dF/dr
          -dhij(0)(i, j), // dF/dtheta
          -dhij(1)(i, j), // dF/dphi / sin(theta)

          0,
          1,
          0,

          0,
          0,
          1,
      };
      if (0)
        std::cout << "  dFdr:" << dFdr << "\n";
      const mat3<T> dFdx([&](int a, int b) {
        return sum3([&](int x) { return dFdr(a, x) * drdx(x, b); });
      });

      // Metric
      const smat3<T> g([&](int a, int b) { return metric.g(a, b)()(i, j); });
      const smat3<T> gu = inv(g);

      // Surface normal
      const vec3<T> grad_F_F{1, 0, 0};

      const vec3<T> grad_F([&](int a) {
        return sum3([&](int x) { return grad_F_F(x) * dFdx(x, a); });
      });
      const vec3<T> grad_F_u([&](int a) {
        return sum3([&](int x) { return gu(a, x) * grad_F(x); });
      });

      using std::sqrt;
      const T len_grad_F =
          sqrt(sum3([&](int x) { return grad_F_u(x) * grad_F(x); }));
      const vec3<T> s = grad_F / len_grad_F;
      const vec3<T> su(
          [&](int a) { return sum3([&](int x) { return gu(a, x) * s(x); }); });
      for (int a = 0; a < 3; ++a)
        (T &)sij(a)()(i, j) = s(a);

      // Auxiliary term rho
      const T rho = 2 * pow2(r) * len_grad_F / sum3([&](int a, int b) {
                      return (gu(a, b) - su(a) * su(b)) *
                             ((a == b) - grad_r(a) * grad_r(b));
                    });
      if (0)
        std::cout << "  rho:" << rho << "\n";
      rhoij()(i, j) = rho;
    }
  }

  const vec3<scalar_alm_t<std::complex<T> > > slm(
      [&](int a) { return expand(sij(a)); });
  const vec3<vector_alm_t<std::complex<T> > > dslm(
      [&](int a) { return gradient(slm(a)); });
  const vec3<vector_aij_t<T> > dsij([&](int a) { return evaluate(dslm(a)); });
  if (0)
    std::cout << "dsij:" << dsij << "\n";

  scalar_aij_t<T> Thetaij(geom); // (9), or (22)
  scalar_aij_t<T> rhoThetaij(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";

      // Coordinates
      const T r = hij()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      if (0)
        std::cout << "  r:" << r << " theta:" << theta << " phi:" << phi
                  << "\n";
      using std::cos, std::sin;
      const vec3<T> x{
          pos(0) + r * sin(theta) * cos(phi),
          pos(1) + r * sin(theta) * sin(phi),
          pos(2) + r * cos(theta),
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const mat3<T> dxdr{
          sin(theta) * cos(phi),
          r * cos(theta) * cos(phi),
          -r * sin(phi), // dx/dphi / sin(theta)

          sin(theta) * sin(phi),
          r * cos(theta) * sin(phi),
          r * cos(phi), // dz/dphi / sin(theta)

          cos(theta),
          -r * sin(theta),
          0, // dz/dphi / sin(theta)
      };
      if (0)
        std::cout << "  dxdr:" << dxdr << "\n";
      const mat3<T> drdx = inv(dxdr);
      if (0)
        std::cout << "  drdx:" << drdx << "\n";

      // Level set function
      // F(r, theta, phi) = r - h(theta, phi)
      // const T F = 0
      const mat3<T> dFdr{
          1,              // dF/dr
          -dhij(0)(i, j), // dF/dtheta
          -dhij(1)(i, j), // dF/dphi / sin(theta)

          0,
          1,
          0,

          0,
          0,
          1,
      };
      if (0)
        std::cout << "  dFdr:" << dFdr << "\n";
      const mat3<T> dFdx([&](int a, int b) {
        return sum3([&](int x) { return dFdr(a, x) * drdx(x, b); });
      });

      // Metric
      const smat3<T> g([&](int a, int b) { return metric.g(a, b)()(i, j); });
      const smat3<T> K([&](int a, int b) { return metric.K(a, b)()(i, j); });
      const smat3<T> gu = inv(g);
      const smat3<vec3<T> > dg([&](int a, int b) {
        return vec3<T>([&](int c) { return metric.dg(a, b)(c)()(i, j); });
      });
      const vec3<smat3<T> > Gamma([&](int a) {
        return smat3<T>([&](int b, int c) {
          return sum3([&](int x) {
            return gu(a, x) * (dg(x, c)(b) + dg(b, x)(c) - dg(b, c)(x)) / 2;
          });
        });
      });

      // Surface normal
      const vec3<T> s([&](int a) { return sij(a)()(i, j); });
      if (0)
        std::cout << "  s:" << s << "\n";
      const mat3<T> dsdF{
          0, dsij(0)(0)(i, j), dsij(0)(1)(i, j),

          0, dsij(1)(0)(i, j), dsij(1)(1)(i, j),

          0, dsij(2)(0)(i, j), dsij(2)(1)(i, j),
      };
      if (0)
        std::cout << "  dsdF:" << dsdF << "\n";
      const mat3<T> dsdx([&](int a, int b) {
        return sum3([&](int x) { return dsdF(a, x) * dFdx(x, b); });
      });
      if (0)
        std::cout << "  dsdx:" << dsdx << "\n";
      const mat3<T> grad_s([&](int a, int b) {
        return dsdx(a, b) + sum3([&](int x) { return Gamma(a)(b, x) * s(x); });
      });
      const vec3<T> su(
          [&](int a) { return sum3([&](int x) { return gu(a, x) * s(x); }); });

      // Expansion Theta_(l) or Theta_(n)
      const int sign = which_Theta == which_Theta_t::Theta_l ? -1 : +1;
      const T Theta = sum3([&](int x, int y) {
        return (gu(x, y) - su(x) * su(y)) * (grad_s(x, y) + sign * K(x, y));
      });
      if (0)
        std::cout << "  Theta:" << Theta << "\n";
      Thetaij()(i, j) = Theta;
      const T rho = rhoij()(i, j);
      rhoThetaij()(i, j) = rho * Theta;
    }
  }

  scalar_alm_t<std::complex<T> > Thetalm = expand(Thetaij);
  scalar_alm_t<std::complex<T> > rhoThetalm = expand(rhoThetaij);

  return Theta_t<T>{which_Theta, std::move(Thetalm), std::move(rhoThetalm)};
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
void update_position(vec3<T> &pos, scalar_alm_t<std::complex<T> > &hlm) {
  using std::sqrt;

  // x = \cos\phi = (\exp i \phi + \exp -i \phi) / 2
  // y = \sin\phi = (\exp i \phi - \exp -i \phi) / 2
  // z = \cos\theta
  const vec3<T> delta{
      real(hlm()(1, -1) - hlm()(1, +1)) / 2 / sqrt(2 * T(M_PI) / 3),
      imag(hlm()(1, -1) + hlm()(1, +1)) / 2 / sqrt(2 * T(M_PI) / 3),
      real(hlm()(1, 0)) / sqrt(4 * T(M_PI) / 3),
  };

  pos += delta;

  for (int m = -1; m <= +1; ++m)
    hlm()(1, m) = 0;
}

template <typename T>
scalar_alm_t<std::complex<T> >
step(const cGH *const cctkGH, const vec3<T> &pos, const T &radius,
     const scalar_alm_t<std::complex<T> > &hlm, const Theta_t<T> &Theta) {
  DECLARE_CCTK_PARAMETERS;

  const geom_t &geom = hlm.geom;
  const auto &rhoThetalm = Theta.rhoThetalm;

  // const T alpha = 1.0;
  // const T beta = 0.5;
  // const T A = alpha / (geom.lmax * (geom.lmax + 1)) + beta;
  // const T B = beta / alpha;

  const T A = fast_flow_A;
  const T B = fast_flow_B;

  scalar_alm_t<std::complex<T> > delta_hlm(geom);
  for (int l = 0; l <= geom.lmax; ++l) {
#pragma omp simd
    for (int m = -l; m <= l; ++m) {
      // const T lambda = A / (1 + B * l * (l + 1));
      const T lambda = l == 0 ? A : A / (B + T(l * (l + 1)));
      delta_hlm()(l, m) = -lambda * rhoThetalm()(l, m);
      // const auto L = geom.lmax;
      // const auto ll1 = l * (l + 1);
      // const auto LL1 = L * (L + 1);
      // const auto Q = (alpha / LL1 + beta) / (1 + beta / alpha * ll1);
      // const auto Q = alpha * (alpha / LL1 + beta) / (alpha + beta * ll1);
      // const auto Q = (alpha / (LL1 * LL1) + beta / LL1) /
      //                (1 / LL1 + beta / alpha * ll1 / LL1);
      // delta_hlm()(l, m) = - Q * rhoThetalm()(l, m);
    }
  }

  // Limit step size to 10% of the current radius
  const T h00 = real(hlm()(0, 0));
  delta_hlm()(0, 0) = clamp(real(delta_hlm()(0, 0)), -0.1 * h00, 0.1 * h00);

  return delta_hlm;
}

template <typename T>
void solve(const cGH *const cctkGH, vec3<T> &pos, T &radius,
           scalar_alm_t<std::complex<T> > hlm) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int iter = 0;
  for (;;) {
    if (iter >= max_iters)
      break;

    ++iter;
    CCTK_VINFO("iter: %d", iter);

    update_position(pos, hlm);
    radius = average(hlm);

    const auto hij = evaluate(hlm);
    CCTK_VINFO("    pos=[%g,%g,%g]", pos(0), pos(1), pos(2));
    CCTK_VINFO("    r_avg=%g   r_min=%g r_max=%g", average(hlm), minimum(hij),
               maximum(hij));
    if (0) {
      const int lmax = hlm.geom.lmax;
      for (int l = 0; l <= min(4, lmax); ++l) {
        using std::abs, std::max, std::min;
        T r = 0.0, rmin = 1.0 / 0.0, rmax = -1.0 / 0.0;
        for (int m = -l; m <= +l; ++m) {
          rmin = min(rmin, real(hlm()(l, m)));
          rmin = min(rmin, imag(hlm()(l, m)));
          rmax = max(rmax, real(hlm()(l, m)));
          rmax = max(rmax, imag(hlm()(l, m)));
          r = max(r, abs(hlm()(l, m)));
        }
        CCTK_VINFO("    |h%dm|=%g   %g   %g", l, r, rmin, rmax);
      }
    }

    const auto coords = coords_from_shape(pos, hij);
    const auto metric = use_Brill_Lindquist_metric
                            ? brill_lindquist_metric(cctkGH, coords)
                            : interpolate_metric(cctkGH, coords);
    const auto Theta = expansion(metric, pos, hlm);
    const auto &Thetalm = Theta.Thetalm;
    const auto Thetaij = evaluate(Thetalm);
    CCTK_VINFO("    Θ_avg=%g   Θ_maxabs=%g", average(Thetalm),
               maxabs(Thetaij()));

    auto delta_hlm = step(cctkGH, pos, radius, hlm, Theta);
    CCTK_VINFO("    Δr_avg=%g", average(delta_hlm));

    if (0) {
      using std::abs, std::max, std::min;
      const int lmax = hlm.geom.lmax;
      for (int l = 0; l <= min(4, lmax); ++l) {
        T r = 0.0, rmin = 1.0 / 0.0, rmax = -1.0 / 0.0;
        for (int m = -l; m <= +l; ++m) {
          rmin = min(rmin, real(delta_hlm()(l, m)));
          rmin = min(rmin, imag(delta_hlm()(l, m)));
          rmax = max(rmax, real(delta_hlm()(l, m)));
          rmax = max(rmax, imag(delta_hlm()(l, m)));
          r = max(r, abs(delta_hlm()(l, m)));
        }
        CCTK_VINFO("    |Δh%dm|=%g   %g   %g", l, r, rmin, rmax);
      }
    }

    hlm = hlm + delta_hlm;
    if (maxabs(Thetaij()) <= max_expansion)
      break;
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void AHFinder_init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_init;
  DECLARE_CCTK_PARAMETERS;

  *ah_pos_x = initial_pos_x;
  *ah_pos_y = initial_pos_y;
  *ah_pos_z = initial_pos_z;

  *ah_radius = initial_radius;
}

extern "C" void AHFinder_find(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_find;
  DECLARE_CCTK_PARAMETERS;

  vec3<CCTK_REAL> pos{*ah_pos_x, *ah_pos_y, *ah_pos_z};
  CCTK_REAL radius{*ah_radius};

  const geom_t geom(npoints);
  scalar_alm_t<CCTK_COMPLEX> hlm =
      scalar_from_const(geom, CCTK_COMPLEX(radius));
  solve(cctkGH, pos, radius, hlm);

  *ah_pos_x = pos(0);
  *ah_pos_y = pos(1);
  *ah_pos_z = pos(2);
  *ah_radius = radius;
}

} // namespace AHFinder
