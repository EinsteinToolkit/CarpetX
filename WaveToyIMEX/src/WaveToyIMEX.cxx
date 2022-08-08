#include <fixmath.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

extern CCTK_REAL ODESolvers_alpha;

namespace WaveToyIMEX {
using namespace std;
using namespace Loop;

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Spline with compact support of radius 1 and volume 1
template <typename T> T spline(T r) {
  if (r >= 1.0)
    return 0.0;
  constexpr CCTK_REAL f = dim == 1   ? 1.0
                          : dim == 2 ? 24.0 / 7.0 / M_PI
                          : dim == 3 ? 4.0 / M_PI
                                     : -1;
  const T r2 = pow(r, 2);
  return f * (r <= 0.5 ? 1 - 2 * r2 : 2 + r * (-4 + 2 * r));
}

// The potential for the spline
template <typename T> T spline_potential(T r) {
  // \Laplace u = 4 \pi \rho
  if (r >= 1.0)
    return -1 / r;
  static_assert(dim == 3, "");
  const T r2 = pow(r, 2);
  return r <= 0.5
             ? -7 / T(3) + r2 * (8 / T(3) - 8 / T(5) * r2)
             : (1 / T(15) +
                r * (-40 / T(15) +
                     r2 * (80 / T(15) + r * (-80 / T(15) + r * 24 / T(15))))) /
                   r;
}

// Time derivative
// TODO: Use dual numbers for derivative
template <typename F, typename T> auto timederiv(const F &f, T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Gradient
template <typename T> auto xderiv(T f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx, y, z) - f(t, x - dx, y, z)) / (2 * dx);
  };
}
template <typename T> auto yderiv(T f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy, z) - f(t, x, y - dy, z)) / (2 * dy);
  };
}
template <typename T> auto zderiv(T f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz) - f(t, x, y, z - dz)) / (2 * dz);
  };
}

////////////////////////////////////////////////////////////////////////////////

// Standing wave
CCTK_REAL standing(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = 2 * M_PI * spatial_frequency_x;
  auto ky = 2 * M_PI * spatial_frequency_y;
  auto kz = 2 * M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

// Periodic Gaussian
CCTK_REAL periodic_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = M_PI * spatial_frequency_x;
  auto ky = M_PI * spatial_frequency_y;
  auto kz = M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return exp(-0.5 * pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
}

// Gaussian
CCTK_REAL gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  // u(t,r) = (f(r-t) - f(r+t)) / r
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  auto f = [&](auto x) { return exp(-0.5 * pow(x / width, 2)); };
  auto fx = [&](auto x) { return -x / pow(width, 2) * f(x); };
  if (r < 1.0e-8)
    // Use L'Hôpital's rule for small r
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}


////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyIMEX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyIMEX_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const array<int, dim> indextype = {0, 0, 0}; // using vertex coordinates
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<CCTK_REAL> gf_zeta(layout, zeta);
  const GF3D2<CCTK_REAL> gf_mu(layout, mu);
  const GF3D2<CCTK_REAL> gf_nu(layout, nu);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = standing(t, p.x, p.y, p.z);
      gf_mu(p.I) = 0;
      gf_zeta(p.I) = standing(t, p.x, p.y, p.z);
      gf_nu(p.I) = 0;
    });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = periodic_gaussian(t, p.x, p.y, p.z);
      gf_mu(p.I) = 0;
      gf_zeta(p.I) = periodic_gaussian(t, p.x, p.y, p.z);
      gf_nu(p.I) = 0;
    });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = gaussian(t, p.x, p.y, p.z);
      gf_mu(p.I) = 0;
      gf_zeta(p.I) = gaussian(t, p.x, p.y, p.z);
      gf_nu(p.I) = 0;
    });

  } else {
    assert(0);
  }
}

extern "C" void WaveToyIMEX_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyIMEX_Sync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

/////////////////////////////////////////////////////////////////////////////
extern "C" void WaveToyIMEX_NonStiffRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyIMEX_NonStiffRHS;
  DECLARE_CCTK_PARAMETERS;


  const array<int, dim> indextype = {0, 0, 0}; // use vertex coordinates
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout,phi);
  const GF3D2<const CCTK_REAL> gf_zeta(layout,zeta);
  const GF3D2<const CCTK_REAL> gf_mu(layout,mu);
  const GF3D2<const CCTK_REAL> gf_nu(layout,nu);

  const GF3D2<CCTK_REAL> gf_phi_NonStiffRHS(layout,phi_NonStiffRHS);
  const GF3D2<CCTK_REAL> gf_zeta_NonStiffRHS(layout,zeta_NonStiffRHS);
  const GF3D2<CCTK_REAL> gf_mu_NonStiffRHS(layout,mu_NonStiffRHS);
  const GF3D2<CCTK_REAL> gf_nu_NonStiffRHS(layout,nu_NonStiffRHS);


  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {

    CCTK_REAL ddx_phi = ( -gf_phi(p.I + 2*p.DI[0]) + 16*gf_phi(p.I + p.DI[0]) - 30*gf_phi(p.I) + 
                        16*gf_phi(p.I - p.DI[0]) - gf_phi(p.I - 2*p.DI[0]) ) / ( 12*pow(p.dx, 2) );

    CCTK_REAL ddy_phi = ( -gf_phi(p.I + 2*p.DI[1]) + 16*gf_phi(p.I + p.DI[1]) - 30*gf_phi(p.I) + 
                        16*gf_phi(p.I - p.DI[1]) - gf_phi(p.I - 2*p.DI[1]) ) / ( 12*pow(p.dx, 2) );

    CCTK_REAL ddz_phi = ( -gf_phi(p.I + 2*p.DI[2]) + 16*gf_phi(p.I + p.DI[2]) - 30*gf_phi(p.I) +
                        16*gf_phi(p.I - p.DI[2]) - gf_phi(p.I - 2*p.DI[2]) ) / ( 12*pow(p.dx, 2) );

    gf_mu_NonStiffRHS(p.I) = pow(wave_speed_phi,2) * (ddx_phi + ddy_phi + ddz_phi);
    gf_phi_NonStiffRHS(p.I) = gf_mu(p.I);


    CCTK_REAL ddx_zeta = ( -gf_zeta(p.I + 2*p.DI[0]) + 16*gf_zeta(p.I + p.DI[0]) - 30*gf_zeta(p.I) + 
			16*gf_zeta(p.I - p.DI[0]) - gf_zeta(p.I - 2*p.DI[0]) ) / ( 12*pow(p.dx, 2) );

    CCTK_REAL ddy_zeta = ( -gf_zeta(p.I + 2*p.DI[1]) + 16*gf_zeta(p.I + p.DI[1]) - 30*gf_zeta(p.I) + 
			16*gf_zeta(p.I - p.DI[1]) - gf_zeta(p.I - 2*p.DI[1]) ) / ( 12*pow(p.dx, 2) );

    CCTK_REAL ddz_zeta = ( -gf_zeta(p.I + 2*p.DI[2]) + 16*gf_zeta(p.I + p.DI[2]) - 30*gf_zeta(p.I) +
			 16*gf_zeta(p.I - p.DI[2]) - gf_zeta(p.I - 2*p.DI[2]) ) / ( 12*pow(p.dx, 2) );

    gf_nu_NonStiffRHS(p.I) = pow(wave_speed_zeta,2) * (ddx_zeta + ddy_zeta + ddz_zeta);
    gf_zeta_NonStiffRHS(p.I) = gf_nu(p.I);
  });
}

extern "C" void WaveToyIMEX_NonStiffRHSSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyIMEX_NonStiffRHSSync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

extern "C" void WaveToyIMEX_UserSolvedFunction_G(CCTK_ARGUMENTS) {
  // beta = y0 + alpha*k1_hat
  // y1 = alpha*g(y1) - beta
  // y1 -> output of 'UserSolvedFunction_G(alpha,beta)'
  // but this routine 'UserSolvedFunction_G()' receives beta as state vector
  // Hence, here state vector = beta
  DECLARE_CCTK_ARGUMENTS_WaveToyIMEX_UserSolvedFunction_G;
  DECLARE_CCTK_PARAMETERS;


  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phi(layout,phi);
  const GF3D2<CCTK_REAL> gf_zeta(layout,zeta);
  const GF3D2<CCTK_REAL> gf_mu(layout,mu);
  const GF3D2<CCTK_REAL> gf_nu(layout,nu);
  
  const GF3D2<CCTK_REAL> gf_phi_1(layout,phi);
  const GF3D2<CCTK_REAL> gf_zeta_1(layout,zeta);
  const GF3D2<CCTK_REAL> gf_mu_1(layout,mu);
  const GF3D2<CCTK_REAL> gf_nu_1(layout,nu);


  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {

    gf_phi_1(p.I) = gf_phi(p.I);
    gf_mu_1(p.I) = gf_mu(p.I);
    gf_zeta_1(p.I) = (gf_zeta(p.I) + ODESolvers_alpha*coupling_factor*gf_phi_1(p.I))/(1 + ODESolvers_alpha*coupling_factor);
    gf_nu_1(p.I) = gf_nu(p.I);

  });
  
  // Now store calculated 'y1' in original state vector
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    gf_phi(p.I) = gf_phi_1(p.I);
    gf_mu(p.I) = gf_mu_1(p.I);
    gf_zeta(p.I) = gf_zeta_1(p.I);
    gf_nu(p.I) = gf_nu_1(p.I);

  });
}

//////////////////////////////////////////////////////////////////////////////


} // namespace WaveToyIMEX