#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "utils.hxx"
#include "reconstruct.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

enum class flux_t { LxF, HLLE };

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  flux_t fluxtype;
  if (CCTK_EQUALS(flux_type, "LxF")) {
    fluxtype = flux_t::LxF;
  } else if (CCTK_EQUALS(flux_type, "HLLE")) {
    fluxtype = flux_t::HLLE;
  } else {
    CCTK_ERROR("Unknown value for parameter \"flux_type\"");
  }

  switch (reconstruction) {
  case reconstruction_t::Godunov:
    assert(cctk_nghostzones[dir] >= 1);
    break;
  case reconstruction_t::minmod:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::monocentral:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::ppm:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  }

  const auto reconstruct_pt =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return reconstruct(var, p, reconstruction, dir);
          };

  const auto eigenvalues =
      [=] CCTK_DEVICE(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                      vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                      vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                      vec<CCTK_REAL, 2> h,
                      vec<CCTK_REAL, 2> bsq) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // computing characteristics for the minus side
        // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0
        vec<CCTK_REAL, 3> a_m{
            (bsq(0) + cs2(0) * h(0) * rho(0)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(0)) * h(0) * rho(0) *
                    pow2(beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            2 * beta_avg * (bsq(0) + cs2(0) * h(0) * rho(0)) -
                2 * (-1 + cs2(0)) * h(0) * rho(0) *
                    (beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            bsq(0) + h(0) * rho(0) *
                         (cs2(0) + pow2(w_lor(0)) - cs2(0) * pow2(w_lor(0)))};

        CCTK_REAL det_m = pow2(a_m(1)) - 4.0 * a_m(2) * a_m(0);
        if (det_m < 0.0)
          det_m = 0.0;

        vec<CCTK_REAL, 4> lambda_m{
            ((-a_m(1) + sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) + sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) - sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg,
            ((-a_m(1) - sqrt(det_m)) / (2.0 * a_m(2))) / alp_avg};

        // computing characteristics for the plus side

        vec<CCTK_REAL, 3> a_p{
            (bsq(1) + cs2(1) * h(1) * rho(1)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(1)) * h(1) * rho(1) *
                    pow2(beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            2 * beta_avg * (bsq(1) + cs2(1) * h(1) * rho(1)) -
                2 * (-1 + cs2(1)) * h(1) * rho(1) *
                    (beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            bsq(1) + h(1) * rho(1) *
                         (cs2(1) + pow2(w_lor(1)) - cs2(1) * pow2(w_lor(1)))};

        CCTK_REAL det_p = pow2(a_p(1)) - 4.0 * a_p(2) * a_p(0);
        if (det_p < 0.0)
          det_p = 0.0;

        vec<CCTK_REAL, 4> lambda_p{
            ((-a_p(1) + sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) + sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) - sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg,
            ((-a_p(1) - sqrt(det_p)) / (2.0 * a_p(2))) / alp_avg};

        // 2D array containing characteristics for left (minus) and right (plus)
        // sides
        vec<vec<CCTK_REAL, 4>, 2> lambda{lambda_m, lambda_p};
        return lambda;
      };

  const auto calcflux =
      [=] CCTK_DEVICE(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
                      vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL flx;
        switch (fluxtype) {
        case flux_t::LxF: {
          const CCTK_REAL charmax =
              max({0.0, fabs(lam(0)(0)), fabs(lam(0)(1)), fabs(lam(0)(2)),
                   fabs(lam(0)(3)), fabs(lam(1)(0)), fabs(lam(1)(1)),
                   fabs(lam(1)(2)), fabs(lam(1)(3))});

          flx = 0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
          break;
        }

        case flux_t::HLLE: {
          const CCTK_REAL charmax =
              max({0.0, lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
                   lam(1)(1), lam(1)(2), lam(1)(3)});

          const CCTK_REAL charmin =
              min({0.0, lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
                   lam(1)(1), lam(1)(2), lam(1)(3)});

          const CCTK_REAL charpm = charmax - charmin;

          flx = (charmax * flux(1) - charmin * flux(0) +
                 charmax * charmin * (var(1) - var(0))) /
                charpm;
          break;
        }

        default:
          assert(0);
        }

        return flx;
      };

  /* grid functions for fluxes */
  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss{fxdens, fydens, fzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBxs{fxBx, fyBx, fzBx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBys{fxBy, fyBy, fzBy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBzs{fxBz, fyBz, fzBz};
  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Reconstruct primitives from the cells on left (indice 0) and right
         * (indice 1) side of this face rc = reconstructed variables or
         * computed from reconstructed variables */
        const vec<CCTK_REAL, 2> rho_rc{reconstruct_pt(rho, p)};
        const vec<vec<CCTK_REAL, 2>, 3> vels_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_vels(i), p)};
        });
        const vec<CCTK_REAL, 2> eps_rc{reconstruct_pt(eps, p)};
        const vec<vec<CCTK_REAL, 2>, 3> Bs_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_Bvecs(i), p)};
        });

        /* Interpolate metric components from vertices to faces */
        const CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
        const vec<CCTK_REAL, 3> betas_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2f(gf_beta(i), p, dir);
        });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2f(gf_g(i, j), p, dir);
        });

        /* determinant of spatial metric */
        const CCTK_REAL detg_avg = calc_det(g_avg);
        const CCTK_REAL sqrtg = sqrt(detg_avg);
        /* co-velocity measured by Euleian observer: v_j */
        const vec<vec<CCTK_REAL, 2>, 3> vlows_rc =
            calc_contraction(g_avg, vels_rc);
        /* vtilde^i = alpha * v^i - beta^i */
        const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return alp_avg * vels_rc(i)(f) - betas_avg(i);
          });
        });
        /* Lorentz factor: W = 1 / sqrt(1 - v^2) */
        const vec<CCTK_REAL, 2> w_lorentz_rc([&](int f) ARITH_INLINE {
          return 1.0 / sqrt(1.0 - calc_contraction(vlows_rc, vels_rc)(f));
        });

        /* alpha * b0 = W * B^i * v_i */
        const vec<CCTK_REAL, 2> alp_b0_rc([&](int f) ARITH_INLINE {
          return w_lorentz_rc(f) * calc_contraction(Bs_rc, vlows_rc)(f);
        });
        /* covariant magnetic field measured by the Eulerian observer */
        const vec<vec<CCTK_REAL, 2>, 3> Blows_rc =
            calc_contraction(g_avg, Bs_rc);
        /* B^2 = B^i * B_i */
        const vec<CCTK_REAL, 2> B2_rc = calc_contraction(Bs_rc, Blows_rc);
        /* covariant magnetic field measured by the comoving observer:
         *  b_i = B_i/W + alpha*b^0*v_i */
        const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return Blows_rc(i)(f) / w_lorentz_rc(f) +
                   alp_b0_rc(f) * vlows_rc(i)(f);
          });
        });
        /* b^2 = b^{\mu} * b_{\mu} */
        const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
          return (B2_rc(f) + pow2(alp_b0_rc(f))) / pow2(w_lorentz_rc(f));
        });

        /* componets correspond to the dir we are considering */
        const CCTK_REAL beta_avg = betas_avg(dir);
        const vec<CCTK_REAL, 2> vel_rc{vels_rc(dir)};
        const vec<CCTK_REAL, 2> B_rc{Bs_rc(dir)};
        const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir)};

        // TODO: Compute pressure based on user-specified EOS.
        // Currently, computing press for classical ideal gas from reconstructed
        // vars

        // Ideal gas case {
        /* pressure for ideal gas EOS */
        const vec<CCTK_REAL, 2> press_rc([&](int f) ARITH_INLINE {
          return eps_rc(f) * rho_rc(f) * (gamma - 1);
        });
        /* cs2 for ideal gas EOS */
        const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
          return (gamma - 1.0) * eps_rc(f) / (eps_rc(f) + 1.0 / gamma);
        });
        /* enthalpy h for ideal gas EOS */
        const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
          return 1.0 + eps_rc(f) + press_rc(f) / rho_rc(f);
        });
        // } Ideal gas case

        /* Computing conservatives from primitives: */

        /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
        const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
          return sqrtg * rho_rc(f) * w_lorentz_rc(f);
        });

        /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
        const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
          return dens_rc(f) * h_rc(f) * w_lorentz_rc(f);
        });
        /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
        const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc =
            dens_h_W_rc + sqrtg * (pow2(alp_b0_rc) + B2_rc);
        /* auxiliary: (pgas + pmag) */
        const vec<CCTK_REAL, 2> press_plus_pmag_rc = press_rc + 0.5 * bsq_rc;

        /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
        const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return dens_h_W_plus_sqrtg_W2b2_rc(f) * vlows_rc(i)(f) -
                   sqrtg * alp_b0_rc(f) * blows_rc(i)(f);
          });
        });

        /* tau = sqrt(g)*t =
         *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
        const vec<CCTK_REAL, 2> tau_rc =
            dens_h_W_rc - dens_rc + sqrtg * (B2_rc - press_plus_pmag_rc);

        /* Btildes^i = sqrt(g) * B^i */
        const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
            [&](int i) ARITH_INLINE { return sqrtg * Bs_rc(i); });

        /* Computing fluxes of conserved variables: */

        /* auxiliary: unit in 'dir' */
        const vec<CCTK_REAL, 3> unit_dir{vec<int, 3>::unit(dir)};
        /* auxiliary: alpha * sqrt(g) */
        const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;
        /* auxiliary: B^i / W */
        const vec<CCTK_REAL, 2> B_over_w_lorentz_rc(
            [&](int f) ARITH_INLINE { return B_rc(f) / w_lorentz_rc(f); });

        /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
        const vec<CCTK_REAL, 2> flux_dens(
            [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

        /* flux(mom_j)^i = sqrt(g)*(
         *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
        const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return moms_rc(j)(f) * vtilde_rc(f) +
                   alp_sqrtg * (press_plus_pmag_rc(f) * unit_dir(j) -
                                blows_rc(j)(f) * B_over_w_lorentz_rc(f));
          });
        });

        /* flux(tau) = sqrt(g)*(
         *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
        const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
          return tau_rc(f) * vtilde_rc(f) +
                 alp_sqrtg * (press_plus_pmag_rc(f) * vel_rc(f) -
                              alp_b0_rc(f) * B_over_w_lorentz_rc(f));
        });

        /* electric field E_i = \tilde\epsilon_{ijk} Btilde_j * vtilde_k */
        const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
            calc_cross_product(Btildes_rc, vtildes_rc);
        /* flux(Btildes) = {{0, Ez, -Ey}, {-Ez, 0, Ex}, {Ey, -Ex, 0}} */
        const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
            calc_cross_product(unit_dir, Es_rc);

        /* Calculate eigenvalues: */

        /* variable for either g^xx, g^yy or g^zz depending on the direction */
        const CCTK_REAL u_avg = calc_inv(g_avg, detg_avg)(dir, dir);
        /* eigenvalues */
        vec<vec<CCTK_REAL, 4>, 2> lambda =
            eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
                        w_lorentz_rc, h_rc, bsq_rc);

        /* Calculate numerical fluxes */
        fluxdenss(dir)(p.I) = calcflux(lambda, dens_rc, flux_dens);
        fluxmomxs(dir)(p.I) = calcflux(lambda, moms_rc(0), flux_moms(0));
        fluxmomys(dir)(p.I) = calcflux(lambda, moms_rc(1), flux_moms(1));
        fluxmomzs(dir)(p.I) = calcflux(lambda, moms_rc(2), flux_moms(2));
        fluxtaus(dir)(p.I) = calcflux(lambda, tau_rc, flux_tau);
        fluxBxs(dir)(p.I) =
            (dir != 0) * calcflux(lambda, Btildes_rc(0), flux_Btildes(0));
        fluxBys(dir)(p.I) =
            (dir != 1) * calcflux(lambda, Btildes_rc(1), flux_Btildes(1));
        fluxBzs(dir)(p.I) =
            (dir != 2) * calcflux(lambda, Btildes_rc(2), flux_Btildes(2));
      });
}

void CalcAuxForAvecPsi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_Avecs{Avec_x, Avec_y, Avec_z};
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* interpolate A to vertices */
        const vec<CCTK_REAL, 3> A_vert([&](int i) ARITH_INLINE {
          return calc_avg_e2v(gf_Avecs(i), p, i);
        });
        const smat<CCTK_REAL, 3> g{gxx(p.I), gxy(p.I), gxz(p.I),
                                   gyy(p.I), gyz(p.I), gzz(p.I)};
        const vec<CCTK_REAL, 3> betas{betax(p.I), betay(p.I), betaz(p.I)};
        const CCTK_REAL detg = calc_det(g);
        const CCTK_REAL sqrtg = sqrt(detg);
        const smat<CCTK_REAL, 3> ug = calc_inv(g, detg);
        const vec<CCTK_REAL, 3> Aup = calc_contraction(ug, A_vert);

        Fx(p.I) = alp(p.I) * sqrtg * Aup(0);
        Fy(p.I) = alp(p.I) * sqrtg * Aup(1);
        Fz(p.I) = alp(p.I) * sqrtg * Aup(2);
        Fbetax(p.I) = betas(0) * Psi(p.I);
        Fbetay(p.I) = betas(1) * Psi(p.I);
        Fbetaz(p.I) = betas(2) * Psi(p.I);
        G(p.I) = alp(p.I) * Psi(p.I) / sqrtg - calc_contraction(betas, A_vert);
      });
}

extern "C" void AsterX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  CalcFlux<0>(cctkGH);
  CalcFlux<1>(cctkGH);
  CalcFlux<2>(cctkGH);

  /* Set auxiliary variables for the rhs of A and Psi  */
  CalcAuxForAvecPsi(cctkGH);
}

} // namespace AsterX
