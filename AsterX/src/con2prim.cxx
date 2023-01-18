#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>


#include <loop_device.hxx>

#include <cmath>

#include "con2primIdealFluid.hxx"
#include "utils.hxx"
#include <boost/math/tools/roots.hpp>

namespace AsterX
{
  using namespace std;
  using namespace Loop;

  /***************************************************************************
  2DNRNoble C2P
  ------------------------------------
  2D-NR Noble scheme for c2p.
  Sources: Noble+2006, Section 3.1 of Siegel+2018,
  NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
  ****************************************************************************/
  template <typename typeEoS>
  CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(
      typeEoS &plasma)
  { // Send con2primFactory object as reference to modify it,
    // and because we can not instantiate abstract class

    /* get Lorentz factor seed, calculated by constructor */
    CCTK_REAL W = plasma.WLorentz_Seed;

    /* get pressure seed */
    plasma.get_Press_Seed();

    /* get Z seed */
    plasma.get_Z_Seed();

    /* initialize unknowns for c2p, Z and vsq: */
    CCTK_REAL x[2];
    CCTK_REAL x_old[2];
    x[0] = fabs(plasma.Z_Seed);
    x[1] = (-1.0 + W * W) / (W * W);

    /* initialize old values */
    x_old[0] = x[0];
    x_old[1] = x[1];

    /* Start Recovery with 2D NR Solver */
    const CCTK_INT n = 2;
    const CCTK_REAL dv = (1. - 1.e-15);
    CCTK_REAL fvec[n];
    CCTK_REAL dx[n];
    CCTK_REAL fjac[n][n];

    CCTK_REAL detjac_inv;
    CCTK_REAL errf;
    plasma.Failed_2DNRNoble = 1;
    CCTK_INT k;
    for (k = 1; k <= plasma.maxIterations; k++)
    {
      fvec[0] = plasma.get_2DNRNoble_f0(x[0], x[1]);
      fvec[1] = plasma.get_2DNRNoble_f1(x[0], x[1]);
      fjac[0][0] = plasma.get_2DNRNoble_df0dZ(x[0], x[1]);
      fjac[0][1] = plasma.get_2DNRNoble_df0dVsq(x[0], x[1]);
      fjac[1][0] = plasma.get_2DNRNoble_df1dZ(x[0], x[1]);
      fjac[1][1] = plasma.get_2DNRNoble_df1dVsq(x[0], x[1]);
      detjac_inv = 1.0 / (fjac[0][0] * fjac[1][1] - fjac[0][1] * fjac[1][0]);
      dx[0] = -detjac_inv * (fjac[1][1] * fvec[0] - fjac[0][1] * fvec[1]);
      dx[1] = -detjac_inv * (-fjac[1][0] * fvec[0] + fjac[0][0] * fvec[1]);

      errf = 0.0;
      for (CCTK_INT i = 0; i < n; i++)
      {
        errf += fabs(fvec[i]);
      }
      if (errf <= plasma.tolerance)
      {
        plasma.Failed_2DNRNoble = 0;
        break;
      }

      /* save old values before calculating the new */
      x_old[0] = x[0];
      x_old[1] = x[1];

      for (CCTK_INT i = 0; i < n; i++)
      {
        x[i] += dx[i];
      }
    }
    plasma.Nit_2DNRNoble = k;

    /* make sure that the new x[] is physical */
    if (x[0] < 0.0)
    {
      x[0] = fabs(x[0]);
    }
    else
    {
      if (x[0] > 1e20)
      {
        x[0] = x_old[0];
      }
    }

    if (x[1] < 0.0)
    {
      x[1] = 0.0;
    }
    else
    {
      if (x[1] >= 1.0)
      {
        x[1] = dv;
      }
    }

    /* Calculate primitives from Z and W */
    plasma.Z_Sol = x[0];
    plasma.vsq_Sol = x[1];
  }

  /***************************************************************************
  1DBrentPalenzuela C2P
  ------------------------------------
  1D-Brent Palenzuela scheme for c2p.
  Sources: Nielsen+2014, Palenzuela+2015, Section 3.4.1 of Siegel+2018,
  https://en.wikipedia.org/wiki/Brent%27s_method
  ****************************************************************************/
  template <typename F, typename typeEoS>
  CCTK_HOST CCTK_DEVICE std::pair<CCTK_REAL, CCTK_REAL> brent(F f, int min_bits, int max_iters,
                                        int &iters, typeEoS &plasma)
  {
    using std::abs, std::min, std::max, std::swap;

    CCTK_REAL qPalenzuela =
        plasma.ConservedVars[TAU] / plasma.ConservedVars[CONS_D];
    // CCTK_REAL rPalenzuela = plasma.Ssq / pow(plasma.ConservedVars[CONS_D], 2);
    CCTK_REAL sPalenzuela = plasma.Bsq / plasma.ConservedVars[CONS_D];
    // CCTK_REAL tPalenzuela = plasma.BiSi / pow(plasma.ConservedVars[CONS_D], 3.
    // / 2.);

    CCTK_REAL xPalenzuela_lowerBound = 1.0 + qPalenzuela - sPalenzuela;
    CCTK_REAL xPalenzuela_upperBound = 2.0 + 2.0 * qPalenzuela - sPalenzuela;

    CCTK_REAL a = xPalenzuela_lowerBound;
    CCTK_REAL b = xPalenzuela_upperBound;
    // CCTK_REAL c = xPalenzuela_upperBound;

    auto tol = boost::math::tools::eps_tolerance<CCTK_REAL>(min_bits);

    iters = 0;
    auto fa = f(a, plasma);
    auto fb = f(b, plasma);
    if (abs(fa) < abs(fb))
    {
      swap(a, b);
      swap(fa, fb);
    }
    if (fb == 0)
      return {b, b};
    if (fa * fb >= 0)
    {
      // Root is not bracketed
      iters = max_iters;
      return {min(a, b), max(a, b)};
    }
    CCTK_REAL c = a;
    auto fc = fa;
    bool mflag = true;
    CCTK_REAL d{};

    while (fb != 0 && !tol(a, b) && iters < max_iters)
    {
      CCTK_REAL s;
      if (fa != fc && fb != fc)
        // inverse quadratic interpolation
        s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
            (b * fa * fc) / ((fb - fa) * (fb - fc)) +
            (c * fa * fb) / ((fc - fa) * (fc - fb));
      else
        // secant method
        s = (a + b) / 2 - (fa + fb) / 2 * (b - a) / (fb - fa);
      CCTK_REAL u = (3 * a + b) / 4;
      CCTK_REAL v = b;
      if (u > v)
        swap(u, v);
      bool cond1 = !(u <= s && s <= v);
      bool cond2 = mflag && abs(s - b) >= abs(b - c) / 2;
      bool cond3 = !mflag && abs(s - b) >= abs(c - d) / 2;
      bool cond4 = mflag && tol(c, b);
      bool cond5 = !mflag && tol(c, d);
      if (cond1 || cond2 || cond3 || cond4 || cond5)
      {
        // bisection
        s = (a + b) / 2;
        mflag = true;
      }
      else
      {
        mflag = false;
      }
      auto fs = f(s, plasma);
      // `d` is assigned for the first time here; it won't be used above on the
      // first iteration because `mflag` is set
      d = c;
      c = b;
      fc = fb;
      if (fa * fs < 0)
      {
        b = s;
        fb = fs;
      }
      else
      {
        a = s;
        fa = fs;
      }
      // CCTK_VINFO("iters=%d mflag=%d   a=%.17g b=%.17g c=%.17g d=%.17g fa=%.17g"
      //            "fb=%.17g fc=%.17g",
      //            iters, int(mflag), double(a), double(b), double(c), double(d),
      //            double(fa), double(fb), double(fc));
      if (fa * fb >= 0)
      {
        return {min(a, b), max(a, b)};
      }
      if (abs(fa) < abs(fb))
      {
        swap(a, b);
        swap(fa, fb);
      }
      ++iters;
    }

    if (fb == 0)
      return {b, b};
    return {min(a, b), max(a, b)};
  }

  template <typename typeEoS>
  CCTK_HOST CCTK_DEVICE CCTK_REAL funcRoot_1DPalenzuela(CCTK_REAL x,
                                                        typeEoS &plasma)
  {

    // computes f(x) from x and q,r,s,t
    const CCTK_REAL qPalenzuela =
        plasma.ConservedVars[TAU] / plasma.ConservedVars[CONS_D];
    const CCTK_REAL rPalenzuela =
        plasma.Ssq / pow(plasma.ConservedVars[CONS_D], 2);
    const CCTK_REAL sPalenzuela = plasma.Bsq / plasma.ConservedVars[CONS_D];
    const CCTK_REAL tPalenzuela =
        plasma.BiSi / pow(plasma.ConservedVars[CONS_D], 3. / 2.);

    // (i)
    CCTK_REAL Wminus2 =
        1.0 - (x * x * rPalenzuela +
               (2.0 * x + sPalenzuela) * tPalenzuela * tPalenzuela) /
                  (x * x * (x + sPalenzuela) * (x + sPalenzuela));
    Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
    const CCTK_REAL W_loc = pow(Wminus2, -0.5);

    // (ii)
    CCTK_REAL rho_loc = plasma.ConservedVars[CONS_D] / W_loc;

    // (iii)
    CCTK_REAL eps_loc = W_loc - 1.0 + (1.0 - W_loc * W_loc) * x / W_loc +
                        W_loc * (qPalenzuela - sPalenzuela +
                                 tPalenzuela * tPalenzuela / (2 * x * x) +
                                 sPalenzuela / (2 * W_loc * W_loc));

    // (iv)
    CCTK_REAL P_loc = plasma.get_Press_funcRhoEps(rho_loc, eps_loc);

    return (x - (1.0 + eps_loc + P_loc / rho_loc) * W_loc);
  }

  template <typename typeEoS>
  CCTK_HOST CCTK_DEVICE void Con2Prim_1DBrentPalenzuela(
      typeEoS &plasma)
  { // Send con2primFactory object as reference to modify it,
    // and because we can not instantiate abstract class

    plasma.Failed_1DBrentPalenzuela = 1;

    // Find x, this is the recovery process
    const int minbits = std::numeric_limits<CCTK_REAL>::digits - 4;
    const int maxiters = plasma.maxIterations;
    int iters;
    std::pair<CCTK_REAL, CCTK_REAL> result =
        brent(*funcRoot_1DPalenzuela<typeEoS>, minbits, maxiters, iters, plasma);

    // Pick best solution
    if (abs(funcRoot_1DPalenzuela(result.first, plasma)) <
        abs(funcRoot_1DPalenzuela(result.second, plasma)))
    {

      plasma.xPalenzuela_Sol = result.first;
    }
    else
    {
      plasma.xPalenzuela_Sol = result.second;
    }

    // Check solution and calculate primitives
    if (iters < maxiters &&
        abs(funcRoot_1DPalenzuela(result.first, plasma)) < plasma.tolerance)
    {
      plasma.Failed_1DBrentPalenzuela = 0;
      plasma.xPalenzuelaToPrim();
    }

    return;
  }

  ///

  template <typename typeEoS>
  void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
    DECLARE_CCTK_PARAMETERS;

    const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
    const vec<GF3D2<CCTK_REAL>, 6> gf_saved_prims{
        saved_rho, saved_velx, saved_vely, saved_velz, saved_eps, press};
    const vec<GF3D2<CCTK_REAL>, 6> gf_prims{rho, velx, vely, velz, eps, press};
    const vec<GF3D2<CCTK_REAL>, 5> gf_cons{dens, momx, momy, momz, tau};

    // Loop over the interior of the grid
    cctk_grid.loop_int_device<
        1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                       const PointDesc
                                           &p) CCTK_ATTRIBUTE_ALWAYS_INLINE
                 {
    CCTK_REAL g_up[4][4];
    CCTK_REAL g_lo[4][4];

    /* Get covariant metric */
    const smat<CCTK_REAL, 3> g3_avg(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });
    g_lo[1][1] = g3_avg(0, 0);
    g_lo[1][2] = g3_avg(0, 1);
    g_lo[1][3] = g3_avg(0, 2);
    g_lo[2][2] = g3_avg(1, 1);
    g_lo[2][3] = g3_avg(1, 2);
    g_lo[3][3] = g3_avg(2, 2);
    CCTK_REAL lapse = calc_avg_v2c(alp, p);
    CCTK_REAL betax_up = calc_avg_v2c(betax, p);
    CCTK_REAL betay_up = calc_avg_v2c(betay, p);
    CCTK_REAL betaz_up = calc_avg_v2c(betaz, p);

    // beta_lo
    g_lo[0][1] =
        g_lo[1][1] * betax_up + g_lo[1][2] * betay_up + g_lo[1][3] * betaz_up;
    g_lo[0][2] =
        g_lo[1][2] * betax_up + g_lo[2][2] * betay_up + g_lo[2][3] * betaz_up;
    g_lo[0][3] =
        g_lo[1][3] * betax_up + g_lo[2][3] * betay_up + g_lo[3][3] * betaz_up;

    g_lo[0][0] = -lapse * lapse + betax_up * g_lo[0][1] +
                 betay_up * g_lo[0][2] + betaz_up * g_lo[0][3];

    // symmetric components
    g_lo[1][0] = g_lo[0][1];
    g_lo[2][0] = g_lo[0][2];
    g_lo[3][0] = g_lo[0][3];
    g_lo[2][1] = g_lo[1][2];
    g_lo[3][1] = g_lo[1][3];
    g_lo[3][2] = g_lo[2][3];

    /* Calculate inverse of 4-dim metric */
    const CCTK_REAL spatial_detg = calc_det(g3_avg);
    const CCTK_REAL sqrtg = sqrt(spatial_detg);
    const smat<CCTK_REAL, 3> g3inv_avg = calc_inv(g3_avg, spatial_detg);
    g_up[0][0] = -1.0 / (lapse * lapse);
    g_up[0][1] = betax_up / (lapse * lapse);
    g_up[0][2] = betay_up / (lapse * lapse);
    g_up[0][3] = betaz_up / (lapse * lapse);
    g_up[1][1] = g3inv_avg(0, 0) - betax_up * betax_up / (lapse * lapse);
    g_up[1][2] = g3inv_avg(0, 1) - betax_up * betay_up / (lapse * lapse);
    g_up[1][3] = g3inv_avg(0, 2) - betax_up * betaz_up / (lapse * lapse);
    g_up[2][2] = g3inv_avg(1, 1) - betay_up * betay_up / (lapse * lapse);
    g_up[2][3] = g3inv_avg(1, 2) - betay_up * betaz_up / (lapse * lapse);
    g_up[3][3] = g3inv_avg(2, 2) - betaz_up * betaz_up / (lapse * lapse);

    g_up[1][0] = g_up[0][1];
    g_up[2][0] = g_up[0][2];
    g_up[3][0] = g_up[0][3];
    g_up[2][1] = g_up[1][2];
    g_up[3][1] = g_up[1][3];
    g_up[3][2] = g_up[2][3];

    const vec<CCTK_REAL, 3> Bup{dBx(p.I) / sqrtg, dBy(p.I) / sqrtg,
                                dBz(p.I) / sqrtg};

    // Set to atmosphere if dens below the threshold
    if (dens(p.I) <= sqrtg * (rho_abs_min * (1.0 + atmo_tol))) {
      const vec<CCTK_REAL, 3> Blow = calc_contraction(g3_avg, Bup);
      const CCTK_REAL Bsq = calc_contraction(Bup, Blow);
      set_to_atmosphere(rho_abs_min, poly_K, gamma, sqrtg, Bsq, gf_saved_prims,
                        gf_cons, p);
    }

    CCTK_REAL cons[NCONS];
    cons[CONS_D] = dens(p.I) / sqrtg;
    cons[S1_COV] = momx(p.I) / sqrtg;
    cons[S2_COV] = momy(p.I) / sqrtg;
    cons[S3_COV] = momz(p.I) / sqrtg;
    cons[TAU] = tau(p.I) / sqrtg;
    cons[B1] = Bup(0);
    cons[B2] = Bup(1);
    cons[B3] = Bup(2);

    CCTK_REAL prims[NPRIMS];
    prims[RHO] = saved_rho(p.I);
    prims[V1_CON] = saved_velx(p.I);
    prims[V2_CON] = saved_vely(p.I);
    prims[V3_CON] = saved_velz(p.I);
    prims[EPS] = saved_eps(p.I);
    prims[B1] = cons[B1];
    prims[B2] = cons[B2];
    prims[B3] = cons[B3];

    // Construct con2primFactory object:
    typeEoS plasma_0(gamma, cons, prims, g_lo, g_up);
    plasma_0.tolerance = c2p_tol;
    plasma_0.maxIterations = max_iter;

    // 1) Try 2DNRNoble
    Con2Prim_2DNRNoble(plasma_0);
    if (plasma_0.Failed_2DNRNoble) {
      // 1) Try 1DBrentPalenzuela
      Con2Prim_1DBrentPalenzuela(plasma_0);
    } else {
      plasma_0.WZ2Prim(); // 2DNRNoble to prims
    }

    if (plasma_0.Failed_2DNRNoble && plasma_0.Failed_1DBrentPalenzuela) {
      //printf("c2p failed.\n");
      /* set flag to failure */
      con2prim_flag(p.I) = 0.0;

      if (debug_mode) {
        printf(
            "WARNING: "
            "2DNRNoble failed. Printing cons and saved prims before set to "
            "atmo: \n"
            "cctk_iteration = %i \n "
            "x, y, z = %26.16e, %26.16e, %26.16e \n "
            "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
            "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
            "dBy = %26.16e \n dBz = %26.16e \n "
            "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e \n "
            "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
            "%26.16e \n "
            "Bvecx = %26.16e \n Bvecy = %26.16e \n "
            "Bvecz = %26.16e \n "
            "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
            cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
            momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), saved_rho(p.I),
            saved_eps(p.I), press(p.I), saved_velx(p.I), saved_vely(p.I),
            saved_velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Avec_x(p.I),
            Avec_y(p.I), Avec_z(p.I));
      }

      rho(p.I) = saved_rho(p.I);
      velx(p.I) = saved_velx(p.I);
      vely(p.I) = saved_vely(p.I);
      velz(p.I) = saved_velz(p.I);
      eps(p.I) = saved_eps(p.I);
      press(p.I) = (gamma - 1.0) * eps(p.I) * rho(p.I);

      // assert(0);
    } else {
      /* set flag to success */
      con2prim_flag(p.I) = 1.0;

      rho(p.I) = plasma_0.PrimitiveVars[RHO];
      velx(p.I) = plasma_0.PrimitiveVars[V1_CON];
      vely(p.I) = plasma_0.PrimitiveVars[V2_CON];
      velz(p.I) = plasma_0.PrimitiveVars[V3_CON];
      eps(p.I) = plasma_0.PrimitiveVars[EPS];
      press(p.I) = (gamma - 1.0) * eps(p.I) * rho(p.I);

      if (rho(p.I) <= rho_abs_min * (1 + atmo_tol)) {
        const vec<CCTK_REAL, 3> Blow = calc_contraction(g3_avg, Bup);
        const CCTK_REAL Bsq = calc_contraction(Bup, Blow);
        set_to_atmosphere(rho_abs_min, poly_K, gamma, sqrtg, Bsq, gf_prims,
                          gf_cons, p);

        if (debug_mode) {
          printf("WARNING: "
                 "Printing cons and prims after set to atmo when dens below "
                 "threshold: \n"
                 "cctk_iteration = %i \n "
                 "x, y, z = %26.16e, %26.16e, %26.16e \n "
                 "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
                 "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
                 "dBy = %26.16e \n dBz = %26.16e \n "
                 "rho = %26.16e \n eps = %26.16e \n press= %26.16e \n "
                 "velx = %26.16e \n vely = %26.16e \n velz = %26.16e \n "
                 "Bvecx = %26.16e \n Bvecy = %26.16e \n "
                 "Bvecz = %26.16e \n "
                 "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
                 cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
                 momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), rho(p.I),
                 eps(p.I), press(p.I), velx(p.I), vely(p.I), velz(p.I),
                 Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Avec_x(p.I), Avec_y(p.I),
                 Avec_z(p.I));
        }
      }
    }

    Bvecx(p.I) = cons[B1];
    Bvecy(p.I) = cons[B2];
    Bvecz(p.I) = cons[B3];
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I); }); // Loop
  }                    // AsterX_Con2Prim_2DNRNoble

  /***************************************************************************
   * AsterX_Con2Prim
   * ---
   *  Routines implemented:
   *   1) 2DNRNoble
   *
   *   Based on con2primFactory (https://github.com/fedelopezar/con2primFactory)
   *   ****************************************************************************/
  extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS)
  {
    if (1)
    { // Use this if for idealFluid/tabeos
      AsterX_Con2Prim_typeEoS<idealFluid>(CCTK_PASS_CTOC);
      // CCTK_PASS_CTOC == cctkGH, and more. Preferred over just cctkGH.
    }
  }

  extern "C" void AsterX_Con2Prim_Interpolate_Failed(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim_Interpolate_Failed;
    DECLARE_CCTK_PARAMETERS;

    const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
    const vec<GF3D2<CCTK_REAL>, 6> gf_prims{rho, velx, vely, velz, eps, press};
    const vec<GF3D2<CCTK_REAL>, 5> gf_cons{dens, momx, momy, momz, tau};

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE
        {
          if (con2prim_flag(p.I) == 0.0)
          {

            const vec<CCTK_REAL, 6> flag_nbs = get_neighbors(con2prim_flag, p);
            const vec<CCTK_REAL, 6> rho_nbs = get_neighbors(rho, p);
            const vec<CCTK_REAL, 6> velx_nbs = get_neighbors(velx, p);
            const vec<CCTK_REAL, 6> vely_nbs = get_neighbors(vely, p);
            const vec<CCTK_REAL, 6> velz_nbs = get_neighbors(velz, p);
            const vec<CCTK_REAL, 6> eps_nbs = get_neighbors(eps, p);
            const vec<CCTK_REAL, 6> saved_rho_nbs = get_neighbors(saved_rho, p);
            const vec<CCTK_REAL, 6> saved_velx_nbs = get_neighbors(saved_velx, p);
            const vec<CCTK_REAL, 6> saved_vely_nbs = get_neighbors(saved_vely, p);
            const vec<CCTK_REAL, 6> saved_velz_nbs = get_neighbors(saved_velz, p);
            const vec<CCTK_REAL, 6> saved_eps_nbs = get_neighbors(saved_eps, p);

            CCTK_REAL sum_nbs =
                sum<6>([&](int i) ARITH_INLINE
                       { return flag_nbs(i); });
            assert(sum_nbs > 0);
            rho(p.I) = calc_avg_neighbors(flag_nbs, rho_nbs, saved_rho_nbs);
            velx(p.I) = calc_avg_neighbors(flag_nbs, velx_nbs, saved_velx_nbs);
            vely(p.I) = calc_avg_neighbors(flag_nbs, vely_nbs, saved_vely_nbs);
            velz(p.I) = calc_avg_neighbors(flag_nbs, velz_nbs, saved_velz_nbs);
            eps(p.I) = calc_avg_neighbors(flag_nbs, eps_nbs, saved_eps_nbs);
            press(p.I) = (gamma - 1.0) * eps(p.I) * rho(p.I);

            /* reset flag */
            con2prim_flag(p.I) = 1.0;

            // set to atmos
            if (rho(p.I) <= rho_abs_min * (1 + atmo_tol))
            {
              const smat<CCTK_REAL, 3> g3_avg([&](int i, int j) ARITH_INLINE
                                              { return calc_avg_v2c(gf_g(i, j), p); });
              const CCTK_REAL sqrtg = sqrt(calc_det(g3_avg));
              const vec<CCTK_REAL, 3> Bup{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
              const vec<CCTK_REAL, 3> Blow = calc_contraction(g3_avg, Bup);
              const CCTK_REAL Bsq = calc_contraction(Bup, Blow);

              set_to_atmosphere(rho_abs_min, poly_K, gamma, sqrtg, Bsq, gf_prims,
                                gf_cons, p);
            };

            saved_rho(p.I) = rho(p.I);
            saved_velx(p.I) = velx(p.I);
            saved_vely(p.I) = vely(p.I);
            saved_velz(p.I) = velz(p.I);
            saved_eps(p.I) = eps(p.I);
          }
        });
  }

} // namespace AsterX
