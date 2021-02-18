/* TwoPunctures:  File  "TwoPunctures.c"*/

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>
#include <loopcontrol.h>
#include "TP_utilities.h"
#include "TwoPunctures.h"

/* Swap two variables */
static inline
void tp_swap (CCTK_REAL * restrict const a, CCTK_REAL * restrict const b)
{
  CCTK_REAL const t = *a; *a=*b; *b=t;
}
#undef SWAP
#define SWAP(a,b) (tp_swap(&(a),&(b)))



static
void set_initial_guess(CCTK_POINTER_TO_CONST cctkGH,
                       derivs v)
{
  DECLARE_CCTK_PARAMETERS;

  int nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  CCTK_REAL *s_x, *s_y, *s_z;
  CCTK_REAL al, A, Am1, be, B, phi, R, r, X;
  CCTK_INT ivar, i, j, k, i3D, indx;
  derivs U;
  FILE *debug_file;

  if (solve_momentum_constraint)
    nvar = 4;

  s_x    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_y    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  s_z    =calloc(n1*n2*n3, sizeof(CCTK_REAL));
  allocate_derivs (&U, nvar);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          i3D = Index(ivar,i,j,k,1,n1,n2,n3);

          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;

          /* Calculation of (X,R)*/
          AB_To_XR (nvar, A, B, &X, &R, U);
          /* Calculation of (x,r)*/
          C_To_c (nvar, X, R, &(s_x[i3D]), &r, U);
          /* Calculation of (y,z)*/
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);
        }
  Set_Initial_Guess_for_u(cctkGH, n1*n2*n3, v.d0, s_x, s_y, s_z);
  for (ivar = 0; ivar < nvar; ivar++)
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
        for (k = 0; k < n3; k++)
        {
          indx = Index(ivar,i,j,k,1,n1,n2,n3);
          v.d0[indx]/=(-cos(Pih * (2 * i + 1) / n1)-1.0);
        }
  Derivatives_AB3 (nvar, n1, n2, n3, v);
  if (do_initial_debug_output && CCTK_MyProc(cctkGH) == 0)
  {
    debug_file=fopen("initial.dat", "w");
    assert(debug_file);
    for (ivar = 0; ivar < nvar; ivar++)
      for (i = 0; i < n1; i++)
        for (j = 0; j < n2; j++)
        {
          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          Am1 = A -1.0;
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 0.0;
          indx = Index(ivar,i,j,0,1,n1,n2,n3);
          U.d0[0] = Am1 * v.d0[indx];        /* U*/
          U.d1[0] = v.d0[indx] + Am1 * v.d1[indx];        /* U_A*/
          U.d2[0] = Am1 * v.d2[indx];        /* U_B*/
          U.d3[0] = Am1 * v.d3[indx];        /* U_3*/
          U.d11[0] = 2 * v.d1[indx] + Am1 * v.d11[indx];        /* U_AA*/
          U.d12[0] = v.d2[indx] + Am1 * v.d12[indx];        /* U_AB*/
          U.d13[0] = v.d3[indx] + Am1 * v.d13[indx];        /* U_AB*/
          U.d22[0] = Am1 * v.d22[indx];        /* U_BB*/
          U.d23[0] = Am1 * v.d23[indx];        /* U_B3*/
          U.d33[0] = Am1 * v.d33[indx];        /* U_33*/
        /* Calculation of (X,R)*/
        AB_To_XR (nvar, A, B, &X, &R, U);
        /* Calculation of (x,r)*/
        C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
        /* Calculation of (y,z)*/
        rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
        fprintf(debug_file,
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g "
                "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                (double)s_x[indx], (double)s_y[indx],
                (double)A,(double)B,
                (double)U.d0[0],
                (double)(-cos(Pih * (2 * i + 1) / n1)-1.0),
                (double)U.d1[0],
                (double)U.d2[0],
                (double)U.d3[0],
                (double)U.d11[0],
                (double)U.d22[0],
                (double)U.d33[0],
                (double)v.d0[indx],
                (double)v.d1[indx],
                (double)v.d2[indx],
                (double)v.d3[indx],
                (double)v.d11[indx],
                (double)v.d22[indx],
                (double)v.d33[indx]
                );
        }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2; i++)
    {
      CCTK_REAL d;
      indx = Index(0,0,i,0,1,n1,n2,n3);
      d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
              s_x[indx], 0.0, 0.0);
      fprintf(debug_file, "%.16g %.16g\n",
                (double)s_x[indx], (double)d);
    }
    fprintf(debug_file, "\n\n");
    for (i=n2-10; i<n2-1; i++)
    {
      CCTK_REAL d;
      int ip= Index(0,0,i+1,0,1,n1,n2,n3);
      indx = Index(0,0,i,0,1,n1,n2,n3);
      for (j=-10; j<10; j++)
      {
        d = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,
                s_x[indx]+(s_x[ip]-s_x[indx])*j/10,
                0.0, 0.0);
        fprintf(debug_file, "%.16g %.16g\n",
                (double)(s_x[indx]+(s_x[ip]-s_x[indx])*j/10), (double)d);
      }
    }
    fprintf(debug_file, "\n\n");
    for (i = 0; i < n1; i++)
      for (j = 0; j < n2; j++)
      {
        X = 2*(2.0*i/n1-1.0);
        R = 2*(1.0*j/n2);
        if (X*X+R*R > 1.0)
        {
          C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[indx]), &(s_z[indx]), U);
          *U.d0  = s_x[indx]*s_x[indx];
          *U.d1  = 2*s_x[indx];
          *U.d2  = 0.0;
          *U.d3 = 0.0;
          *U.d11 = 2.0;
          *U.d22 = 0.0;
          *U.d33 = *U.d12 = *U.d23 = *U.d13 = 0.0;
          C_To_c (nvar, X, R, &(s_x[indx]), &r, U);
          fprintf(debug_file,
                  "%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n",
                  (double)s_x[indx], (double)r, (double)X, (double)R, (double)U.d0[0],
                  (double)U.d1[0],
                  (double)U.d2[0],
                  (double)U.d3[0],
                  (double)U.d11[0],
                  (double)U.d22[0],
                  (double)U.d33[0]);
        }
      }
    fclose(debug_file);
  }
  free(s_z);
  free(s_y);
  free(s_x);
  free_derivs (&U, nvar);
  /*exit(0);*/
}

/* -------------------------------------------------------------------*/
void TwoPunctures(CCTK_ARGUMENTS);
void
TwoPunctures (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TwoPunctures;
  DECLARE_CCTK_PARAMETERS;

  * mp = par_m_plus;
  * mm = par_m_minus;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  int const nvar = 1, n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int const ntotal = n1 * n2 * n3 * nvar;
#if 0
  int percent10 = 0;
#endif
  static bool did_setup = false;
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  static CCTK_REAL mp_saved, mm_saved, mp_adm_saved, mm_adm_saved, E_saved, J1_saved, J2_saved, J3_saved;
  CCTK_REAL admMass;

  if (! did_setup) {
#pragma omp critical
    if (! did_setup) {
    CCTK_REAL up, um;
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);
    allocate_derivs (&cf_v, ntotal);

    if (use_sources) {
      CCTK_INFO ("Solving puncture equation for BH-NS/NS-NS system");
    } else {
      CCTK_INFO ("Solving puncture equation for BH-BH system");
    }
    CCTK_VINFO ("b = %g", par_b);
    
    /* initialise to 0 */
    for (int j = 0; j < ntotal; j++)
    {
      cf_v.d0[j] = 0.0;
      cf_v.d1[j] = 0.0;
      cf_v.d2[j] = 0.0;
      cf_v.d3[j] = 0.0;
      cf_v.d11[j] = 0.0;
      cf_v.d12[j] = 0.0;
      cf_v.d13[j] = 0.0;
      cf_v.d22[j] = 0.0;
      cf_v.d23[j] = 0.0;
      cf_v.d33[j] = 0.0;
      v.d0[j] = 0.0;
      v.d1[j] = 0.0;
      v.d2[j] = 0.0;
      v.d3[j] = 0.0;
      v.d11[j] = 0.0;
      v.d12[j] = 0.0;
      v.d13[j] = 0.0;
      v.d22[j] = 0.0;
      v.d23[j] = 0.0;
      v.d33[j] = 0.0;
    }
    /* call for external initial guess */
    if (use_external_initial_guess)
    {
      set_initial_guess(cctkGH, v);
    }

    /* If bare masses are not given, iteratively solve for them given the 
       target ADM masses target_M_plus and target_M_minus and with initial 
       guesses given by par_m_plus and par_m_minus. */
    if(!(give_bare_mass)) {
      CCTK_REAL tmp, mp_adm_err, mm_adm_err;
      char valbuf[100];

      CCTK_REAL M_p = target_M_plus;
      CCTK_REAL M_m = target_M_minus;

      CCTK_VINFO ("Attempting to find bare masses.");
      CCTK_VINFO ("Target ADM masses: M_p=%g and M_m=%g",
                  (double) M_p, (double) M_m);
      CCTK_VINFO ("ADM mass tolerance: %g", (double) adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        CCTK_VINFO ("Bare masses: mp=%.15g, mm=%.15g",
                    (double)*mp, (double)*mm);
        Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, 1);

        F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

        up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * par_b);
        *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p-*mp_adm);
        mm_adm_err = fabs(M_m-*mm_adm);
        CCTK_VINFO ("ADM mass error: M_p_err=%.15g, M_m_err=%.15g",
                    (double)mp_adm_err, (double)mm_adm_err);
        
        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*par_b*( 1 + um + up + um*up ) + 
                sqrt(16*par_b*M_m*(1 + um)*(1 + up) + 
                  pow(-M_m + M_p + 4*par_b*(1 + um)*(1 + up),2));
        *mp = (tmp + M_p - M_m)/(2.*(1 + up));
        *mm = (tmp - M_p + M_m)/(2.*(1 + um));
        
        /* Set the par_m_plus and par_m_minus parameters */
        sprintf (valbuf,"%.17g", (double) *mp);
        CCTK_ParameterSet ("par_m_plus", "TwoPunctures", valbuf);
        
        sprintf (valbuf,"%.17g", (double) *mm);
        CCTK_ParameterSet ("par_m_minus", "TwoPunctures", valbuf);
        
      } while ( (mp_adm_err > adm_tol) ||
                (mm_adm_err > adm_tol) );
                
      CCTK_VINFO ("Found bare masses.");
    }

    Newton (cctkGH, nvar, n1, n2, n3, v, Newton_tol, Newton_maxit);
 
    F_of_v (cctkGH, nvar, n1, n2, n3, v, F, u);

    SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

    CCTK_VINFO ("The two puncture masses are mp=%.17g and mm=%.17g",
                (double) *mp, (double) *mm);

    up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, par_b, 0., 0.);
    um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * par_b);
    *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * par_b);

    CCTK_VINFO ("Puncture 1 ADM mass is %g", (double) *mp_adm);
    CCTK_VINFO ("Puncture 2 ADM mass is %g", (double) *mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (*mp + *mm
               - 4*par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    CCTK_VINFO ("The total ADM mass is %g", (double) admMass);
    *E = admMass;

    /*
      Run this in Mathematica (version 8 or later) with
        math -script <file>

      Needs["SymbolicC`"];
      co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
      r1 = co + {"par_b", 0, 0};
      r2 = co + {-"par_b", 0, 0};
      {p1, p2} = Table["par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];
      {s1, s2} = Table["par_S_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

      J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

      JVar = Table["*J" <> ToString[i], {i, 1, 3}];
      Print[OutputForm@StringReplace[
        ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}], 
        "\"" -> ""]];
     */

    *J1 = -(center_offset[2]*par_P_minus[1]) + center_offset[1]*par_P_minus[2] - center_offset[2]*par_P_plus[1] + center_offset[1]*par_P_plus[2] + par_S_minus[0] + par_S_plus[0];
    *J2 = center_offset[2]*par_P_minus[0] - center_offset[0]*par_P_minus[2] + par_b*par_P_minus[2] + center_offset[2]*par_P_plus[0] - center_offset[0]*par_P_plus[2] - par_b*par_P_plus[2] + par_S_minus[1] + par_S_plus[1];
    *J3 = -(center_offset[1]*par_P_minus[0]) + center_offset[0]*par_P_minus[1] - par_b*par_P_minus[1] - center_offset[1]*par_P_plus[0] + center_offset[0]*par_P_plus[1] + par_b*par_P_plus[1] + par_S_minus[2] + par_S_plus[2];

    // store these in local variables so that we can restore them once CarpetX
    // can wipes the grid scalars
    mp_saved = *mp;
    mm_saved = *mm;
    mp_adm_saved = *mp_adm;
    mm_adm_saved = *mm_adm;
    E_saved = *E;
    J1_saved = *J1;
    J2_saved = *J2;
    J3_saved = *J3;

    did_setup = true;
    }
  }

  // before each call CarpetX wipes the grid scalars so I need to restore them
  *mp = mp_saved;
  *mm = mm_saved;
  *mp_adm = mp_adm_saved;
  *mm_adm = mm_adm_saved;
  *E = E_saved;
  *J1 = J1_saved;
  *J2 = J2_saved;
  *J3 = J3_saved;

  if (CCTK_EQUALS(grid_setup_method, "Taylor expansion"))
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (CCTK_EQUALS(grid_setup_method, "evaluation"))
  {
    gsm = GSM_evaluation;
  }
  else
  {
    CCTK_ERROR ("internal error");
  }

  antisymmetric_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-antisymmetric");
  averaged_lapse = CCTK_EQUALS(initial_lapse, "twopunctures-averaged");
	pmn_lapse = CCTK_EQUALS(initial_lapse, "psi^n");
  if (pmn_lapse)
		CCTK_VINFO("Setting initial lapse to psi^%f profile.",
               (double)initial_lapse_psi_exponent);
  brownsville_lapse = CCTK_EQUALS(initial_lapse, "brownsville");
  if (brownsville_lapse)
    CCTK_VINFO("Setting initial lapse to a Brownsville-style profile "
               "with exp %f.",
               (double)initial_lapse_psi_exponent);

  CCTK_INFO ("Interpolating result");

  const int di = 1;
  const int dj = di * (cctk_ash[0] + 1); // one extra grid point for vertex centering
  const int dk = dj * (cctk_ash[1] + 1); // one extra grid point for vertex centering
  const int np = dk * (cctk_ash[2] + 1); // one extra grid point for vertex centering

  for(int i=0; i < cctk_lsh[0]+1; i++)
  for(int j=0; j < cctk_lsh[1]+1; j++)
  for(int k=0; k < cctk_lsh[2]+1; k++)
//  CCTK_LOOP3_ALL(TwoPunctures, cctkGH, i,j,k)
      {

        const int ind = i*di + j*dj + k*dk;

        CCTK_REAL xx, yy, zz;
        xx = vcoordx[ind] - center_offset[0];
        yy = vcoordy[ind] - center_offset[1];
        zz = vcoordz[ind] - center_offset[2];

        /* We implement swapping the x and z coordinates as follows.
           The bulk of the code that performs the actual calculations
           is unchanged.  This code looks only at local variables.
           Before the bulk --i.e., here-- we swap all x and z tensor
           components, and after the code --i.e., at the end of this
           main loop-- we swap everything back.  */
        if (swap_xz) {
          /* Swap the x and z coordinates */
          SWAP (xx, zz);
        }

        CCTK_REAL r_plus
          = sqrt(pow(xx - par_b, 2) + pow(yy, 2) + pow(zz, 2));
        CCTK_REAL r_minus
          = sqrt(pow(xx + par_b, 2) + pow(yy, 2) + pow(zz, 2));

        CCTK_REAL U;
        switch (gsm)
        {
        case GSM_Taylor_expansion:
          U = PunctTaylorExpandAtArbitPosition
            (0, nvar, n1, n2, n3, v, xx, yy, zz);
          break;
        case GSM_evaluation:
          U = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, xx, yy, zz);
          break;
        default:
          assert (0);
        }
        r_plus = pow (pow (r_plus, 4) + pow (TP_epsilon, 4), 0.25);
        r_minus = pow (pow (r_minus, 4) + pow (TP_epsilon, 4), 0.25);
        if (r_plus < TP_Tiny)
            r_plus = TP_Tiny;
        if (r_minus < TP_Tiny)
            r_minus = TP_Tiny;
        CCTK_REAL psi1 = 1
          + 0.5 * *mp / r_plus
          + 0.5 * *mm / r_minus + U;
#define EXTEND(M,r) \
          ( M * (3./8 * pow(r, 4) / pow(TP_Extend_Radius, 5) - \
                 5./4 * pow(r, 2) / pow(TP_Extend_Radius, 3) + \
                 15./8 / TP_Extend_Radius))
        if (r_plus < TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mp,r_plus)
             + 0.5 * *mm / r_minus + U;
        }
        if (r_minus < TP_Extend_Radius) {
          psi1 = 1
             + 0.5 * EXTEND(*mm,r_minus)
             + 0.5 * *mp / r_plus + U;
        }
        CCTK_REAL static_psi = 1;
        
        CCTK_REAL Aij[3][3];
        BY_Aijofxyz (xx, yy, zz, Aij);

        CCTK_REAL old_alp=1.0;
        if (multiply_old_lapse)
            old_alp = alp[ind];

        if ((pmn_lapse) || (brownsville_lapse)) {

          CCTK_REAL xp, yp, zp, rp, ir;
          CCTK_REAL s1, s3, s5;
          CCTK_REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
          p = 1.0;
          px = py = pz = 0.0;
          pxx = pxy = pxz = 0.0;
          pyy = pyz = pzz = 0.0;

          /* first puncture */
          xp = xx - par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mp *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          /* second puncture */
          xp = xx + par_b;
          yp = yy;
          zp = zz;
          rp = sqrt (xp*xp + yp*yp + zp*zp);
          rp = pow (pow (rp, 4) + pow (TP_epsilon, 4), 0.25);
          if (rp < TP_Tiny)
              rp = TP_Tiny;
          ir = 1.0/rp;

          if (rp < TP_Extend_Radius) {
            ir = EXTEND(1., rp);
          }

          s1 = 0.5* *mm *ir;
          s3 = -s1*ir*ir;
          s5 = -3.0*s3*ir*ir;

          p += s1;

          px += xp*s3;
          py += yp*s3;
          pz += zp*s3;

          pxx += xp*xp*s5 + s3;
          pxy += xp*yp*s5;
          pxz += xp*zp*s5;
          pyy += yp*yp*s5 + s3;
          pyz += yp*zp*s5;
          pzz += zp*zp*s5 + s3;

          if (pmn_lapse)
            alp[ind] = pow(p, initial_lapse_psi_exponent);
          if (brownsville_lapse)
            alp[ind] = 2.0/(1.0+pow(p, initial_lapse_psi_exponent));

        } /* if brownsville_lapse */

        puncture_u[ind] = U;

        gxx[ind] = pow (psi1 / static_psi, 4);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow (psi1 / static_psi, 4);
        gyz[ind] = 0;
        gzz[ind] = pow (psi1 / static_psi, 4);

        kxx[ind] = Aij[0][0] / pow(psi1, 2);
        kxy[ind] = Aij[0][1] / pow(psi1, 2);
        kxz[ind] = Aij[0][2] / pow(psi1, 2);
        kyy[ind] = Aij[1][1] / pow(psi1, 2);
        kyz[ind] = Aij[1][2] / pow(psi1, 2);
        kzz[ind] = Aij[2][2] / pow(psi1, 2);

        if (antisymmetric_lapse || averaged_lapse) {
          alp[ind] =
            ((1.0 -0.5* *mp /r_plus -0.5* *mm/r_minus)
            /(1.0 +0.5* *mp /r_plus +0.5* *mm/r_minus));

          if (r_plus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(*mp, r_plus) -0.5* *mm/r_minus)
              /(1.0 +0.5*EXTEND(*mp, r_plus) +0.5* *mm/r_minus));
          }
          if (r_minus < TP_Extend_Radius) {
            alp[ind] =
              ((1.0 -0.5*EXTEND(*mm, r_minus) -0.5* *mp/r_plus)
              /(1.0 +0.5*EXTEND(*mp, r_minus) +0.5* *mp/r_plus));
          }
          
          if (averaged_lapse) {
            alp[ind] = 0.5 * (1.0 + alp[ind]);
          }
        }
        if (multiply_old_lapse)
          alp[ind] *= old_alp;

        if (swap_xz) {
          /* Swap the x and z components of all tensors */
          SWAP (gxx[ind], gzz[ind]);
          SWAP (gxy[ind], gyz[ind]);
          SWAP (kxx[ind], kzz[ind]);
          SWAP (kxy[ind], kyz[ind]);
        } /* if swap_xz */

      }// CCTK_ENDLOOP3_ALL(TwoPunctures);

  if (use_sources && rescale_sources)
  {
    assert(0);                  // TODO: Implement via critical region
#pragma omp single
    Rescale_Sources(cctkGH,
                    np,
                    vcoordx, vcoordy, vcoordz,
                    NULL,
                    gxx, gyy, gzz,
                    gxy, gxz, gyz);
  }
}

