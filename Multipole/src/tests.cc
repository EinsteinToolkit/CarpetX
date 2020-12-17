#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"
#include "integrate.hh"
#include "interpolate.hh"
#include "sphericalharmonic.hh"

static const int max_l_modes = 10;
static const int max_m_modes = 2 * max_l_modes + 1;

static CCTK_REAL test_integral(int n, CCTK_REAL (*integration_fn) (const CCTK_REAL *, int, int, CCTK_REAL, CCTK_REAL),const int is_midpoint)
{
  const int nx = n;
  const int ny = n;
  const int array_size=(nx+1)*(ny+1);

  CCTK_REAL *f = new CCTK_REAL[array_size];

  const CCTK_REAL dx = 1./(nx + is_midpoint);
  const CCTK_REAL dy = 1./(ny + is_midpoint);

  for (int ix = 0; ix <= nx; ix++)
  {
    for (int iy = 0; iy <= ny; iy++)
    {
      const int i = Multipole_Index(ix, iy, nx);

      const CCTK_REAL x = ix*dx + 0.5*dx*is_midpoint;
      const CCTK_REAL y = iy*dy + 0.5*dy*is_midpoint;
      const CCTK_REAL PI = acos(-1.0);

      f[i] = x*pow(y,2)*pow(cos(2*PI*y),2)*pow(sin(2*PI*x),2);
    }
  }

  const CCTK_REAL result = integration_fn(f, nx, ny, dx, dy);
  delete [] f;
  return result;
}

static CCTK_REAL test_pi_symmetric_sphere_integral(CCTK_REAL (*integration_fn) (const CCTK_REAL *, int, int, CCTK_REAL, CCTK_REAL),const int is_midpoint)
{
  const int n = 100;
  const int nth = n;
  const int nph = n;
  const int array_size=(nth+1)*(nph+1);
  const CCTK_REAL PI = acos(-1.0);

  CCTK_REAL *f = new CCTK_REAL[array_size];

  const CCTK_REAL dth = PI/(nth + is_midpoint);
  const CCTK_REAL dph = 2*PI/(nph + is_midpoint);

  for (int ith = 0; ith <= nth; ith++)
  {
    for (int iph = 0; iph <= nph; iph++)
    {
      const int i = Multipole_Index(ith, iph, nth);

      const CCTK_REAL th = ith*dth + 0.5*dth*is_midpoint;
      const CCTK_REAL ph = iph*dph + 0.5*dph*is_midpoint;

      f[i] = -(cos(ph)*sqrt(5/PI)*pow(cos(th/2.),3)*sin(th/2.));
    }
  }

  const CCTK_REAL result = integration_fn(f, nth, nph, dth, dph);
  delete [] f;
  return result;
}

CCTK_REAL integration_convergence_order(CCTK_REAL (*integration_fn) (const CCTK_REAL *, int, int, CCTK_REAL, CCTK_REAL), CCTK_REAL *store_low, CCTK_REAL *store_high, const int is_midpoint)
{
  const int n1 = 100;
  const int n2 = 200;
  const CCTK_REAL PI = acos(-1.0);
  const CCTK_REAL result1 = test_integral(100, integration_fn, is_midpoint);
  *store_low = result1;
  const CCTK_REAL result2 = test_integral(200, integration_fn, is_midpoint);
  *store_high = result2;
  const CCTK_REAL exact = 1./24 + 1./(64 * pow(PI,2));
  const CCTK_REAL error1 = fabs(result1 - exact);
  const CCTK_REAL error2 = fabs(result2 - exact);
  return log10(error1/error2) / log10((CCTK_REAL) n2/n1);
}

void Multipole_TestIntegrationConvergence(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  *test_simpson_convergence_order = integration_convergence_order(
    &Simpson2DIntegral, test_simpson_result_low, test_simpson_result_high,0);
  *test_trapezoidal_convergence_order = integration_convergence_order(
    &Trapezoidal2DIntegral, test_trapezoidal_result_low, test_trapezoidal_result_high,0);
  *test_midpoint_convergence_order = integration_convergence_order(
    &Midpoint2DIntegral, test_midpoint_result_low, test_midpoint_result_high,1);
}

void Multipole_TestIntegrationSymmetry(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  *test_simpson_pi_symmetry = test_pi_symmetric_sphere_integral(&Simpson2DIntegral,0);
  *test_midpoint_pi_symmetry = test_pi_symmetric_sphere_integral(&Midpoint2DIntegral,1);
  *test_trapezoidal_pi_symmetry = test_pi_symmetric_sphere_integral(&Trapezoidal2DIntegral,0);
  *test_driscollhealy_pi_symmetry = test_pi_symmetric_sphere_integral(&DriscollHealy2DIntegral,0);
  printf("Pi symmetry Simpson integral: %.19g\n", *test_simpson_pi_symmetry);
  printf("Pi symmetry midpoint integral: %.19g\n", *test_midpoint_pi_symmetry);
  printf("Pi symmetry trapezoidal integral: %.19g\n", *test_trapezoidal_pi_symmetry);
  printf("Pi symmetry Driscoll and Healy integral: %.19g\n", *test_driscollhealy_pi_symmetry);
}

// void Multipole_TestIntegrate(CCTK_ARGUMENTS)
// {
//   const int n = 100;

//   const int nth = n;
//   const int nph = n;
//   const int array_size=(nth+1)*(nph+1);

//   CCTK_REAL *f = new CCTK_REAL[array_size];

//   const CCTK_REAL dth = 1/nx;
//   const CCTK_REAL dph = 1/ny;

//   for (int ix = 0; ix <= nx; ix++)
//   {
//     for (int iy = 0; iy <= ny; iy++)
//     {
//       const int i = Multipole_Index(ix, iy, nx);

//       const CCTK_REAL x = ix*dx;
//       const CCTK_REAL y = iy*dy;

//       f[i] = sin(2*PI*x)*cos(2*PI*y);
//     }
//   }




//   CCTK_REAL result = Multipole_Integrate(int array_size, int nthetap,
//                                          CCTK_REAL const array1r[], CCTK_REAL const array1i[],
//                                          CCTK_REAL const array2r[], CCTK_REAL const array2i[],
//                                          CCTK_REAL const th[], CCTK_REAL const ph[], 
//                                          CCTK_REAL *outre, CCTK_REAL *outim)



//   printf("Integration result: %.19g\n", result);
// }

void Multipole_TestOrthonormality(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Campute Cartesian coordinates of points on the sphere */
  int array_size=(ntheta+1)*(nphi+1);

  CCTK_REAL *th = new CCTK_REAL[array_size];
  CCTK_REAL *ph = new CCTK_REAL[array_size];
  CCTK_REAL *xhat = new CCTK_REAL[array_size];
  CCTK_REAL *yhat = new CCTK_REAL[array_size];
  CCTK_REAL *zhat = new CCTK_REAL[array_size];

  Multipole_CoordSetup(xhat, yhat, zhat, th, ph);

  /* Populate spherical-harmonic array */
  CCTK_REAL *reY[1][max_l_modes][max_m_modes];
  CCTK_REAL *imY[1][max_l_modes][max_m_modes];

  for (int sw = 0; sw <= 0; sw++)
  {
    for (int l = 0; l < max_l_modes; l++)
    {
      for (int m = -l; m <= l; m++)
      {
        reY[sw][l][m+l] = new CCTK_REAL[array_size];
        imY[sw][l][m+l] = new CCTK_REAL[array_size];

        Multipole_HarmonicSetup(sw, l, m, array_size, th, ph,
                                reY[sw][l][m+l], imY[sw][l][m+l]);
      }
    }
  }

  /* Loop over l and m, assign Ylm to (rel,imag), and compute the scalar
     product with all spherical harmonics (loop over li, mi) */
  // [0..max_l_modes) has N=max_l_modes^2
  // comparing each mode with each other but skipping the duplicates
  // gives N*(N+1)/2
  // only 1 spin-weight mode is tested
  const int N = max_l_modes*max_l_modes;
  int idx = 0;
  for (int sw = 0; sw <= 0; sw++)
  {
    for (int l = 0; l < max_l_modes; l++)
    {
      for (int m = -l; m <= l; m++)
      {

        /* Compute scalar product of (real,imag) and all the Ylimi */
        for (int li = 0; li < max_l_modes; li++)
        {
          for (int mi = -li; mi <= li; mi++)
          {
            // only handle lower triangle in ((l,m),(li,mi)) space
            if(l*l+l+m < li*li+li+mi)
              continue;

            CCTK_REAL real_lm = 0.0, imag_lm = 0.0;
            Multipole_Integrate(array_size, ntheta,
                                reY[sw][li][mi+li], imY[sw][li][mi+li],
                                reY[sw][l][m+l], imY[sw][l][m+l], th, ph,
                                &real_lm, &imag_lm);

            assert(idx < 1*N*(N+1)/2);
            test_orthonormality[idx++] = sqrt(real_lm*real_lm + imag_lm*imag_lm);
          }
        }

      }
    }
  }
  assert(idx == 1*N*(N+1)/2);

  for (int sw = 0; sw <= 0; sw++)
  {
    for (int l = 0; l < max_l_modes; l++)
    {
      for (int m = -l; m <= l; m++)
      {
        delete[] imY[0][l][m+l];
        delete[] reY[0][l][m+l];
      }
    }
  }
  delete[] zhat;
  delete[] yhat;
  delete[] xhat;
  delete[] ph;
  delete[] th;

  return;
}
