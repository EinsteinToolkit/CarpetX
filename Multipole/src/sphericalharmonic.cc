#include <math.h>
#include <assert.h>
#include <iostream>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

static const CCTK_REAL PI = acos(-1.0);

static double factorial(int n)
{
  double returnval = 1;
  for (int i = n; i >= 1; i--)
  {
    returnval *= i;
  }
  return returnval;
}

static inline double combination(int n, int m)
{
  // Binomial coefficient is undefined if these conditions do not hold
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);

  return factorial(n) / (factorial(m) * factorial(n-m));
}

static inline int imin(int a, int b)
{
  return a < b ? a : b;
}

static inline int imax(int a, int b)
{
  return a > b ? a : b;
}

void Multipole_SphericalHarmonic(int s, int l, int m, 
                                 CCTK_REAL th, CCTK_REAL ph,
                                 CCTK_REAL *reY, CCTK_REAL *imY)
{
//  assert(s == -2 && l == 2 && m == 2);
//  *reY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * cos(2*ph);
//  *imY = 1.0/2.0 * sqrt(5/PI) * pow(cos(th/2), 4) * sin(2*ph);
  double all_coeff = 0, sum = 0;
  all_coeff = pow(-1.0, m);
  all_coeff *= sqrt(factorial(l+m)*factorial(l-m)*(2*l+1) / (4.*PI*factorial(l+s)*factorial(l-s)));
  sum = 0.;
  for (int i = imax(m - s, 0); i <= imin(l + m, l - s); i++)
  {
    double sum_coeff = combination(l-s, i) * combination(l+s, i+s-m);
    sum += sum_coeff * pow(-1.0, l-i-s) * pow(cos(th/2.), 2 * i + s - m) * 
      pow(sin(th/2.), 2*(l-i)+m-s);
  }
  *reY = all_coeff*sum*cos(m*ph);
  *imY = all_coeff*sum*sin(m*ph);
}

void Multipole_HarmonicSetup(int s, int l, int m,
                             int array_size, CCTK_REAL const th[], CCTK_REAL const ph[],
                             CCTK_REAL reY[], CCTK_REAL imY[])
{
  for (int i = 0; i < array_size; i++)
  {
    Multipole_SphericalHarmonic(s,l,m,th[i],ph[i],&reY[i], &imY[i]);
  }
}


// Fill a grid function with a given spherical harmonic
extern "C" void Multipole_SetHarmonic(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_Multipole_SetHarmonic;
  DECLARE_CCTK_PARAMETERS;

  for (int k = 0; k < cctk_lsh[2]; k++)
  {
    for (int j = 0; j < cctk_lsh[1]; j++)
    {
      for (int i = 0; i < cctk_lsh[0]; i++)
      {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k) ;

        CCTK_REAL theta = acos(z[index]/r[index]);
        CCTK_REAL phi = atan2(y[index],x[index]);

        CCTK_REAL re = 0;
        CCTK_REAL im = 0;

        Multipole_SphericalHarmonic(test_sw,test_l,test_m,theta,phi,
                                    &re, &im);

        CCTK_REAL fac = test_mode_proportional_to_r ? r[index] : 1.0;

        harmonic_re[index] = re * fac;
        harmonic_im[index] = im * fac;
      }
    }
  }
  return;
}
