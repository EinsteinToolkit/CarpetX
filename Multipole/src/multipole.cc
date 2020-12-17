#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "multipole.hh"
#include "interpolate.hh"
#include "utils.hh"
#include "sphericalharmonic.hh"

using namespace std;

static const int max_vars = 10;

// since each variable can have at most one spin weight, max_vars is a good
// upper limit for the number of spin weights to expect
static const int max_spin_weights = max_vars;

// modes are 0...max_l_modes-1
static const int max_l_modes = 10;
static const int max_m_modes = 2 * max_l_modes + 1;

typedef struct
{
  int n_vars;
  Multipole::variable_desc *vars;
}
variables_desc;

static void fill_variable(int idx, const char *optstring, void *callback_arg)
{
  assert(idx >= 0);
  assert(callback_arg != 0);

  variables_desc *vs = (variables_desc * ) callback_arg;

  assert(vs->n_vars < max_vars); // Too many variables in the variables list
  Multipole::variable_desc *v = &vs->vars[vs->n_vars];

  v->index = idx;

  // Default values if there is no option string or if the options are
  // not present
  v->imag_index = -1;
  v->spin_weight = 0;
  v->name = string(CCTK_VarName(v->index));

  if (optstring != 0)
  {
    int table = Util_TableCreateFromString(optstring);

    if (table >= 0)
    {
      const int buffer_length = 256;
      char buffer[buffer_length];
      
      Util_TableGetInt(table, &v->spin_weight , "sw");
      if (Util_TableGetString(table, buffer_length, buffer , "cmplx") >= 0)
      {
        v->imag_index = CCTK_VarIndex(buffer);
      }
      if (Util_TableGetString(table, buffer_length, buffer , "name") >= 0)
      {
        v->name = string(buffer);
      }

      const int ierr = Util_TableDestroy(table);
      if (ierr)
      {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not destroy table: %d", ierr);
      }
    }
  }
  vs->n_vars++;
}

static void parse_variables_string(const string &var_string, Multipole::variable_desc v[max_vars], int *n_variables)
{
  variables_desc  vars;

  vars.n_vars = 0;
  vars.vars = v;

  int ierr = CCTK_TraverseString(var_string.c_str(), fill_variable, &vars, CCTK_GROUP_OR_VAR);
  assert(ierr >= 0);

  *n_variables = vars.n_vars;
}

static void output_modes(CCTK_ARGUMENTS, const Multipole::variable_desc vars[], const CCTK_REAL radii[],
                     const Multipole::mode_array& modes)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_ascii)
  {
    if (CCTK_MyProc(cctkGH) == 0)
    {
      for (int v = 0; v < modes.get_nvars(); v++)
      {
        for(int i=0; i<modes.get_nradii(); i++)
        {
          const CCTK_REAL rad = radii[i];
          for (int l=0; l <= modes.get_lmax(); l++)
          {
            for (int m=-l; m <= l; m++)
            {
              ostringstream name;
              name << "mp_" << vars[v].name << "_l" << l << "_m" << m <<
                "_r" << setiosflags(ios::fixed) << setprecision(2) <<  rad << ".asc";
              Multipole_OutputComplexToFile(CCTK_PASS_CTOC, name.str(),
                                            modes(v, i, l, m, 0), modes(v, i, l, m, 1));
            }
          }
        }
      }
    }
  }
  if (output_hdf5)
  {
    if (CCTK_MyProc(cctkGH) == 0)
    {
      Multipole_OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
    }
  }
}

static void output_1D(CCTK_ARGUMENTS, const Multipole::variable_desc *v, CCTK_REAL rad,
                      CCTK_REAL *th, CCTK_REAL *ph,
                      CCTK_REAL *real, CCTK_REAL *imag, int array_size)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (CCTK_MyProc(cctkGH) == 0 && output_ascii)
  {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0)
    {
      ostringstream real_base;
      real_base << "mp_" << string(CCTK_VarName(v->index)) << "_r" << setiosflags(ios::fixed) << setprecision(2) << rad;
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str()+string(".th.asc"), array_size, th, ph, mp_theta, real);
      Multipole_Output1D(CCTK_PASS_CTOC, real_base.str()+string(".ph.asc"), array_size, th, ph, mp_phi, real);

      if (v->imag_index != -1)
      {
        ostringstream imag_base;
        imag_base << "mp_" << string(CCTK_VarName(v->imag_index)) << "_r" << setiosflags(ios::fixed) << setprecision(2) << rad;
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str()+string(".th.asc"), array_size, th, ph, mp_theta, imag);
        Multipole_Output1D(CCTK_PASS_CTOC, imag_base.str()+string(".ph.asc"), array_size, th, ph, mp_phi, imag);
      }
    }
  }
}

bool int_in_array(int a, const int array[], int len)
{
  for (int i = 0; i < len; i++)
  {
    if (array[i] == a)
      return true;
  }
  return false;
}

int find_int_in_array(int a, const int array[], int len)
{
  for (int i = 0; i < len; i++)
  {
    if (array[i] == a)
      return i;
  }
  return -1;
}


static
void get_spin_weights(Multipole::variable_desc vars[], int n_vars, int spin_weights[max_spin_weights], int *n_weights)
{
  int n_spin_weights = 0;

  for (int i = 0; i < n_vars; i++)
  {
    if (!int_in_array(vars[i].spin_weight, spin_weights, n_spin_weights))
    {
      assert(n_spin_weights < max_spin_weights);
      spin_weights[n_spin_weights] = vars[i].spin_weight;
      n_spin_weights++;
    }
  }
  *n_weights = n_spin_weights;
}

// For backward compatibility we allow the user to set l_mode instead
// of l_max, but if it is left at the default of -1, l_max is used.
static int get_l_max()
{
  DECLARE_CCTK_PARAMETERS;
  return l_mode == -1 ? l_max : l_mode;
}

static
void setup_harmonics(const int spin_weights[max_spin_weights],
                     int n_spin_weights, int lmax,
                     CCTK_REAL th[], CCTK_REAL ph[], int array_size,
                     CCTK_REAL *reY[max_spin_weights][max_l_modes][max_m_modes],
                     CCTK_REAL *imY[max_spin_weights][max_l_modes][max_m_modes])
{
  for (int si = 0; si < max_spin_weights; si++)
  {
    for (int l = 0; l < max_l_modes; l++)
    {
      for (int mi = 0; mi < max_m_modes; mi++)
      {
        reY[si][l][mi] = 0;
        imY[si][l][mi] = 0;
      }
    }
  }
  for (int si = 0; si < n_spin_weights; si++)
  {
    int sw = spin_weights[si];

    for (int l=0; l <= lmax; l++)
    {
      for (int m=-l; m <= l; m++)
      {
        reY[si][l][m+l] = new CCTK_REAL[array_size];
        imY[si][l][m+l] = new CCTK_REAL[array_size];

        Multipole_HarmonicSetup(sw, l, m, array_size, th, ph,
                                reY[si][l][m+l], imY[si][l][m+l]);
      }
    }
  } 
}

extern "C" void Multipole_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (l_mode != -1)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter l_mode is deprecated. Use l_max instead.  For compatibility, l_max = l_mode is being used.");
  }

  if (!CCTK_Equals(mode_type, "deprecated"))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter mode_type is deprecated and is no longer used.  All modes will be computed.");
  }

  if (l_min != -1)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter l_min is deprecated and is no longer used.  Modes from l = 0 will be computed.");
  }

  if (m_mode != -100)
  {
    CCTK_WARN(CCTK_WARN_ALERT, "The parameter m_mode is deprecated. All m modes will be computed.");
  }
}

extern "C" void Multipole_Calc(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  static CCTK_REAL *xs, *ys, *zs;
  static CCTK_REAL *xhat, *yhat, *zhat;
  static CCTK_REAL *th, *ph;
  static CCTK_REAL *real = 0, *imag = 0;
  static CCTK_REAL *reY[max_spin_weights][max_l_modes][max_m_modes];
  static CCTK_REAL *imY[max_spin_weights][max_l_modes][max_m_modes];
  static Multipole::variable_desc   vars[max_vars];
  static int n_variables = 0;
  static int spin_weights[max_spin_weights];
  static int n_spin_weights = 0;

  static bool initialized = false;

  const int array_size=(ntheta+1)*(nphi+1);

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  int lmax = get_l_max();

  assert(lmax + 1 <= max_l_modes);

  if (!initialized)
  {
    real = new CCTK_REAL[array_size];
    imag = new CCTK_REAL[array_size];
    th = new CCTK_REAL[array_size];  
    ph = new CCTK_REAL[array_size];
    xs = new CCTK_REAL[array_size];
    ys = new CCTK_REAL[array_size];
    zs = new CCTK_REAL[array_size];
    xhat = new CCTK_REAL[array_size];
    yhat = new CCTK_REAL[array_size];
    zhat = new CCTK_REAL[array_size];
  
    parse_variables_string(string(variables), vars, &n_variables);
    get_spin_weights(vars, n_variables, spin_weights, &n_spin_weights);
    Multipole_CoordSetup(xhat, yhat, zhat, th, ph);
    setup_harmonics(spin_weights, n_spin_weights, lmax, th, ph,
                    array_size, reY, imY);
    initialized = true;
  }

  Multipole::mode_array modes(n_variables, nradii, lmax);
  for (int v = 0; v < n_variables; v++)
  {
    //assert(vars[v].spin_weight == -2);

    int si = find_int_in_array(vars[v].spin_weight, spin_weights, n_spin_weights);
    assert(si != -1);

    for(int i=0; i<nradii; i++)
    {
      // Compute x^i = r * \hat x^i
      Multipole_ScaleCartesian(ntheta, nphi, radius[i], xhat, yhat, zhat, xs, ys, zs);
      
      // Interpolate Psi4r and Psi4i
      Multipole_Interp(CCTK_PASS_CTOC, xs, ys, zs, vars[v].index, vars[v].imag_index, 
                       real, imag);
      for (int l=0; l <= lmax; l++)
      {
        for (int m=-l; m <= l; m++)
        {
          // Integrate sYlm (real + i imag) over the sphere at radius r
          Multipole_Integrate(array_size, ntheta,
                              reY[si][l][m+l], imY[si][l][m+l],
                              real, imag, th, ph,
                              &modes(v, i, l, m, 0), &modes(v, i, l, m, 1));
    
        }//loop over m
      }//loop over l
      output_1D(CCTK_PASS_CTOC, &vars[v], radius[i], th, ph, real, imag, array_size);
    }//loop over radii
  }//loop over variables
  output_modes(CCTK_PASS_CTOC, vars, radius, modes);
}
