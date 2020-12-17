#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

// Multipole_Interp:
//      This function interpolates psi4 onto the sphere in cartesian 
//      coordinates as created by Multipole_CoordSetup.
void Multipole_Interp(CCTK_ARGUMENTS,
                      CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[], 
                      int real_idx, int imag_idx,
                      CCTK_REAL psi4r[], CCTK_REAL psi4i[]);

void Multipole_InterpVar(CCTK_ARGUMENTS,
                         CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[], const char *var_name,
                         CCTK_REAL sphere_var[]);
