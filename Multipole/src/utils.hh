
#ifndef __utils_h
#define __utils_h

#include "cctk.h"
#include <string>

using namespace std;

enum mp_coord {mp_theta, mp_phi};

namespace Multipole {
  class mode_array;
  struct variable_desc;
}

void Multipole_OutputArrayToFile(CCTK_ARGUMENTS, const string &name, int array_size,
                                 CCTK_REAL const th[], CCTK_REAL const ph[],
                                 CCTK_REAL const x[], CCTK_REAL const y[], CCTK_REAL const z[],
                                 CCTK_REAL const data[]);

void Multipole_Output1D(CCTK_ARGUMENTS, const string &name, int array_size,
                        CCTK_REAL const th[], CCTK_REAL const ph[], mp_coord coord,
                        CCTK_REAL const data[]);

void Multipole_OutputComplexToFile(CCTK_ARGUMENTS, const string &name, CCTK_REAL redata, CCTK_REAL imdata);

void Multipole_OutputComplexToH5File(CCTK_ARGUMENTS, const Multipole::variable_desc vars[], const CCTK_REAL radii[],
                                     const Multipole::mode_array& modes);

void Multipole_CoordSetup(CCTK_REAL xhat[], CCTK_REAL yhat[],
                          CCTK_REAL zhat[], CCTK_REAL th[], 
                          CCTK_REAL ph[]);

void Multipole_ScaleCartesian(int ntheta, int nphi, CCTK_REAL r,
                              CCTK_REAL const xhat[], CCTK_REAL const yhat[], CCTK_REAL const zhat[],
                              CCTK_REAL x[], CCTK_REAL y[], CCTK_REAL z[]);

static inline int Multipole_Index(int it, int ip, int ntheta)
{
  return it + (ntheta+1)*ip;
}

void Multipole_Integrate(int array_size, int ntheta,
    CCTK_REAL const array1r[], CCTK_REAL const array1i[],
    CCTK_REAL const array2r[], CCTK_REAL const array2i[],
    CCTK_REAL const th[], CCTK_REAL const pph[], 
    CCTK_REAL out_arrayr[], CCTK_REAL out_arrayi[]);

#endif
