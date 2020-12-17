
#ifndef __integrate_h
#define __integrate_h
#include "cctk.h"

CCTK_REAL Midpoint2DIntegral(CCTK_REAL const *f, int nx, int ny,
                             CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL Trapezoidal2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, 
                            CCTK_REAL hx, CCTK_REAL hy);

CCTK_REAL DriscollHealy2DIntegral(CCTK_REAL const *f, int nx, int ny,
                                  CCTK_REAL hx, CCTK_REAL hy);

#endif
