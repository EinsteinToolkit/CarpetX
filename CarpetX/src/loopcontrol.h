#ifdef HAVE_CAPABILITY_CarpetX

#ifndef CARPETX_LOOPCONTROL_H
#define CARPETX_LOOPCONTROL_H

#include <cctk_Loop.h>

#ifdef __cplusplus
extern "C" {
#endif

// This file must mimic the GridDescBase in loop.hxx.

#define LC_DIM 3

typedef struct {
  int gsh[LC_DIM];
  int lbnd[LC_DIM], ubnd[LC_DIM];
  int lsh[LC_DIM];
  int ash[LC_DIM];
  int bbox[2 * LC_DIM];
  int nghostzones[LC_DIM];
  int tmin[LC_DIM], tmax[LC_DIM];
  CCTK_REAL x0[LC_DIM];
  CCTK_REAL dx[LC_DIM];
} GridDescBase_t;

GridDescBase_t LC_CreateGridDesc(const cGH *cctkGH);

// Replicate looping constructs in loop.hxx as good as we can, given that there
// is only one CCTK_LOOP macro that cannot distinguish between cell- and
// vertex-centred variables.
//
// This is not quite the way LoopControl does it, which would involve changing
// LC_LOOP3STROFF_NORMAL, which no longer knows what it needs to loop over and
// is only given the imin/imax ranges that Cactus thinks are required, but has
// no idea about tiles.

#define LC_LOOP3_BOX(name, i, j, k, grid, imin, imax, inormal)                 \
  do {                                                                         \
    for (int k = imin[2]; k < imax[2]; ++k) {                                  \
      for (int j = imin[1]; j < imax[1]; ++j) {                                \
        for (int i = imin[0]; i < imax[0]; ++i) {

#define LC_ENDLOOP3_BOX(name)                                                  \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  }                                                                            \
  while (0)

// Replace CCTK_LOOP macros
#if !defined CCTK_LOOP3STROFF_NORMAL || !defined CCTK_ENDLOOP3STROFF_NORMAL
#error "internal error"
#endif
#undef CCTK_LOOP3STROFF_NORMAL
#undef CCTK_ENDLOOP3STROFF_NORMAL

#undef CCTK_LOOP3_ALL
#undef CCTK_ENDLOOP3_ALL

#define CCTK_LOOP3_ALL(name, cctki3_cctkGH_, i, j, k)                          \
  do {                                                                         \
    const GridDescBase_t lc_grid = LC_CreateGridDesc(cctki3_cctkGH_);          \
    int imin[LC_DIM] CCTK_ATTRIBUTE_UNUSED,                                    \
        imax[LC_DIM] CCTK_ATTRIBUTE_UNUSED,                                    \
        inormal[LC_DIM] CCTK_ATTRIBUTE_UNUSED;                                 \
    const int lc_offset[LC_DIM] = {0, 0, 0};                                   \
    for (int d = 0; d < LC_DIM; ++d) {                                         \
      imin[d] = lc_grid.tmin[d];                                               \
      imax[d] = lc_grid.tmax[d];                                               \
      if (imax[d] > lc_grid.lsh[d] - lc_offset[d])                             \
        imax[d] = lc_grid.lsh[d] - lc_offset[d];                               \
      inormal[d] = 0;                                                          \
    }                                                                          \
    LC_LOOP3_BOX(name, i, j, k, lc_grid, imin, imax, inormal) {

#define CCTK_ENDLOOP3_ALL(name)                                                \
  }                                                                            \
  LC_ENDLOOP3_BOX(name);                                                       \
  }                                                                            \
  while (0)

#ifdef __cplusplus
}
#endif

#endif // CARPETX_LOOPCONTROL_H

#endif // HAVE_CAPABILITY_CarpetX
