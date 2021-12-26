// a bit of a hack since we must include fixmath.hxx before cctk.h but need to
// find out about capabilities to decide if fixmath (for CarpetX) is required
#include "cctk_Capabilities.h"
#ifdef HAVE_CAPABILITY_CarpetX
#include <fixmath.hxx>
#endif
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace SyncTestX {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

extern "C" void SyncTestX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_SyncTestX_Initialize;
  DECLARE_CCTK_PARAMETERS;

// Can Carpet and CarpetX coexist in a single compilation?
// If not, then the double if will be unnecessary. If so,
// Then the code will eventually need to be changed to
// allow for this.
#ifdef HAVE_CAPABILITY_Carpet
#ifdef HAVE_CAPABILITY_CarpetX
#error "Can't have Carpet and CarpetX"
#endif
  CCTK_LOOP3_ALL(SyncTest_Init,cctkGH,i,j,k) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    phi[index] = 0;
    t1[index] = 0;
    t2[index] = 0;
    t3[index] = 0;
    t4[index] = 0;
    t5[index] = 0;
  }
  CCTK_ENDLOOP3_ALL(SyncTest_Init);
#endif

#ifdef HAVE_CAPABILITY_CarpetX
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    phi[p.idx] = 0;
    t1[p.idx] = 0;
    t2[p.idx] = 0;
    t3[p.idx] = 0;
    t4[p.idx] = 0;
    t5[p.idx] = 0;
  });

  Loop::loop_bnd<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    phi[p.idx] = 0;
    t1[p.idx] = 0;
    t2[p.idx] = 0;
    t3[p.idx] = 0;
    t4[p.idx] = 0;
    t5[p.idx] = 0;
  });
#endif
}

extern "C" void SyncTestX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_SyncTestX_Evolve;
  DECLARE_CCTK_PARAMETERS;

#ifdef HAVE_CAPABILITY_Carpet
#ifdef HAVE_CAPABILITY_CarpetX
#error "Can't have Carpet and CarpetX"
#endif
  CCTK_LOOP3_ALL(SyncTest_Init,cctkGH,i,j,k) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    phi[index] = 0;
    t1[index] = 0;
    t2[index] = 0;
    t3[index] = 0;
    t4[index] = 0;
    t5[index] = 0;
  }
  CCTK_ENDLOOP3_ALL(SyncTest_Init);
#endif

#ifdef HAVE_CAPABILITY_CarpetX
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    phi[p.idx] = 0;
    t1[p.idx] = 0;
    t2[p.idx] = 0;
    t3[p.idx] = 0;
    t4[p.idx] = 0;
    t5[p.idx] = 0;
  });

  Loop::loop_bnd<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    phi[p.idx] = 0;
    t1[p.idx] = 0;
    t2[p.idx] = 0;
    t3[p.idx] = 0;
    t4[p.idx] = 0;
    t5[p.idx] = 0;
  });
#endif
}

extern "C" void SyncTestX_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

extern "C" void SyncTestX_Sync2(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace SyncTestX
