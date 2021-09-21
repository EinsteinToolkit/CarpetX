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

#  CCTK_LOOP3_ALL(SyncTest_Init,cctkGH,i,j,k) {
#    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
#    phi[index] = 0;
#    t1[index] = 0;
#    t2[index] = 0;
#    t3[index] = 0;
#    t4[index] = 0;
#    t5[index] = 0;
#  }
#  CCTK_ENDLOOP3_ALL(SyncTest_Init);

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
}

extern "C" void SyncTestX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_SyncTestX_Evolve;
  DECLARE_CCTK_PARAMETERS;

#  CCTK_LOOP3_ALL(SyncTest_Init,cctkGH,i,j,k) {
#    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
#    phi[index] = 0;
#    t1[index] = 0;
#    t2[index] = 0;
#    t3[index] = 0;
#    t4[index] = 0;
#    t5[index] = 0;
#  }
#  CCTK_ENDLOOP3_ALL(SyncTest_Init);

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
}

extern "C" void SyncTestX_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

extern "C" void SyncTestX_Sync2(CCTK_ARGUMENTS) {
  // do nothing
}

} // namespace SyncTestX
