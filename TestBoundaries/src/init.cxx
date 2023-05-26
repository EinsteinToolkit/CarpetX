#include <cctk.h>
#include <cctk_Arguments.h>

#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#endif

namespace TestBoundaries {

int bbox_any[6];

extern "C" void TestBoundaries_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_Init;

  for (int d = 0; d < 6; ++d)
    bbox_any[d] = 0;
}

extern "C" void TestBoundaries_ReduceLocal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_ReduceLocal;

  for (int d = 0; d < 6; ++d)
#pragma omp atomic
    bbox_any[d] |= cctk_bbox[d];
}

extern "C" void TestBoundaries_ReduceGlobal(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_ReduceGlobal;

#ifdef HAVE_CAPABILITY_MPI
  MPI_Allreduce(MPI_IN_PLACE, bbox_any, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
#endif
}

extern "C" void TestBoundaries_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TestBoundaries_Check;

  CCTK_VINFO("bbox=[%d,%d,%d,%d,%d,%d]", bbox_any[0], bbox_any[1], bbox_any[2],
             bbox_any[3], bbox_any[4], bbox_any[5]);
  for (int d = 0; d < 6; ++d)
    if (!bbox_any[d])
      CCTK_VERROR("bbox[%d] is nowhere set", d);
}

} // namespace TestBoundaries
