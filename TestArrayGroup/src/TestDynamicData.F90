#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

subroutine TestArrayGroup_DynamicDataF(CCTK_ARGUMENTS)
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  ! Validate grid array dynamic data
  if(RANK(test1) /= 3) then ! note rank 3 b/c of vector of rank=2 arrays
      call CCTK_ERROR("incorrect dimension in test1 array dynamic data")
  endif
  if(SIZE(test1, 1) /= 5 .or. SIZE(test1, 2) /= 6 .or. SIZE(test1, 3) /= 4) then
      call CCTK_ERROR("incorrect size in test1 array dynamic data")
  endif

  if(RANK(test2) /= 3) then ! note rank 3 b/c of vector of rank=2 arrays
      call CCTK_ERROR("incorrect dimension in test2 array dynamic data")
  endif
  if(SIZE(test2, 1) /= 5 .or. SIZE(test2, 2) /= 6 .or. SIZE(test2, 3) /= 4) then
      call CCTK_ERROR("incorrect size in test2 array dynamic data")
  endif

  if(RANK(test3) /= 3) then ! note rank 3 b/c of vector of rank=2 arrays
      call CCTK_ERROR("incorrect dimension in test3 array dynamic data")
  endif
  if(SIZE(test3, 1) /= 5 .or. SIZE(test3, 2) /= 6 .or. SIZE(test3, 3) /= 4) then
      call CCTK_ERROR("incorrect size in test3 array dynamic data")
  endif

end subroutine TestArrayGroup_DynamicDataF
