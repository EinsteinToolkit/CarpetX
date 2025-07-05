#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

subroutine TestArrayGroup_DynamicDataF(CCTK_ARGUMENTS)
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  integer, dimension(3) :: sizes1, sizes2, sizes3

  ! check that grid variable is of rank 3. This fails to compiler otherwise.
  sizes1 = SHAPE(test1)
  if(sizes1(1) /= 5 .or. sizes1(2) /= 6 .or. sizes1(3) /= 4) then
      call CCTK_ERROR("incorrect size in test1 array dynamic data")
  endif

  sizes2 = SHAPE(test2)
  if(sizes2(1) /= 5 .or. sizes2(2) /= 6 .or. sizes2(3) /= 4) then
      call CCTK_ERROR("incorrect size in test2 array dynamic data")
  endif

  sizes3 = SHAPE(test3)
  if(sizes3(1) /= 5 .or. sizes3(2) /= 6 .or. sizes3(3) /= 4) then
      call CCTK_ERROR("incorrect size in test3 array dynamic data")
  endif

end subroutine TestArrayGroup_DynamicDataF
