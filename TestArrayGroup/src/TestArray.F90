#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>

subroutine TestArrayGroup_CompareF(CCTK_ARGUMENTS)
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

  !// test array group
  integer, parameter :: imax = 5
  integer, parameter :: jmax = 6
  integer, parameter :: nmax = 4
  integer, parameter :: size = imax * jmax * nmax

  integer :: array_error_count(3)
  integer :: scalar_error_count

  integer :: i,j,n
  character*256 :: warnline

  array_error_count = (/0,0,0/)
  do n = 1,nmax
    do j = 1,jmax
      do i = 1,imax
        ! -1 is due to data being filled by C code counting from 0
        if (test1(i,j,n) /= 1 + (i-1) * (j-1) * (n-1)) then
          array_error_count(1) = array_error_count(1) + 1
        end if
        if (test2(i,j,n) /= 1 + 7 * (i-1) * (j-1) * (n-1)) then
          array_error_count(2) = array_error_count(2) + 1
        end if
        if (test3(i,j,n) /= 1 + 13 * (i-1) * (j-1) * (n-1)) then
          array_error_count(3) = array_error_count(3) + 1
        end if
      end do 
    end do
  end do

  if (array_error_count(1) > 0) then
    write(warnline, "('TestArrayGroup: Fortran grid array test1 failed in ', i3, ' of ', i3, ' elements')") array_error_count(1), size
    call CCTK_ERROR(warnline)
  end if

  if (array_error_count(2) > 0) then
    write(warnline, "('TestArrayGroup: Fortran grid array test2 failed in ', i3, ' of ', i3, ' elements')") array_error_count(2), size
    call CCTK_ERROR(warnline)
  end if

  if (array_error_count(3) > 0) then
    write(warnline, "('TestArrayGroup: Fortran grid array test3 failed in ', i3, ' of ', i3, ' elements')") array_error_count(3), size
    call CCTK_ERROR(warnline)
  end if

  ! test grid scalar
  scalar_error_count = 0

  if(test_scalar /= 1) then
    scalar_error_count = scalar_error_count + 1
  end if

  if (scalar_error_count > 0) then
    write(warnline, "('TestArrayGroup: Fortran grid scalr test_scalar failed in ', i3, ' of ', i3, ' elements')") scalar_error_count, 1
    call CCTK_ERROR(warnline)
  end if
end subroutine TestArrayGroup_CompareF
