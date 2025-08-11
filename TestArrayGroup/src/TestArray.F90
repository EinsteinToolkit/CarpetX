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

  integer :: array_real_error_count(3)
  integer :: array_int_error_count(3)
  integer :: array_complex_error_count(3)

  integer :: i,j,n
  character*256 :: warnline

  array_real_error_count = (/0,0,0/)
  array_int_error_count = (/0,0,0/)
  array_complex_error_count = (/0,0,0/)
  do n = 1,nmax
    do j = 1,jmax
      do i = 1,imax
        ! -1 is due to data being filled by C code counting from 0
        if (test1(i,j,n) /= 1 + i * j * n) then
          array_real_error_count(1) = array_real_error_count(1) + 1
        end if
        if (test_int1(i,j,n) /= 2 + i * j * n) then
          array_int_error_count(1) = array_int_error_count(1) + 1
        end if
        if (test_complex1(i,j,n) /= CMPLX(3 + i * j * n, 1, KIND(test_complex1(1,1,1)))) then
          array_complex_error_count(1) = array_complex_error_count(1) + 1
        end if
        if (test2(i,j,n) /= 1 + 7 * i * j * n) then
          array_real_error_count(2) = array_real_error_count(2) + 1
        end if
        if (test_int2(i,j,n) /= 2 + 7 * i * j * n) then
          array_int_error_count(2) = array_int_error_count(2) + 1
        end if
        if (test_complex2(i,j,n) /= CMPLX(3 + 7 * i * j * n, 1, KIND(test_complex2(1,1,1)))) then
          array_complex_error_count(2) = array_complex_error_count(2) + 1
        end if
        if (test3(i,j,n) /= 1 + 13 * i * j * n) then
          array_real_error_count(3) = array_real_error_count(3) + 1
        end if
        if (test_int3(i,j,n) /= 2 + 13 * i * j * n) then
          array_int_error_count(3) = array_int_error_count(3) + 1
        end if
        if (test_complex3(i,j,n) /= CMPLX(3 + 13 * i * j * n, 1, KIND(test_complex3(1,1,1)))) then
          array_complex_error_count(3) = array_complex_error_count(3) + 1
        end if
      end do 
    end do
  end do

  do i = 1,3
    if (array_real_error_count(i) > 0) then
      write(warnline, "('TestArrayGroup: Fortran CCTK_REAL grid array test', i1, ' failed in ', i3, ' of ', i3, ' elements')") i, array_real_error_count(i), size
      call CCTK_WARN(CCTK_WARN_ALERT, warnline)
      call TestArrayGroup_FoundError()
    end if
  end do
  do i = 1,3
    if (array_int_error_count(i) > 0) then
      write(warnline, "('TestArrayGroup: Fortran CCTK_INT grid array test', i1, ' failed in ', i3, ' of ', i3, ' elements')") i, array_int_error_count(i), size
      call CCTK_WARN(CCTK_WARN_ALERT, warnline)
      call TestArrayGroup_FoundError()
    end if
  end do
  do i = 1,3
    if (array_complex_error_count(i) > 0) then
      write(warnline, "('TestArrayGroup: Fortran CCTK_COMPLEX grid array test', i1, ' failed in ', i3, ' of ', i3, ' elements')") i, array_complex_error_count(i), size
      call CCTK_WARN(CCTK_WARN_ALERT, warnline)
      call TestArrayGroup_FoundError()
    end if
  end do

  ! test grid scalar
  if(test_scalar /= 1) then
    call CCTK_WARN(CCTK_WARN_ALERT, "TestArrayGroup: CCTK_REAL grid scalar failed")
    call TestArrayGroup_FoundError()
  end if
  if(test_scalar_int /= 2) then
    call CCTK_WARN(CCTK_WARN_ALERT, "TestArrayGroup: CCTK_INT grid scalar failed")
    call TestArrayGroup_FoundError()
  end if
  if(test_scalar_complex /= 3) then
    call CCTK_WARN(CCTK_WARN_ALERT, "TestArrayGroup: CCTK_COMPLEX grid scalar failed")
    call TestArrayGroup_FoundError()
  end if
end subroutine TestArrayGroup_CompareF
