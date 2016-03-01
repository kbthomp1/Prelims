module test_helper

  use kinddefs,  only : dp
  use gridtools, only : gridtype, preprocess_grid

  implicit none
  public

  type(gridtype) :: grid

  character(len=20) :: gridfile = "testgrid"

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp

  integer :: ndof

  interface wrt
    module procedure :: print_vector
    module procedure :: print_matrix
  end interface wrt

contains

  subroutine setup_cube
    integer :: nnode = 3
  continue

  ! Read the grid
    call preprocess_grid(grid,gridfile,nnode)

  end subroutine

  subroutine print_vector(v,title)
    real(dp), dimension(:), intent(in) :: v
    character(len=*), optional, intent(in) :: title
    integer :: d1, i
  continue
    d1 = size(v,1)
    if (present(title)) then
      do i = 1, d1
        write(*,30) trim(adjustl(title)), v(i)
      end do
    else
      do i = 1, d1
        write(*,30) "V", v(i)
      end do
    end if
30 format(A,": ",g11.4)
  end subroutine print_vector

  subroutine print_matrix(m,title)
    real(dp), dimension(:,:), intent(in) :: m
    character(len=*), optional, intent(in) :: title
    integer :: d1, d2, i, j
  continue
    d1 = size(m,1); d2 = size(m,2)
    if (present(title)) then
      do i = 1, d1
        write(*,20) trim(adjustl(title)), (m(i,j),j=1,d2)
      end do
    else
      do i = 1, d1
        write(*,20) "M", (m(i,j),j=1,d2)
      end do
    end if
20 format(A,": ",300g11.4)
  end subroutine print_matrix


end module test_helper
