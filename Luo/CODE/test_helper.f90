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

contains

  subroutine setup_cube
    integer :: nnode = 3
  continue

  ! Read the grid
    call preprocess_grid(grid,gridfile,nnode)

  end subroutine

end module test_helper
