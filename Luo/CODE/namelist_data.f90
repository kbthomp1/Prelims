module namelist_data

!============================================================================80
! This module is a common block intended to hold all data from the namelist
!============================================================================80

  use kinddefs, only : dp

  implicit none
  private

  public :: uinf, vinf, tolerance
  public :: nsteps, nnode, dt, cfl, rk_order
  public :: gridfile, tec_dataname, lin_solver
  public :: restart_file, read_restart

  real(dp) :: uinf
  real(dp) :: vinf
  real(dp) :: tolerance
  real(dp) :: dt
  real(dp) :: cfl

  integer :: nsteps
  integer :: nnode
  integer :: rk_order

  character(len=100) :: gridfile
  character(len=100) :: tec_dataname
  character(len=100) :: lin_solver
  character(len=100) :: restart_file

  logical :: read_restart

end module namelist_data
