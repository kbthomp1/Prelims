module namelist_data

!============================================================================80
! This module is a common block intended to hold all data from the namelist
!============================================================================80

  use kinddefs, only : dp

  implicit none
  private

  public :: uinf, vinf, tolerance
  public :: nsteps, nnode, nvars, neqns
  public :: gridfile, tec_dataname, lin_solver

  real(dp) :: uinf
  real(dp) :: vinf
  real(dp) :: tolerance

  integer :: nsteps
  integer :: nnode
  integer :: nvars !Number of unknowns to solve
  integer :: neqns !Number of equations

  character(len=100) :: gridfile
  character(len=100) :: tec_dataname
  character(len=100) :: lin_solver

end module namelist_data
