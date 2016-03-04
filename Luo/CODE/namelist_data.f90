module namelist_data

!============================================================================80
! This module is a common block intended to hold all data from the namelist
!============================================================================80

  use kinddefs, only : dp

  implicit none
  private

  public :: uinf, vinf, tolerance
  public :: relative_tol, absolute_tol
  public :: nsteps, nnode, dt, cfl, rk_order
  public :: gridfile, tec_dataname, lin_solver
  public :: restart_file, read_restart, output_freq
  public :: read_tec_restart, tec_output_freq
  public :: write_restart_freq

  real(dp) :: uinf
  real(dp) :: vinf
  real(dp) :: tolerance
  real(dp) :: relative_tol
  real(dp) :: absolute_tol
  real(dp) :: dt
  real(dp) :: cfl

  integer :: nsteps
  integer :: nnode
  integer :: rk_order
  integer :: output_freq
  integer :: tec_output_freq
  integer :: write_restart_freq

  character(len=100) :: gridfile
  character(len=100) :: tec_dataname
  character(len=100) :: lin_solver
  character(len=100) :: restart_file

  logical :: read_restart
  logical :: read_tec_restart

end module namelist_data
