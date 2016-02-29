program main

  use kinddefs,      only : dp
  use gridtools,     only : gridtype, preprocess_grid
  use solver,        only : iterate, get_soln, init_freestream
  use io_helpers,    only : write_tec_volume, write_tec_surface, read_namelist
  use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, nsteps,   &
                            nvars, tec_dataname, lin_solver, tolerance
  implicit none
  
  integer :: ndof
  
  real(dp), dimension(:),   allocatable :: Vx, Vy, Vt
  
  real(dp), dimension(:),   allocatable :: residual, phi
  
  type(gridtype) :: grid

  real(dp) :: dt

  integer :: i

continue

  call read_namelist
  
  call preprocess_grid(grid,gridfile,nnode)

! Number of degrees of freedom to solve in global system
  ndof = grid%npoin
  write(*,*) "CHECK: ndof = ",ndof
  
! Allocate the work arrays
  allocate(residual(ndof))
  allocate(phi(ndof))
  allocate(Vx(grid%npoin))
  allocate(Vy(grid%npoin))
  allocate(Vt(grid%npoin))

  phi = init_freestream(ndof,grid,uinf,vinf)

  call iterate(phi,grid,ndof,uinf,vinf)
 
  do i = 1, ndof
    write(*,*) "CHECK: phi = ",phi(i)
  end do

  call get_soln(Vx,Vy,Vt,phi,grid)
  
  call write_tec_volume(tec_dataname,grid,phi,Vx,Vy,Vt)
  call write_tec_surface(grid,Vt)

end program main
