program main

  use kinddefs,      only : dp
  use gridtools,     only : gridtype, preprocess_grid
  use solver,        only : iterate, get_soln, init_freestream
  use io_helpers,    only : write_tec_volume, write_tec_surface, read_namelist
  use namelist_data, only : gridfile, nnode, tec_dataname, uinf, vinf, tolerance

  implicit none
  
  integer :: ndof
  
  real(dp), dimension(:),   allocatable :: Vx, Vy, Vt
  
  real(dp), dimension(:),   allocatable :: residual, phi, nodal_phi
  
  type(gridtype) :: grid

continue

  call read_namelist
  
  call preprocess_grid(grid,gridfile,nnode)

! Number of degrees of freedom to solve in global system
  ndof = grid%nnode*grid%nelem
  
! Allocate the work arrays
  allocate(residual(ndof))
  allocate(phi(ndof))
  allocate(nodal_phi(grid%npoin))
  allocate(Vx(grid%npoin))
  allocate(Vy(grid%npoin))
  allocate(Vt(grid%npoin))

  !phi = init_freestream(ndof,grid,uinf,vinf)
  phi = 0._dp

  call iterate(phi,grid,ndof,tolerance)

  call get_soln(Vx,Vy,Vt,phi,nodal_phi,grid)
  
  call write_tec_volume(tec_dataname,grid,phi,Vx,Vy,Vt)
  call write_tec_surface(grid,Vt)

end program main
