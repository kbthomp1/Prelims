program main

  use kinddefs,      only : dp
  use gridtools,     only : gridtype, preprocess_grid
  use solver,        only : get_lhspo, get_rhspo, get_soln, solve
  use io_helpers,    only : write_tec_volume, write_tec_surface, read_namelist
  use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, nsteps,   &
                            nvars, tec_dataname, lin_solver, tolerance
  implicit none
  
  integer :: ndof
  
  real(dp), dimension(:),   allocatable :: Vx, Vy, Vt
  
  real(dp), dimension(:),   allocatable :: rhspo, phi
  real(dp), dimension(:,:), allocatable :: lhspo
  
  type(gridtype) :: grid

  integer :: i

continue

  call read_namelist
  
  call preprocess_grid(grid,gridfile,nnode)

! Number of degrees of freedom to solve
  ndof = grid%npoin + grid%numfac
  
! Allocate the work arrays
  allocate(lhspo(ndof,ndof))
  allocate(rhspo(ndof))
  allocate(phi(ndof))
  allocate(Vx(grid%npoin))
  allocate(Vy(grid%npoin))
  allocate(Vt(grid%npoin))
  
! Formulate the load vector (RHS)
  rhspo = get_rhspo(grid,ndof,uinf,vinf)
  
! Formulate the stiffness matrix (LHS)
  lhspo = get_lhspo(grid,ndof)
   
  call solve(lin_solver,lhspo,rhspo,phi,ndof,nsteps,tolerance)

  do i = 1, ndof
    write(*,*) "CHECK: phi = ",phi(i)
  end do

  call get_soln(Vx,Vy,Vt,phi,grid)
  
  call write_tec_volume(tec_dataname,grid,phi,Vx,Vy,Vt)
  call write_tec_surface(grid,Vt)

end program main
