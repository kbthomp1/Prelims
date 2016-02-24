program main

  use kinddefs,      only : dp
  use gridtools,     only : readgrid_lou, basis_function, face_norm
  use solver,        only : get_lhspo, get_rhspo, get_soln, solve
  use io_helpers,    only : write_tec_volume, write_tec_surface, read_namelist
  use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, nsteps,   &
                            nvars, neqns, tec_dataname, lin_solver, tolerance
  implicit none
  
  integer :: ndimn, ntype, nelem, npoin, nface
  
  real(dp), dimension(:),   allocatable :: Vx, Vy, Vt
  
  real(dp), dimension(:,:), allocatable :: coord, geoel, rface
  real(dp), dimension(:),   allocatable :: rhspo, phi
  real(dp), dimension(:),   allocatable :: lhspo
  
  integer, dimension(:,:),  allocatable :: inpoel, bcface

continue

  call read_namelist
  
! Read the grid
  call readgrid_lou(ndimn,ntype,nelem,npoin,nface,nnode,inpoel,gridfile,       &
                    coord,bcface)
  
! Allocate the work arrays
  allocate(lhspo(nvars,neqns,npoin,npoin))
  allocate(rhspo(nvars,npoin))
  allocate(phi(nvars,npoin))
  allocate(Vx(npoin))
  allocate(Vy(npoin))
  allocate(Vt(npoin))
  
! assemble geometry related matricies
  call basis_function(nelem,geoel,inpoel,coord)
  call face_norm(rface,coord,bcface,nface,ndimn)
  
! Formulate the load vector (RHS)
  rhspo = get_rhspo(bcface,rface,nface,npoin,neqns,uinf,vinf)
  
! Formulate the stiffness matrix (LHS)
  lhspo = get_lhspo(npoin,nelem,nnode,nvars,neqns,inpoel,geoel)
   
  call solve(lin_solver,lhspo,rhspo,phi,npoin,nvars,neqns,nsteps,tolerance)
  
  call get_soln(Vx,Vy,Vt,phi,geoel,inpoel,npoin,nelem)
  
  call write_tec_volume(tec_dataname,npoin,nelem,coord,inpoel,phi,Vx,Vy,Vt)
  call write_tec_surface(nface,bcface,coord,Vt)

end program main
