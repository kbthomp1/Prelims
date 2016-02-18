program main

  use kinddefs,      only : dp
  use gridtools,     only : readgrid_lou, basis_function, face_norm
  use solver,        only : get_lhspo, get_rhspo, set_bc, get_soln, solve
  use io_helpers,    only : write_tec_volume, write_tec_surface, read_namelist
  use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, nsteps,   &
                            tec_dataname, lin_solver

  implicit none
  
  integer      :: ndimn, ntype, nelem, npoin, nface, nsteps, nnode
  
  real(dp)      :: uinf, vinf, tolerance
  
  real(dp), dimension(:),   allocatable :: rhspo, phi
  real(dp), dimension(:),   allocatable :: Vx, Vy, Vt
  
  real(dp), dimension(:,:), allocatable :: coord, geoel, lhspo, rface
  
  integer, dimension(:,:), allocatable :: inpoel, bcface
  
  character(len=100) :: gridfile, bc_case, tec_dataname, lin_solver
  namelist /fe_input/ gridfile, nnode, uinf, vinf, bc_case, nsteps,           &
                      tec_dataname, lin_solver

continue

  call read_namelist

! Read the namelist
  open(11,file="fe_input.nml",status="old")
  read(11,nml=fe_input)
  close(11)
  
! Read the grid
  call readgrid_lou(ndimn,ntype,nelem,npoin,nface,nnode,inpoel,gridfile,       &
                    coord,bcface)
  
! Allocate the work arrays
  allocate(lhspo(npoin,npoin))
  allocate(rhspo(npoin))
  allocate(phi(npoin))
  allocate(Vx(npoin))
  allocate(Vy(npoin))
  allocate(Vt(npoin))
  
! assemble geometry related matricies
  call basis_function(nelem,geoel,inpoel,coord)
  call face_norm(rface,coord,bcface,nface,ndimn,npoin)
  
! Formulate the load vector (RHS)
  rhspo = get_rhspo(bcface,rface,nface,npoin,uinf,vinf)
  
! Formulate the stiffness matrix (LHS)
  lhspo = get_lhspo(npoin,nelem,nnode,inpoel,geoel)
   
  call set_bc(phi,lhspo,rhspo,npoin,bcface)
  
  tolerance = 1.E-1_dp
  call solve(lin_solver,lhspo,rhspo,phi,npoin,nsteps,tolerance)
  
  call get_soln(Vx,Vy,Vt,phi,geoel,inpoel,npoin,nelem)
  
  call write_tec_volume(tec_dataname,npoin,nelem,coord,inpoel,phi,Vx,Vy,Vt)
  call write_tec_surface(npoin,bcface,coord,Vt)

end program main
