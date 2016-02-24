module test_helper

  use kinddefs, only : dp

  implicit none
  public

  integer :: ndimn, ntype, nelem, npoin, nface, nnode

  real(dp), dimension(:,:),     allocatable :: coord, geoel, rface
  
  integer, dimension(:,:),  allocatable :: inpoel, bcface

  character(len=20) :: gridfile = "testgrid"

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp

contains

  subroutine setup_cube
    use gridtools, only : readgrid_lou, basis_function, face_norm
  continue
    nnode = 3
  ! Read the grid
    call readgrid_lou(ndimn,ntype,nelem,npoin,nface,nnode,inpoel,gridfile,       &
                    coord,bcface)
  ! assemble geometry related matricies
    call basis_function(nelem,geoel,inpoel,coord)
    call face_norm(rface,coord,bcface,nface,ndimn)
  end subroutine

end module test_helper
