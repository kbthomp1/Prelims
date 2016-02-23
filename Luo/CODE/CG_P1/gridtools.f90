module gridtools

  use kinddefs, only : dp

  implicit none
  private

  public :: readgrid_lou
  public :: basis_function
  public :: face_norm

contains

!=============================== READGRID_LOU ===============================80
! Read in grid file that is in Dr. Lou's format
!============================================================================80
  subroutine readgrid_lou(ndimn,ntype,nelem,npoin,nface,nvertices,inpoel,     &
                          gridfile,coord,bcface)
  
    integer      :: titleline,ndimn,ntype,nelem,npoin,nface,nvertices
    integer      :: i,j,elem_dummy,poin_dummy,face_dummy
    
    real(dp), dimension(:,:),allocatable  :: coord
    
    integer, dimension(:,:), allocatable :: inpoel,bcface
    
    character(len=100) :: gridfile

  continue
    
    open(12,file=gridfile,status="old")
    read(12,*) titleline
    
    do i=1,titleline
      read(12,*)
    end do
    
    read(12,*) 
    read(12,*) ndimn,ntype
    read(12,*)
    read(12,*) nelem,npoin,nface
    read(12,*)
    
    allocate(inpoel(3,nelem))
    allocate(coord(ndimn,npoin))
    allocate(bcface(ndimn+1,nface))
    
    do j=1,nelem
      read(12,*) elem_dummy,(inpoel(i,j),i=1,nvertices)
    end do
    
    read(12,*)
    
    do j=1,npoin
      read(12,*) poin_dummy,(coord(i,j),i=1,ndimn)
    end do
    
    do j=1,npoin+2
      read(12,*)
    end do
    
    do j=1,nface
      read(12,*) face_dummy,(bcface(i,j),i=1,ndimn+1)
    end do
  
  end subroutine

!=========================== BASIS_FUNCTION =================================80
! Compute the 2D finite element basis functions for triangles
!
! B_i = (a_i*x + b_i*y + c_i)/D
!
! a_i =  (y_j - y_k)
! b_i = -(x_j - x_k)
! c_i =  (x_j*y_k - x_k*y_j)
!
! D = 2*(area of element)   !NOTE: this is twice the area of the element
!   = a_1*b_2 - a_2*b_1
!   = a_2*b_3 - a_3*b_2   All of these
!   = a_3*b_1 - a_1*b_3   are equivalent
!   = c_1 + c_2 + c_3
!
!============================================================================80
  subroutine basis_function(nelem,geoel,inpoel,coord)
 
    integer,                 intent(in) :: nelem
    integer, dimension(:,:), intent(in) :: inpoel

    real(dp), dimension(:,:),              intent(in)  :: coord
    real(dp), dimension(:,:), allocatable, intent(out) :: geoel
    
    real(dp), dimension(3) :: a,b

    real(dp) :: D,dd
    integer  :: ielem,ip1,ip2,ip3
    
    allocate(geoel(5,nelem))
    
    do ielem=1,nelem
      
      ip1  = inpoel(1,ielem)
      ip2  = inpoel(2,ielem)
      ip3  = inpoel(3,ielem)
    
      a(1) = coord(2,ip2)-coord(2,ip3)
      a(2) = coord(2,ip3)-coord(2,ip1)
    
      b(1) = -(coord(1,ip2)-coord(1,ip3))
      b(2) = -(coord(1,ip3)-coord(1,ip1))
    
      D    = a(1)*b(2)-a(2)*b(1)
      dd   = 1/D

      !NOTE: this is twice the area of the element
      geoel(5,ielem) = D
    
      a(1) = a(1)*dd
      a(2) = a(2)*dd
      b(1) = b(1)*dd
      b(2) = b(2)*dd
   
      !Coefficients of the basis function
      geoel(1,ielem) = a(1)
      geoel(2,ielem) = a(2)
      geoel(3,ielem) = b(1)
      geoel(4,ielem) = b(2)
    
    end do
    
    end subroutine
  
!=========================== SET_BC =========================================80
! Calculate normal vector at each interface
!============================================================================80
  subroutine face_norm(rface,coord,bcface,nface,ndimn)

    integer,                 intent(in) :: nface, ndimn
    integer, dimension(:,:), intent(in) :: bcface

    real(dp), dimension(:,:),               intent(in) :: coord
    real(dp), dimension(:,:), allocatable, intent(out) :: rface
    
    integer :: ip1,ip2,iface

    real(dp) :: x1,x2,y1,y2

  continue
    
    allocate(rface(ndimn,nface))
    
    do iface=1,nface
      ip1  = bcface(1,iface)
      ip2  = bcface(2,iface)
      
      x1   = coord(1,ip1)
      x2   = coord(1,ip2)
    
      y1   = coord(2,ip1)
      y2   = coord(2,ip2)
      
      rface(1,iface) = y2-y1
      rface(2,iface) = x1-x2
    end do
  
  end subroutine

end module gridtools
