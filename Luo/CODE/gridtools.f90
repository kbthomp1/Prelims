module gridtools

  use kinddefs, only : dp

  implicit none
  private

  public :: readgrid_lou
  public :: basis_function
  public :: set_bc
  public :: getesuel
  public :: getesup
  public :: getnormal
  public :: getintfac
  public :: preprocess_grid
  public :: gridtype
  public :: get_global_dof

  type :: gridtype

    integer :: ndimn   ! Number of spatial dimensions
    integer :: ntype   ! Type of element selector
    integer :: nelem   ! Number of elements
    integer :: npoin   ! Number of points
    integer :: nnode   ! Number of nodes per element
    integer :: nfael   ! Number of faces per element
    integer :: numfac  ! Number of interior faces
    integer :: nbcface ! Number of boundary faces
    integer :: nface   ! Number of faces total

    real(dp), dimension(:,:), allocatable :: coord  ! (x,y) coordinates
    real(dp), dimension(:,:), allocatable :: geoel  ! basis function info
    real(dp), dimension(:,:), allocatable :: del    ! interior face normal
  
    integer, dimension(:),    allocatable :: esup1, esup2 ! helper arrays
    integer, dimension(:,:),  allocatable :: gdof   ! dofs on a face
    integer, dimension(:,:),  allocatable :: inpoel ! points on an element
    integer, dimension(:,:),  allocatable :: bcface ! boundary face nodes
    integer, dimension(:,:),  allocatable :: esuel  ! helper array
    integer, dimension(:,:),  allocatable :: intfac ! interior face nodes

  end type gridtype

contains

  subroutine preprocess_grid(grid,gridfile,nnode)

    type(gridtype), intent(inout) :: grid

    integer, intent(in) :: nnode

    character(len=*), intent(in) :: gridfile

  continue

    grid%nnode = nnode
    grid%nfael = nnode

    call readgrid_lou(grid,gridfile)

    call setup_grid_arrays(grid)

  ! assemble geometry related matricies
    call basis_function(grid%nelem, grid%geoel, grid%inpoel, grid%coord)

    call getesup(grid%inpoel, grid%npoin, grid%nelem,  &
                 grid%nnode,  grid%esup1, grid%esup2)

    call getesuel(grid%inpoel, grid%esup1, grid%esup2, &
                  grid%nfael,  grid%nelem, grid%npoin, &
                  grid%esuel)

    call getintfac(grid%esuel,  grid%inpoel, grid%nelem, &
                   grid%nfael,  grid%numfac, grid%nbcface, &
                   grid%intfac, grid%nface)

    call set_bc(grid%intfac, grid%bcface, grid%nbcface)

    call getnormal(grid%intfac, grid%coord, grid%nface, &
                   grid%del)

    call map_gdof(grid)

  end subroutine preprocess_grid


!=============================== READGRID_LOU ===============================80
! Read in grid file that is in Dr. Lou's format
!============================================================================80
  subroutine readgrid_lou(grid,gridfile)

    type(gridtype), intent(inout) :: grid

    integer :: titleline
    integer :: i,j,elem_dummy,poin_dummy,face_dummy
    
    character(len=*) :: gridfile

  continue
    
    open(12,file=gridfile,status="old")
    read(12,*) titleline
    
    do i=1,titleline
      read(12,*)
    end do
    
    read(12,*) 
    read(12,*) grid%ndimn, grid%ntype
    read(12,*)
    read(12,*) grid%nelem, grid%npoin, grid%nbcface
    read(12,*)
    
    allocate(grid%inpoel(3,grid%nelem))
    allocate(grid%coord(grid%ndimn,grid%npoin))
    allocate(grid%bcface(grid%ndimn+1,grid%nbcface))
    
    do j=1,grid%nelem
      read(12,*) elem_dummy,(grid%inpoel(i,j),i=1,grid%nnode)
    end do
    
    read(12,*)
    
    do j=1,grid%npoin
      read(12,*) poin_dummy,(grid%coord(i,j),i=1,grid%ndimn)
    end do
    
    do j=1,grid%npoin+2
      read(12,*)
    end do
    
    do j=1,grid%nbcface
      read(12,*) face_dummy,(grid%bcface(i,j),i=1,grid%ndimn+1)
    end do
 
    close(12)
  end subroutine

!======================= SETUP_GRID_ARRAYS ==================================80
! Allocate storage to grid arrays
!============================================================================80
  subroutine setup_grid_arrays(grid)

    type(gridtype), intent(inout) :: grid

  continue

    allocate(grid%geoel(5,grid%nelem))
    allocate(grid%esup1(grid%npoin*6))
    allocate(grid%esup2(grid%npoin+1))
    allocate(grid%esuel(grid%nfael,grid%nelem))
    
  end subroutine setup_grid_arrays

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

    real(dp), dimension(:,:), intent(in)  :: coord
    real(dp), dimension(:,:), intent(out) :: geoel
    
    real(dp), dimension(3) :: a,b

    real(dp) :: D,dd
    integer  :: ielem,ip1,ip2,ip3

  continue

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
! Set the boundary condition flag in the intfac array
!============================================================================80
  subroutine set_bc(intfac, bcface, nbcface)

    integer,                 intent(in)    :: nbcface
    integer, dimension(:,:), intent(in)    :: bcface
    integer, dimension(:,:), intent(inout) :: intfac
    
    integer :: bc1, bc2, f1, f2, iface, ibc

  continue
    
    loop_bcfaces: do ibc = 1, nbcface
      bc1 = bcface(1,ibc)
      bc2 = bcface(2,ibc)
      loop_intfac: do iface = 1, nbcface
        f1  = intfac(3,iface)
        f2  = intfac(4,iface)
        if ( (bc1 == f1 .and. bc2 == f2) .or. (bc1 == f2 .and. bc2 == f1) ) then
          intfac(5,iface) = bcface(3,ibc)
          exit
        end if
      end do loop_intfac
    end do loop_bcfaces
  
  end subroutine

!=============================================================================

  subroutine getesup(inpoel,npoin,nelem,nnode,esup1,esup2)
 
    integer, intent(in) :: npoin, nelem, nnode

    integer, dimension(:,:),     intent(in)  :: inpoel
    integer, dimension(npoin+1), intent(out) :: esup2
    integer, dimension(npoin*6), intent(out) :: esup1
    
    integer :: ipoin, ielem, inode, istor

  continue
  
    !initialize esup2
    esup2 = 0
    
    !1st pass
    do ielem=1,nelem
      do inode=1,nnode
        esup2(inpoel(inode,ielem)+1)=esup2(inpoel(inode,ielem)+1)+1
      end do
    end do
    
    do ipoin=2,npoin+1
      esup2(ipoin)=esup2(ipoin)+esup2(ipoin-1)
    end do
    
    !2nd pass
    do ielem=1,nelem
      do inode=1,nnode
        ipoin=inpoel(inode,ielem)
        istor=esup2(ipoin)+1
        esup2(ipoin)=istor
        esup1(istor)=ielem
      end do
    end do
    
    do ipoin=npoin+1,2,-1
      esup2(ipoin)=esup2(ipoin-1)
    end do
    
    esup2(1)=0
  
  end subroutine
  
  !===========================================================================
  
  subroutine getesuel(inpoel,esup1,esup2,nfael,nelem,npoin,esuel)

    integer, intent(in) :: nfael, nelem, npoin

    integer, dimension(:,:), intent(in)  :: inpoel
    integer, dimension(:),   intent(in)  :: esup1,esup2
    
    integer, dimension(nfael,nelem), intent(out) :: esuel

    integer, dimension(2)     :: lhelp, lhelp2
    integer, dimension(npoin) :: lpoin
    
    integer :: ifael,jfael,ielem,jelem,ipoin,jpoin
    integer :: jnofa,istor,icoun

  continue
  
    lpoin(1:npoin)=0
    esuel(1:nfael,1:nelem)=0
    
    do ielem=1,nelem
      do ifael=1,nfael
        if(ifael<nfael) then
          lhelp(1:2)=inpoel(ifael:ifael+1,ielem)
        else
          lhelp(1)=inpoel(ifael,ielem)
          lhelp(2)=inpoel(1,ielem)
        end if
        lpoin(lhelp(1:2))=1
        ipoin=lhelp(1)
    
        do istor=esup2(ipoin)+1,esup2(ipoin+1)
          jelem=esup1(istor)
          if(jelem/=ielem) then
            do jfael=1,nfael
              if(jfael<nfael) then
                lhelp2(1:2)=inpoel(jfael:jfael+1,jelem)
              else
                lhelp2(1)=inpoel(jfael,jelem)
                lhelp2(2)=inpoel(1,jelem)
              end if
    
              icoun=0
              do jnofa=1,2
                jpoin=lhelp2(jnofa)
                if(lpoin(jpoin)==1) then
                  icoun=icoun+1
                end if
              end do
              if(icoun==2) then
                esuel(ifael,ielem)=jelem
                esuel(jfael,jelem)=ielem
              end if
            end do
          end if
        end do
        lpoin(lhelp(1:2))=0
      end do
    end do
    
  end subroutine
  
  !=====================================================================
  
  subroutine getintfac(esuel,inpoel,nelem,nfael,numfac,nbcface,intfac,  &
                       nface)
 
    integer, intent(in)  :: nelem, nfael, nbcface
    integer, intent(out) :: numfac, nface

    integer, dimension(:,:),              intent(in)  :: esuel,inpoel
    integer, dimension(:,:), allocatable, intent(out) :: intfac

    integer, dimension(4,nbcface) :: intbcfac
  
    integer :: ielem, jelem, ifael, iface
    integer :: ier, iel, ip1, ip2

  continue
  
    numfac=0
    iface=1
    
    do ielem=1,nelem
      do ifael=1,nfael
        iel=ielem
        ier=esuel(ifael,ielem)
        if(iel<ier) numfac=numfac+1
      end do
    end do
    
    do ielem=1,nelem
      do ifael=1,nfael
        jelem=esuel(ifael,ielem)
    
        if(ifael<nfael) then
          ip1=inpoel(ifael,ielem)
          ip2=inpoel(ifael+1,ielem)
        else
          ip1=inpoel(ifael,ielem)
          ip2=inpoel(1,ielem)
        end if
    
        if(jelem==0) then
          intbcfac(1,iface)=nelem+iface
          intbcfac(2,iface)=ielem
          intbcfac(3,iface)=ip1
          intbcfac(4,iface)=ip2
          iface=iface+1
        end if
     
      end do
    end do

    nface = nbcface+numfac
    
    allocate(intfac(5,nface))
    intfac(1:4,1:nbcface)=intbcfac(1:4,1:nbcface)
    intfac(5,1:(nface))=0
    
    do ielem=1,nelem
      do ifael=1,nfael
        iel=ielem
        ier=esuel(ifael,ielem)
        if(ifael<nfael) then
          ip1=inpoel(ifael,ielem)
          ip2=inpoel(ifael+1,ielem)
        else
          ip1=inpoel(ifael,ielem)
          ip2=inpoel(1,ielem)
        end if
        if(iel<ier) then
          intfac(1,iface)=iel
          intfac(2,iface)=ier
          intfac(3,iface)=ip1
          intfac(4,iface)=ip2
          iface=iface+1
        end if
      end do
    end do
       
  end subroutine
       
  !========================================================================
  
  subroutine getnormal(intfac,coord,nface,del)

    integer, intent(in) :: nface

    integer, dimension(:,:), intent(in) :: intfac

    real(dp), dimension(:,:), intent(in) :: coord
  
    real(dp),dimension(:,:), allocatable, intent(out) :: del
  
    real(dp) :: x1,x2,y1,y2

    integer :: iface,ip1,ip2

  continue

    allocate(del(3,nface))
  
    do iface=1,nface

      ip1=intfac(3,iface)
      ip2=intfac(4,iface)
      
      x1=coord(1,ip1)
      x2=coord(1,ip2)
      
      y1=coord(2,ip1)
      y2=coord(2,ip2)
    
      del(1,iface)=y2-y1
      del(2,iface)=-(x2-x1)
      del(3,iface)=sqrt(del(1,iface)**2+del(2,iface)**2)
    
    end do
  
  end subroutine

  !========================================================================
  
  subroutine map_gdof(grid)
    type(gridtype), intent(inout) :: grid
    integer :: ielem, ip, i
  continue
    allocate(grid%gdof(grid%npoin,grid%nelem))
    grid%gdof = 0
    do ielem = 1, grid%nelem
      do i = 1, grid%nnode
        ip = grid%inpoel(i,ielem)
        grid%gdof(ip,ielem) = (ielem-1)*3 + i
      end do
    end do
  end subroutine map_gdof

!============================ GET_GLOBAL_DOF ================================80
! Get the global dof number in solution
!============================================================================80
  function get_global_dof(node,elem,grid) result(global_dof)
    type(gridtype), intent(in) :: grid
    integer,        intent(in) :: node,elem
    integer :: global_dof
  continue
    global_dof = grid%gdof(node,elem)
    if ( global_dof == 0 ) then
      write(*,20) "Error: could not match node", node," on element", elem
      stop
    end if
20 format(2(A,x,i5))
  end function

end module gridtools
