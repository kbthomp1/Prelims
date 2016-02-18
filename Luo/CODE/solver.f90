module solver

  use kinddefs, only : dp

  implicit none
  private

  public :: get_lhspo
  public :: set_bc

contains

  function get_lhspo(npoin,nelem,nnode,inpoel,geoel) result (lhspo)
  
    integer, intent(in) :: npoin, nelem, nnode
    integer,  dimension(:,:), intent(in) :: inpoel
    real(dp), dimension(:,:), intent(in) :: geoel
  
    real(dp), dimension(npoin,npoin) :: lhspo
  
    real(dp), dimension(3) :: bx, by
  
    real(dp) :: rjac
  
    integer :: ielem, ip1, ip2,ip3
    integer :: i, ip, j, jp
  
  continue
  
    lhspo(1:npoin,1:npoin) = 0.0
     
    do ielem=1,nelem
    
      ip1   = inpoel(1,ielem)
      ip2   = inpoel(2,ielem)
      ip3   = inpoel(3,ielem)
      
      bx(1) = geoel(1,ielem)
      bx(2) = geoel(2,ielem)
      bx(3) = -(bx(1)+bx(2))
      by(1) = geoel(3,ielem)
      by(2) = geoel(4,ielem)
      by(3) = -(by(1)+by(2))
    
      rjac = 0.5*geoel(5,ielem)
    
      do i=1,nnode
        ip = inpoel(i,ielem)
        do j=1,nnode
          jp = inpoel(j,ielem)
          lhspo(ip,jp) = lhspo(ip,jp) + (bx(i)*bx(j) + by(i)*by(j))*rjac
        end do
      end do
    
    end do
  
  end function get_lhspo
  
!=========================== SET_BC =========================================80
!                        IMPOSE DIRCHLET BC
!============================================================================80
  subroutine set_bc(phi,lhspo,rhspo,npoin,bcface)
 
    integer,                 intent(in) :: npoin
    integer, dimension(:,:), intent(in) :: bcface
  
    real(dp), dimension(npoin), intent(out)   :: phi
    real(dp), dimension(:,:),   intent(inout) :: lhspo
    real(dp), dimension(:),     intent(inout) :: rhspo
  
    real(dp), parameter :: phi_ib = 1.0_dp

    integer :: ib
  
  continue
  
    phi(1:npoin) = 0.0_dp
    
    ib = bcface(1,1)
    
    lhspo(ib,ib) = lhspo(ib,ib)*1.0e+20_dp
    rhspo(ib)   = lhspo(ib,ib)*phi_ib
  
  end subroutine

end module solver
