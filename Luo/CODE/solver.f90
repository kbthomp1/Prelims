module solver

  implicit none
  private

  public :: get_lhspo

contains

function get_lhspo(npoin,nelem,inpoel,geoel) result (lhspo)

  integer, intent(in) :: npoin, nelem
  integer,  dimension(:,:), intent(in) :: inpoel
  real(dp), dimension(:,:), intent(in) :: geoel

  real(dp), dimension(npoin,npoin), intent(out) :: lhspo

  integer :: ielem, ip1, ip2,ip3
  integer :: 

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

end module solver
