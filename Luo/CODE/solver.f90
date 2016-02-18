module solver

  use kinddefs, only : dp

  implicit none
  private

  public :: get_lhspo
  public :: set_bc
  public :: get_soln

contains

!================================ GET_LHSPO =================================80
! Form the global stiffness matrix
!============================================================================80

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
! Set the Dirchlet boundary condition
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

!============================ GET_SOLN ======================================80
! Calculate velocities Vx, Vy, and Vt by interpolating solution to nodes
!============================================================================80

  subroutine get_soln(Vx,Vy,Vt,phi,geoel,inpoel,npoin,nelem)

    integer,                     intent(in) :: npoin,nelem
    integer, dimension(3,npoin), intent(in) :: inpoel

    real(dp), dimension(5,nelem), intent(in)  :: geoel
    real(dp), dimension(npoin),   intent(in)  :: phi
    real(dp), dimension(npoin),   intent(out) :: Vx, Vy, Vt

    real(dp), dimension(nelem) :: Vx_local, Vy_local
    real(dp), dimension(npoin) :: Vxarea, Vyarea, area

    integer :: ielem, ipoin, ip1, ip2, ip3

    real(dp), parameter :: zero = 0.0_dp

  continue

    Vx_local(1:nelem) = zero
    Vy_local(1:nelem) = zero
    Vxarea(1:npoin)   = zero
    Vyarea(1:npoin)   = zero
    area(1:npoin)     = zero

    do ielem=1,nelem

      ip1=inpoel(1,ielem)
      ip2=inpoel(2,ielem)
      ip3=inpoel(3,ielem)

      Vx_local(ielem) = (geoel(1,ielem)*(phi(ip1)-phi(ip3))       &
                       + geoel(2,ielem)*(phi(ip2)-phi(ip3)))*geoel(5,ielem)

      Vy_local(ielem) = (geoel(3,ielem)*(phi(ip1)-phi(ip3))       &
                       + geoel(4,ielem)*(phi(ip2)-phi(ip3)))*geoel(5,ielem)

      Vxarea(ip1) = Vxarea(ip1)+Vx_local(ielem)
      Vxarea(ip2) = Vxarea(ip2)+Vx_local(ielem)
      Vxarea(ip3) = Vxarea(ip3)+Vx_local(ielem)

      Vyarea(ip1) = Vyarea(ip1)+Vy_local(ielem)
      Vyarea(ip2) = Vyarea(ip2)+Vy_local(ielem)
      Vyarea(ip3) = Vyarea(ip3)+Vy_local(ielem)

      area(ip1)  = area(ip1)+geoel(5,ielem)
      area(ip2)  = area(ip2)+geoel(5,ielem)
      area(ip3)  = area(ip3)+geoel(5,ielem)

    end do

    do ipoin=1,npoin

      Vx(ipoin) = Vxarea(ipoin)/area(ipoin)
      Vy(ipoin) = Vyarea(ipoin)/area(ipoin)

      Vt(ipoin) = sqrt(Vx(ipoin)**2 + Vy(ipoin)**2)

    end do
  end subroutine get_soln

end module solver
