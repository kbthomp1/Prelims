module solver

  use kinddefs, only : dp

  implicit none
  private

  public :: get_lhspo
  public :: get_rhspo
  public :: set_bc
  public :: get_soln
  public :: solve

contains

!============================= GET_RHSPO ====================================80
! Formulate the RHS load vector
!============================================================================80
  function get_rhspo(bcface,rface,nface,npoin,uinf,vinf) result(rhspo)
  
    integer,                 intent(in) :: nface,npoin
    integer, dimension(:,:), intent(in) :: bcface
  
    real(dp), dimension(:,:), intent(in) :: rface
    real(dp),                 intent(in) :: uinf,vinf
  
    real(dp), dimension(npoin) :: rhspo
    real(dp) :: cface

    integer :: ip1, ip2, iface
    
  continue 
  
    rhspo(1:npoin) = 0.0
    
    do iface=1,nface
      if(bcface(3,iface) == 4) then
        
        ip1   = bcface(1,iface)
        ip2   = bcface(2,iface)
        
    !    itype = bcface(4,iface)
    !
    !    roinf = uchar(1,itype)
    !     uinf = uchar(2,itype)
    !     vinf = uchar(3,itype)
    !     pinf = uchar(4,itype)
         
        cface = 0.5*(uinf*rface(1,iface) + vinf*rface(2,iface))
    
        rhspo(ip1) = rhspo(ip1) + cface
        rhspo(ip2) = rhspo(ip2) + cface
    
      end if
    end do
  
  end function get_rhspo

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

    integer,                  intent(in) :: npoin,nelem
    integer, dimension(:,:),  intent(in) :: inpoel

    real(dp), dimension(:,:), intent(in)  :: geoel
    real(dp), dimension(:),   intent(in)  :: phi
    real(dp), dimension(:),   intent(out) :: Vx, Vy, Vt

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

!============================== SOLVE =======================================80
! Solve the linear system per user's specified method
!============================================================================80
  subroutine solve(lin_solver,lhspo,rhspo,phi,npoin,nsteps,tolerance)

    use linalg, only : gauss_seidel,jacobi,conjgrad

    character(len=*), intent(in) :: lin_solver

    integer, intent(in) :: npoin, nsteps

    real(dp), dimension(npoin),       intent(in)    :: rhspo
    real(dp), dimension(npoin),       intent(inout) :: phi
    real(dp), dimension(npoin,npoin), intent(in)    :: lhspo

    real(dp), intent(in) :: tolerance

  continue

    select case(lin_solver)
    case("conjugate gradient")
      call conjgrad(lhspo,rhspo,phi,npoin,nsteps)
    case("point jacobi")
      call jacobi(lhspo,rhspo,phi,npoin,nsteps)
    case("gauss seidel")
      call gauss_seidel(lhspo,rhspo,phi,npoin,nsteps,tolerance)
    case default
      write(*,*) "Error: linear solver specified not available:",lin_solver
      stop 1
    end select

  end subroutine solve

end module solver
