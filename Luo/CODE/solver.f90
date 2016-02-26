module solver

  use kinddefs,  only : dp
  use gridtools, only : gridtype


  implicit none
  private

  public :: get_lhspo
  public :: get_rhspo
  public :: get_soln
  public :: solve

  !unit tests
  public :: add_flux_contributions
  public :: add_lift_primal_domain
  public :: add_lift_second_domain
  public :: add_jump_second_face
  public :: add_flux_primal_face

  real(dp), parameter :: my_4th = 0.25_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 1.0_dp

contains

!============================= GET_RHSPO ====================================80
! Formulate the RHS load vector with neumann BC's
!============================================================================80
  function get_rhspo(grid,ndof,uinf,vinf) result(rhspo)

    type(gridtype), intent(in) :: grid

    integer, intent(in) :: ndof

    real(dp), intent(in) :: uinf,vinf
  
    real(dp), dimension(ndof) :: rhspo

    real(dp) :: cface

    integer :: ip1, ip2, iface
    
  continue 
  
    rhspo = zero
    
    do iface=1,grid%nface
      if(grid%bcface(3,iface) == 4) then
        
        ip1 = grid%bcface(1,iface)
        ip2 = grid%bcface(2,iface)
        
        cface = 0.5_dp*(uinf*grid%rface(1,iface) + vinf*grid%rface(2,iface))
    
        rhspo(ip1) = rhspo(ip1) + cface
        rhspo(ip2) = rhspo(ip2) + cface
    
      end if
    end do
  
  end function get_rhspo

!================================ GET_LHSPO =================================80
! Form the global stiffness matrix
!============================================================================80

  function get_lhspo(grid,ndof) result (lhspo)

    integer, intent(in) :: ndof

    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof,ndof) :: lhspo

    real(dp), dimension(3) :: bx, by

    real(dp) :: area, D

    integer :: ielem, ip1, ip2,ip3
    integer :: i, ip, j, jp

    real(dp), parameter :: fact = one/12._dp

  continue

  lhspo(1:ndof,1:ndof) = 0.0

    do ielem=1,grid%nelem

      ip1   = grid%inpoel(1,ielem)
      ip2   = grid%inpoel(2,ielem)
      ip3   = grid%inpoel(3,ielem)

      bx(1) = grid%geoel(1,ielem)
      bx(2) = grid%geoel(2,ielem)
      bx(3) = -(bx(1)+bx(2))
      by(1) = grid%geoel(3,ielem)
      by(2) = grid%geoel(4,ielem)
      by(3) = -(by(1)+by(2))

      D = grid%geoel(5,ielem)
      area = half*D

      do i=1,grid%nnode
        ip = grid%inpoel(i,ielem)
        do j=1,grid%nnode
          jp = grid%inpoel(j,ielem)
          lhspo(ip,jp) = lhspo(ip,jp) + (bx(i)*bx(j) + by(i)*by(j))*area
        end do
      end do

    end do

    call add_flux_contributions(lhspo,grid,ndof)

  end function get_lhspo

!===================== ADD_LIFT_PRIMAL_DOMAIN ===============================80
! Add local lift functions to primal domain integral
!============================================================================80
  subroutine add_lift_primal_domain(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid
    real(dp), dimension(:,:), intent(inout) :: lhspo

    real(dp), dimension(3) :: bx, by
    integer,  dimension(2) :: cells
    integer  :: ilift, ip, i, cell, l
    real(dp) :: cell_area

  continue

    ilift = global_lift_coord(iface,grid)
    cells(1) = icell
    cells(2) = jcell

    do l = 1,2
      cell = cells(l)
      call get_basis(bx,by,grid,cell)
      cell_area = half*grid%geoel(5,cell)
      do i=1,grid%nnode
        ip = grid%inpoel(i,cell)
        lhspo(ip,ilift) = lhspo(ip,ilift) - (bx(i) + by(i))*cell_area
      end do
    end do

  end subroutine add_lift_primal_domain

!===================== ADD_LIFT_SECOND_DOMAIN ===============================80
! Add local lift to lift equation domain integral
!============================================================================80
  subroutine add_lift_second_domain(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    integer :: ilift

    real(dp) :: icell_area, jcell_area

  continue

    ilift = global_lift_coord(iface,grid)
    icell_area = half*grid%geoel(5,icell)
    jcell_area = half*grid%geoel(5,jcell)
    lhspo(ilift,ilift) = lhspo(ilift,ilift) + two*(icell_area+jcell_area)

  end subroutine add_lift_second_domain

!===================== ADD_LIFT_SECOND_DOMAIN ===============================80
! Add jump condition to lift equation
!============================================================================80
  subroutine add_jump_second_face(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    integer :: ilift, ip, jp, i

    real(dp) :: nx, ny, face_area, integrated_const
    real(dp) :: eta = 4._dp

  continue

    ilift = global_lift_coord(iface,grid)
    call get_face_normal(nx,ny,face_area,grid,iface)

    ! result of integral of jump function at iface
    integrated_const = my_4th*eta*(nx+ny)*face_area

    do i=1,grid%nnode
      ip = grid%inpoel(i,icell)
      jp = grid%inpoel(i,jcell)
      lhspo(ilift,jp) = lhspo(ilift,jp) + integrated_const
      lhspo(ilift,ip) = lhspo(ilift,ip) - integrated_const
    end do

  end subroutine add_jump_second_face

!======================= ADD_FLUX_PRIMAL_FACE ===============================80
! Add flux to primal equation at iface
!============================================================================80
  subroutine add_flux_primal_face(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    real(dp), dimension(3) :: bx, by, grad
    real(dp), dimension(3) :: phi_coefs
    real(dp), dimension(2) :: sgn

    integer, dimension(2) :: cells
    integer :: ilift, ip, jp, i, j, l, s
    integer :: cell, cell_idx

    real(dp) :: nx, ny, face_area, lift_coef

  continue

    ilift = global_lift_coord(iface,grid)
    cells(1) = icell
    cells(2) = jcell
    sgn(1)   = one   ! add flux to left cell
    sgn(2)   = -one  ! subtract flux from right cell
    call get_face_normal(nx,ny,face_area,grid,iface)

    flip_sign: do s = 1,2
      cell_idx = cells(s) !index of cell to add/subtract flux
      left_right_cells: do l = 1,2

        cell = cells(l) !index of cell to form contibutions
        ! Basis functions -> dB/dx, dB/dy
        call get_basis(bx,by,grid,cell)

        ! half of averaged gradient
        grad(1:3) = half*(bx(1:3)*nx + by(1:3)*ny)

        ! integrated coefficients of phi
        phi_coefs(1:3) = grad(1:3)*half*face_area

        ! integrated coefficient of local lift
        lift_coef = my_4th*face_area

        do i=1,grid%nnode
          ip = grid%inpoel(i,cell_idx)
          lhspo(ip,ilift) = lhspo(ip,ilift) + sgn(s)*lift_coef
          do j=1,grid%nnode
            jp = grid%inpoel(j,cell)
            lhspo(ip,jp) = lhspo(ip,jp) - sgn(s)*phi_coefs(j)
          end do
        end do

      end do left_right_cells
    end do flip_sign

  end subroutine add_flux_primal_face

!===================== ADD_FLUX_CONTRIBUTIONS ===============================80
! Add the flux contributions to the LHS matrix
!============================================================================80
  subroutine add_flux_contributions(lhspo,grid,ndof)

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof,ndof), intent(inout) :: lhspo

    integer :: iface, icell, jcell

  continue

    interior_faces: do iface = grid%nface+1, grid%nface+grid%numfac

      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      call add_lift_primal_domain(icell,jcell,iface,grid,lhspo)
      call add_lift_second_domain(icell,jcell,iface,grid,lhspo)
      call add_jump_second_face(icell,jcell,iface,grid,lhspo)
      call add_flux_primal_face(icell,jcell,iface,grid,lhspo)

    end do interior_faces

  end subroutine add_flux_contributions

!============================ GET_SOLN ======================================80
! Calculate velocities Vx, Vy, and Vt by interpolating solution to nodes
!
! Note that the gradients of the linear basis functions are constant across
! each element, thus we can use an area weighted average to recover the 
! the gradient at a node.  Since phi = sum(phi_i*B_i)
! 
!   phi_x = sum(phi_i*Bx_i) = phi_1*a_1 + phi_2*a_2 + phi_3*a_3
!   phi_y = sum(phi_i*By_i) = phi_1*b_1 + phi_2*b_2 + phi_3*b_3
! 
! NOTE: since B_1 + B_2 + B_3 = 1, a/b_3 = 1 - a/b_2 - a/b_1 
!============================================================================80

  subroutine get_soln(Vx,Vy,Vt,phi,grid)

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in)  :: phi
    real(dp), dimension(:),   intent(out) :: Vx, Vy, Vt

    real(dp), dimension(grid%npoin) :: Vxarea, Vyarea, area

    integer :: ielem, ipoin, ip1, ip2, ip3

    real(dp) :: local_area, Vx_local, Vy_local

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: half = 0.5_dp

  continue

    Vx_local = zero
    Vy_local = zero
    Vxarea   = zero
    Vyarea   = zero
    area     = zero

    do ielem=1,grid%nelem

      ip1=grid%inpoel(1,ielem)
      ip2=grid%inpoel(2,ielem)
      ip3=grid%inpoel(3,ielem)

      local_area = half*grid%geoel(5,ielem)

      Vx_local = (grid%geoel(1,ielem)*(phi(ip1)-phi(ip3))       &
                + grid%geoel(2,ielem)*(phi(ip2)-phi(ip3)))      &
                * local_area

      Vy_local = (grid%geoel(3,ielem)*(phi(ip1)-phi(ip3))       &
                + grid%geoel(4,ielem)*(phi(ip2)-phi(ip3)))      &
                * local_area

      Vxarea(ip1) = Vxarea(ip1)+Vx_local
      Vxarea(ip2) = Vxarea(ip2)+Vx_local
      Vxarea(ip3) = Vxarea(ip3)+Vx_local

      Vyarea(ip1) = Vyarea(ip1)+Vy_local
      Vyarea(ip2) = Vyarea(ip2)+Vy_local
      Vyarea(ip3) = Vyarea(ip3)+Vy_local

      area(ip1)  = area(ip1)+local_area
      area(ip2)  = area(ip2)+local_area
      area(ip3)  = area(ip3)+local_area

    end do

    do ipoin=1,grid%npoin

      Vx(ipoin) = Vxarea(ipoin)/area(ipoin)
      Vy(ipoin) = Vyarea(ipoin)/area(ipoin)

      Vt(ipoin) = sqrt(Vx(ipoin)**2 + Vy(ipoin)**2)

    end do
  end subroutine get_soln

!============================== SOLVE =======================================80
! Solve the linear system per user's specified method
!============================================================================80
  subroutine solve(lin_solver,lhspo,rhspo,phi,ndof,nsteps,tolerance)

    use linalg, only : gauss_seidel,jacobi,conjgrad

    character(len=*), intent(in) :: lin_solver

    integer, intent(in) :: nsteps, ndof

    real(dp), dimension(ndof),      intent(in)    :: rhspo
    real(dp), dimension(ndof),      intent(inout) :: phi
    real(dp), dimension(ndof,ndof), intent(in)    :: lhspo

    real(dp), intent(in) :: tolerance

  continue

    select case(lin_solver)
    case("conjugate gradient")
      call conjgrad(lhspo,rhspo,phi,ndof,nsteps)
    case("point jacobi")
      call jacobi(lhspo,rhspo,phi,ndof,nsteps,tolerance)
    case("gauss seidel")
      call gauss_seidel(lhspo,rhspo,phi,ndof,nsteps,tolerance)
    case default
      write(*,*) "Error: linear solver specified not available:",lin_solver
      stop 1
    end select

    !phi(1,:) = x

  end subroutine solve

  subroutine get_basis(bx,by,grid,icell)
    integer,                intent(in)  :: icell
    type(gridtype),         intent(in)  :: grid
    real(dp), dimension(3), intent(out) :: bx, by
  continue
    bx(1:2) = grid%geoel(1:2,icell)
    bx(3)   = -(bx(1) - bx(2))
    by(1:2) = grid%geoel(3:4,icell)
    by(3)   = -(by(1) - by(2))
  end subroutine get_basis

  subroutine get_face_normal(nx,ny,face_area,grid,iface)
    integer,        intent(in) :: iface
    type(gridtype), intent(in) :: grid
    real(dp),      intent(out) :: nx, ny, face_area
  continue
    !compute unit normal vector
    nx = grid%del(1,iface)/grid%del(3,iface)
    ny = grid%del(2,iface)/grid%del(3,iface)
    face_area = grid%del(3,iface)
  end subroutine get_face_normal

  function global_lift_coord(iface,grid) result(ilift)
    integer,        intent(in)  :: iface
    type(gridtype), intent(in)  :: grid
    integer :: ilift
  continue
    ilift = grid%npoin + iface - grid%nface
  end function global_lift_coord

end module solver
